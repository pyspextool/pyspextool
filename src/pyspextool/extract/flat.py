import numpy as np
import numpy.typing as npt
from scipy import ndimage, interpolate
from scipy.signal import medfilt2d
from astropy.io import fits
import os
import matplotlib.pyplot as pl

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils.add_entry import add_entry
from pyspextool.utils.arrays import idl_rotate, idl_unrotate
from pyspextool.utils.math import moments
from pyspextool.utils.loop_progress import loop_progress

from pyspextool.fit.polyfit import polyfit_1d, poly_1d
from pyspextool.plot.plot_image import plot_image
from pyspextool.extract.wavecal import make_interp_indices_1d
from pyspextool.extract.images import rectify_order
from pyspextool.fit.fiterpolate import fiterpolate
from pyspextool.plot.limits import get_image_range


def find_top_bot(
    fcol: npt.ArrayLike,
    rownum: int,
    imgcol: npt.ArrayLike,
    sobcol: npt.ArrayLike,
    yguess,
    imgguess,
    frac: int | float,
    halfwin: int,
    debug: bool = False,
):
    """
    To identify the top and bottom of the slit in a given column.

    Called from within locate_orders as a subroutine

    Parameters
    ----------
    fcol :

    rownum :

    imgcol

    """

    # Moving up from the guess position, find the index of the first pixel in
    # the column that falls below the flux threshold, i.e. the edge of the slit

    ztop = np.where((imgcol < frac * imgguess) & (rownum > yguess))
    if len(ztop[0]) > 0:

        # Compute the index of the bottom and top of the center-of-mass window
        # being careful not to fall off of the array

        bidx = max(0, ztop[0][0] - halfwin)
        tidx = min(ztop[0][0] + halfwin, max(rownum))

        # Yank out the pixels to compute the center of mass

        yvals = rownum[bidx:tidx]
        zvals = sobcol[bidx:tidx]

        # Compute the center of mass

        com_top = np.sum(yvals * zvals) / np.sum(zvals)

    else:
        com_top = np.nan

    # Moving down from the guess position, find the index of the first pixel in
    # the column that falls below the flux threshold, i.e. the edge of the slit

    zbot = np.where((imgcol < frac * imgguess) & (rownum < yguess))
    if len(zbot[0]) > 0:

        # Compute the index of the bottom and top of the center-of-mass window
        # being careful not to fall off of the array

        bidx = max(0, zbot[0][-1] - halfwin)
        tidx = min(zbot[0][-1] + halfwin, max(rownum))

        # Yank out the pixels to compute the center of mass

        yvals = rownum[bidx:tidx]
        zvals = sobcol[bidx:tidx]

        # Compute the center of mass

        com_bot = np.sum(yvals * zvals) / np.sum(zvals)

    else:
        com_bot = np.nan

    if debug is True:

        plotfig2 = pl.figure(2)
        pl.plot(rownum, sobcol)

        pl.axvline(yguess, color="g")
        pl.axvline(com_top, color="r")
        pl.axvline(com_bot, color="r")
        pl.title(fcol)

        plotfig3 = pl.figure(3)
        pl.plot(rownum, imgcol)

        pl.axvline(yguess, color="g")
        pl.axvline(com_top, color="r")
        pl.axvline(com_bot, color="r")
        pl.title(fcol)
        plotfig2.clear()
        plotfig3.clear()

    return com_bot, com_top


def locate_orders(
    img: npt.ArrayLike,
    guess_positions: npt.ArrayLike,
    search_ranges: npt.ArrayLike,
    step_size: int,
    slit_height_range: npt.ArrayLike,
    poly_degree: int,
    ybuffer: int,
    intensity_fraction: int or float,
    com_width: int,
    qa_plotnumber: int = None,
    qa_figuresize: tuple = (7, 7),
    qa_fontsize: int = 12,
    qa_show: bool = None,
    qa_showscale: int | float = 1,
    qa_showblock: bool = False,
    qa_fullpath: str = None,
    debug: bool = False,
):
    """
    Locates orders in a (cross-dispersed) spectral image


    Parameters
    ----------
    img : ndarray
        An (nrows, ncols) image with (cross-dispersed) spectral orders.  It is
        assumed that the dispersion direction is roughly aligned with the
        rows of `img` and the spatial axis is roughly aligned with the
        columns of `img.  That is, orders go left-right and not up-down.

    guess_positions : ndarray
        An (norders,2) int array giving the positions to start the search.
        guesspos[0,0] gives the column number for the order closest to the
        bottom of `img` and guesspos[0,1] gives the row number for said
        order.  Typically, the positions are near the center of the image
        and center of the slit.

    search_ranges : ndarray
        An (norders,2) int array giving the column numbers over which to
        search.  search_ranges[0,0] gives the starting column number for
        the order closest to the bottom of `img` and search_ranges[0,1]
        gives the end column number for said order.

    step_size : int
        The step size in the dispersion or column direction

    slit_height_range: ndarray
        (2, ) array giving the minimum and maxmimum possible slit heights in
        units of pixels.

    poly_degree : int
        Polynomial fit degree for the edges of the orders

    ybuffer : int
        Number of pixels to buffer from the top and bottom of `img`
        during the search

    intensity_fraction : float
        The fraction of the flux of the center of the slit used to
        identify the location of the edge of the order

    com_width : int
        The window in units of pixels used to compute the center-of-mass (COM).

    qa_figuresize : tuple, default=(7,7)
        A (2,) tuple giving the plot size.

    qa_fullpath : str, optional
        The fullpath to write the QA file to disk.

    debug : bool, default=False
        If True, the function will plot the results of the order

    Returns
    -------
    (edgecoeffs, xranges) : tuple
    edgecoeffs : ndarray
        (norders, 2, ncoeffs) array giving the polynomial coefficients
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for
        the top of said order.

    xranges : ndarray
        An (norders, 2) array giving the columns where the orders fall entirely.
        xranges[0, 0] is the left column of the order nearest the bottom of the
        array and xranges[0, 1] is the right column of the same order.

    Notes
    -----
        The IDL sobel function is used to enhance the edges of the orders.
        Within an order, and at a particular column, the positions above
        and below the guess position at which the flux falls to `frac` the
        value at the guess position is determined.  Centroids are computed
        at these positions on the Sobeled image to determine the location
        of the edge of the order.  As the program steps away from the guess
        position, the new guess position is estimated using the previously
        determined edges.  After the entire array has been checked, the
        results are fit with a robust least-squares polynomial of degree
        `deg`.

    """

    #
    # Check parameters
    #

    check_parameter("locate_orders", "img", img, "ndarray", 2)

    check_parameter("locate_orders", "guess_positions", guess_positions, "ndarray", 2)

    check_parameter("locate_orders", "search_ranges", search_ranges, "ndarray", 2)

    check_parameter("locate_orders", "step_size", step_size, "int")

    check_parameter("locate_orders", "slit_height_range", slit_height_range, "ndarray")

    check_parameter("locate_orders", "poly_degree", poly_degree, "int")

    check_parameter("locate_orders", "ybuffer", ybuffer, "int")

    check_parameter("locate_orders", "intensity_fraction", intensity_fraction, "float")

    check_parameter("locate_orders", "com_width", com_width, "int")

    check_parameter("locate_orders", "qa_figuresize", qa_figuresize, "tuple")

    qa = check_qakeywords(show=qa_show, showblock=qa_showblock, showscale=qa_showscale)

    # Get set up to collect plotting info

    plguesspos = []
    plcols = []
    pledges = []
    pledgecoeffs = []
    plgoodbad = []

    # Get basic info and do basic things

    nrows, ncols = img.shape
    norders = len(guess_positions)

    rownum = np.arange(nrows)

    halfwin = int(com_width / 2.0)

    edgecoeffs = np.empty((norders, 2, poly_degree + 1))
    xranges = np.empty((norders, 2), dtype=int)

    # Sobel the image

    scl = np.max(img)
    simg1 = ndimage.sobel(img / scl, axis=0)
    simg2 = ndimage.sobel(img / scl, axis=1)
    simg = np.sqrt(simg1**2 + simg2**2)

    if debug is True:

        print('hi')
        pl.ion()
        # plotfig2 = pl.figure(2, figsize=(6, 6))  # sobel profile
        # plotfig3 = pl.figure(3, figsize=(6, 6))  # image profile

        minmax = get_image_range(simg, "zscale")
        if minmax[0] > minmax[1]:
            minmax = (np.min(simg), np.max(simg))

        # plotfig3 = pl.figure(4, figsize=(10, 10))  # sobel image

        cmap = pl.cm.gray
        pl.imshow(simg, vmin=minmax[0], vmax=minmax[1], cmap=cmap, origin="lower")
        pl.pause(1)
        return

    # Start looping over each order

    for i in range(0, norders):

        # Collect the guess positions for plotting

        plguesspos.append(list(guess_positions[i, :]))

        # Determine the start and stop column

        start = search_ranges[i, 0] + step_size - 1
        stop = search_ranges[i, 1] - step_size + 1

        # Generate an array of columns where you will find the slit

        fcols = np.arange(start, stop, step_size)  # "f"ind columns
        nfcols = len(fcols)

        # Create some empty arrays to fill

        edges = np.full((2, nfcols), np.nan)  # edges of slit
        cens = np.full(nfcols, np.nan)  # center point of slit

        # Fill the first few points in the cens array with the guess position

        delta = fcols - guess_positions[i, 0]
        gpidx = np.ndarray.item(
            np.where(np.absolute(delta) == np.absolute(delta).min())[0]
        )

        cens[(gpidx - poly_degree) : (gpidx + poly_degree)] = guess_positions[i, 1]
        cens_mask = ~np.isnan(cens)

        if debug is True:
            # plotfig3 = pl.figure(4)
            pl.plot(fcols, cens, "ro", markersize=0.5)

        #
        # Now move left from the guess position
        #

        for j in range(gpidx, 0, -1):

            imgcol = img[:, fcols[j]].ravel()  # the image column
            sobcol = simg[:, fcols[j]].ravel()  # the sobel image column

            # Fit the centers, so you can project away from the guess position

            r = polyfit_1d(
                fcols, cens, max(1, poly_degree - 2), justfit=True, silent=True
            )

            # Find the new guess position yguess

            yguess = np.polynomial.polynomial.polyval(fcols[j], r["coeffs"])

            if debug is True:

                # plotfig2 = pl.figure(3)
                pl.plot(fcols[j], yguess, "ro", markersize=0.5)

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess, ybuffer, (nrows - ybuffer - 1)))
            iguess = imgcol[yguess]

            if debug is True:
                print(f"yguess: {yguess}")

            bot, top = find_top_bot(
                fcols[j],
                rownum,
                imgcol,
                sobcol,
                yguess,
                iguess,
                intensity_fraction,
                halfwin,
                debug=debug,
            )

            if debug is True:
                print("Bot/Top", bot, top)

            # Confirm that both COMs were computed

            if np.isfinite(bot):

                # Now check to make sure calculated slit height falls within the
                # parameter slith_pix

                if debug is True:
                    print(f"fcols: {fcols}, fcols[j]: {fcols[j]}")
                    print("Slit height range", slit_height_range)
                    print("Top-Bot", top - bot)

                if slit_height_range[0] <= top - bot <= slit_height_range[1]:

                    # Store the results, update the cens array
                    edges[:, j] = np.array([bot, top])

                    z = np.where((fcols <= fcols[j]) & (cens_mask == 1))[0]
                    if z.size == 0:

                        cens[j] = (bot + top) / 2

                    else:

                        cens[z] = (bot + top) / 2

        #
        # Now move right from the guess position
        #

        for j in range(gpidx, nfcols, 1):

            imgcol = img[:, fcols[j]].ravel()  # the image column
            sobcol = simg[:, fcols[j]].ravel()  # the sobel image column

            # Fit the centers, so you can project away from the guess position

#            print(f"debug: {debug}")
            r = polyfit_1d(fcols, cens, max(1, poly_degree - 2), silent=not debug)

            # Find the new guess position yguess

            yguess = np.polynomial.polynomial.polyval(fcols[j], r["coeffs"])

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess, ybuffer, (nrows - ybuffer - 1)))
            iguess = imgcol[yguess]

            bot, top = find_top_bot(
                fcols[j],
                rownum,
                imgcol,
                sobcol,
                yguess,
                iguess,
                intensity_fraction,
                halfwin,
                debug=debug,
            )

            # Confirm that both COMs were computed

            if np.isfinite(bot + top):

                # Now check to make sure calculated slit height falls within the
                # parameter slith_pix

                if slit_height_range[0] <= top - bot <= slit_height_range[1]:
                    # Store the results, update the cens array

                    edges[:, j] = np.array([bot, top])
                    cens[j] = (bot + top) / 2

        # Now fit the results

        tmp = np.empty([2, poly_degree + 1])
        for j in range(0, 2):
            fit = polyfit_1d(
                fcols,
                edges[j, :],
                poly_degree,
                robust={"thresh": 4, "eps": 0.1},
                justfit=True,
                silent=not debug,
            )
            tmp[j, :] = fit["coeffs"]

            # Store the results for possible plotting

            plcols.append(fcols)
            pledges.append(edges[j, :])
            pledgecoeffs.append(fit["coeffs"])
            plgoodbad.append(fit["goodbad"])

        edgecoeffs[i, :, :] = tmp

        #
        # Confirm the orders fall on the arrays within search_ranges
        #

        xs = np.arange(search_ranges[i, 0], search_ranges[i, 1] + 1)
        top = poly_1d(xs, edgecoeffs[i, 1, :])
        z = top <= nrows - 1
        xranges[i, :] = [np.min(xs[z]), np.max(xs[z])]

    # Make the plotinfo dictionary

    plotinfo = {
        "guess_positions": plguesspos,
        "x": plcols,
        "y": pledges,
        "goodbad": plgoodbad,
        "edgecoeffs": pledgecoeffs,
    }

    if qa["show"] is True:

        plot_image(
            img,
            plot_number=qa_plotnumber,
            figure_size=(
                qa_figuresize[0] * qa_showscale,
                qa_figuresize[1] * qa_showscale,
            ),
            font_size=qa_fontsize * qa_showscale,
            showblock=qa_showblock,
            locateorders_plotinfo=plotinfo,
        )

    if qa_fullpath is not None:

        plot_image(
            img,
            figure_size=qa_figuresize,
            font_size=qa_fontsize,
            output_fullpath=qa_fullpath,
            locateorders_plotinfo=plotinfo,
        )

    return edgecoeffs, xranges


def normalize_flat(
    img: npt.ArrayLike,
    edgecoeffs: npt.ArrayLike,
    xranges: npt.ArrayLike,
    slith_arc: int | float,
    nxgrid: int,
    nygrid: int,
    var: npt.ArrayLike = None,
    oversamp: int = 1,
    ybuffer: int = 0,
    verbose: bool = False,
):
    """
    Normalize spectral flat field image using Tonry's fiterpolate routine


    Parameters
    ----------
    img : ndarray
        An (nrows, ncols) image with (cross-dispersed) spectral orders.
        It is assumed that the dispersion direction is roughly aligned
        with the rows of `img` and the spatial axis is roughly aligned
        with the columns of `img.  That is, orders go left-right and
        not up-down.

    edgecoeffs : ndarray
        (norders, ncoeffs+1, 2) array giving the polynomial coefficients
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for
        the top of said order.

    xranges : ndarray
        An (norders,2) array giving the column numbers over which to
        operate.  xranges[0,0] gives the starting column number for the
        order nearest the bottom of `img` and xranges[0,1] gives the end
        column number for said order.

    slith_arc : float
        The slit height in arcseconds

    nxgrid : int

    nygrid : int

    var : ndarray

    ybuffer : int, default 1, optional
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels.

    Returns
    -------
    ndarray, ndarray, ndarray

        nimg :  An image of the same size as `img` with the spectral ordered
        normalized to unity, and all other pixels set to unity.


        nvar :  A variance image of the same size as `img` with the spectral
            ordered normalized by the same model as `nimg`.  All other pixels
            are set to NaN.


        rms : (norders,) array giving the RMS value of the pixels


    Notes
    -----
    Resampling each order onto a rectangular grid.  Uses fiterpolate
    to create a surface model and then resamples the model back onto the
    original pixels.

    """

    #
    # Check parameters
    #

    check_parameter("normalize_flat", "img", img, "ndarray", 2)

    check_parameter("normalize_flat", "edgecoeffs", edgecoeffs, "ndarray")

    check_parameter("normalize_flat", "xranges", xranges, "ndarray")

    check_parameter("normalize_flat", "slith_arc", slith_arc, ["int", "float"])

    check_parameter("normalize_flat", "nxgrid", nxgrid, "int")

    check_parameter("normalize_flat", "nygrid", nygrid, "int")

    check_parameter("normalize_flat", "ybuffer", ybuffer, "int")

    #
    # Get basic info and do basic things
    #

    nrows, ncols = img.shape

    ndimen = edgecoeffs.ndim
    if ndimen == 2:
        norders = 1
    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    nimg = np.full((nrows, ncols), 1.0)
    nvar = np.full((nrows, ncols), np.nan) if var is not None else None
    rms = np.empty(norders)

    #
    # Start the loop
    #

    for i in range(norders):

        if verbose is True and i == 0:
            loop_progress(0, 0, norders, message="")

        #
        # Rectify the order
        #

        # Get the rectification indices

        xidx, yidx, wavemap, spatmap = make_interp_indices_1d(
            edgecoeffs[i, :, :], xranges[i, :], slith_arc
        )

        # Do the rectification

        order = rectify_order(img, xidx, yidx, ybuffer=ybuffer)

        # Fiterpolate the results after median smoothing to minimize bad pixels

        model = fiterpolate(
            medfilt2d(order["image"], kernel_size=(5, 5)), nxgrid, nygrid
        )

        # Normalize the raw data using the fiterpolated model

        # Get useful things and set up

        startcol = xranges[i, 0]
        stopcol = xranges[i, 1]

        order_ncols = stopcol - startcol + 1
        order_cols = np.arange(order_ncols, dtype=int) + startcol

        botedge = np.polynomial.polynomial.polyval(order_cols, edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(order_cols, edgecoeffs[i, 1, :])
        dif = topedge - botedge

        # Loop over each column, interpolate the model onto the data
        # sampling, and divide.

        # NOTE:  This has to be redone at some point...

        mask = np.full((nrows, ncols), 0)
        for j in range(order_ncols):

            # get linear conversion from pixels to arcseconds

            b = slith_arc / dif[j]
            m = -1 * b * botedge[j]

            # get range over which the slit falls and create y values

            yrange = [
                np.ceil(botedge[j]).astype("int") + ybuffer,
                np.floor(topedge[j]).astype("int") - ybuffer,
            ]

            ypix_slit = np.arange(yrange[0], yrange[1] + 1)

            # Do the linterpolation

            f = interpolate.interp1d(spatmap[:, 0], model[:, j])

            slit_model = f(np.polynomial.polynomial.polyval(ypix_slit, [m, b]))

            # divide the data by the interpolated model

            nimg[yrange[0] : yrange[1] + 1, j + startcol] = (
                img[yrange[0] : yrange[1] + 1, j + startcol] / slit_model
            )

            # divide the variance by the interpolated model if need be

            if var is not None:
                nvar[yrange[0] : yrange[1] + 1, j + startcol] = (
                    var[yrange[0] : yrange[1] + 1, j + startcol] / slit_model
                )

                # fill the mask

            mask[yrange[0] : yrange[1] + 1, j + startcol] = 1

        # Compute the RMS of the result

        z = np.where(mask == 1)
        m = moments(nimg[z], robust=4.5)
        rms[i] = m["stddev"]

        if verbose:
            loop_progress(i, 0, norders)

    return nimg, nvar, rms


def read_flat_fits(filename: str):
    """
    To read a pySpextool FITS flat image into memory.

    Parameters
    ----------
    filename : str
        The full path to a pySpextool FITS flat file.

    Returns
    -------
    dict
        'flat' : ndarray
            An (nrows, ncols) flat image array

        'var' : ndarray
            An (nrows, ncols) variance image array

        'bitmask' : ndarray of int
            An (nrows, ncols) bitmask image array.  First bit is set if a
            pixel in the original stack was over the linmax limit.

        'ordermask' : ndarray of int
            An (nrows, ncols) order mask image array.  A pixel is set to its
            order number.  Inter-order pixels are set to zero.

        'ncols' : float
            The number of columns

        'nrows' : float
            The number of rows

        'mode' : str
            The instrument mode name

        'norders' : int
            The number of orders on the image.

        'orders' : ndarray of int
            The orders numbers.  orders[0] is the order number of the order
            nearest the bottom of the image after rotation.

        'edgedeg' : int
            The polynomial degree for the fits to the edge of the order.

        'ps' : float
            The plate scale (arcseconds per pixel).

        'ybuffer' : int
            The number of native pixels from the top and bottom of the slit to
            avoid during the operation.  Useful to account for the fact that
            the drop-off in intensity at the edge of the slit is not a
            heaviside function but rather occurs over a few pixels.

        'slith_pix' : float
            The nominal slit height (pixels).

        'slith_arc' : float
            The slit height (arcseconds).

        'slitw_arc' : float
            The slit width (arcseconds).

        'slitw_pix' : float
            The slit width (pixels).

        'rp' : float
            The nominal resolving power.

        'rotation' : {0, 1, 2, 3, 4, 5, 6, 7}
            Direction to rotate a raw image so that the dispersion direction
            roughly aligns with the columns of the array and wavelength
            increases to the right, and the spatial axis roughly aligns with
            the rows of the array.

            Direction  Transpose?  Rotation Counterclockwise
            -------------------------------------------------

            0          No          None
            1          No          90 deg
            2          No          180 deg
            3          No          270 deg
            4          Yes         None
            5          Yes         90 deg
            6          Yes         180 deg
            7          Yes         270 deg

            The directions follow the IDL rotate function convention.
            Choices in brackets, default first when optional.

        'edgecoeffs' : ndarray
            An (norders, `edgedeg`+1, 2) array giving the polynomial
            coefficients delineating the top and bottom of each order.
            edgecoeffs[0,0,:] gives the coefficients for the bottom of the
            order closest to the bottom of the image and edgecoeffs[0,1,:]
            gives the coefficients for the top of said order.

        'xranges' : ndarray
            An (norders, 2) array giving the column numbers over which to
            operate.  xranges[0,0] gives the starting column number for the
            order nearest the bottom of the image and xranges[0,1] gives the
            end column number for said order.

        'rms' : list
            An (norders,) list of RMS values for each order.

    """

    #
    # Check parameters
    #

    check_parameter("read_flat_fits", "filename", filename, "str")

    #
    # Read the data
    #

    hdul = fits.open(filename)
    hdul[0].verify("silentfix")

    hdr = hdul[0].header

    flat = idl_rotate(hdul[1].data, hdr["ROTATION"])
    var = idl_rotate(hdul[2].data, hdr["ROTATION"])
    mask = idl_rotate(hdul[3].data, hdr["ROTATION"])
    ordermask = idl_rotate(hdul[4].data, hdr["ROTATION"])

    hdul.close()

    #
    # create flatinfo dictionary and add values
    #

    flatinfo = {"flat": flat}
    flatinfo.update({"var": var})
    flatinfo.update({"bitmask": np.uint8(mask)})

    shape = np.shape(flat)
    flatinfo.update({"ncols": shape[1]})
    flatinfo.update({"nrows": shape[0]})

    flatinfo.update({"mode": hdr["MODE"]})

    norders = hdr["NORDERS"]
    flatinfo.update({"norders": norders})

    orders = hdr["ORDERS"].split(",")
    orders = [int(o) for o in orders]
    flatinfo.update({"orders": np.array(orders, dtype=int)})

    edgedeg = hdr["EDGEDEG"]
    flatinfo.update({"edgedeg": edgedeg})

    flatinfo.update({"ps": hdr["PLTSCALE"]})
    flatinfo.update({"ybuffer": hdr["YBUFFER"]})
    flatinfo.update({"slith_arc": hdr["SLTH_ARC"]})
    flatinfo.update({"slith_pix": hdr["SLTH_PIX"]})
    flatinfo.update({"slitw_arc": hdr["SLTW_ARC"]})
    flatinfo.update({"slitw_pix": hdr["SLTW_PIX"]})
    flatinfo.update({"rp": round(hdr["RP"])})
    flatinfo.update({"rotation": hdr["ROTATION"]})

    # Grab the edge coeffiecients, xranges, and rms values

    edgecoeffs = np.empty([norders, 2, edgedeg])
    xranges = np.empty([norders, 2], dtype=int)
    rms = np.empty([norders])

    for i in range(norders):

        root = "OR" + str(orders[i]).zfill(3)

        for j in range(edgedeg):
            edgecoeffs[i, 0, j] = hdr[root + "_B*"][j]
            edgecoeffs[i, 1, j] = hdr[root + "_T*"][j]

        xranges[i, :] = [int(x) for x in hdr[root + "_XR"].split(",")]

        #  May not have an RMS if it wasn't normalized

        try:

            rms[i] = hdr[root + "RMS"]

        except KeyError:

            rms[i] = np.nan

    flatinfo.update({"edgecoeffs": edgecoeffs})
    flatinfo.update({"xranges": xranges})
    flatinfo.update({"rms": rms})

    flatinfo = add_entry(flatinfo, "bitmask", "after", "ordermask", ordermask)

    return flatinfo


def read_flatcal_file(filename: str):
    """
    Reads a Spextool flatinfo calibration file.


    Parameters
    ----------------
    file : str
        A Spextool flatinfo file, e.g. ShortXD_flatinfo.fits

    Returns
    --------
    dict
       A dictionary with the following keywords:

       'rotation' : int
           IDL rotation command for the order numbers to increase
           upwards and wavelength increase to the right

           Direction  Transpose?  Rotation Counterclockwise
           -------------------------------------------------

           0          No          None
           1          No          90 deg
           2          No          180 deg
           3          No          270 deg
           4          Yes         None
           5          Yes         90 deg
           6          Yes         180 deg
           7          Yes         270 deg

       'slith_arc' : float
           The slit height in arcseconds.

       'slith_pix' : float
           The nominal slit height in pixels.

       'slith_range'  : list of int
           (2, ) list range of slit heights in pixels

       'orders' : ndarray of int
           The order numbers.

       'rppix' : float
           resolving power per pixel

       'ps' : float
           plate scale in arcseconds per pixel

       'step' : int
           step size in pixels for tracing

       'flatfrac' : float
           see findorders.py

       'comwidth' : int
           The window in units of pixels used to compute the
           center-of-mass (COM) (see findorders.py)

       'edgedeg' : int
           Polynomial fit degree for the edges of the orders
           (see findorders.py)

       'norm_nxg' : int
           See normspecflat.py and fiterpolate.py

       'norm_nyg' : int
           See normspecflat.py and fiterpolate.py

       'oversamp' : float
           See normspecflat.py

       'ybuffer'  : int
           See normspecflat.py

       'ycororder' : int
           See adjustguesspos.py

       'xranges' : ndarray of int
           An (norders,2) array giving the column numbers over which to
           search.  sranges[0,0] gives the starting column number for
           the first order and sranges[0,1] gives the end column number
           for the first order.

       'edgecoeffs' : ndarray
           (norders,2,ncoeffs) array giving the polynomial
           coefficients delineating the top and bottom of each order.
           edgecoeffs[0,0,:] gives the coefficients for the bottom of
           the order closests to the bottom of `img` and
           edgecoeffs[0,1,:] gives the coefficients for the top of said
           order.

       'guesspos' : ndarray
           An (norders,2) array giving the positions to start the
           search.  guesspos[0,0] gives the column number for the
           first order and guesspos[0,1] gives the row number for the
           first order.  Typically, the positions are near the center
           of the image and center of the slit.

    Example
    --------
    result = readflatinfo('ShortXD_flatinfo.fits')

    """

    #
    # Check parameters
    #

    check_parameter("read_flatcal_file", "filename", filename, "str")

    # Open the file, grab the mask

    hdul = fits.open(filename)

    # Clean the header and grab important keywords

    hdul[0].verify("silentfix")  # this was needed for to correct hdr problems

    val = hdul[0].header["ROTATION"]
    result = {"rotation": val}

    val = hdul[0].header["SLTH_ARC"]
    result.update({"slith_arc": val})

    val = hdul[0].header["SLTH_PIX"]
    result.update({"slith_pix": val})

    val = hdul[0].header["SLTH_RNG"].split(",")
    val = [int(x) for x in val]
    result.update({"slith_range": np.array(val)})

    val = hdul[0].header["ORDERS"].split(",")
    orders = [int(x) for x in val]
    norders = len(orders)
    result.update({"orders": np.array(orders)})

    val = hdul[0].header["RPPIX"]
    result.update({"rpppix": val})

    val = hdul[0].header["PLTSCALE"]
    result.update({"ps": val})

    val = hdul[0].header["STEP"]
    result.update({"step": val})

    val = hdul[0].header["FLATFRAC"]
    result.update({"flatfrac": val})

    val = hdul[0].header["COMWIN"]
    result.update({"comwidth": val})

    deg = int(hdul[0].header["EDGEDEG"])
    result.update({"edgedeg": deg})

    val = hdul[0].header["NORM_NXG"]
    result.update({"nxgrid": int(val)})

    val = hdul[0].header["NORM_NYG"]
    result.update({"nygrid": int(val)})

    val = hdul[0].header["OVERSAMP"]
    result.update({"oversamp": val})

    val = hdul[0].header["YBUFFER"]
    result.update({"ybuffer": val})

    val = hdul[0].header["YCORORDR"]
    result.update({"ycororder": val})

    # Get the edge coefficients and the xranges

    xranges = np.empty((norders, 2), dtype=int)
    edgecoeffs = np.empty((norders, 2, deg + 1))
    guesspos = np.empty((norders, 2), dtype=int)

    for i in range(0, norders):

        # Get the xrange and guess position x position

        val = hdul[0].header["OR" + str(orders[i]).zfill(3) + "_XR"].split(",")
        val = [int(x) for x in val]
        xranges[i, :] = val

        guesspos[i, 0] = sum(xranges[i, :]) / 2

        # Now grab the edgecoefficients

        for j in range(0, deg + 1):
            keyt = "OR" + str(orders[i]).zfill(3) + "_T" + str(j + 1).zfill(1)
            keyb = "OR" + str(orders[i]).zfill(3) + "_B" + str(j + 1).zfill(1)

            edgecoeffs[i, 0, j] = hdul[0].header[keyb]
            edgecoeffs[i, 1, j] = hdul[0].header[keyt]

        # Now determine the guess position y position

        bot = np.polynomial.polynomial.polyval(guesspos[i, 0], edgecoeffs[i, 0, :])
        top = np.polynomial.polynomial.polyval(guesspos[i, 0], edgecoeffs[i, 1, :])

        guesspos[i, 1] = int((bot + top) / 2)

    result.update({"xranges": xranges})
    result.update({"edgecoeffs": edgecoeffs})
    result.update({"guesspos": guesspos})
    hdul.close()

    return result


def write_flat(
    flat: npt.ArrayLike,
    var: npt.ArrayLike,
    flag: npt.ArrayLike,
    order_mask: npt.ArrayLike,
    hdrinfo: dict,
    rotate: int,
    orders: npt.ArrayLike,
    edgecoeffs: npt.ArrayLike,
    xranges: npt.ArrayLike,
    ybuffer: int,
    ps: int | float,
    slith_pix: int | float,
    slith_arc: int | float,
    slitw_pix: int | float,
    slitw_arc: int | float,
    modename: str,
    rms: float,
    rp: int | float,
    version: str,
    history: list,
    oname: str,
    linmax: int | float = None,
    overwrite: bool = True,
):
    """
    To write a Spextool flat FITS file to disk.


    Parameters
    ----------
    flat : numpy.ndarray
        (nrows, ncols) flat field image.

    var : numpy.ndarray
        (nrows, ncols) variance image.

    flag : numpy.ndarray
        (nrows, ncols) int bitset image.

    hdrinfo : dict

    rotate : {0, 1, 2, 3, 4, 5, 6, 7}, optional

        Direction to rotate a raw image so that the dispersion direction
        roughly aligns with the columns of the array and wavelength
        increases to the right, and the spatial axis roughly aligns with
        the rows of the array.

        Direction  Transpose?  Rotation Counterclockwise
        -------------------------------------------------

        0          No          None
        1          No          90 deg
        2          No          180 deg
        3          No          270 deg
        4          Yes         None
        5          Yes         90 deg
        6          Yes         180 deg
        7          Yes         270 deg

        The directions follow the IDL rotate function convention.
        Choices in brackets, default first when optional.

    orders : list of int
        The order numbers.  orders[0] is the order number of the
        order nearest the bottom of the image after rotation.

    edgecoeffs : array_like of float
        (norders,ncoeffs+1,2) array giving the polynomial coefficients
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients
        for the top of said order.

    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to
        operate.  xranges[0,0] gives the starting column number for the
        order nearest the bottom of the image and xranges[0,1] gives the
        end column number for said order.

    ybuffer : int
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels.

    ps : float
        The plate scale (arcseconds per pixel).

    slith_pix : float
        The slit height (pixels).

    slith_arc : float
        The nominal slit height (arcseconds).

    slitw_arc : float
        The slit width (arcseconds).

    slitw_pix : float
        The slit width (pixels).

    modename : str
        The name of the instrument mode.

    rms : list of float
        (norders,) list of RMS values for each order.

    rp: float
        The nominal resolving power.

    version : str
        The version of the pySpextool software used to create the flat.

    oname : str
        The filename of the flat field image to written to disk.

    linmax : int, optional
        The linearity maximum used to identified pixels beyond the
        linearity limit.

    history : list of str, optional
        The history string spliced to fit in a FITS file.

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.

    Returns
    -------
    None
        Writes a FITS file to PROCPATH.

    Examples
    --------
    TODO: add example

    """

    # Get basic things

    norders = len(orders)
    edgedeg = edgecoeffs.shape[-1]

    # Create the primary HDU

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the hdrlist keywords and values

    keys = hdrinfo.keys()
    for key in keys:

        if key == "COMMENT":

            comments = hdrinfo["COMMENT"]

        elif key != "HISTORY":  # probably not necessary, just in case

            hdr[key] = tuple(hdrinfo[key])

        # Now add new ones

        hdr["FILENAME"] = (os.path.basename(oname), " Filename")
        hdr["MODE"] = (modename, " Instrument Mode")
        hdr["NORDERS"] = (norders, " Number of orders identified")
        hdr["ORDERS"] = (",".join(str(o) for o in orders), " Orders identified")
        hdr["PLTSCALE"] = (ps, " Plate scale (arcseconds per pixel)")
        hdr["YBUFFER"] = (ybuffer, " y buffer (pixels)")
        hdr["SLTH_PIX"] = (slith_pix, " Nominal slit length (pixels)")
        hdr["SLTH_ARC"] = (slith_arc, " Slit length (arcseconds)")
        hdr["SLTW_PIX"] = (slitw_pix, " Slit width (pixels)")
        hdr["SLTW_ARC"] = (slitw_arc, " Slit width (arcseconds)")
        hdr["RP"] = (rp, " Nominal resolving power")
        hdr["ROTATION"] = (rotate, " IDL rotate value")
        hdr["VERSION"] = (version, " Spextool version")

    # Record linearity maximum if given

    if linmax is not None:
        hdr["LINMAX"] = (linmax, " Linearity maximum")

    # Add the RMS values.  Check to make sure not NaN

    if np.sum(np.isnan(rms)) == 0:

        for i in range(norders):
            name = "OR" + str(orders[i]).zfill(3) + "RMS"
            comment = " RMS of normalized order " + str(orders[i]).zfill(3)
            hdr[name] = (rms[i], comment)

    # Add the xranges

    for i in range(norders):
        name = "OR" + str(orders[i]).zfill(3) + "_XR"
        comment = " Extraction range for order " + str(orders[i]).zfill(3)
        hdr[name] = (",".join(str(x) for x in xranges[i, :]), comment)

    # Add the edgecoeffs

    hdr["EDGEDEG"] = (edgedeg, " Degree of the polynomial fit to order edges")
    for i in range(norders):

        for j in range(2):

            for k in range(edgedeg):

                if j == 0:
                    name = "OR" + str(orders[i]).zfill(3) + "_B" + str(k + 1)
                    comment = (
                        " a"
                        + str(k)
                        + " edge coefficient for bottom of order "
                        + str(orders[i]).zfill(3)
                    )

                    hdr[name] = (edgecoeffs[i, j, k], comment)

                if j == 1:
                    name = "OR" + str(orders[i]).zfill(3) + "_T" + str(k + 1)
                    comment = (
                        " a"
                        + str(k)
                        + " edge coefficient for top of order "
                        + str(orders[i]).zfill(3)
                    )

                    hdr[name] = (edgecoeffs[i, j, k], comment)

    # Now add the comments

    if "comments" in locals():

        for com in comments:
            hdr["COMMENT"] = com

    # and then the history

    label = "         ============ Spextool History ============"
    hdr["HISTORY"] = label
    for hist in history:
        hdr["HISTORY"] = hist

    # Write the results

    flat = np.float32(flat)
    var = np.float32(var)
    flag = np.int8(flag)
    order_mask = np.int8(order_mask)

    img_hdu = fits.ImageHDU(idl_unrotate(flat, rotate))
    var_hdu = fits.ImageHDU(idl_unrotate(var, rotate))
    flg_hdu = fits.ImageHDU(idl_unrotate(flag, rotate))
    ord_hdu = fits.ImageHDU(idl_unrotate(order_mask, rotate))

    hdu = fits.HDUList([phdu, img_hdu, var_hdu, flg_hdu, ord_hdu])
    hdu.writeto(oname, overwrite=overwrite)
