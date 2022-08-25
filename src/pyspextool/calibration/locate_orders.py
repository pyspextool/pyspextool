import numpy as np
import matplotlib.pyplot as pl
import os
from scipy import ndimage
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.plot.limits import get_image_range

np.set_printoptions(threshold=np.inf)


def locate_orders(img, guess_positions, search_ranges, step_size,
                  slit_height_range, poly_degree, pixel_buffer,
                  intensity_fraction, com_width, qafileinfo=None):
    """
    Locates orders in a (cross-dispersed) spectral image


    Parameters
    ----------
    img : array_like of float
        An (nrows, ncols) image with (cross-dispersed) spectral orders.  It is
        assumed that the dispersion direction is roughly aligned with the
        rows of `img` and the spatial axis is roughly aligned with the
        columns of `img.  That is, orders go left-right and not up-down.

    guess_positions : array_like
        An (norders,2) int array giving the positions to start the search.
        guesspos[0,0] gives the column number for the order closest to the
        bottom of `img` and guesspos[0,1] gives the row number for said
        order.  Typically, the positions are near the center of the image
        and center of the slit.

    search_ranges : array_like
        An (norders,2) int array giving the column numbers over which to
        search.  search_ranges[0,0] gives the starting column number for
        the order closest to the bottom of `img` and search_ranges[0,1]
        gives the end column number for said order.

    step_size : int
        The step size in the dispersion or column direction

    slit_height_range: array_like
        (2,) array giving the minimum and maxmimum possible slit heights in
        units of pixels.

    poly_degree : int
        Polynomial fit degree for the edges of the orders

    pixel_buffer : int
        Number of pixels to buffer from the top and bottom of `img`
        during the search

    intensity_fraction : float
        The fraction of the flux of the center of the slit used to
        identify the location of the edge of the order

    com_width : int
        The window in units of pixels used to compute the center-of-mass (COM)

    qafileinfo : dict, optional
        `"figsize"` : tuple
            (2,) tuple of the figure size (inches).

        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

    Returns
    -------
    edgecoeffs : numpy.ndarray
        (norders,2,ncoeffs) array giving the polynomial coefficients
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for
        the top of said order.

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

    Examples
    --------
    later?

    """

    debug = False

    # Check for qaplotting and set up empty lists

    if qafileinfo is not None:
        plcols = []
        pledges = []
        pledgecoeffs = []
        plogoodbad = []

    # Get basic info and do basic things

    nrows, ncols = img.shape
    norders = len(guess_positions)

    rownum = np.arange(nrows)

    halfwin = int(com_width / 2.)

    edgecoeffs = np.empty((norders, 2, poly_degree + 1))

    # Sobel the image

    scl = np.max(img)
    simg1 = ndimage.sobel(img / scl, axis=0)
    simg2 = ndimage.sobel(img / scl, axis=1)
    simg = np.sqrt(simg1 ** 2 + simg2 ** 2)

    if debug is not False:
        pl.ion()
        plotfig = pl.figure(2, figsize=(6, 6))

    # Start looping over each order

    for i in range(0, norders):

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

        gpidx = np.asscalar(np.where(np.absolute(fcols - guess_positions[i, 0]) ==
                                     np.absolute(fcols - guess_positions[i, 0]).min())[0])

        cens[(gpidx - poly_degree):(gpidx + poly_degree)] = \
            guess_positions[i, 1]

        #
        # --------------Now move left from the guess position-------------------
        #

        for j in range(gpidx, 0, -1):

            imgcol = img[:, fcols[j]].ravel()  # the image column
            sobcol = simg[:, fcols[j]].ravel()  # the sobel image column

            # Fit the centers, so you can project away from the guess position

            r = poly_fit_1d(fcols, cens, max(1, poly_degree - 2), justfit=True,
                            silent=True)

            # Find the new guess position yguess 

            yguess = np.polynomial.polynomial.polyval(fcols[j], r['coeffs'])

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess, pixel_buffer, (nrows - pixel_buffer - 1)))
            iguess = imgcol[yguess]

            bot, top = find_top_bot(fcols[j], rownum, imgcol, sobcol, yguess,
                                    iguess, intensity_fraction, halfwin,
                                    debug=debug)

            # Confirm that both COMs were computed

            if np.isfinite(bot):

                # Now check to make sure calculated slit height falls within the
                # parameter slith_pix

                if slit_height_range[0] <= top - bot <= slit_height_range[1]:

                    # Store the results, update the cens array
                    edges[:, j] = np.array([bot, top])
                    cens[j] = (bot + top) / 2

        #
        # --------------Now move right from the guess position-------------------
        #

        for j in range(gpidx, nfcols, 1):

            imgcol = img[:, fcols[j]].ravel()  # the image column
            sobcol = simg[:, fcols[j]].ravel()  # the sobel image column

            # Fit the centers, so you can project away from the guess position

            r = poly_fit_1d(fcols, cens, max(1, poly_degree - 2), silent=True)

            # Find the new guess position yguess

            yguess = np.polynomial.polynomial.polyval(fcols[j], r['coeffs'])

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess, pixel_buffer, (nrows - pixel_buffer - 1)))
            iguess = imgcol[yguess]

            bot, top = find_top_bot(fcols[j], rownum, imgcol, sobcol, yguess,
                                    iguess, intensity_fraction, halfwin,
                                    debug=debug)

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

            fit = poly_fit_1d(fcols, edges[j, :], poly_degree,
                              robust={'thresh': 4, 'eps': 0.1},
                              justfit=True, silent=True)
            tmp[j, :] = fit['coeffs']

            # Store the results for possible plotting         

            if qafileinfo is not None:
                plcols.append(fcols)
                pledges.append(edges[j, :])
                pledgecoeffs.append(fit['coeffs'])
                plogoodbad.append(fit['goodbad'])

        edgecoeffs[i, :, :] = tmp

    if qafileinfo is not None:

        minmax = get_image_range(img, 'zscale')
        fig = pl.figure(figsize=(7, 7))
        pl.imshow(img, vmin=minmax[0], vmax=minmax[1], cmap='gray', origin='lower')
        pl.xlabel('Columns (pixels)')
        pl.ylabel('Rows (pixels)')

        # Overplot everything

        for i in range(norders):
            pl.plot(guess_positions[i, 0], guess_positions[i, 1], 'go')

        for i in range(len(pledges)):

            pl.plot(plcols[i], pledges[i], 'go', markersize=1)

            z = (np.where(plogoodbad[i] == 0))[0]
            if z.size != 0:
                pl.plot(plcols[i][z], pledges[i][z], 'bo', markersize=1)

            pl.plot(plcols[i],
                    np.polynomial.polynomial.polyval(plcols[i], pledgecoeffs[i]),
                    'r-', linewidth=0.3)

        pl.savefig(os.path.join(qafileinfo['filepath'],
                                qafileinfo['filename'] + '_locateorders'
                                + qafileinfo['extension']))

    return edgecoeffs


def find_top_bot(fcol, rownum, imgcol, sobcol, yguess, imgguess, frac, halfwin,
                 debug=False):

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

    if debug is not False:
        plotfig = pl.figure(2)
        pl.plot(rownum, sobcol)
        pl.xlim([com_bot - 10, com_top + 10])
        pl.axvline(com_top, color='r')
        pl.axvline(com_bot, color='r')
        pl.axvline(yguess, color='g')
        pl.title(fcol)
        pl.pause(0.1)
        #        val = input("-->")
        plotfig.clear()

    return com_bot, com_top
