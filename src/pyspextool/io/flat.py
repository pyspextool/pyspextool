import numpy as np
import logging
from astropy.io import fits
import os

from pyspextool.utils.add_entry import add_entry
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils.arrays import idl_unrotate


def read_flat_fits(flatfile):
    """
    To read a pySpextool FITS flat image into memory.


    Parameters
    ----------
    flatfile : str
        The full path to a pySpextool FITS flat file.


    Returns
    -------
    dict
        ``"flat"``
            (numpy.ndarray) -> (nrows,cols) flat image

        ``"var"``
            (numpy.ndarray) -> (nrows,cols) variance image

        ``"bitmask"``
            (numpy.ndarray of int) -> (nrows,cols) bitmask image.  First
            bit is set if a pixel in the original stack was over the
            linmax limit.

        ``"ordermask"``
            (numpy.ndarray of int) -> (nrows,cols) order mask image.
            A pixel is set to its order number.

        ``"ncols"``
            float -> The number of columns.

        ``"nrows"``
            float -> The number of rows.

        ``"mode"``
            str -> The instrument mode name.

        ``"norders"``
            int -> The number of orders on the image.

        ``"orders"``
            numpy.ndarray of int -> The orders numbers.  orders[0] is
            the order number of the order nearest the bottom of the
            image after rotation.

        ``"edgedeg"``
            int -> The polynomial degree for the fits to the edge of
            the order.

        ``"ps"``
            float -> The plate scale (arcseconds per pixel).

        ``"ybuffer"``
            int -> The number of native pixels from the top and bottom of the
                   slit to avoid during the operation.  Useful to account for
                   the fact that the drop-off in intensity at the edge of the
                   slit is not a heaviside function but rather occurs over a
                   few pixels.

        ``"slith_pix"``
            float -> The slit height (pixels).

        ``"slith_arc"``
            float -> The nominal slit height (arcseconds).

        ``"slitw_arc"``
            float -> The slit width (arcseconds).

        ``"slitw_pix"``
            float -> The slit width (pixels).

        ``"rp"``
            float -> The nominal resolving power.

        ``"rotation"``
            {0, 1, 2, 3, 4, 5, 6, 7} -> Direction to rotate a raw image
            so that the dispersion direction roughly aligns with the
            columns of the array and wavelength increases to the right,
            and the spatial axis roughly aligns with the rows of the array.
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

        ``"edgecoeffs"``
            array_like of float -> (norders,`edgedeg`+1,2) array giving
            the polynomial coefficients delineating the top and bottom
            of each order.  edgecoeffs[0,0,:] gives the coefficients for
            the bottom of the order closest to the bottom of the image
            and edgecoeffs[0,1,:] gives the coefficients for the top of
            said order.

        ``"xranges"``
            array_like of float -> An (norders,2) array giving the
            column numbers over which to operate.  xranges[0,0] gives
            the starting column number for the order nearest the bottom
            of the image and xranges[0,1] gives the end column number
            for said order.

        ``"rms"``
            list of float -> (norders,) list of RMS values for each order.

    """

    # Read the data

    hdul = fits.open(flatfile)
    hdul[0].verify("silentfix")

    hdr = hdul[0].header

    flat = idl_rotate(hdul[1].data, hdr["ROTATION"])
    var = idl_rotate(hdul[2].data, hdr["ROTATION"])
    mask = idl_rotate(hdul[3].data, hdr["ROTATION"])

    hdul.close()

    # create flatinfo dictionary

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
    flatinfo.update({"rp": hdr["RP"]})
    flatinfo.update({"rotation": hdr["ROTATION"]})

    # Grab the edge coeffiecients, xranges, and rms values

    ordermask = np.zeros([flatinfo["nrows"], flatinfo["ncols"]], dtype=int)
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

        except KeyError as err:
            rms[i] = np.nan

        # Create order mask

        x = np.arange(xranges[i, 0], xranges[i, 1] + 1, 1, dtype=int)
        bot = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 0, :])
        top = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 1, :])

        for j in range(len(x)):
            ordermask[
                np.floor(bot[j]).astype("int") : np.ceil(top[j]).astype("int"), x[j]
            ] = orders[i]

    flatinfo.update({"edgecoeffs": edgecoeffs})
    flatinfo.update({"xranges": xranges})
    flatinfo.update({"rms": rms})
    flatinfo = add_entry(flatinfo, "bitmask", "after", "ordermask", ordermask)

    return flatinfo


def read_flatcal_file(file):
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

       rotation : int
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

       slith_arc : int
           The slit height in arcseconds.

       slith_pix : int
           The nominal slit height in pixels.

       slith_range  : list of int
           (2,) list range of slit heights in pixels

       orders : list of int
           The order numbers.

       rppix : float
           resolving power per pixel

       ps : float
           plate scale in arcseconds per pixel

       step : int
           step size in pixels for tracing

       flatfrac : float
           see findorders.py

       comwidth : int
           The window in units of pixels used to compute the
           center-of-mass (COM) (see findorders.py)

       edgedeg : int
           Polynomial fit degree for the edges of the orders
           (see findorders.py)

       norm_nxg : int
           See normspecflat.py and fiterpolate.py

       norm_nyg : int
           See normspecflat.py and fiterpolate.py

       oversamp : float
           See normspecflat.py

       ybuffer  : int
           See normspecflat.py

       ycororder : int
           See adjustguesspos.py

       xranges : array_like of int, [norders,2]
           An (norders,2) array giving the column numbers over which to
           search.  sranges[0,0] gives the starting column number for
           the first order and sranges[0,1] gives the end column number
           for the first order.

       edgecoeffs : array_like of float
           (norders,2,ncoeffs) array giving the polynomial
           coefficients delineating the top and bottom of each order.
           edgecoeffs[0,0,:] gives the coefficients for the bottom of
           the order closests to the bottom of `img` and
           edgecoeffs[0,1,:] gives the coefficients for the top of said
           order.

       guesspos : array_like of int
           An (norders,2) array giving the positions to start the
           search.  guesspos[0,0] gives the column number for the
           first order and guesspos[0,1] gives the row number for the
           first order.  Typically, the positions are near the center
           of the image and center of the slit.

    Procedure
    ---------
    Just reads FITS header information, and calculates the guess
    positions based on the edgecoeffs and the xranges

    Example
    --------
    result = readflatinfo('ShortXD_flatinfo.fits')

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """

    # Open the file, grab the mask
    try:
        hdul = fits.open(file)
        logging.debug(f"Opened flatcal file: {file}")
    except OSError as e:
        msg = f"""
        Could not open the flatcal file: {file}.
        Please check that the file downloaded from Git LFS properly. 
        It should be approximately 2 MB.
        Download directly from <URL TBD>"""
        logging.error(msg)
        raise OSError(e)

    omask = hdul[0].data

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
    result.update({"slith_range": val})

    val = hdul[0].header["ORDERS"].split(",")
    orders = [int(x) for x in val]
    norders = len(orders)
    result.update({"orders": orders})

    val = hdul[0].header["RPPIX"]
    result.update({"rpppix": val})

    val = hdul[0].header["PLTSCALE"]
    result.update({"ps": val})

    #    val = hdul[0].header['FIXED']
    #    result.update({'fixed':val})

    #    if not val:

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

    # Now get the edge coefficients and the xranges

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
    flat,
    var,
    flag,
    hdrinfo,
    rotate,
    orders,
    edgecoeffs,
    xranges,
    ybuffer,
    ps,
    slith_pix,
    slith_arc,
    slitw_pix,
    slitw_arc,
    modename,
    rms,
    rp,
    version,
    history,
    oname,
    linmax=None,
    overwrite=True,
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
    Later

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

    img_hdu = fits.ImageHDU(idl_unrotate(flat, rotate))
    var_hdu = fits.ImageHDU(idl_unrotate(var, rotate))
    flg_hdu = fits.ImageHDU(idl_unrotate(flag, rotate))

    hdu = fits.HDUList([phdu, img_hdu, var_hdu, flg_hdu])
    hdu.writeto(oname, overwrite=overwrite)
