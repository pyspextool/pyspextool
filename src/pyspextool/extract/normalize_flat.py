import numpy as np
from scipy import interpolate
from scipy.signal import medfilt2d
from scipy import interpolate

from pyspextool.extract.make_interp_indices_1d import make_interp_indices_1d
from pyspextool.extract.rectify_order import rectify_order
from pyspextool.fit.fiterpolate import *
from pyspextool.utils.math import moments
from pyspextool.utils.loop_progress import loop_progress


def normalize_flat(img, edgecoeffs, xranges, slith_arc, nxgrid, nygrid,
                   var=None, oversamp=1, ybuffer=0, verbose=False):
    """
    Normalize spectral flat field image


    Parameters
    ----------
    img : array_like of float
        An (nrows,ncols) image with (cross-dispersed) spectral orders.
        It is assumed that the dispersion direction is roughly aligned
        with the rows of `img` and the spatial axis is roughly aligned
        with the columns of `img.  That is, orders go left-right and
        not up-down.

    edgecoeffs : array_like of float
        (norders,ncoeffs+1,2) array giving the polynomial coefficients
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for
        the top of said order.

    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to
        operate.  xranges[0,0] gives the starting column number for the
        order nearest the bottom of `img` and xranges[0,1] gives the end
        column number for said order.

    slith_arc : float
        The slit height in arcseconds

    nxgrid : int

    nygrid : int

    var : numpy.ndarray

    ybuffer : int, default 1, optional
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels.

    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray

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


    Examples
    --------
    later

    """

    # Get basic info and do basic things

    nrows, ncols = img.shape

    ndimen = edgecoeffs.ndim
    if ndimen == 2:
        norders = 1
    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    nimg = np.full((nrows, ncols), 1.0)
    nvar = np.full((nrows, ncols), np.nan) if var is not None else None
    rms = np.empty(norders)

    # Start the loop

    for i in range(norders):

        #
        # Rectify the order
        #

        # Get the rectification indices

        xidx,yidx,wavemap,spatmap = make_interp_indices_1d(edgecoeffs[i, :, :],
                                                            xranges[i,:],
                                                            slith_arc)

        # Do the rectification

        order = rectify_order(img, xidx, yidx, ybuffer=ybuffer)
                
        # Fiterpolate the results after median smoothing to minimize bad pixels

        model = fiterpolate(medfilt2d(order['img'], kernel_size=(5, 5)), nxgrid,
                            nygrid)

        # Now normalize the raw data using the fiterpolated model

        # Get useful things and set up
        
        startcol = xranges[i, 0]
        stopcol = xranges[i, 1]

        order_ncols = stopcol - startcol + 1
        order_cols = np.arange(order_ncols, dtype=int) + startcol

        botedge = np.polynomial.polynomial.polyval(order_cols,
                                                   edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(order_cols,
                                                   edgecoeffs[i, 1, :])
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

            yrange = [np.ceil(botedge[j]).astype('int') + ybuffer,
                      np.floor(topedge[j]).astype('int') - ybuffer]

            ypix_slit = np.arange(yrange[0], yrange[1] + 1)

            # Do the linterpolation

            f = interpolate.interp1d(spatmap[:,0], model[:, j])
            tmp = np.polynomial.polynomial.polyval(ypix_slit, [m, b])

            slit_model = f(np.polynomial.polynomial.polyval(ypix_slit, [m, b]))

            # divide the data by the interpolated model

            nimg[yrange[0]:yrange[1] + 1, j + startcol] = \
                img[yrange[0]:yrange[1] + 1, j + startcol] / slit_model

            # divide the variance by the interpolated model if need be

            if var is not None:
                nvar[yrange[0]:yrange[1] + 1, j + startcol] = \
                    var[yrange[0]:yrange[1] + 1, j + startcol] / slit_model

                # fill the mask

            mask[yrange[0]:yrange[1] + 1, j + startcol] = 1

        # Compute the RMS of the result

        z = np.where(mask == 1)
        m = moments(nimg[z], robust=4.5)
        rms[i] = m['stddev']

        if verbose:
            loop_progress(i, 0, norders, message='Normalizing the flat...')

    return nimg, nvar, rms
