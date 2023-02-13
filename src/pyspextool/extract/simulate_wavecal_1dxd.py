import numpy as np


def simulate_wavecal_1dxd(ncols, nrows, edgecoeffs, xranges, slith_arc):

    """
    To simulate Spextool wavecal and spatcal arrays.

    Will generate wavecal and spatcal files in the 1DXD case with the
    wavelengths replaced with the column numbers.


    Input Parameters
    ----------------
    ncols : int
        The number of columns of the image.

    nrows : int
        The number of rows of the image.

    edgecoeffs : array_like of float
        (norders,`edgedeg`+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  

    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.

    slith_arc : float
        The nominal slit height (arcseconds).

    Returns
    -------
    wavecal, spatcal : numpy.ndarray, numpy.ndarray
        - wavecal (nrows,ncols) array where each pixel is set to its 
          wavelength which in this case is the column number.
        - spatcal (nrows,ncols) array where each pixel is set to its 
          angular position on the sky (in arcseconds).

    Notes
    -----
    None

    Examples
    --------
    later

    Modification History
    --------------------
    2022-06-28 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_simwavecal2d.pro.

    """

    # Get basic info

    ndimen = edgecoeffs.ndim
    if ndimen == 2:
        norders = 1

    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    # Make some NaN arrays

    wavecal = np.full([nrows, ncols], np.nan)
    spatcal = np.full_like(wavecal, np.nan)

    # start the loop over order and column number

    y = np.arange(nrows)

    for i in range(norders):

        start = xranges[i, 0]
        stop = xranges[i, 1]

        x = np.arange(stop - start + 1) + start

        # Get the top and bottom positions of the slit

        botedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 1, :])

        # Creat the pixel to arcsecond transformation

        pixtoarc = np.empty([2, stop - start + 1])
        pixtoarc[1, :] = slith_arc / (topedge - botedge)
        pixtoarc[0, :] = -1 * pixtoarc[1, :] * botedge

        # Fill things in

        for j in range(stop - start + 1):

            wavecal[np.floor(botedge[j]).astype('int'):
                    np.ceil(topedge[j]).astype('int'), x[j]] = x[j]

            # Create ysub to make things readable...

            ysub = y[np.floor(botedge[j]).astype('int'):
                     np.ceil(topedge[j]).astype('int')]

            spatcal[np.floor(botedge[j]).astype('int'):
                    np.ceil(topedge[j]).astype('int'), x[j]] = \
                np.polynomial.polynomial.polyval(ysub, pixtoarc[:, j])

    return wavecal, spatcal
