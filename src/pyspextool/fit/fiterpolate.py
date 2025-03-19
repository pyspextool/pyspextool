"""Functions for interpolation."""

import numpy as np
from pyspextool.utils.arrays import make_image_indices


def bicucof(z, dz1, dz2, dz12, nd1, nd2):

    """
    Returns the coefficients used in bicuvals.py.

    Parameters
    ----------------
    z : array_like
        The grid points values starting at the lower left and moving 
        counterclockwise.

    dz1 : float
        The gradient in dimension 1 evaluated at z.

    dz2 : float
        The gradient in dimension 2 evaluated at z.

    dz12 : float
        The cross derivative evaluated at z.

    nd1 : int
        Length of the grid cell in dimension 1.

    nd2 : int
        Length of the grid cell in dimension 2.

    Returns
    --------
    numpy.ndarray
        Returns a (4,4) set of coefficients used by bicuval.py.

    Notes
    ---------
    See bicucof in Numerical Recipes

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_bicucoeffs.pro.

    """

    wt = np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 0., 1.,
                    0., 0., 0., 0., 0., 0., 0.],
                   [-3., 0., 0., 3., 0., 0., 0., 0.,
                    -2., 0., 0., -1., 0., 0., 0., 0.],
                   [2., 0., 0., -2., 0., 0., 0., 0.,
                    1., 0., 0., 1., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 1., 0., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 0.,
                    0., 0., 0., 0., 1., 0., 0., 0.],
                   [0., 0., 0., 0., -3., 0., 0., 3.,
                    0., 0., 0., 0., -2., 0., 0., -1.],
                   [0., 0., 0., 0., 2., 0., 0., -2.,
                    0., 0., 0., 0., 1., 0., 0., 1.],
                   [-3., 3., 0., 0., -2., -1., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 0.,
                    -3., 3., 0., 0., -2., -1., 0., 0.],
                   [9., -9., 9., -9., 6., 3., -3., -6.,
                    6., -6., -3., 3., 4., 2., 1., 2.],
                   [-6., 6., -6., 6., -4., -2., 2., 4.,
                    -3., 3., 3., -3., -2., -1., -1., -2.],
                   [2., -2., 0., 0., 1., 1., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 0., 0., 0., 2.,
                    -2., 0., 0., 1., 1., 0., 0.],
                   [-6., 6., -6., 6., -3., -3., 3., 3.,
                    -4., 4., 2., -2., -2., -2., -1., -1.],
                   [4., -4., 4., -4., 2., 2., -2., -2.,
                    2., -2., -2., 2., 1., 1., 1., 1.]])

    x = np.concatenate((z, dz1 * nd1, dz2 * nd2, dz12 * nd1 * nd2))

    return np.reshape(np.matmul(wt, x), (4, 4))


def bicuval(z, dz1, dz2, dz12, xl, xu, yl, yu, x, y):

    """
    Evaluates a bicubic interpolation

    Parameters
    ----------------
    z : array_like
        The grid points values starting at the lower left and moving 
        counterclockwise.

    dz1 : float
        The gradient in dimension 1 evaluated at z.

    dz2 : float
        The gradient in dimension 2 evaluated at z.

    dz12 : float
        The cross derivative evaluated at z.

    xl : int
        The lower coordinate of the grid in the x direction.

    xu : int
        The upper coordinate of the grid in the x direction.

    yl : int
        The lower coordinate of the grid in the y direction.

    yu : int
        The upper coordinate of the grid in the y direction.

    x : array_like
        The x coordinates of the desired interpolation point(s)

    y : array_like
        The y values of the desired interpolation point(s)

    Returns
    --------
    numpy.ndarray
        The bicubic interpolation values at x and y.

    Notes
    ---------
    See Numerical Recipes bcucof

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_bicuval.pro.

    """

    c = bicucof(z, dz1, dz2, dz12, xu - xl, yu - yl)

    t = (x - xl) / (xu - xl)
    u = (y - yl) / (yu - yl)

    nz = np.zeros(x.shape)

    for i in range(3, -1, -1):
        nz = t * nz + c[i, 0] + u * (c[i, 1] + u * (c[i, 2] + u * c[i, 3]))

    return nz


def fiterpolate(img, ncg, nrg):

    """
    Not sure how to explain!

    Parameters
    ----------------
    img : array_like
        The image to be fiterpolate-d.  

    ncg : int
        The number of columnn grid cells
       
    nrg : int
        The number of row grid cells

    Returns
    --------
    numpy.ndarray
        The fitted image

    Notes
    -----
    Based on John Tonry's fiterpolate program.  Going to require
    me going into this again and remembering exactly what it does.
    

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_fiterpolate.pro.

    """

    # Get basic infor and create basic things

    nrows, ncols = img.shape

    fitimg = np.empty((nrows, ncols))

    # Determine the number of grid points

    nxgrid = ncg + 1
    nygrid = nrg + 1

    # Compute the actual grid points

    gx = np.rint(np.arange(nxgrid) * (ncols - 1) / (nxgrid - 1)).astype(int)
    gy = np.rint(np.arange(nygrid) * (nrows - 1) / (nygrid - 1)).astype(int)

    # Set up the grid information for ease of use later

    gridinfo = []

    for i in range(ncg):

        for j in range(nrg):
            dict = {}
            dict.update({'xrng': [gx[i], gx[i + 1]]})
            dict.update({'yrng': [gy[j], gy[j + 1]]})
            dict.update({'xsize': gx[i + 1] - gx[i] + 1})
            dict.update({'ysize': gy[j + 1] - gy[j] + 1})
            gridinfo.append(dict)

    # Now determine the values of z, dy1, dy2, dy12 at the grid points.
    # They will be stored in the val array

    idx = 0
    vals = np.empty((nxgrid * nygrid, 4))

    for i in range(nxgrid):

        x1 = np.rint((gx[np.max([i - 1, 0])] + gx[i]) / 2).astype(int)
        x2 = np.rint((gx[i] + gx[np.min([nxgrid - 1, i + 1])]) / 2).astype(int)

        for j in range(nygrid):
            y1 = np.rint((gy[np.max([j - 1, 0])] + gy[j]) / 2).astype(int)
            y2 = np.rint((gy[j] + gy[np.min([nygrid - 1, j + 1])]) / 2). \
                astype(int)

            c = imgquadfit(img[y1:y2 + 1, x1:x2 + 1])

            x = gx[i] - x1  # offset to account for zero-based fit
            y = gy[j] - y1

            vals[idx, 0] = c[0] + c[1] * x + c[2] * y + c[3] * x ** 2 + c[4] * y ** 2 + c[5] * x * y
            vals[idx, 1] = c[1] + 2. * c[3] * x + c[5] * y
            vals[idx, 2] = c[2] + 2. * c[4] * y + c[5] * x
            vals[idx, 3] = c[5]

            idx += 1

            # Now perform the bicubic interpolation and reconstruct the image

    k = 0
    for i in range(ncg):

        for j in range(nrg):
            idx = [i * nygrid + j, (i + 1) * nygrid + j, (i + 1) * nygrid + j + 1,
                   (1 + i * nygrid) + j]

            # Create coordinate images

            ximg, yimg = make_image_indices(int(gridinfo[k]['ysize']),
                                            int(gridinfo[k]['xsize']))

            nimg = bicuval(vals[idx, 0], vals[idx, 1], vals[idx, 2],
                           vals[idx, 3], 0.0, gridinfo[k]['xsize'], 0.0,
                           gridinfo[k]['ysize'], ximg, yimg)

            fitimg[gridinfo[k]['yrng'][0]:gridinfo[k]['yrng'][1] + 1,
            gridinfo[k]['xrng'][0]:gridinfo[k]['xrng'][1] + 1] = nimg

            k += 1

    return fitimg


def imgquadfit(img, imgunc=None, doalpha=False):
    """
    Fits a quadratic surface to an image

    Parameters
    ----------------
    img : array_like
        The image to fit.

    imgunc : array_like, optional
        The uncertainty image.

    doalpha : {False, True}, optional
        If true, then A**T * A is computed directly instead of A.

    Returns
    --------
    numpy.ndarray
        The coefficients of the fit

    Notes
    -----
    Standard least squares fitting (e.g. Numerical Recipes).  
    The program fits the following function:

    z = c[0] + c[1]*x + c[2]*y + c[3]x^2 + c[4]y^2 + c[5]*x*y

    where x is the column number and y is the row number.

    By default, the design matrix A is constructed.  The alpha 
    (A^T ## A) and beta (A^T ## b) arrays of the normal equations 
    are then constructed and then passed to np.linalg.solve.  

    If the dataset to be fit is large and/or the number of coefficients 
    is large, the design matrix A can be large (ncoeff,ndat) and may 
    take too much memory and/or too much time to construct.  In this case,
    set `doalpha` to True and the alpha and beta arrays are constructed 
    directly.  There is a trade-off since constructing the alpha and beta
    arrays directly takes longer for smaller surfaces.

    Unlike polyfit1d and polyfit2d, this function ONLY computes the 
    coefficients.

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_quadfit.pro

    """

    # Do basic things and create basic things

    nrows, ncols = img.shape
    if imgunc is None:
        imgunc = np.full((nrows, ncols), 1.0)

    # Create coordinate images

    ximg = np.tile(np.arange(ncols, dtype=float), (nrows, 1))
    yimg = np.tile(np.reshape(np.arange(nrows, dtype=float),
                              (nrows, 1)), (1, ncols))

    # Ravel the images

    x = np.ravel(ximg)
    y = np.ravel(yimg)
    z = np.ravel(img)
    zunc = np.ravel(imgunc)

    # Get rid of NaNs 

    znonan = ~np.isnan(z)
    nnan = np.size(znonan) - np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    zz = z[znonan]
    zzunc = zunc[znonan]

    # Set up to construct the equations

    ndat = len(xx)
    ncoeffs = 6
    exp = np.array([[0, 0], [1, 0], [0, 1], [2, 0], [0, 2], [1, 1]])

    b = zz / zzunc

    # Now start the linear algebra construction

    if doalpha is True:

        alpha = np.empty((ncoeffs, ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i, ncoeffs):
                at = (xx ** exp[i, 0] * yy ** exp[i, 1]) / zzunc
                a = (xx ** exp[j, 0] * yy ** exp[j, 1]) / zzunc

                val = np.sum(at * a)

                alpha[j, i] = val
                alpha[i, j] = val
                beta[i] = np.sum(at * b)

    elif doalpha is False:

        at = np.empty((ncoeffs, ndat))
        for i in range(ncoeffs):
            at[i, :] = xx ** exp[i, 0] * yy ** exp[i, 1] / zzunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)

    # Solve things 

    coeffs = np.linalg.solve(alpha, beta)

    return coeffs
