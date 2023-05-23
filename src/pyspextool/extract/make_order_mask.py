import numpy as np

def make_order_mask(ncols, nrows, edgecoeffs, xranges, orders, ybuffer=0):

    """
    To create an order mask.

    Each pixel in an order is set to its order number.


    Parameters
    ----------
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

    orders : list of int
        The order numbers.  orders[0] is the order number of the 
        order nearest the bottom of the image after rotation. 

    ybuffer : int, default=0
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels. 


    Returns
    -------
    omask : numpy.ndarray
        An (nrows,ncols) int array where each pixel is set to its order number.
        Inter-order pixels are set to zero. 

    Notes
    -----
    None

    Examples
    --------
    later

    """

    # Get basic info

    ndimen = edgecoeffs.ndim
    if ndimen == 2:
        norders = 1

    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    # Make some NaN arrays

    order_mask = np.zeros([nrows, ncols])

    # start the loop over order and column number

    y = np.arange(nrows)

    for i in range(norders):

        start = xranges[i, 0]
        stop = xranges[i, 1]

        x = np.arange(stop - start + 1) + start

        # Get the top and bottom positions of the slit

        botedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 1, :])

        # Fill things in

        for j in range(stop - start + 1):

            order_mask[np.floor(botedge[j]+ybuffer).astype('int'):
                    np.ceil(topedge[j]-ybuffer).astype('int'), x[j]] = orders[i]

    return order_mask
