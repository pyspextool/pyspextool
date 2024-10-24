import numpy as np
import numpy.typing as npt

from pyspextool.io.check import check_parameter


def median_1dxd(image:npt.ArrayLike,
                edgecoeffs:npt.ArrayLike,
                xranges:npt.ArrayLike,
                var:npt.ArrayLike=None,
                ybuffer:int=1):
    
    """
    To background subtraction a 1DXD image.

    The median value is subtracted off each column of each order.

    Parameters
    ----------

    image : ndarray
        (nrows, ncols) image.

    edgecoeffs : ndarray
        (norders,2,ncoeffs) array giving the polynomial 
        coefficients delineating the top and bottom of each order.  
        edgecoeffs[0,0,:] gives the coefficients for the bottom of 
        the order closests to the bottom of `img` and 
        edgecoeffs[0,1,:] gives the coefficients for the top of said 
        order.  

    xranges : ndarray of int
        An (norders,2) array giving the column numbers over which to 
        search.  sranges[0,0] gives the starting column number for 
        the first order and sranges[0,1] gives the end column number 
        for the first order.

    var : ndarray, default None
        (nrows, ncols) variance image.    

    ybuffer : int, default 1
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels.

    Returns
    -------
    ndarray, optional ndarray
    The background subtracted image where the median intensity is subtracted
    from each column in each row.

    If var is not None, then (1.4826 * MAD) ** 2, where MAD is the median
    absolute deviation given by median(abs(med-data)), is added to each
    column.  
    
    """

    #
    # Check parameters
    #

    check_parameter('median_1dxd', 'image', image, 'ndarray')

    check_parameter('median_1dxd', 'edgecoeffs', edgecoeffs, 'ndarray')

    check_parameter('median_1dxd', 'xranges', xranges, 'ndarray')

    check_parameter('median_1dxd', 'ybuffer', ybuffer, 'int')

    check_parameter('median_1dxd', 'var', var, ['ndarray', 'NoneType'])

    # Do basic things

    norders = np.shape(xranges)[0]

    # Star the loop over the orders

    for i in range(norders):

        # Determine the start and stop column numbers

        startcol = xranges[i, 0]
        stopcol = xranges[i, 1]

        order_ncols = stopcol - startcol + 1
        order_cols = np.arange(order_ncols, dtype=int) + startcol

        # Determine the top and the bottom of the slit

        botedge = np.polynomial.polynomial.polyval(order_cols,
                                                   edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(order_cols,
                                                   edgecoeffs[i, 1, :])

        # Start the loop over the columns

        for j in range(order_ncols):

            # Get the indices

            bot_idx = int(np.ceil(botedge[j]) + ybuffer)
            top_idx = int(np.floor(topedge[j]) - ybuffer)

            # Get the median and do the subtraction

            med = np.median(image[bot_idx:top_idx, j])

            # Do the subtraction

            image[bot_idx:top_idx, j] -= med

            # Do the error propagation if need be

            if var is not None:
                mad = np.median(np.abs(image[bot_idx:top_idx+1, j] - med))
                var[bot_idx:top_idx+1, j] += (1.4826 * mad) ** 2

    # Return the results

    if var is None:

        return image

    else:

        return image, var
