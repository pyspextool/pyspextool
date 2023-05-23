import numpy as np

from pyspextool.io.check import check_parameter
from pyspextool.extract.make_order_mask import make_order_mask


def scale_orders(stack, orders, edgecoeffs, xranges, var_stack=None, ybuffer=1):
    """
    Scale orders to a common a flux level across a data stack.

    Parameters
    ----------
    stack : ndarray
        (nrows, ncols, nimages) stack of spectral images.

    orders : ndarray of int
        (norders,) array of order numbers.  orders[0] is the order number of 
        the order nearest the bottom of the image after rotation. 

    edgecoeffs : array_like of float
        (norders,2,ncoeffs) array giving the polynomial 
        coefficients delineating the top and bottom of each order.  
        edgecoeffs[0,0,:] gives the coefficients for the bottom of 
        the order closests to the bottom of `img` and 
        edgecoeffs[0,1,:] gives the coefficients for the top of said 
        order.  

    xranges : array_like of int, [norders,2] 
        An (norders,2) array giving the column numbers over which to 
        search.  sranges[0,0] gives the starting column number for 
        the first order and sranges[0,1] gives the end column number 
        for the first order.

    var_stack : ndarray, optional
        (nrows, ncols) stack of variance images.

    ybuffer : int, default=1
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels.

    Returns
    -------
    Later

    """

    #
    # Check parameters
    #

    check_parameter('scale_orders', 'stack', stack, 'ndarray')

    check_parameter('scale_orders', 'orders', orders, 'ndarray')

    check_parameter('scale_orders', 'edgecoeffs', edgecoeffs, 'ndarray')

    check_parameter('scale_orders', 'xranges', xranges, 'ndarray')

    check_parameter('scale_orders', 'ybuffer', ybuffer, 'int')

    check_parameter('scale_orders', 'var_stack', var_stack,
                    ['ndarray', 'NoneType'])

    # Do basic things

    nimages, nrows, ncols = stack.shape

    norders = len(orders)

    # Create the order mask

    order_mask = make_order_mask(ncols, nrows, edgecoeffs, xranges, orders,
                                 ybuffer=ybuffer)

    # Star the loop over order number

    for i in range(norders):

        # Get the median values of the orders for each image

        meds = np.zeros(nimages)

        for j in range(nimages):
            # Grab the image and protect against gettting a median of zero.

            image = stack[j, :, :].copy()  # need a deep copy

            z = np.where(image == 0.0)
            image[z] = np.nan

            # Find the order pixels and compute the median

            z = np.where(order_mask == orders[i])
            meds[j] = np.nanmedian(image[z])

        # Compute the median of the median value

        med_all = np.median(meds)

        # Compute the scale factor

        scales = med_all / meds

        for j in range(nimages):

            # Grab the image 

            image = stack[j, :, :]

            # Find the order pixels

            z = np.where(order_mask == orders[i])

            # Scale the values and store

            image[z] *= scales[j]

            # Don't need to store because I did a shallow copy of the stack

            if var_stack is not None:
                image = var_stack[j, :, :]
                image *= scales[j] ** 2

                # Don't need to store because I did a shallow copy of the
                # stack                

    if var_stack is None:

        return stack

    else:

        return stack, var_stack
