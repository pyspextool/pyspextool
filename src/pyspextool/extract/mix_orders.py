import numpy as np
import numpy.typing as npt

from pyspextool.io.check import check_parameter
from pyspextool.pyspextoolerror import pySpextoolError


def mix_orders(image1:npt.ArrayLike,
               image2:npt.ArrayLike,
               order_mask:npt.ArrayLike,
               orders:npt.ArrayLike | list,
               orders_from_image1:npt.ArrayLike | list,
               orders_from_image2:npt.ArrayLike | list,
               interorder_value=np.nan):

    """
    To create an image with orders from two different images.

    Parameters
    ----------
    image1 : ndarray
        An (nrows, ncols) array.

    image2 : ndarray
        An (nrows, ncols) array.

    order_mask : ndarray
        An (nrows, ncols) array where each pixel is set to its order number.

    orders : ndarray
        An (norders,) array of the orders numbers in order_mask.

    orders_from_image1: ndarray
        A array of orders numbers to grab from image 1.

    orders_from_image2: ndarray
        A array of orders numbers to grab from image 2.

    interorder_value : float, optional
        The value to use for inter-order pixels.  The default is np.nan.
             
    Returns
    -------
    ndarray

       
    """

    #
    # Check parameters
    #

    check_parameter('mix_orders', 'image1', image1, 'ndarray')

    check_parameter('mix_orders', 'image2', image2, 'ndarray')

    check_parameter('mix_orders', 'order_mask', order_mask, 'ndarray')

    check_parameter('mix_orders', 'orders', orders, ['ndarray','list'])

    check_parameter('mix_orders', 'orders_from_image1', orders_from_image1,
                    ['ndarray','list'])                

    check_parameter('mix_orders', 'orders_from_image2', orders_from_image2,
                    ['ndarray','list'])

    check_parameter('mix_orders', 'interorder_value', interorder_value,
                    ['int','float'])

    #
    # Do checking to make sure the requested orders exist.
    #

    # Check the orders from image 1

    if np.sum(np.in1d(orders_from_image1, orders)) != len(orders_from_image1):

        message = 'Order number in parameter `orders_from_image1` is not in ' \
            '`orders`.'
        raise pySpextoolError(message)
    
    # Check the orders from image 2

    if np.sum(np.in1d(orders_from_image2, orders)) != len(orders_from_image2):

        message = 'Order number in parameter `orders_from_image2` is not in ' \
            '`orders`.'
        raise pySpextoolError(message)

    # Check whether threre is overlap in orders_from_image1 and
    # orders_from_image2

    if np.sum(np.in1d(orders_from_image2, orders_from_image1)) != 0:

        message = 'At least one order number is in both parameters ' \
            '`orders_from_image1` and `orders_from_image2`.'
        raise pySpextoolError(message)
       
    #
    # Loop over each request
    #

    output_image = np.full_like(image1, interorder_value, dtype=float)
    mask = np.full_like(image1, 1, dtype=float)    

    for order in orders_from_image1:

        z = np.where(order_mask == order)
        output_image[z] = image1[z]

    for order in orders_from_image2:
       
        z = np.where(order_mask == order)
        output_image[z] = image2[z]

    #
    # Return the new array
    #

    return output_image
