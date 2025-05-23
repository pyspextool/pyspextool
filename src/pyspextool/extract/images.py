import numpy as np
import numpy.typing as npt
from scipy import interpolate

from pyspextool.io.check import check_parameter
from pyspextool.utils.math import bit_set
from pyspextool.utils.loop_progress import loop_progress

def make_ordermask(ncols:int,
                   nrows:int,
                   edgecoeffs:npt.ArrayLike,
                   xranges:npt.ArrayLike,
                   orders:npt.ArrayLike,
                   ybuffer:int=0):

    """
    To create an order mask.

    Parameters
    ----------
    ncols : int
        The number of columns of the image.

    nrows : int
        The number of rows of the image.

    edgecoeffs : ndarray
        (norders,`edgedeg`+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  

    xranges : ndarray
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.

    orders : list of int
        The order numbers.  orders[0] is the order number of the 
        order nearest the bottom of the image after rotation. 

    ybuffer : int, default 0
        The number of native pixels from the top and bottom of the slit to
        avoid during the operation.  Useful to account for the fact that
        the drop-off in intensity at the edge of the slit is not a
        heaviside function but rather occurs over a few pixels. 


    Returns
    -------
    ndarray
        An (nrows, ncols) int array where each pixel is set to its order number.
        Inter-order pixels are set to zero. 

    """

    #
    # Check parameters
    #

    check_parameter('make_ordermask', 'ncols', ncols, 'int')

    check_parameter('make_ordermask', 'nrows', nrows, 'int')

    check_parameter('make_ordermask', 'edgecoeffs', edgecoeffs, 'ndarray')

    check_parameter('make_ordermask', 'xranges', xranges, 'ndarray')

    check_parameter('make_ordermask', 'orders', orders, 'ndarray')

    check_parameter('make_ordermask', 'ybuffer', ybuffer, 'int')            
    
    
    #
    # Get basic info and get set up
    #
    
    ndimen = edgecoeffs.ndim
    if ndimen == 2:
        norders = 1

    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    order_mask = np.zeros([nrows, ncols])

    y = np.arange(nrows)

    # start the loop over order and column number
    
    for i in range(norders):

        start = xranges[i, 0]
        stop = xranges[i, 1]

        x = np.arange(stop - start + 1) + start

        # Get the top and bottom positions of the slit

        botedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(x, edgecoeffs[i, 1, :])

        # Fill things in

        for j in range(stop - start + 1):

            bottom = np.floor(botedge[j]+ybuffer).astype('int')
            top = np.ceil(topedge[j]-ybuffer).astype('int')
            
            order_mask[bottom:top, x[j]] = orders[i]

    return order_mask




def rectify_order(image:npt.ArrayLike,
                  xidx:npt.ArrayLike,
                  yidx:npt.ArrayLike,
                  variance:npt.ArrayLike=None,
                  bad_pixel_mask:npt.ArrayLike=None,
                  flag_mask:npt.ArrayLike=None,
                  ybuffer:int=0,
                  nbits:int=8):

    """
    To rectify a spectral order 

    The function "straightens" a spectral order onto a uniform rectangular 
    grid

    Parameters
    ----------

    image : ndarray 
        An (nrows, ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    xidx : ndarray
        An (nwave, nspat) array giving x coordinates in `image` 
    
    yidx: ndarray

    ybuffer : int, optional
        The number of native pixels from the top and bottom of the slit to 
        avoid during the operation.  Useful to account for the fact that 
        the drop-off in intensity at the edge of the slit is not a 
        heaviside function but rather occurs over a few pixels.  

    Returns
    -------
    dict

    """

    #
    # Check the parameters
    #

    check_parameter('rectify_order','image', image, 'ndarray', 2)

    check_parameter('rectify_order','variance', variance,
                    ['NoneType', 'ndarray'])

    check_parameter('rectify_order','bad_pixel_mask', bad_pixel_mask,
                    ['NoneType', 'ndarray'])    

    check_parameter('rectify_order','flag_mask', flag_mask,
                    ['NoneType', 'ndarray'])    
    
    check_parameter('rectify_order','xidx', xidx, 'ndarray', 2)

    check_parameter('rectify_order','yidx', yidx, 'ndarray', 2)
    
    # Get basic info and create basic things

    nrows, ncols = image.shape

    points = (np.arange(nrows), np.arange(ncols))

    # Do the resampling
    
    f = interpolate.RegularGridInterpolator(points, image,method='cubic')
    image = f((yidx, xidx))

    ny, nx = image.shape

    # Do the buffering if requested

    if ybuffer > 0:

        image[0:ybuffer,:] = np.tile(image[ybuffer,:],(ybuffer,1))
        image[ny-ybuffer:,:] = np.tile(image[ny-ybuffer-1,:],(ybuffer,1))
        
    # Pack things up

    order = {'image':image}

    if variance is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, variance,method='cubic')
        var = f((yidx, xidx))

        order.update({'variance':var})

    if bad_pixel_mask is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, bad_pixel_mask,
                                                fill_value=1)

        # The interpolation alone will give values between [0,1] because
        # the original mask has just zeros or ones.  So we floor the values
        # to convert any number that isn't zero to zero.
        
        tmp = np.floor(f((yidx, xidx))).astype('uint8')

        order.update({'bpmask':tmp})        

    if flag_mask is not None:

        output = np.zeros((ny,nx),dtype=np.uint8)
        
        for i in range(nbits):

            set = bit_set(flag_mask, i)

            # Do the resampling
    
            f = interpolate.RegularGridInterpolator(points, set, fill_value=0)
            set = np.floor(f((yidx, xidx))).astype(np.uint8)
                                   
            mask = set > 0

            output += np.multiply(mask, (2**i), dtype='uint8')
        
        order.update({'flagmask':output})        

    return order



def scale_order_background(stack:npt.ArrayLike,
                           orders:npt.ArrayLike,
                           edgecoeffs:npt.ArrayLike,
                           xranges:npt.ArrayLike,
                           var_stack:npt.ArrayLike=None,
                           ybuffer:int=1,
                           verbose:bool=False):

    """
    Scale orders to a common a flux level across a data stack.

    Parameters
    ----------
    stack : ndarray
        (nrows, ncols, nimages) stack of spectral images.

    orders : ndarray of int
        (norders,) array of order numbers.  orders[0] is the order number of 
        the order nearest the bottom of the image after rotation. 

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

    var_stack : ndarray, optional
        (nrows, ncols) stack of variance images.

    ybuffer : int, default 1
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

    check_parameter('scale_order_background', 'stack', stack, 'ndarray')

    check_parameter('scale_order_background', 'orders', orders, 'ndarray')

    check_parameter('scale_order_background', 'edgecoeffs', edgecoeffs,
                    'ndarray')

    check_parameter('scale_order_background', 'xranges', xranges, 'ndarray')

    check_parameter('scale_order_background', 'ybuffer', ybuffer, 'int')

    check_parameter('scale_order_background', 'var_stack', var_stack,
                    ['ndarray', 'NoneType'])

    check_parameter('scale_order_background', 'verbose', verbose, 'bool')    

    # Do basic things

    nimages, nrows, ncols = stack.shape

    norders = len(orders)

    # Create the order mask

    order_mask = make_ordermask(ncols,
                                nrows,
                                edgecoeffs,
                                xranges,
                                orders,
                                ybuffer=ybuffer)

    # Star the loop over order number

    for i in range(norders):

        if verbose is True:

            loop_progress(i, 0, norders)

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

