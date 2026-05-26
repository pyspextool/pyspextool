import numpy as np
import numpy.typing as npt
from scipy import interpolate

from pyspextool.io.check import check_parameter
from pyspextool.utils.math import bit_set
from pyspextool.utils.loop_progress import loop_progress


def rectify_orders(    
    image:npt.ArrayLike,
    indices:list,
    interpolation_method:str='cubic',
    variance:npt.ArrayLike=None,
    badpixel_mask:npt.ArrayLike=None,
    flag_mask:npt.ArrayLike=None,
    ybuffer:int=0,
    nbits:int=8):

    """
    To rectify spectral orders.

    The function "straightens" a spectral order onto a uniform rectangular 
    ndarray of size (nwavelengths, nangles).

    Parameters
    ----------

    image : ndarray 
        An (nrows, ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    indices : list
        An (norders,) list of dictionaries with the following keys:

        'order' : int
            The order number.
    
        'w' : ndarray
            An (nwavelengths,) array of wavelengths on which the order 
            is rectified.

        'a' : ndarray
            An (nangles,) array of angles on which the rectified 
            order.

        'xidx': ndarray
            An (nwavelengths, nangles) array of zero-based x positions in `image`
            at which to interpolate.

        'yidx': ndarray
            An (nwavelengths, nangles) array of zero-based y positions in `image`
            at which to interpolate.

    interpolation_method : {'cubic', 'linear'}
        A string giving the interpolation method passed to 
        sci.interpolate.RegularGridInterpolator.  

    badpixel_mask : ndarray
        An (nrows, ncols) bad pixel mask.  Good=1, bad=0.

    flag_mask : ndarray
        An (nrows, ncols) flag array.  This is a bit-set array with values 
        up to `nbit`.

    ybuffer : int
        Number of pixels at the top and bottom of the rectified image set 

    nbits : int, default=8
        The number of bits used in each pixel of `flag_mask`.

    Returns
    -------
    list
        An (norders,) list of dictionaries with the following keys:

        'wavelengths' : ndarray
            An (nwavelengths,) array of wavelengths on which the order 
            is rectified.

        'angles' : ndarray
            An (nangles,) array of angles on which the rectified 
            order.

        'image' : ndarray
            An (nangles, nwavelengths) array of interpolated values from `image` 
            at positions indices['xidx'] and indices['yidx'].

        'variance' : ndarray, None
            An (nangles, nwavelengths) array of interpolated values from `variance` 
            at positions indices['xidx'] and indices['yidx'].

        'badpixel_mask' : ndarray, None
            An (nangles, nwavelengths) array of interpolated values from 
            `badpixel_mask` at positions indices['xidx'] and indices['yidx'].  

        'flag_mask' : ndarray, None
            An (nangles, nwavelengths) array of interpolated values from 
            `flagmask_mask` at positions indices['xidx'] and indices['yidx'].  

    """

    #
    # Check the parameters
    #

    check_parameter('rectify_orders','image', 
                    image, 'ndarray', 2)

    check_parameter('rectify_orders','indices', 
                    indices, 'list', 2)

    check_parameter('rectify_orders','interpolation_method', 
                    interpolation_method, 'str', possible_values=['linear', 'cubic'])

    check_parameter('rectify_orders','variance', 
                    variance, ['NoneType', 'ndarray'])

    check_parameter('rectify_orders','badpixel_mask', 
                    badpixel_mask, ['NoneType', 'ndarray'], 2)    

    check_parameter('rectify_orders','flag_mask', 
                    flag_mask, ['NoneType', 'ndarray'], 2)    
    
    check_parameter('rectify_orders','ybuffer', 
                    ybuffer, 'int')

    check_parameter('rectify_orders','nbits', 
                    nbits, 'int')

    
    # Get basic info and create basic things

    nrows, ncols = image.shape

    points = (np.arange(nrows), np.arange(ncols))

    #
    # Get the functions defined first
    #


    image_function = interpolate.RegularGridInterpolator(
        points, 
        image,
        method=interpolation_method)

    if variance is not None:

        variance_function = interpolate.RegularGridInterpolator(
            points, 
            variance,
            method=interpolation_method)

    if badpixel_mask is not None:

        badpixel_function = interpolate.RegularGridInterpolator(
            points, 
            badpixel_mask,
            fill_value=1)

    if flag_mask is not None:

        
        flag_function = []
        for i in range(nbits):

            set = bit_set(flag_mask, i)

            # Do the resampling
    
            f = interpolate.RegularGridInterpolator(
                points, 
                set, fill_value=0)
            
            flag_function.append(f)

    #
    # Now do the rectifications
    #

    rectorders = []
    for order in indices:

        # Do the image first

        rimg = image_function((order['yidx'],order['xidx']))
    
        ny, nx = rimg.shape

        # Do the buffering if requested

        if ybuffer > 0:

            rimg[0:ybuffer,:] = np.tile(rimg[ybuffer,:],(ybuffer,1))
            rimg[ny-ybuffer:,:] = np.tile(rimg[ny-ybuffer-1,:],(ybuffer,1))

        # do the variance if requested.

        rvar=None
        if variance is not None:

            rvar = variance_function((order['yidx'],order['xidx']))

            if ybuffer > 0:
                
                rvar[0:ybuffer,:] = np.tile(rvar[ybuffer,:],(ybuffer,1))
                rvar[ny-ybuffer:,:] = np.tile(rvar[ny-ybuffer-1,:],(ybuffer,1))
                
        # Do the bad pixel mask if requested.
                
        rbp = None
        if badpixel_mask is not None:

            # The interpolation alone will give values between [0,1] because
            # the original mask has just zeros or ones.  So we floor the values
            # to convert any number that isn't zero to zero.
        
            rbp = np.floor(
                badpixel_function(
                    (order['yidx'],order['xidx']))).astype('uint8')

        # Do the flag mask if requested.

        rfl = None
        if flag_mask is not None:

            rfl = np.zeros((ny,nx),dtype=np.uint8)
            for func in flag_function:

                set = np.floor(
                    func(
                        (order['yidx'],order['xidx']))).astype(np.uint8)
                                   
                mask = set > 0
                
                rfl += np.multiply(mask, (2**i), dtype='uint8')

        # Store the results for return

        rectorders.append(
            {'wavelengths':order['w'],
             'angles':order['a'],
             'image':rimg,
             'variance':rvar,
             'badpixel_mask':rbp,
             'flag_mask':rfl})

    return rectorders

        



        

                

            

    




