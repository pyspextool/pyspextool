import numpy as np
from scipy import interpolate

from pyspextool.io.check import check_parameter

def rectify_order(img, xidx, yidx, var=None, bpmask=None, bsmask=None,
                  ybuffer=0):

    """
    To rectify a spectral order 

    The function "straightens" a spectral order onto a uniform rectangular 
    grid

    Parameters
    ----------

    img : numpy.ndarray 
        An (nrows, ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    xidx : numpy.ndarray

    yidx: numpy.ndarray

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
    check_parameter('rectify_order1d','img', img, 'ndarray', 2)
    
    check_parameter('rectify_order1d','xidx', xidx, 'ndarray', 2)

    check_parameter('rectify_order1d','yidx', yidx, 'ndarray', 2)
    
    # Get basic info and create basic things

    nrows, ncols = img.shape

    points = (np.arange(nrows), np.arange(ncols))

    # Do the resampling
    
    f = interpolate.RegularGridInterpolator(points, img)
    img = f((yidx, xidx))

    ny, nx = img.shape

    # Do the buffering if requested

    if ybuffer > 0:

        img[0:ybuffer,:] = np.tile(img[ybuffer,:],(ybuffer,1))
        img[ny-ybuffer:,:] = np.tile(img[ny-ybuffer-1,:],(ybuffer,1))
        
    # Pack things up

    order = {'image':img}

    if var is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, var)
        var = f((yidx, xidx))

        order.update({'variance':var})

    if bpmask is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, bpmask, fill_value=1)
        bpmask = f((yidx, xidx))

        order.update({'bpmask':bpmask})        

    if bsmask is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, bsmask, fill_value=0)
        bsmask = f((yidx, xidx))

        order.update({'bsmask':bsmask})        

    return order    
