import numpy as np
from scipy import interpolate

from pyspextool.io.check_parameter import check_parameter

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
    spatmap, rimg : numpy.ndarray of float, numpy.ndarray of float
    
        spatmap - each pixel is set to its angular position on the sky

        rimg - the resampled order

    Notes
    -----
    Later

    Examples
    --------
    later?

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

    order = {'img':img}

    if var is not None:

        # Do the resampling
    
        f = interpolate.RegularGridInterpolator(points, var)
        var = f((yidx, xidx))

        order.update({'var':var})

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
    

#    startcol = xranges[0]
#    stopcol = xranges[1]
#
#    order_ncols = stopcol-startcol+1
#    cols = np.arange(order_ncols, dtype=int)+startcol
#    rows = np.arange(nrows, dtype=int)
#
#    # Get the bottom and top of the slit 
#
#    botedge = np.polynomial.polynomial.polyval(cols, edgecoeffs[0, :])
#    topedge = np.polynomial.polynomial.polyval(cols, edgecoeffs[1, :])
#
#    # Get the conversion between pixels and arcseconds
#
#    dif = topedge-botedge
#    arctopix = dif/slith_arc
#    
#    # Create a uniform grid of arcsecond values for the new order
#
#    mindif = np.floor(min(dif))
#
#    nrslit = round(oversamp*mindif)  # number of pixels in resampled slit
#    
#    rslit_arc = np.arange(nrslit)*slith_arc/(nrslit-1)
#
#    # Set the pixels ybuffer from the ends to the values at ybuffer from the edge
#
#    if ybuffer > 0:
#    
#        ybuffer = round(ybuffer*oversamp)
#    
#        rslit_arc[0:ybuffer-1] = rslit_arc[ybuffer-1]
#        rslit_arc[(nrslit-ybuffer):nrslit] = rslit_arc[nrslit-ybuffer]
#        
#    # Now do the interpolation one column at a time
#
#    # NOTE: this is not the fast way, but I don't see an IDL equivalent of
#    # interpolate in python.  Looking through SOFIA's repo they seem to have
#    # written a bunch of their own, so it may not be possible with standard
#    # scipy packages.
#
#    order = np.empty((nrslit, order_ncols))
#    
#    for i in range(0, order_ncols):
#
#        # Get the bottom and top of the slit positions at cols[i].   Use floor and
#        # ceil to ensure the interpolation has no extrapolation.
#        
#        slitbot_pix = np.floor(botedge[i]).astype(int)
#        slittop_pix = np.ceil(topedge[i]).astype(int)
#
#        # Yank out the row values and convert to arcseconds
#        
#        slit_pix = rows[slitbot_pix:(slittop_pix+1)]-botedge[i]
#        slit_arc = slit_pix/arctopix[i]
#
#        # Do the interpolation and store
#        
#        f = interpolate.interp1d(slit_arc, img[slitbot_pix:slittop_pix+1,
#                                 i+startcol])
#        order[:, i] = f(rslit_arc)
#
#    # Create the spatial map using rslit_arc
#
#    spatmap = np.rot90(np.tile(rslit_arc, (order_ncols, 1)), k=3)
#        
#    return spatmap, order
