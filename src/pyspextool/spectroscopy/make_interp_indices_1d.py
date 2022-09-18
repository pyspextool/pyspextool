import numpy as np

from pyspextool.fit.polyfit import poly_1d
from pyspextool.io.check_parameter import check_parameter

def make_interp_indices_1d(edgecoeffs, xranges, slith_arc, array_output=False):

    """
    To generate 1D indices for resampling of a spectral order.

    The program generates indices that can be used with 
    scipy.interpolate.RegularGridInterpolator to "straighten" a spectral order.
    The "1D" refers to the fact that it is assumed that each column of the 
    detector represents a single wavelength.  


    Parameters
    ----------
    edgecoeffs : numpy.ndarray
        A (`edgedeg`+1,2) float array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,:]
        gives the coefficients for the bottom of the order and 
        edgecoeffs[1,:] gives the coefficients for the top of the order.  

    xranges : numpy.ndarray
        An (2, ) float array giving the column numbers over which to 
        operate.  xranges[0] gives the starting column number for the 
        order and xranges[1] gives the end column number for the order.

    slith_arc : float
        The nominal slit height (arcseconds).

    array_output : {False, True}, optional
        Set to pack the results into a single array (see Returns).

    Returns
    -------
    array_output=False
        wavemap : numpy.ndarray
        The (ny, nx) array of columns associated with the resampled order
        spatmap : numpy.ndarray
        The (ny, nx) spatial positions associated with the resampled order 
        in arcseconds.
        xidx : numpy.ndarray
        The (ny, nx) array of x coordinates.
        yidx : numpy.ndarray
        The (ny, nx) array of y coordinates.

    numpy.ndarray 
        A (2, nx+1, ny+1) float array.          
        [0,0,1:] = x_pix
        [0,1:,0] = y_arc
        [0,:,:] = the array of x coordinates.
        [1,0,1:] = x_pix
        [1,1:,0] = y_arc
        [1,:,:] = the array of y coordinates.    

    Examples
    --------
    later

    """

    # Check the parameters

    check_parameter('make_interp_indices_1d', 'edgecoeffs', edgecoeffs,
                    'ndarray', 2)
    check_parameter('make_interp_indices_1d', 'xranges', xranges, 'ndarray', 1)
    check_parameter('make_interp_indices_1d', 'slith_arc', slith_arc,
                    ['int', 'float'])
        
    # Generate the the positions at the top and bottom of the slit
    
    x_pix = np.arange(xranges[0], xranges[1]+1)
    top = poly_1d(x_pix, edgecoeffs[1,:])
    bot = poly_1d(x_pix, edgecoeffs[0,:])    

    # Find the minimum number of pixels spanned by the slit
    
    dif = top-bot
    ny = np.floor(np.min(dif)).astype(int)

    #
    # Generate the rectification grid
    # 

    # Build the xidx array
    
    xidx = np.tile(x_pix, (ny, 1))

    # Now do the yidx array.

    y_pix = np.arange(ny)

    slope = dif/(ny-1)
    nx = len(x_pix)

    yimg = np.tile(np.reshape(y_pix,(ny, 1)), (1, nx))
    
    scale = np.tile(slope, (ny, 1))
    zpt = np.tile(bot, (ny, 1))    

    yidx = yimg*scale+zpt

    y_arc = y_pix/y_pix[-1]*slith_arc

    #
    # Now return the results
    #
    
    if array_output is False:

        wavemap = np.tile(x_pix,(ny,1))
        spatmap = np.rot90(np.tile(y_arc, (nx, 1)), k=3)
        
        return xidx, yidx, wavemap, spatmap

    else:
        
        indices = np.full((2, ny+1, nx+1),np.nan)
        
        indices[0,1:, 1:] = xidx
        indices[0,0,1:] = x_pix
        indices[0,1:,0] = y_arc
        
        indices[1,1:, 1:] = yidx
        indices[1,0,1:] = x_pix
        indices[1,1:,0] = y_arc
        
        return indices
