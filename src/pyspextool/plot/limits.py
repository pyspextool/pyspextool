import numpy as np
import numpy.typing as npt
from astropy.visualization import PercentileInterval, ZScaleInterval, MinMaxInterval

from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.io.check import check_parameter


def buffer_range(range:tuple,
                 frac:float=0.1):

    """
    To expand a numerical range by a given fraction.

    Typically used to make a nice y-axis plot range.

    Parameters
    ----------
    range : tuple
        (2,) range to be expanded.

    frac : float, default 0.1
        The fraction by which to expand the range if desired.

    Returns
    -------
    tuple
         (2,) tuple of the expanded range.

    Examples
    --------
    > bufrange((1,11), frac=0.1)
      (0.0, 12.0)

    """
    
    # Do the calculation
    
    delta = range[1]-range[0]

    return range[0]-delta*frac, range[1]+delta*frac



def get_image_range(arr, info):

    """
    a wrapper for astropy to get image ranges

    Parameters
    ----------
    arr : numpy.ndarray
        an array (typically a 2D image)

    info : float or str
        float - the fraction of pixels to keep in the range (0,1)
        str - 'zscale' or 'minmax'

    Returns
    -------
    tuple
         interval based on user request

    Procedure
    ---------
    calls astropy.visualization routines and uses all default settings

    Examples
    --------
    later

    """
    
    if isinstance(info, float) is True:

        # This is if the users requests a fraction

        interval = PercentileInterval(info, n_samples=int(np.size(arr)/500))
        range = interval.get_limits(arr)
        return range

    elif info == 'zscale':

        # This is if the users requests the zscale
        
        interval = ZScaleInterval()
        range = interval.get_limits(arr)
        return range

    elif info == 'minmax':

        # This is if the users requests the zscale

        interval = MinMaxInterval()
        range = interval.get_limits(arr)
        return range

    else:

        return None

def get_spectra_range(*args,
                      robust=False,
                      savgol_window:int=11,
                      frac:float=0.0):

    """
    To return a nice y-axis plot range

    Parameters
    ----------
    *args : array-like
        Array(s) of y values to be plotted.

    robust : bool, default False
        Set to True to deal with outliers.  Each spectrum is smoothed
        with a robust Savitsky-Golay filter to "fix" bad before identifying
        the minimum and maximum.
        Set to False to keep outliers in the spectra.

    savgol_window : int, default 11
        The number of 
    
    frac : float, default 0.0
        The fraction by which to expand the range if desired.

    Returns
    -------
    tuple
         (2,) tuple of the range.

    Notes
    -----
    The minimum and maximum values of `*args` are first determined.
    If `frac` is not 0, then the range is expanded by `frac`*(max-min).


    """

    #
    # Check parameters
    #

    check_parameter('get_spectra_range', 'args', args, 'tuple')

    check_parameter('get_spectra_range', 'robust', robust, 'bool')

    check_parameter('get_spectra_range', 'savgol_window', savgol_window, 'int')
    
    check_parameter('get_spectra_range', 'frac', frac, ['int','float'])
    
    # Convert each element into a numpy.ndarray

    longlist = []
    for arg in args:

        if robust is True:

            array = np.ravel(arg)
            x_values = np.arange(len(array))
            tmp = list(robust_savgol(x_values, array, savgol_window)['fit'])
            
        else:

            tmp = list(np.ravel(arg))
        
        longlist = longlist+tmp
    
    # Get the min and max value, ignoring NaNs

    minval = np.nanmin(longlist)

    maxval = np.nanmax(longlist)

    # Expand the range if asked
    
    range = buffer_range((minval, maxval), frac=frac)

    return range


def get_stack_range(stack:npt.ArrayLike,
                    savgol:bool=False,
                    savgol_window:int=11,
                    frac:float=0.0):

    """
    To obtain the plotting range of a stack of data.

    Parameters
    ----------
    stack : ndarray

    sg_window : int, default=11

    frac : int or float, default=0.0

    Returns
    -------
    tuple

    """

    #
    # Check parameters
    #

    check_parameter('get_stack_range', 'stack', stack, 'ndarray')

    check_parameter('get_stack_range', 'savgol', savgol, 'bool')

    check_parameter('get_stack_range', 'savgol_window', savgol_window, 'int')

    check_parameter('get_stack_range', 'frac', frac, ['int', 'float'])    

    #
    # Determine useful things and make useful things.
    #

    
    norders, npixels = np.shape(stack)
    x_values = np.arange(npixels)
    duplicate = np.copy(stack)

    #
    # Smooth the spectra if requested
    #
    
    if savgol is True:

        for i in range(norders):

            result = robust_savgol(x_values, stack[i,:], savgol_window)
            duplicate[i,:] = result['fit']

    #
    # Find the range
    #
            
    minval = np.nanmin(duplicate)
    maxval = np.nanmax(duplicate)    

    #
    # Expand the range if asked
    #
    
    rng = buffer_range((minval, maxval), frac=frac) 

    return rng


