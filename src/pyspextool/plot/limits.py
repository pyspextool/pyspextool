import numpy as np
from astropy.visualization import PercentileInterval, ZScaleInterval, MinMaxInterval


def buffer_range(range, frac=0.1):

    """
    To expand a numerical range by a given fraction.

    Typically used to make a nice y-axis plot range.

    Input Parameters
    ----------------
    range : tuple
        (2,) range to be expanded.

    frac : float, default 0.1
        The fraction by which to expand the range if desired.

    Returns
    -------
    tuple
         (2,) tuple of the expanded range.

    Notes
    -----
    None

    Examples
    --------
    > bufrange((1,11), frac=0.1)
      (0.0, 12.0)

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_bufrange.pro.

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

    Modification History
    --------------------
    2022-06-07 - Written by M. Cushing, University of Toledo.

    """
    
    if isinstance(info, float) is True:

        # This is if the users requests a fraction

        interval = PercentileInterval(info)
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

def get_spec_range(*args, frac=0.0):

    """
    To return a nice y-axis plot range

    Parameters
    ----------------
    *args : array-like
        Array(s) of y values to be plotted.

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

    Examples
    --------
    later

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_burange.pro.

    """

    # Convert each element into a numpy.ndarray

    longlist = []
    for arg in args:

        longlist = longlist+list(np.ravel(arg))
    
    # Get the min and max value, ignoring NaNs

    minval = np.nanmin(longlist)

    maxval = np.nanmax(longlist)

    # Expand the range if asked
    
    range = buffer_range((minval, maxval), frac=frac)

    return range    


    
