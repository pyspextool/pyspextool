import numpy as np
from bufrange import bufrange

def getyrange(*args,frac=0.1):

    '''
    To return a nice y-axis plot range


    Input Parameters
    ----------------
    *args : array-like
        Array(s) of y values to be plotted.

    frac : float, default 0.1
        The fraction by which to expand the range if desired.

    Returns
    -------
    tuple
         (2,) tuple of the range.

    Notes
    -----
    The minumum and maximum values are first determined.  If `frac` is not 0,
    then the range is expanded by `frac`*(max-min).  

    Examples
    --------
    later

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_burange.pro.

    '''

# Convert each element into a numpy.ndarray

    longlist = []
    for arg in args:

        longlist = longlist+list(np.ravel(arg))
    
# Get the min and max value, ignoring NaNs

    minval = np.nanmin(longlist)

    maxval = np.nanmax(longlist)

# Expand the range if asked
    
    range = bufrange((minval,maxval),frac=frac)

    return range
