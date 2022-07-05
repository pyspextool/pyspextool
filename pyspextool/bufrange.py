import numpy as np

def bufrange(range,frac=0.1):

    '''
    To expand a range by a given fraction.

    Typically used to make a nice y-axis plot range.

    Input Parameters
    ----------------
    range : tuple
        (,2) range to be expanded.

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
    > bufrange((1,11),frac=0.1)
      (0.0, 12.0)

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_burange.pro.

    '''
    
    # Do the calculation
    
    delta = range[1]-range[0]

    return (range[0]-delta*frac,range[1]+delta*frac)

    
    
    
