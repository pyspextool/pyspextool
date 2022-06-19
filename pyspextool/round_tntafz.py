import numpy as np

def round(x):

    '''
    To round numbers in a numpy array

   
    Input Parameters
    ----------------
    x : numpy.ndarray


    Returns
    --------
    numpy.ndarray like `x`
        The numbers are rounded to the nearest integer, and tied away \
        from zero. 

    Notes
    -----

    Based on:  https://stackoverflow.com/questions/53201470/how-to-always-round-up-a-xx-5-in-numpy

    Examples
    --------
    > import numpy as np
    > round(np.array([-3.5,-2.5,2.5,3.5]))
      [-4. -3.  3.  4.]

    Modification History
    --------------------
    2022-06-19 - Written by M. Cushing, University of Toledo.

    '''
    
    mask = (x >= 0)
    out = np.empty_like(x)
    out[mask] = np.floor(x[mask] + 0.5)
    out[~mask] = np.ceil(x[~mask] - 0.5)

    return out
