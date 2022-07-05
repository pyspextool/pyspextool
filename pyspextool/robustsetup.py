import numpy as np

def robustsetup(*args,goodbad=None):

    '''
    To set up the goodbad array for a robust fit.


    Input Parameters
    ----------------
    *args : array_like or None
        (ndat,) arrays to be used in a robust fit.  Typically the are the
        dependent variable array and its associated uncertainty.  

    goodbad : array_like of int, default=None, optional
         An (ndat,) array of values where 0=bad, 1=good, 2=NaN.

    Returns
    -------
    numpy.ndarray
         A goodbad day consistent with the properties of `goodbad` and
         `*args`.

    Notes
    -----
    The function identifies NaN values in *args and sets those pixels 
    in either a user-passed goodbad array or a function-created goodbad array
    to 2.  

    Examples
    --------
    > x= [1,2,np.nan,5,6] 
    > y = [0,1,2,3,np.nan]

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program XS

    '''
    
    if goodbad == None:

# Create a good bad array
        
        goodbad = np.full(len(args[0]),1,dtype=int)

    else:

# Convert to a numpy array if need be
        
        goodbad = np.array(goodbad,dtype=int)        

# Find NaNs if need be

    for arg in args:

        if arg is not None:
        
            znan = np.isnan(arg)
            if np.sum(znan) > 0:
            
                goodbad[znan] = 2

    return goodbad
