import numpy as np

def nantrim(arr,flag=0,trim=False):

    '''
    To clip NaNs from an array.


    Input Parameters
    ----------------
    arr : array_like
        An (ndat,) array

    flag : {0,1,2,3}, optional
        0 - The trailing NaNs are removed (default).
        1 - The leading NaNs are removed.
        2 - Both leading and trailing NaNs are removed.
        3 -  All NaNs are removed.

    trim : {False, True}, optional
        If True, the program will return the trimmed array instead of a bool
        array.

    Returns
    -------
    numpy.ndarray
         A bool array where True indicates the number is not a NaN and False 
         indicates the number is a NaN. 

    Notes
    -----
    None

    Examples
    --------
    > x = np.array([np.nan,2,3,4,np.nan,7,89,90,np.nan])
    > nantrim(x,flag=0,trim=True)
      [nan  2.  3.  4. nan  7. 89. 90.]
    > nantrim(x,flag=1,trim=True)
      [ 2.  3.  4. nan  7. 89. 90. nan]
    > nantrim(x,flag=2,trim=True)
      [ 2.  3.  4. nan  7. 89. 90.]
    > nantrim(x,flag=3,trim=True)
      [ 2.  3.  4.  7. 89. 90.]

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_nantrim.pro.

    '''
    
# Test for correct flags

    if flag not in [0,1,2,3]:

        raise ValueError('unknown flag.')

# Convert to numpy array
    
    arr = np.array(arr)

# Test for 1D
    
    if arr.ndim != 1:

        raise ValueError('must be 1D.')        

    ndat = len(arr)
    
# test the easy cases of non NaNs

    test = np.isnan(arr)

    if flag == 3 or sum(test) == 0:

        z = np.invert(test)
        
    if flag == 0:

        cumsum = np.nancumsum(np.flip(arr))
        z = np.flip(cumsum) != 0
        
    elif flag == 1:

        cumsum = np.nancumsum(arr)
        z = cumsum !=0

    elif flag == 2:

        cumsum = np.nancumsum(arr)
        zleading = cumsum != 0

        cumsum = np.nancumsum(np.flip(arr))
        ztrailing = np.flip(cumsum) != 0

        z = zleading*ztrailing

    if trim is True:

        return(arr[z])

    else:

        return(z)
