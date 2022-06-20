import numpy as np
from math import nan as mnan

def medcomb(data,mask=None,stderr=True):

    """
    Median a spectral or image stack with optional mask


    Input Parameters
    ----------------
    data : numpy.ndarray
        either a stack of spectra [nspec,npoints] or a stack of images 
        [nimgs,nrows,ncols].  

    mask : numpy.ndarray, optional
        a mask array with the same shape as `data`.  
        0 = bad, 1=good

    stderr : {True, False}, optional
        Set to return 1.4826*MAD/sqrt(n) instead of just 1.4826*MAD 
        (see Procedure)

    Returns
    --------
    list
        list[0] : numpy.ndarray 
            the median of the spectral or image stack

        list[1] : numpy.ndarray
            the uncertainty of the spectral or image stack (see Procedure)

    Procedure
    ---------
    Spectral stack:

        in this case, the data have the shape [nspec,npoints].  The 
        median of the stack is computed producing an array of size 
        [npoints]. At each spectral point, the median absolute deviation 
        (MAD=median(|data-med}) is computed.  

    Image stack:

        in this case, the data have the shape [nspec,nrows,ncols].  The 
        median of the stack is computed producing an array of size 
        [nrows,ncols]. At each image point, the median absolute deviation 
        (MAD=median(|data-med}) is computed.  


    The estimate of the standard deviation assuming the data arises from 
    a gaussian is given by 1.4826*MAD.  Finally, if stderr is set, then the 
    standard error is computed as 1.4826*MAD/root(n), where n is the number 
    of data points at a given spectral point.

    Note:  Points excluded by the mask are ignore in all calculations.

    Examples
    --------

    > import numpy as np
    > ss = np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
    > msk = np.ones((3,4),dtype=int)
    > msk[0,0] = 0
    > print(ss)
    > print(msk)
    > med,unc=medcomb(ss,mask=msk)
    > print(med)
    > print(unc)

      [[1 2 3 7]
       [0 5 8 2]
       [2 9 4 6]]
      [[0 1 1 1]
       [1 1 1 1]
       [1 1 1 1]]
      [1. 5. 4. 6.]
      [1.04835651 2.56793853 0.85597951 0.85597951]

    > istack = np.array([[[1,2,3],[4,5,6],[7,8,9]],\
                         [[6,3,1],[9,2,4],[1,5,0]],\
                         [[3,4,9],[5,7,7],[3,9,1]],\
                         [[1,6,5],[2,1,9],[5,2,7]]])              
    > msk = np.ones((4,3,3),dtype=int)
    > msk[0,0,0] = 0
    > print('Image Stack Test')
    > print(istack)
    > print(msk)
    > med,unc=medcomb(istack,mask=msk)
    > print(med)
    > print(unc)

      [[[1 2 3]
        [4 5 6]
        [7 8 9]]

       [[6 3 1]
        [9 2 4]
        [1 5 0]]

       [[3 4 9]
        [5 7 7]
        [3 9 1]]

       [[1 6 5]
        [2 1 9]
        [5 2 7]]]
      [[[0 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]]
      [[3.  3.5 4. ]
       [4.5 3.5 6.5]
       [4.  6.5 4. ]]
      [[1.71195902 0.7413     1.4826    ]
       [1.11195    1.4826     1.11195   ]
       [1.4826     1.4826     2.59455   ]]

    Modification History
    --------------------
    2022-06-01 - Written by M. Cushing, University of Toledo.
        Based on the Spextool IDL program mc_medcomb.pro.

    """

# Get array dimensions

    ndimen = np.ndim(data)
    shape = np.shape(data)

# If no mask passed, create one.

    if mask is None:

        mask = np.ones(shape,dtype=int)

# Now search and replace any masked pixels with NaNs
        
    data = np.where(mask != 0, data,mnan)

# Spectral or image stack?

    if ndimen == 2:

        tileshape = (shape[0],1) # spectral stack

    elif ndimen == 3:

        tileshape = (shape[0],1,1) # image stack

    else:

        print('Unknown data shape.')
        return
    
# Compute the median and median absolute deviation

    med = np.nanmedian(data,axis=0)

    mad = np.nanmedian(np.abs(data-np.tile(med,tileshape)),axis=0)

    if stderr is not None:
            
        mad *= 1.4826   # assume gaussian distribution
        unc = mad/np.sqrt(np.sum(mask, axis=0))
            
    else:

        unc = mad / 1.4826   # assume gaussian distribution            

    return[med,unc]
