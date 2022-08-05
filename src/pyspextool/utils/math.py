"""Functions for math utlitiy management."""


import numpy as np
from scipy.stats import describe


def find_index(x_array, x):

    """
    Finds the effective index of a function value in an ordered array.

    Parameters
    ----------------
    x_array : array_like
        (N,) The array of floats or integers to be searched,
        must be monotonically increasing or decreasing.
    x : array_like or float or int
        (M,) The value or values whose indices are required.

    Returns
    -------
    numpy.ndarray or float
        Effective indices in `x_array` for `x`.

    Raises
    ------
    ValueError
        If input arrays are not array_like or
        if input array is not monotonic

    Notes
    -----
    Near line-by-line copy of IDL Astronomy User's Library routine tabinv.pro. 

    From its documentation:  A binary search is used to find the values XARR(I) 
    and XARR(I+1) where XARR(I) < X < XARR(I+1).  IEFF is then computed using 
    linear interpolation between I and I+1.
    IEFF = I + (X-XARR(I)) / (XARR(I+1)-XARR(I)).

    IEFF values outside of `xarr` are set to NaN.  

    Examples
    --------
    >>> x_array = [1, 2.5, 3, 5.5]
    >>> x = [1.5, 3.1, 6]
    >>> findidx(x_array, x)
    array([0.33333333, 2.04      ,        nan])

    Modification History
    --------------------
    2022-06-30 - Written by M. Cushing, University of Toledo.
    """

    # Convert to numpy arrays and get basic things
    try:
        x_array = np.asarray(x_array, dtype='float')
    except ValueError:
        raise

    ndat = len(x_array)

    # Deal with array_like versus int/float
    if isinstance(x, int) or isinstance(x, float) is True:
        single = True
        x = np.asarray([x], dtype='float')
    else:
        single = False
        x = np.asarray(x, dtype='float')

    # Initialize binary search area and compute number of divisions needed
    ileft = np.zeros(len(x), dtype=int)
    iright = np.zeros(len(x), dtype=int)

    n_divisions = int(np.log10(ndat)/np.log10(2) + 1)

    # Test array for monotonicity. If not, raise ValueError
    i = x_array - np.roll(x_array, 1)
    i = i[1:]
    # check if array is increasing or decreasing
    a = i >= 0
    if np.sum(a) == ndat-1:
        iright += ndat-1
    else:
        a = i <= 0
        if np.sum(a) == ndat-1:
            ileft += ndat-1
        else:
            raise ValueError('find_index:  array is not monotonic.')

    # Perform binary search by dividing search interval in half NDIVISIONS times
    for i in range(n_divisions):
        idiv = (ileft + iright)//2   # Split interval in half
        x_val = x_array[idiv]            # Find function values at center
        greater = x > x_val           # Determine which side of X is on
        less = x <= x_val
        np.multiply(ileft, less, out=ileft)
        np.add(ileft, idiv*greater, out=ileft)
        np.multiply(iright, greater, out=iright)
        np.add(iright, idiv*less, out=iright)

    # linearly interpolate
    x_left = x_array[ileft]
    x_right = x_array[iright]

    mask = x_left == x_right
    ieff = (x_right-x) * ileft
    np.add(ieff, (x-x_left) * iright, out=ieff)
    np.add(ieff, mask * ileft, out=ieff)
    np.divide(ieff, (x_right - x_left + mask), out=ieff)

    # Clip points beyond interpolation to NaN
    z = ieff < 0
    ieff[z] = np.nan

    z = ieff > ndat-1
    ieff[z] = np.nan

    # Return numpy.ndarray of float or float
    return ieff[0] if single is True else ieff
    

def loopprogress(idx,bot,top,message=None):

    '''
    To print a loop update message on the command line

    Parameters
    ----------------
    idx : int
          the current loop index number, e.g. i

    bot : int
          the bottom of the loop, e.g. range(bot,top)

    top : int
          the bottom of the loop, e.g. range(bot,top)

    message : str
              A string to print before the update begin, optional.

    Returns
    --------
    None

    Procedure
    ---------
    Just prints a updating message on the command line

    Example
    --------
    >>> looprogress(i,0,100)
     90% |*************************************                       |


    Modification History
    --------------------
    2022-05-25 - Written by M. Cushing, University of Toledo.
                Based on the Spextool mc_loopprogress.pro IDL program.
    '''
    
    # Print the message if necessary
    
    if idx == 0 and message is not None:

        print(message)

    # Run counter        
            
    frac = (idx+1-bot)/(top-bot)
    stars = "{:<70}".format('*'*round(frac*70))
    print(str(round(frac*100)).rjust(3),'% ', '|', stars,'|',\
          sep='',end='\r')

    # Now get your prompt back

    if idx == top-1:

        print()


def mkimgidxs(nrows,ncols,dtype=int):

    '''

    To generate indice grids for an image

    Input Parameters
    ----------------
    nrows : int
        The number of rows in the image

    ncols : int
        The number of columns in the image
        
    dtype : dtype, default, `int`, optional
        The type of the output grids


    Returns
    --------
    ximg, yimg : numpy.ndnarray, numpy.ndnarray of `dtype`
         A list of integers giving the individual file numbers


    Notes
    -----
    https://stackoverflow.com/questions/1550130/\
    cloning-row-or-column-vectors

    Examples
    --------
    > ximg, yimg = mkimgidxs(3,5)
    > print(ximg)
      [[0 1 2 3 4]
       [0 1 2 3 4]
       [0 1 2 3 4]]

    > print(yimg)
      [[0 0 0 0 0]
      [1 1 1 1 1]
      [2 2 2 2 2]]

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Just the old reform/rebin trick from IDL.

    '''
        
    ximg = np.tile(np.arange(ncols,dtype=dtype),(nrows,1))
    yimg = np.tile(np.reshape(np.arange(nrows,dtype=dtype),\
                   (nrows,1)),(1,ncols))

    return (ximg, yimg)


def moments(data,goodbad=False,robust=None,silent=True):

    '''
    (Robustly) computes various statistics

    Input Parameters
    ----------------
    data : numpy.ndarray
        
    goodbad : numpy.ndarray, optional
        An array with the same shape as `data` that identifies good and 
        bad data points.  0=bad, 1=good, 2=NaN

    robust : float, optional
        If given, outliers are identified before computing the stats.  
        See notes for details.
        
    silent : {True, False}, optional
        If False, the result will be written to the command line.

    Returns
    --------
    dict
        A dict where with the following entries:

        ndat : int
            The numger of data points used in the calculation.

        mean : float
            The mean of the (non-NaN, good) data points.

        variance : float
            Estimate of the variance of the (non-NaN, good) data points.  
            That is, the denominator is 1/(ndat-1).

        stddev : float
            Estimate of the standard deviation of the (non-NaN, good) 
            data points.  That is, the denominator is 1/(ndat-1).

        stderr : float
            The standard error of the (non-NaN, good) data points.  
            `stddev`/sqrt(ndat)

        skewness : float
            The skewness of the (non-NaN, good) data points.

        kurtosis : float
            The kurtosis of the (non-NaN, good) data points.
        
        goodbad : numpy.ndarray of int
            An array with the same shape as `data` that identifies good and 
            bad data points.  0=bad, 1=good, 2=NaN        

    Notes
    -----
    If goodbad is passed, only points with values of 1 are used.  If 
    robust is passed, the median and median absolute deviation and 
    points are idetentified as an outlier if:

    |x_i - MED|/(1.4826*MAD) > robust

    where MAD = median(|x_i - MED|) and 1.4826*MAD is a robust estimate 
    of the standard deviation for a gaussian distribution. Outliers are 
    labelled `bad` in the goodbad array.  Finally, the statistics are 
    computed using scipy.stats.describe.  

    NOTE:  The variance and standard deviations are *estimates* of the 
    variance and standard deviation of the parent population and so 
    have 1/(ndat-1) in the denominator. 

    Examples
    --------
    > import numpy as np
    > data = np.array([[np.nan,1.2,20],[2.,1.2,2.5]])
    > m = moments(data,robust=4,silent=False)
      Moments results:

      Total number of input points =  6
              Number of good points =  4
                     Number of NaNs =  1

                    Mean =  1.725
                Variance =  0.4091666666666667
       Standar Deviation =  0.6396613687465162
          Standard Error =  0.3198306843732581
                Skewness =  0.28952649685958215
                Kurtosis =  -1.6237779003737334
      [[2 1 0]
       [1 1 1]]


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_moments.pro.

    '''
    
    # Set up goodbad array if need be
    
    if goodbad is False: goodbad = np.full_like(data,1,dtype=int)

    # Check for NaNs and update goodbad array

    nanbool = np.isnan(data)
    goodbad[nanbool] = 2

    nnan = np.sum((goodbad == 2))
    
    # Now find the sample you are working with
    
    zsample = np.where(goodbad == 1)

    sample = data[zsample]

    # Robust calculation?

    if robust is not None:

        # Compute the median and 1.4826*median absolute deviation (gmad).

        med = np.median(sample)
        gmad = 1.4826*np.median(np.abs(sample-med))

        #  Generate the mask
        
        mask = ((sample-med)/gmad) <= robust

    elif robust is None:

        mask = np.full_like(sample,True, dtype=bool)
        
    # update the goodbad array

    goodbad[zsample] = np.array(mask)
        
    # Do the stats
        
    stats = describe(sample[mask])

    ngood = stats.nobs
    mean = stats.mean
    var = stats.variance
    stddev = np.sqrt(stats.variance)
    stderr = np.sqrt(var/ngood)
    skewness = stats.skewness
    kurtosis = stats.kurtosis

    # Report the results if asked

    if silent is False:

        print('Moments results:')        
        print()
        print(' Total number of input points = ',np.size(data))
        print('        Number of good points = ',ngood)
        print('               Number of NaNs = ',nnan)
        print()
        print('              Mean = ',mean)
        print('          Variance = ',var)
        print(' Standar Deviation = ',stddev)
        print('    Standard Error = ',stderr)
        print('          Skewness = ', skewness)
        print('          Kurtosis = ', kurtosis)
        
    return {'ndat':ngood,'mean':mean,'var':var,'stddev':stddev,\
            'stderr':stderr,'skewness':skewness,'kurtosis':kurtosis,\
            'goodbad':goodbad}
    

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