import numpy as np
import numpy.typing as npt
import warnings
from scipy.stats import describe

from pyspextool.io.check import check_parameter

def bit_set(array:npt.ArrayLike,
            bits:int | list | npt.ArrayLike):

    """
    To determine if the given bits are set in an array.

    Parameters
    ----------
    array : ndarray
        numpy array to search

    bits : int, list, ndarray
        bit values to search in `array`.

        bit=0 -> the first bit
        bit=1 -> the second bit

    Returns
    --------
    ndarray
        Returns a byte array of the same size as `array`.  An element
        is set if any of the bits requested are set in the same element
        of `array`.

    Example
    -------

    > import numpy as np
    > bit_set(np.array([3,4,1]),0)

    [1, 0, 1]

    > import numpy as np
    > bit_set(np.array([3,4,1]),[0,3])

    [1, 0, 1]

    > import numpy as np
    > bit_set(np.array([3,4,1]),[2,3])

    [0, 1, 0]

    """

    #
    # Check parameters
    #

    check_parameter('bit_set', 'array', array, 'ndarray')

    check_parameter('bit_set', 'bits', bits, ['int', 'list', 'ndarray'])    

    #  Define empty mask

    mask = np.zeros_like(array, dtype=np.uint8)

    #  test to see if bits is iterable

    try:

        iter(bits)

    except TypeError:

        bits = [bits]

    #  Loop over every bit requested and identify those pixels for
    #  which that bit is set.

    for val in bits:

        tmp = (array >> val) & 1
        mask = (mask | tmp)

    return mask



def combine_flag_stack(stack:npt.ArrayLike,
                       nbits:int=8):

    """
    To combine bit-set flag arrays.

    Parameters
    ----------
    stack : ndarray
        The stack of bit-set flag arrays to combine.  The stack
        can either be a stack of spectra [nspec,ndat] or a stack
        of images [nimg,nx,ny].

    nbits : int, default=8
        The number of bits that can potentially be set.  This
        routine assumes the bits are set sequentially, starting
        with the zeroth bit.  So if nbits is 2, then it will
        check the 0th and 1st bit.  The default is to check all
        eight bits

    Returns
    -------
    ndarray
        A bit-set flag array that reflects the bit-set flags from all
        the spectra or images.

    Example
    -------
    Consider a two spectra masks

    > spec1 = np.array([0,4,4,2])
    > spec2 = np.array([0,0,3,1])
    > stack = np.stack((spec1,spec2))
    > combflagstack(stack)

    [0 4 7 3]

    Consider two image masks

    > img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
    > img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
    > stack = np.stack((img1,img2))
    > combine_flag_stack(stack)

    [[1 2 0]
     [3 0 4]
     [0 0 0]]

    """

    #
    # Check parameters
    #

    check_parameter('combine_flag_stack', 'stack', stack, 'ndarray')

    check_parameter('combine_flag_stack', 'nbits', nbits, 'int')    

    # Determine whether this is a spectral stack or image stack

    ndim = stack.ndim
    shape = stack.shape

    # Set up the output array
    # Replace with match case statement when you upgrade to 3.10.

    if ndim == 2:
        comb = np.zeros(shape[1], dtype=np.uint8)

    if ndim == 3:
        comb = np.zeros(shape[1:3], dtype=np.uint8)

        
    # Now just loop over each bit requested.

    for i in range(0, nbits):
        #  Identify the pixels with the particular bit set

#        print(i, 'in', stack.dtype)
        set = bit_set(stack, i)
#        print(i, 'out', set.dtype)
        #  Collapse everything down one dimension

        sum = np.sum(set, axis=0)

        #  Identify which pixels are set

        mask = sum > 0

        #  Set them to the proper bit value and add to the comb

        comb +=  np.multiply(mask,(2 ** i), dtype='uint8')

    return comb



def find_outliers(data, thresh, leave_nans=True, silent=False):

    '''
    To identify outliers in a distribution of data.

    Parameters
    ----------
    data : ndarray
        A array of data with potential outliers and possibly NaNs.  

    thresh : int or float
        The value over which data is identified as an outlier (see Notes).

    leave_nans : {True, False} optional
        Set to True to not identify NaNs as an outlier. Set to False to 
        identify NaNs as outliers.

    silent : {False, True} optional
        Set to True to report to the command line that no identification 
        of outliers can be performed because the median absolute deviation is 
        zero.

    Returns
    -------
        ndarray of bool

    Notes
    -----

    Computes the median and median absolute deviation
    mad=median(abs(median-data))
    of the data excluding NaNs and identifies pixels as outliers if 
    abs( (data-med)/mad ) > thresh).

    Examples
    --------
    > x = np.array([1,2,np.nan,3,4.1,6.7,8.2,10.2,1.4])
    > mask = find_outliers(x, 5)
      [ True  True  True  True  True  True  True  True  True]

    > x = np.array([1,2,np.nan,3,4.1,6.7,8.2,10.2,1.4])
    > mask = find_outliers(x, 5, leave_nans=False)
      [ True  True False  True  True  True  True  True  True]

    > x = np.array([1,2,np.nan,3,4.1,6.7,8.2,10.2,1.4,100])
    > mask = find_outliers(x, 5)
      [ True  True  True  True  True  True  True  True  True False]

    > x = np.array([1,2,np.nan,3,4.1,6.7,8.2,10.2,1.4,100])
    > mask = find_outliers(x, 5, leave_nans=False)
      No outliers identified because MAD=0.0.
      [ True  True  True]

    '''

    #
    # Check parameters
    #

    check_parameter('find_outliers', 'data', data, 'ndarray')

    check_parameter('find_outliers', 'thresh', thresh, ['int', 'float'])

    check_parameter('find_outliers', 'silent', silent, 'bool')        

    #
    # Deal with NaNs
    #
    
    # Create a mask to avoid NaNs

    mask = np.full_like(data, 1, dtype=int)
    mask[np.isnan(data)] = 2

    # Check to see whether there are at least two data points that aren't NaN

    z = mask != 2
    if np.sum(z) < 2:

        if silent is False:
            print('Fewer than 2 not-NaN points found.')

        
        if leave_nans is True:

            return mask != 0

        else:
        
            return mask == 1
        
    #
    # Do the search
    #
    
    # Compute the median and median absolute deviation
    
    med = np.median(data[z])
    mad = np.median( np.abs(data[z]-med))
    
    # If the MAD is zero, just return the non NaN data
    
    if mad == 0.0:
        z = mask == 1
        if silent is False:
            print('No outliers identified because MAD=0.0.')

        return z

    # Now do the real search

    z = np.where(np.abs( (data-med)/mad ) > thresh)
    mask[z] = 0

    if leave_nans is True:

        return mask != 0

    else:
        
        return mask == 1


    
def mean_data_stack(data,
                    weights=None,
                    goodbad=None,
                    robust=None,
                    stderr=True):

    '''
    (Robustly) Compute the mean of a spectral or image stack with optional mask
    
    Parameters
    ----------
    data : ndarray
        either a stack of spectra (nspec, npoints) or a stack of images
        (nimgs, nrows, ncols).

    weights : ndarray , optional
        either a stack for spectra (nspec, npoints) or a stack of 
        for images (nimgs, nrows, ncols).  For a standard weighted mean, 
        pass 1/var.

    goodbad : numpy.ndarray, optional
        an int mask array with the same shape as `data`.  
        0 = bad, 1=good

    robust : float, optional
        The threshold for identifying outliers.  See find_outliers.

    stderr : {True, False}, optional
        If True, and `weights` is None, then the output `mvar` will be the 
        standard error**2.

    Returns
    -------
    tuple

        tuple(0): nndarray
            The mean of the data stack

        tuple(1) : ndarray
            The sample variance of the stack.  That is, the result is divided 
              by (n-1).  
            If `stderr` is True, and `weights` is None, then the 
            `mvar` is the standard error which is the sample standard 
            deviation divided by sqrt(n).  

        tuple(2) : ndarray
        

    Examples
    --------
    later
       

    '''

    #
    # Check the parameters
    #

    check_parameter('mean_data_stack', 'data', data, 'ndarray', [2, 3])

    check_parameter('mean_data_stack', 'weights', weights,
                    ['NoneType','ndarray'], [2, 3])

    check_parameter('mean_data_stack', 'goodbad', goodbad,
                    ['NoneType','ndarray'], [2, 3])

    check_parameter('mean_data_stack', 'robust', robust,
                    ['NoneType','int', 'float'])    

    check_parameter('mean_data_stack', 'stderr', stderr, 'bool')        

    # Get array dimensions and get set up

    ndimen = np.ndim(data)
    shape = np.shape(data)

    #
    # Get the mask set up
    #
    
    if goodbad is None:

        # Create a blank mask where all pixels are good
        
        goodbad = np.full_like(data, 1,dtype=int)
        
    else:

        # Convert any pixel identified as a NaN to a bad pixel 
        
        z = goodbad == 2
        goodbad[z] = 0
    
    # Check for NaNs in the data and weights just to make sure.
        
    z = np.isnan(data)
    goodbad[z] = 0

    if weights is not None:

        z = np.isnan(weights)
        goodbad[z] = 0

    #
    # Identify outliers if requested
    #

    if robust is not None:

        if ndimen == 2:  #  This is a spectral stack

            for i in range(shape[1]):
                        
                submask = find_outliers(data[:,i], robust, leave_nans=False,
                                        silent=True)
                goodbad[~submask,i] = 0                

        else:  #  This is an image stack

            for i in range(shape[1]):

                for j in range(shape[2]):
                        
                    submask = find_outliers(data[:,i,j], robust,
                                            leave_nans=False, silent=True)
                    goodbad[~submask,i,j] = 0            

    #
    # Now set bad pixels and/or NaN pixels to zero.
    #
    
    z = goodbad == 0
    data[z] = 0.0

    if weights is not None:

        weights[z] = 0.0
                        
    #    
    # Now do the combining
    #

    if weights is not None:

        sum_weights = np.nansum(weights, dtype=float, axis=0)
        z_allnan = sum_weights == 0
        sum_weights[z_allnan] = np.nan
 
        mvar = 1/sum_weights
        mean = np.nansum(np.multiply(data, weights),axis=0)*mvar
        
    else:

        ndat = np.nansum(goodbad, axis=0)

        # Find the pixels that were all NaNs
        
        z_allnan = ndat == 0
        
        # set those column in ndat to two to avoid a division error

        ndat[z_allnan] = 2

        # Do the calculation

        mean = np.sum(data, axis=0)/ndat
        mvar = (np.sum(data**2,axis=0)-ndat*mean**2)/(ndat-1)

        if stderr is True:
            np.divide(mvar,ndat,out=mvar)

        # Set the all-NaN columns to NaN

        mean[z_allnan] = np.nan
        mvar[z_allnan] = np.nan

    return (mean, mvar, goodbad)
    


def median_data_stack(data,
                      mask=None,
                      stderr=True):

    """
    Median a spectral or image stack with optional mask


    Parameters
    ----------
    data : numpy.ndarray
        either a stack of spectra (nspec, npoints) or a stack of images
        (nimgs, nrows, ncols).

    mask : numpy.ndarray, optional
        a mask array with the same shape as `data`.  
        0 = bad, 1=good

    stderr : {True, False}, optional
        Set to return 1.4826*MAD/sqrt(n) instead of just 1.4826*MAD 
        (see Procedure)

    Returns
    --------
    tuple
        tuple[0] : numpy.ndarray 
            the median of the spectral or image stack

        tuple[1] : numpy.ndarray
            the uncertainty of the spectral or image stack (see Procedure)

    Procedure
    ---------
    Spectral stack:

        in this case, the data have the shape (nspec, npoints).  The
        median of the stack is computed producing an array of size 
        (npoints,). At each spectral point, the median absolute deviation
        (MAD=median(|data-med}) is computed.  

    Image stack:

        in this case, the data have the shape (nspec, nrows, ncols).  The
        median of the stack is computed producing an array of size 
        (nrows, ncols). At each image point, the median absolute deviation
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

    """

    # Get array dimensions

    ndimen = np.ndim(data)
    shape = np.shape(data)

    # If no mask passed, create one.

    if mask is None:
        mask = np.ones(shape, dtype=int)

    # Now search and replace any masked pixels with NaNs

    data = np.where(mask != 0, data, np.nan)

    # Spectral or image stack?

    if ndimen == 2:

        tileshape = (shape[0], 1)  # spectral stack

    elif ndimen == 3:

        tileshape = (shape[0], 1, 1)  # image stack

    else:

        print('Unknown data shape.')
        return

    # Compute the median and median absolute deviation.  We know there can be
    # all-NaN columns, so shut down the warning of all-NaN slices pre-emptively.
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
    
        med = np.nanmedian(data, axis=0)
        mad = np.nanmedian(np.abs(data - np.tile(med, tileshape)), axis=0)

    if stderr is not None:

        mad *= 1.4826  # assume gaussian distribution
        unc = mad / np.sqrt(np.sum(mask, axis=0))

    else:

        unc = mad / 1.4826  # assume gaussian distribution

    return (med, unc)


def moments(data:npt.ArrayLike,
            goodbad:npt.ArrayLike=False,
            robust:int | float=None,
            silent:bool=True):

    """
    (Robustly) computes basic statistics

    Parameters
    ----------
    data : numpy
        An array of numbers.

    goodbad : numpy, optional
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
    points are identified as an outlier if:
    |x_i - MED|/(MAD) > robust where MAD = median(|x_i - MED|) (1.4826*MAD 
    is a robust estimate of the standard deviation for a gaussian 
    distribution.). Outliers are labelled `bad` in the goodbad array.  
    Finally, the statistics are computed using scipy.stats.describe.
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

    """

    # Set up goodbad array if need be

    if goodbad is False:
        goodbad = np.full_like(data, 1, dtype=int)

    # Check for NaNs and update goodbad array

    nanbool = np.isnan(data)
    goodbad[nanbool] = 2

    nnan = np.sum((goodbad == 2))

    # Now find the sample you are working with

    zsample = np.where(goodbad == 1)

    sample = data[zsample]

    # Robust calculation?

    if robust is not None:

        # Compute the median and median absolute deviation (gmad).

        med = np.median(sample)
        gmad = np.median(np.abs(sample - med))

        #  Generate the mask

        mask = ((sample - med) / gmad) <= robust

    elif robust is None:

        mask = np.full_like(sample, True, dtype=bool)

    # update the goodbad array

    goodbad[zsample] = np.array(mask)

    # Do the stats

    stats = describe(sample[mask])
    ngood = stats.nobs
    mean = stats.mean
    var = stats.variance
    stddev = np.sqrt(stats.variance)
    stderr = np.sqrt(var / ngood)
    skewness = stats.skewness
    kurtosis = stats.kurtosis

    # Report the results if asked

    if silent is False:
        print('Moments results:')
        print()
        print(' Total number of input points = ', np.size(data))
        print('        Number of good points = ', ngood)
        print('               Number of NaNs = ', nnan)
        print()
        print('              Mean = ', mean)
        print('          Variance = ', var)
        print(' Standar Deviation = ', stddev)
        print('    Standard Error = ', stderr)
        print('          Skewness = ', skewness)
        print('          Kurtosis = ', kurtosis)

    return {'ndat': ngood, 'mean': mean, 'var': var, 'stddev': stddev,
            'stderr': stderr, 'skewness': skewness, 'kurtosis': kurtosis,
            'goodbad': goodbad}

def round(x):

    """
    To round numbers in a numpy array


    Parameters
    ----------
    x : ndarray


    Returns
    --------
    ndarray like `x`
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

    """
    
    mask = (x >= 0)
    out = np.empty_like(x)
    out[mask] = np.floor(x[mask] + 0.5)
    out[~mask] = np.ceil(x[~mask] - 0.5)

    return out



def scale_data_stack(stack:npt.ArrayLike,
                     var:npt.ArrayLike,
                     mask:npt.ArrayLike=None,
                     index=None):

    """
    Scales a stack of spectra or images to a common intensity level.

    Parameters
    ----------
    stack : ndarray
        An (nspec, nwave) array of spectra or an (nimgs, nrows, cols)
        array of images.

    var : ndarray or None
        If not None, the variances associated with the spectra or images.

    mask : ndarray, optional
        A mask array with the same shape as `stack`.  0 = bad, 1=good.
        
    index : int, optional
        The index into the first dimension of `stack` giving the
        spectrum/image that the rest of the spectra/images should be scaled to.

    Returns
    --------
    sstack : ndarray
        An (nspec, nwave) or (nimgs, nrows, ncols) array of the scaled stack.

    svar : ndarray or None
        An (nspec, nwave) or (nimgs, nrows, ncols) array of the scaled 
        variance, if `var` is not None.

    scales : ndarray
        An (nspec, ) array of scale factors.

    Notes
    -----
    Computes the median of each spectrum or image and then determines
    scale factors to scale each spectrum or image to either the median
    value of the median values or the median value associated
    with a particular spectrum or image.

    Examples
    --------
    Spectral stack, scale to median:

        > spec1 = np.array([1,1,1,1])
        > spec2 = np.array([2,2,2,2])
        > spec3 = np.array([4,4,4,4])
        > specstack = np.stack((spec1,spec2,spec3))
        > scaledstack, var, scales = scalestack(specstack,None)
        > print(scaledstack)
        > print(scales)
          [[2. 2. 2. 2.]
           [2. 2. 2. 2.]
           [2. 2. 2. 2.]]
          [2.  1.  0.5]

    Spectral stack, scale to first spectrum:

        > spec1 = np.array([1,1,1,1])
        > spec2 = np.array([2,2,2,2])
        > spec3 = np.array([4,4,4,4])
        > specstack = np.stack((spec1,spec2,spec3))
        > scaledstack, var, scales = scalestack(specstack,None,idx=0)
        > print(scaledstack)
        > print(scales)
          [[1. 1. 1. 1.]
           [1. 1. 1. 1.]
           [1. 1. 1. 1.]]
          [1.   0.5  0.25]

    Image stack, scale to median:

        > img1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
        > img2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
        > img3 = np.array([[4,4,4],[4,4,4],[4,4,4]])
        > imgstack = np.stack((img1,img2,img3))
        > scaledstack, var, scales = scalestack(imgstack, None)
        > print(scaledstack)
        > print(scales)
          [[[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]

           [[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]

           [[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]]
          [2.  1.  0.5]

    Image stack, scale to first image:

        > img1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
        > img2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
        > img3 = np.array([[4,4,4],[4,4,4],[4,4,4]])
        > imgstack = np.stack((img1,img2,img3))
        > scaledstack, var, scales = scalestack(imgstack, None)
        > print(scaledstack)
        > print(scales)

          [[[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]

           [[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]

           [[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]]
          [1.   0.5  0.25]


    """

    #
    #  Check parameters
    #

    check_parameter('scale_data_stack', 'stack', stack, 'ndarray')

    check_parameter('scale_data_stack', 'var', var, ['ndarray', 'NoneType'])

    check_parameter('scale_data_stack', 'var', mask, ['ndarray', 'NoneType'])

    check_parameter('scale_data_stack', 'index', index, ['int', 'NoneType'])            
    
    # Get array dimensions

    ndimen = np.ndim(stack)
    shape = np.shape(stack)

    # If no mask passed, create one.

    if mask is None:
        mask = np.ones(shape, dtype=int)

    # Now search and replace any masked pixels with NaNs

    stack = np.where(mask != 0, stack, np.nan)

    # Set up median, reshape, and axis variables

    if ndimen == 2:

        # Spectral stack

        axis_info = 1
        reshape_info = (shape[0], 1)
        tile_info = (1, shape[1])

    elif ndimen == 3:

        # Image stack

        axis_info = (1, 2)
        reshape_info = (shape[0], 1, 1)
        tile_info = (1, shape[1], shape[2])

    else:

        print('scalestack:  Unknown data shape.')
        return -1

    # Do the calculation

    medvals = np.nanmedian(stack, axis=axis_info)

    # Check whether you are scaling to the median or spectrum/image 

    if index is None:

        # Scale to the median 

        scales = np.median(medvals) / medvals

    else:

        # Scale to a spectrum or image

        scales = medvals[index] / medvals

    # Build the final scale array and do the math

    sclarr = np.tile(np.reshape(scales, reshape_info), tile_info)

    sstack = stack * sclarr

    # Now return the scaled stack and potentially the scaled variance

    if var is not None:

        return sstack, var * sclarr**2, scales

    else:
        return sstack, None, scales



