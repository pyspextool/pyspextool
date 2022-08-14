import numpy as np
from scipy.stats import describe


def moments(data, goodbad=False, robust=None, silent=True):
    """
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
    points are identified as an outlier if:

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

        # Compute the median and 1.4826*median absolute deviation (gmad).

        med = np.median(sample)
        gmad = 1.4826 * np.median(np.abs(sample - med))

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


def median_data_stack(data, mask=None, stderr=True):
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

    # Compute the median and median absolute deviation

    med = np.nanmedian(data, axis=0)

    mad = np.nanmedian(np.abs(data - np.tile(med, tileshape)), axis=0)

    if stderr is not None:

        mad *= 1.4826  # assume gaussian distribution
        unc = mad / np.sqrt(np.sum(mask, axis=0))

    else:

        unc = mad / 1.4826  # assume gaussian distribution

    return [med, unc]


def scale_data_stack(stack, *var, mask=None, idx=None):

    """
    Scales a stack of spectra or images to a common intensity level.

    Input Parameters
    ----------------
    stack : array_like
        (nspec,nwave) - a stack of spectra.

        (nimgs,nrows,cols) - a stack of images.


    *var : array_like, optional
        The variances associated with the spectra or images.


    mask : array_like, optional
        A mask array with the same shape as `stack`.
        0 = bad, 1=good


    idx : int, optional
        A mask array with the same shape as `data`.
        0 = bad, 1=good


    Returns
    --------
    sstack : numpy.ndarray
        The scaled stack.


    svar : numpy.ndarray or None
        The scaled variance, if *var is passed.


    scales : numpy.ndarray
        The scale factors.


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

    Modification History
    --------------------
    2022-06-18 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_getspecscales.pro

    """

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

    if idx is None:

        # Scale to the median 

        scales = np.median(medvals) / medvals

    else:

        # Scale to a spectrum or image

        scales = medvals[idx] / medvals

    # Build the final scale array and do the math

    sclarr = np.tile(np.reshape(scales, reshape_info), tile_info)

    sstack = stack * sclarr

    # Now return the scaled stack and potentially the scaled variance

    if len(var) == 1:

        return sstack, var * sclarr**2, scales

    else:
        return sstack, None, scales


def combine_flag_stack(stack, nbits=8):

    """
    To combine bit-set flag arrays.

    Input Parameters
    ----------------
    stack : numpy.ndarray
        The stack of bit-set flag arrays to combine.  The stack
        can either be a stack of spectra [nspec,ndat] or a stack
        of images [nimg,nx,ny].

    nbits : int, optional
        The number of bits that can potentially be set.  This
        routine assumes the bits are set sequentially, starting
        with the zeroth bit.  So if nbits is 2, then it will
        check the 0th and 1st bit.  The default is to check all
        eight bits

    Output Parameters
    ------------------
    numpy.ndarray
        A bit-set flag array that reflects the bit-set flags from all
        the spectra or images.

    Procedure
    ---------
    Just some basic math.

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
    > combflagstack(stack)

    [[1 2 0]
     [3 0 4]
     [0 0 0]]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.
    Based on the mc_combflagstack.pro IDL program.

    """

    # Determine whether this is a spectral stack or image stack

    ndim = stack.ndim
    shape = stack.shape

    # Set up the output array
    # Replace with match case statement when you upgrade to 3.10.

    if ndim == 2:
        comb = np.zeros(shape[1], dtype=np.int8)

    if ndim == 3:
        comb = np.zeros(shape[1:2], dtype=np.int8)

    # Now just loop over each bit requested.

    for i in range(0, nbits):
        #  Identify the pixels with the particular bit set

        set = bit_set(stack, i)

        #  Collapse everything down one dimension

        sum = np.sum(set, axis=0)

        #  Identify which pixels are set

        mask = sum > 0

        #  Set them to the proper bit value and add to the comb

        comb = comb + mask * 2 ** i

    return comb


def bit_set(array, bits):

    """
    To determine if the given bits are set in an array.

    Input Parameters
    ----------------
    array : numpy.ndarray
        numpy array to search

    bits : int, list, numpy.ndarray
        bit values to search in `array`

    Returns
    --------
    numpy.ndarray
        Returns a byte array of the same size as `array`.  An element
        is set if any of the bits requested are set in the same element
        of `array`.

    Procedure
    ---------
    Uses the Gumley IDL ishft technique.  Note that the "first" bit 
    is denoted as zero, while the "second" bit is denoted as 1.


    Example
    --------

    > import numpy as np
    > bitset(np.array([3,4,1]),0)

    [1, 0, 1]

    > import numpy as np
    > bitset(np.array([3,4,1]),[0,3])

    [1, 0, 1]

    > import numpy as np
    > bitset(np.array([3,4,1]),[2,3])

    [0, 1, 0]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.
                 Based on the Spextool IDL mc_bitset.pro program.
    """

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
        mask = mask | tmp

    return mask


def round(x):

    """
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

    """
    
    mask = (x >= 0)
    out = np.empty_like(x)
    out[mask] = np.floor(x[mask] + 0.5)
    out[~mask] = np.ceil(x[~mask] - 0.5)

    return out
