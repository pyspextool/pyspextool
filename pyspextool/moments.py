import numpy as np
from scipy.stats import describe

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
    
