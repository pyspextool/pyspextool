"""Functions for robust procedures."""

import numpy as np
import scipy
from scipy import signal
from scipy import interpolate


def robust_savgol(x, y, window_length, polyorder=2, thresh=4, eps=0.1,
                  goodbad=None):

    """
    To robustly perform a savizky-golay smoothing.


    Parameters
    ----------
    x : array_like
        (ndat,) array of independent values.

    y : array_like
        (ndat,) array of dependent values.

    window_length : int
        The length of the filter window (i.e., the number of coefficients).
        `width` must be less than or equal to the size of `x`
        (from scipy.signal.savlgol_filter).

    polyorder : int, default 2
        The order of the polynomial used to fit the samples. `polyorder`
        must be less than window_length.
        (from scipy.signal.savlgol_filter).

    thresh : float, default 4
        The standard deviation threshold used identify outliers.

    eps : float, default 0.1
        The limit of the fractional change in the standar deviation of the fit.

    goodbad : array-like of int, optional
        An (ndat,) array identifying good and bad pixels. 0=bad, 1=good, 2=NaN.

    Returns
    -------
    dict
        fit : numpy.ndarray of float
            (ndat,) array of smoothed values.

        goodbad : numpy.ndarray of int
            (ndat,) array identifying good and bad pixels. 0=bad, 1=good, 2=NaN.

    Notes
    -----
    While `x` is not required for savitzky-golay smoothing, it is required
    to do the smoothing robustly because we have to interpolate over pixels
    identified as bad.

    The function uses scipy.signal.savgol_filter.  All savgol_filters take
    on their default values.


    Examples
    --------

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program XS

    """

    # Convert to numpy arrays and do basic things

    goodbad = robustsetup(y, goodbad=goodbad)

    # Grab the good points

    initgood = goodbad == 1

    initgood_cnt = np.sum(initgood)
    xx = x[initgood]
    yy = y[initgood]

    # Do an initial fit

    newyy = scipy.signal.savgol_filter(yy, window_length, polyorder)

    # Compute the residuals and the mean and sigma of the residuals

    residual = yy - newyy

    mean = np.mean(residual)
    stddev = np.std(residual)

    # Now check for residual points that are outside the threshhold

    good = np.abs((residual - mean) / stddev) <= thresh
    good_cnt = np.sum(good)

    # If it doesn't find anything, move on.  If not, get to the looping

    if good_cnt != initgood_cnt:

        for i in range(0, 11):

            oldstddev = stddev

            smthyy = scipy.signal.savgol_filter(yy[good], window_length,
                                                polyorder)
            residual = yy[good] - smthyy

            mean = np.mean(residual)
            stddev = np.std(residual)

            # Check to see if the new stddev isn't much of a change from the old one

            if (oldstddev - stddev) / oldstddev < eps:
                break

            # Now generate a new full residual array

            f = scipy.interpolate.interp1d(xx[good], smthyy)
            newyy = f(xx)
            residual = yy - newyy

            good = np.abs((residual - mean) / stddev) <= thresh
            good_cnt = np.sum(good)

    tmp = np.full(initgood_cnt, 0, dtype=int)
    tmp[good] = 1
    goodbad[initgood] = tmp

    return {'fit': newyy, 'goodbad': goodbad}


def robustsetup(*args, goodbad=None):

    """
    To set up the goodbad array for a robust fit.

    Parameters
    ----------
    *args : array_like or None
        (ndat,) arrays to be used in a robust fit.  Typically, the
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

    """

    if goodbad is None:

        # Create a good bad array

        goodbad = np.full(len(args[0]), 1, dtype=int)

    else:

        # Convert to a numpy array if need be

        goodbad = np.array(goodbad, dtype=int)

        # Find NaNs if need be

    for arg in args:

        if arg is not None:

            znan = np.isnan(arg)
            if np.sum(znan) > 0:
                goodbad[znan] = 2

    return goodbad
