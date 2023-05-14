"""Functions for robust procedures."""

import numpy as np
import scipy
from scipy import signal
from scipy import interpolate

from pyspextool.fit.polyfit import goodbad_init


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

    """

    # Convert to numpy arrays and do basic things

    goodbad = goodbad_init(y, goodbad=goodbad)
    ndat = np.size(y)
        
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

            # Check to see if the new stddev isn't much of a change
            # from the old one

            if (oldstddev - stddev) / oldstddev < eps:
                break

            # Now generate a new full residual array

            f = scipy.interpolate.interp1d(xx[good], smthyy,                                                               fill_value='extrapolate')
            newyy = f(xx)
            residual = yy - newyy

            good = np.abs((residual - mean) / stddev) <= thresh
            good_cnt = np.sum(good)


    newy = np.full(ndat, np.nan)
    newy[initgood] = newyy

            
    tmp = np.full(initgood_cnt, 0, dtype=int)
    tmp[good] = 1
    goodbad[initgood] = tmp

    return {'fit': newy, 'goodbad': goodbad}



