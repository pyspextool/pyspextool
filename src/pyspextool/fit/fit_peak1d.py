import numpy as np
import scipy

import warnings

def fit_peak1d(x, y, type='gaussian', nparms=4, p0=None, positive=False,
              negative=False, nan=True, ignore_optimizewarning=False):

    """
    To fit 1D data with either a Gaussian/Lorentzian+polynomial function.


    At the moment, very basic.  Will add complexity as necessary.
    
    Parameters
    ----------
    x : ndarray
        (nx,) array of independent variables.

    y : ndarray
        (ny,) array of dependent variables.

    type: {'gaussian','lorentzian'}, optional
        The type of fit.

    nparms : {4,3,5,6,7,8}, optional 
        The number of parameters.  See Notes.

    p0 : array-like, optional
        (nparms,) array of parameters estimates.  If given, the function 
        will use these to start the fit.  Otherwise, it will estimate them
        on its own.  

        See Notes for description of each element.
        
    positive : {False, True}, optional
        if `p0` is not given, set to search for a positive peak.  If 
        `p0` is not given, and `positive` or `negative` is not given, 
        the function determines it automatically. 

    negative : {False, True}, optional
        if `p0` is not given, set to search for a negative peak.  If 
        `p0` is not given, and `positive` or `negative` is not given, 
        the function determines it automatically. 

    nan : {True, False}, optional
        Set to True to ignore NaN values.

    ignore_optimizewarning : {False, True}, optional
        Set to True to suppress any OptimizeWarning messages.

    Returns
    -------
    dict
        `"fit"` : numpy.ndarray of float
            (nx,) array of fitted values to `y`.  

        `"parms"` : numpy.ndarray of float
            (ncoeffs,) array of fitted parameters.

        `"goodfit"` : bool
            Set to True if the fit is good.

    Notes
    -----

    Let u = (x - p[1])/p[2].  The two possible models are then:

                  GAUSSIAN                Lorentzian

    Model     p[0]*exp(-0.5*u^2)        p[0]/(u^2 + 1) 

    P[0]         Peak Value               Peak Value    
    p[1]        Peak Centroid           Peak Centroid  
    p[2]       Gaussian Sigma               HWHM%      
    p[3]         + p[3]                    + p[3]      
    p[4]         + p[4]*x                  + p[4]*x    
    p[5]         + p[5]*x**2               + p[5]*x**2
    p[6]           ....                      ....
                                         
    Notes: 
          % Half-width at half maximum
          * Optional depending on NTERMS
    

    """

    # Convert inputs to numpy arrays

    x = np.array(x)
    y = np.array(y)
    
    if p0 is None:

        # Estimate the start parameters
        
        p0 = cmest(x, y, positive=positive, negative=negative, nan=nan)

        if nparms == 3:

            p0 = p0[0:3]
            
        elif nparms == 4:

            p0 = p0

        else:  
        
            p0 = np.append(p0, np.full(nparms-4,0.0))
                    
    else:

        p0 = np.array(p0)
        
        # Check to make sure the number of user parameter estimates equals
        # the number of parameters reqeusted.

    if len(p0) != nparms:

        exception = 'fitpeak1d: number parameter estimates must ' + \
          'equal the number of parameters.'

        raise Exception(exception)

    # Select the function

    if type.lower() == 'gaussian':

        f = gauss1d    
        
    elif type.lower() == 'lorentzian':

        f = lorentz1d
        
    else:

        exception = 'fitpeak1d: Unknown type.'
        raise Exception(exception)        

    # Now do the call
    
    try:

        if ignore_optimizewarning is True:

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", scipy.optimize.OptimizeWarning)

                popt, pcov = scipy.optimize.curve_fit(f, x, y, p0=p0)

        else:

            popt, pcov = scipy.optimize.curve_fit(f, x, y, p0=p0)            
                
        fit = {'parms': popt}
        fit['fit'] = f(x, *popt)
        fit['goodfit'] = True
            
    except (RuntimeError, scipy.optimize.OptimizeWarning):

        fit = {'parms': np.full_like(p0, np.nan)}
        fit['fit'] = np.full_like(y, np.nan)
        fit['goodfit']=False
        
    return fit
    
        
def lorentz1d(x, amp, mean, hwhm, *args):

    """
    Computes values of a lorentzian function given user parameters.

    The function is given by,

    u = (x - mean)/hwhm
    y =  amp/(u**2 + 1) + aergs[0] + args[1]*x

    Parameters
    ----------
    x : array-like
        (nx,) array of independent variables.

    amp: float
        The amplitude

    mean: float
        The mean value

    hwhm: float
        The Half Width at Half Maximum.

    args[0] : float, optional
        A constant background level

    args[1] : float, optional
        The slope of the backgrond level given by args[0] + args[1]*x

    Returns
    -------
    numpy.ndarray
        The lorentzian function evaluated at `x`. 

    Notes
    -----
    Based on the mpfitpeak_lorentz code written by Craig Markwardt.

    $Id: mpfitpeak.pro,v 1.19 2011/12/08 17:51:33 cmarkwar Exp $

    Copyright (C) 1997-2001, 2003, 2005, 2007, 2008, 2009, 2010,
    Craig Markwardt.  This software is provided as is without any
    warranty whatsoever.  Permission to use, copy, modify, and
    distribute modified or unmodified copies is granted, provided
    this copyright and disclaimer are included unchanged.


    """

    u = (x-mean)/hwhm

    # Deal with the baseline
        
    nargs = len(args)

    bl = 0
    for i, coeff in enumerate(args):
        bl += coeff * x**i

    

#    bl = args[0] if nargs >= 1 else 0
#    if nargs == 2:
#        bl = bl + args[1]*x

    #   Compute the values and return 

    return bl + amp/(u**2 + 1)


def gauss1d(x, amp, mean, sigma, *args):

    """
    Computes values of a gaussian function given user parameters.

    The function is given by,

    u = (x - mean)/hwhm
    y =  amp*exp(-0.5*u**2)

    Parameters
    ----------
    x : array-like
        (nx,) array of independent variables.

    amp: float
        The amplitude

    mean: float
        The mean value

    sigma: float
        The gaussian standard deviation 

    args[0] : float, optional
        A constant background level

    args[1] : float, optional
        The slope of the backgrond level given by args[0] + args[1]*x

    Returns
    -------
    numpy.ndarray
        The gaussian function evaluated at `x`. 

    Notes
    -----
    Based on the mpfitpeak_gauss code written by Craig Markwardt.

    $Id: mpfitpeak.pro,v 1.19 2011/12/08 17:51:33 cmarkwar Exp $

    Copyright (C) 1997-2001, 2003, 2005, 2007, 2008, 2009, 2010,
    Craig Markwardt.  This software is provided as is without any
    warranty whatsoever.  Permission to use, copy, modify, and
    distribute modified or unmodified copies is granted, provided
    this copyright and disclaimer are included unchanged.

    """

    u = (x-mean)/sigma
    uz = np.exp(-0.5*u**2)

    # Deal with the baseline
        
    nargs = len(args)
    bl = args[0] if nargs >= 1 else 0
    if nargs == 2:
        bl = bl + args[1]*x

    #   Compute the values and return

    return bl + amp*uz


def cmest(x, y, nan=False, positive=False, negative=False):

    """
    To estimate the parameters for peak fitting.


    Parameters
    ----------
    x : arrqy-like
        (nx,) array of independent variables.

    y : arrqy-like
        (nx,) array of dependent variables.

    nan : {False, True}, optional
        Set to True to ignore NaN values.

    positive : {False, True}, optional
        Set to True if the peak is positive.  If neither `positive` nor
        `negative` is set to True, the function decides automatically.

    negative : {False, True}, optional
        Set to True if the peak is negative.  If neither `positive` nor
        `negative` is set to True, the function decides automatically.

    Returns
    -------
    numpy.ndarray of float
        array[0] : the peak value 

        array[1] : the centroid location 

        array[2] : the width of the peak

        array[3] : the average height

    Notes
    -----
    A line-by-line copy of the IDL mpfitpeak_est program written by 
    Craig Markwardt.

    $Id: mpfitpeak.pro,v 1.19 2011/12/08 17:51:33 cmarkwar Exp $

    Copyright (C) 1997-2001, 2003, 2005, 2007, 2008, 2009, 2010,
    Craig Markwardt.  This software is provided as is without any
    warranty whatsoever.  Permission to use, copy, modify, and
    distribute modified or unmodified copies is granted, provided
    this copyright and disclaimer are included unchanged.

    If neither positive nor negative are True, then the largest
    magnitude peak is identified.


    """

    # Here is the secret - the width is estimated based on the area
    # above/below the average.  Thus, as the signal becomes more
    # noisy the width automatically broadens as it should.

    # Determine which sum function to use based on the nan keyword
    
    sumfunc = np.sum if nan is False else np.nansum
    
    nx = len(x)
    
    idx = np.argsort(x)
    xs = x[idx]
    ys = y[idx]

    minx = np.min(xs)
    maxx = np.max(xs)

    miny = np.min(ys)
    maxy = np.max(ys)    

    dx = 0.5*np.concatenate(([xs[1]-xs[0]], xs[2:]-xs[0:-2],
                             [xs[-1]-xs[-2]]))
    totarea = sumfunc(dx*ys)  # Total area under the curve
    av = totarea/(maxx-minx)  # Average height

    # Check degenerate case:  all flat with no noise    

    if miny == maxy:  

        return np.array([0, xs[nx//2], (xs[-1]-xs[0])//2, ys[0]])

    # Compute the spread in values above and below average... we
    # take the narrowest one as the one with the peak
    
    wh1 = y >= av
    ct1 = np.sum(wh1)

    wh2 = y <= av
    ct2 = np.sum(wh2)

    if ct1 == 0 or ct2 == 0:

        raise Exception('Average y value should fall within the range' + \
                        'of y values but does not.')
                         
    sd1 = np.sum(x[wh1]**2)/ct1 - (np.sum(x[wh1])/ct1)**2
    sd2 = np.sum(x[wh2]**2)/ct2 - (np.sum(x[wh2])/ct2)**2

    # Compute area above/below average
    
    if sd1 < sd2 or positive is True:  # Positive peak
        cent = x[y == maxy]
        cent = cent[0]
        peak = maxy-av
        
    if sd1 >= sd2 or negative is True:  # Negative peak
        cent = x[y == miny]
        cent = cent[0]
        peak = miny-av

    peakarea = totarea - sumfunc(dx*np.where(ys < av, ys, av))
    if peakarea == 0:
        peak = 0.5*peakarea
    width = peakarea/(2*np.abs(peak))
    if width == 0 or width == np.nan:
        width = np.median(dx)

    est = np.array([peak, cent, width, av])

    return est
        
