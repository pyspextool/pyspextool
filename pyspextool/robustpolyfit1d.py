from numpy import isnan as npisnan
from numpy import size as npsize
from numpy import sum as npsum
from numpy import std as npstd
from numpy import full as npfull
from numpy import mean as npmean
from numpy import std as npstd
from numpy import abs as npabs
from numpy import sqrt as npsqrt
from numpy.polynomial.polynomial import polyval as nppolyval
from polyfit1d import polyfit1d

def robustpolyfit1d(x,y,order,thresh,eps,yunc=None,goodbad=None,silent=True):

    '''
    Fits a "robust" polynomial of a given order to a set of 1-D data.

    Input Parameters
    ----------------
    x : numpy.ndarray
        an array of independent values

    y : numpy.ndarray
        an array of dependent values

    order : int
        the order of the polynomial

    thresh : float
        sigma threshold to identify outliers.

    eps : float
        limit of the fractional change in the standar deviation of the fit.
      
    yunc : numpy.ndarray, optional
        an array of uncertainties on the dependent values

    goodbad : np.ndarray, optional 


    silent : {True, False}, optional
        If False, the result of the fit will be written to the command line

    Output Parameters
    -----------------
    dict
        a dict where with the following entries:

        coeffs : numpy.ndarray
            the polynomial coefficients

        var : numpy.ndarray
            the variances of the coefficients

        covar : numpy.ndarray
            the covariance matrix

        yfit : numpy.ndarray 
            the polynomial evaluated at `x`

        nparm : int
            the number of parameters of the fit

        ndof : int
            the number of degrees of freedom

        chi2 : float
            the chi^2 value of the fit

        rchi2 : float
            the reduced chi^2 value of the fit

        rms : flat
            the rms of the fit (1/N, not 1/(N-1))

    Procedure
    ---------
 

    Examples
    --------
    > import numpy as np
    > import math
    > N = 100
    > x = np.linspace(0, 4, N)
    > y = x**3 - 6*x**2 + 12*x - 9 
    > y += np.random.normal(0, 2, size=len(y))
    > sigma = 1.5
    > yerr = np.ones(N)*sigma 
    > y[50] = 30
    > y[30] = math.nan
    > 
    > result = robustpolyfit1d(x,y,3,3,0.1,yunc=yerr,silent=False)
    > yfit = result['yfit']
    > pl.plot(x,y,'o')
    > pl.plot(x,yfit,'r-')
    > z = np.where(result['ogoodbad'] == 0)
    > pl.plot(x[z],y[z],'ro')
    > pl.show()

    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_robustpoly1d.pro
    '''

# Check to see if the user passed an uncertainty array.
    
    if yunc is None: yunc = npfull(len(x),1.0)

# Set up the ogoodbad array

    ogoodbad = npfull(len(x),2)
        
# Get rid of NaNs

    znotnan = ~npisnan(y)
    init_goodcnt = npsum(znotnan)
    nnan = npsum(~znotnan)
    xx = x[znotnan]
    yy = y[znotnan]
    yyunc = yunc[znotnan]

    ogoodbad[znotnan] = 1

# Do a first pass of the data

    fit = polyfit1d(xx,yy,order,yunc=yyunc,silent=True)

    
# Compute the residuals and the mean and sigma of the residuals
    
    residual = yy-fit['yfit']

    mean = npmean(residual)
    stddev  = npstd(residual)

# Now check for residual points that are outside the threshhold
    
    good = npabs((residual-mean)/stddev) <= thresh
    goodcnt = npsum(good)

    if goodcnt != init_goodcnt:

        itter = 0
        for i in range(0,11):

            itter +=1

            oldstddev = stddev
            fit = polyfit1d(xx[good],yy[good],order,yunc=yyunc[good],silent=True)
            residual = yy[good]-fit['yfit']
            mean = npmean(residual)
            stddev  = npstd(residual)        

# Check to see if the new stddev isn't much of a change from the old one

            if (oldstddev-stddev)/oldstddev < eps: break

# Now generate a new full residual array

            residual = yy-nppolyval(xx,fit['coeffs'])

            good = npabs((residual-mean)/stddev) <= thresh
            goodcnt = npsum(good)

# Let's reconstuct the goodbad array

    tmp = npfull(init_goodcnt,0)
    tmp[good] = 1
    ogoodbad[znotnan] = tmp

# Create a full yfit

    yfit = nppolyval(x,fit['coeffs'])

# Report results if requested

    if silent is not True:

        nbad = npsum(good == 0)

        print(' ')
        print('    Number of original points = ',len(x))
        print('              Sigma threshold = ',thresh)
        print('              Number of loops = ',itter)
        print('           Number of outliers = ',nbad)
        print('          Number of NaNs in y = ',nnan)
        print('         Number of parameters = ',order+1)
        print(' Number of degrees of freedom = ',fit["ndof"])
        print('                  Chi-Squared = ',fit["chi2"])
        print('          Reduced Chi-Squared = ',fit["rchi2"])
        print('         RMS deviation of fit = ',fit["rms"])
        print(' ')
        print('Coefficients:')
        print(' ')

        for i in range(0,order+1):
        
            print('Coeff #',str(i).zfill(2),': ',fit['coeffs'][i],'+-',\
                  npsqrt(fit['var'][i]),sep='')

        print(' ')
        print('Covariance Matrix:')
        print(fit['covar'])
        print(' ')

    return({"coeffs":fit['coeffs'],"var":fit['var'],"covar":fit['covar'],
            "ogoodbad":ogoodbad,"yfit":yfit,"nparm":order+1,"ndof":fit['ndof'],
            "chi2":fit['chi2'],"rchi2":fit['rchi2'],"rms":fit['rms']})    

