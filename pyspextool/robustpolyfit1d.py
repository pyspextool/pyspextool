import numpy as np
from polyfit1d import polyfit1d

def robustpolyfit1d(x,y,order,thresh,eps,yunc=None,goodbad=None,justfit=False,\
                    silent=True):

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
        an array identifying good and bad pixels.  0=bad, 1=good, 2=NaN

    justfit : {False, True}, optional
        set to only compute the coefficients

    silent : {True, False}, optional
        if False, the result of the fit will be written to the command line

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
    
    if yunc is None: yunc = np.full(len(x),1.0)

# Set up the ogoodbad array

    ogoodbad = np.full(len(x),2)
        
# Get rid of NaNs

    znotnan = ~np.isnan(y)
    init_goodcnt = np.sum(znotnan)
    nnan = np.sum(~znotnan)
    xx = x[znotnan]
    yy = y[znotnan]
    yyunc = yunc[znotnan]

    ogoodbad[znotnan] = 1

# Do a first pass of the data

    fit = polyfit1d(xx,yy,order,yunc=yyunc,justfit=justfit,silent=True)
    
# Compute the residuals and the mean and sigma of the residuals
    
    residual = yy-np.polynomial.polynomial.polyval(xx,fit['coeffs'])

    mean = np.mean(residual)
    stddev  = np.std(residual)

# Now check for residual points that are outside the threshhold
    
    good = np.abs((residual-mean)/stddev) <= thresh
    goodcnt = np.sum(good)

    if goodcnt != init_goodcnt:

        itter = 0
        for i in range(0,11):

            itter +=1

            oldstddev = stddev
            fit = polyfit1d(xx[good],yy[good],order,yunc=yyunc[good],\
                            justfit=justfit,silent=True)
            residual = yy[good]-\
              np.polynomial.polynomial.polyval(xx[good],fit['coeffs'])
            mean = np.mean(residual)
            stddev  = np.std(residual)        

# Check to see if the new stddev isn't much of a change from the old one

            if (oldstddev-stddev)/oldstddev < eps: break

# Now generate a new full residual array

            residual = yy-np.polynomial.polynomial.polyval(xx,fit['coeffs'])

            good = np.abs((residual-mean)/stddev) <= thresh
            goodcnt = np.sum(good)

# Let's reconstuct the goodbad array

    tmp = np.full(init_goodcnt,0)
    tmp[good] = 1
    ogoodbad[znotnan] = tmp

# Create a full yfit

    yfit = np.polynomial.polynomial.polyval(x,fit['coeffs'])

# Report results if requested

    if silent is not True:

        nbad = np.sum(good == 0)

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
                  np.sqrt(fit['var'][i]),sep='')

        print(' ')
        print('Covariance Matrix:')
        print(fit['covar'])
        print(' ')

    return({"coeffs":fit['coeffs'],"var":fit['var'],"covar":fit['covar'],
            "ogoodbad":ogoodbad,"yfit":yfit,"nparm":order+1,"ndof":fit['ndof'],
            "chi2":fit['chi2'],"rchi2":fit['rchi2'],"rms":fit['rms']})    

