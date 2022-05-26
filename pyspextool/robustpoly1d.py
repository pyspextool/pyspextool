import numpy as np
from mc_polyfit1d import mc_polyfit1d
import matplotlib.pyplot as pl

def robustpoly1d(x,y,order,thresh,eps,yunc=None,goodbad=None,justfit=None,silent=None):

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

    fit = mc_polyfit1d(xx,yy,order,yunc=yyunc,silent=True)

# Compute the residuals and the mean and sigma of the residuals
    
    residual = yy-fit['yfit']

    mean = np.mean(residual)
    stddev  = np.std(residual)

# Now check for residual points that are outside the threshhold
    
    good = np.abs((residual-mean)/stddev) <= thresh
    good[5] = ~good[5]
    goodcnt = np.sum(good)

    if goodcnt != init_goodcnt:

        itter = 0
        for i in range(0,11):

            itter +=1

            oldstddev = stddev
            fit = mc_polyfit1d(xx[good],yy[good],order,yunc=yyunc[good],
                               silent=True)
            residual = yy[good]-fit['yfit']
            mean = np.mean(residual)
            stddev  = np.std(residual)        

# Check to see if the new stddev isn't much of a change from the old one

            if (oldstddev-stddev)/oldstddev < eps: break

# Now generate a new full residual array

            residual = yy-np.polynomial.polynomial.polyval(xx,fit['coeffs'])

            good = np.abs((residual-mean)/stddev) <= thresh
            goodcnt = np.sum(good)

# Let's reconstuct the goodbad array

    tmp = np.zeros(init_goodcnt)
    tmp[good] = 1
    ogoodbad[znotnan] = tmp

# Report results if requested

    if silent is not True:

        nbad = np.sum(good == 0)

        print(' ')
        print('    Number of original points = ',len(x))
        print('              Sigma threshold = ',thresh)
        print('              Number of loops = ',itter)
        print('           Number of outliers = ',nbad)
        print('               Number of Nans = ',nnan)
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
            "yfit":fit['yfit'],"nparm":order+1,"ndof":fit['ndof'],
            "chi2":fit['chi2'],"rchi2":fit['rchi2'],"rms":fit['rms']})
