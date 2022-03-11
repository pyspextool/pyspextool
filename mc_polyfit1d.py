# the np.std function does 1/N not 1/N-1

import numpy as np

def mc_polyfit1d(x,y,order,yunc=None,silent=None):

    '''
    Fits a polynomial of a given order to a set of 1-D data.

    Input Parameters:
        x - A numpy array of independent values.
        y - A numpy array of dependent values.
        order - The order of the polynomial.

    Optional Parameters:
        yunc   - A numpy array of uncertainties on the dependent values.
        silent - Set to note report the results at the command line. 

    Output Parameters:
        A dict where with the following keys:
        coeffs - the polynomial coefficients
        var    - the variances of the coefficients
        covar  - the covariance matrix
        yfit   - the polynomial evaluated at x.
        nparm  - the number of parameters of the fit
        ndof   - the number of degrees of freedom
        chi2   - the chi^2 value of the fit
        rchi2  - the reduced chi^2 value of the 
        rms    - the rms of the fit

    Procedure:
        This program is based on solving the equation A ## coeff = b.  The 
        alpha (A^T ## A) and beta (A^T ## b) arrays of the normal equations 
        are constructed and then solved with np.linalg.solve.

    Example:
        NA

    Modification History:
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the gethdrinfo.pro IDL program.
    '''

    
    if yunc is None: yunc = np.full(len(x),1.0)

# Get rid of NaNs

    znonan = ~np.isnan(y)
    nnan = np.size(znonan)-np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    yyunc = yunc[znonan]

# Get basic parameters
    
    ndat = len(xx)      # number of idependent variables
    ncoeffs = order+1   # number of coefficients
    ndof = ndat-ncoeffs # number of degrees of freedom

#  Construct the alpha and beta matrix of the normal equations.  
#  Build only the upper triangle of the alpha array since it is symmetric.  
    
    exp = np.arange(0,order+1)
    alpha = np.zeros((ncoeffs,ncoeffs))
    beta = np.zeros(ncoeffs)

    b = yy/yyunc

    for i in range(0,ncoeffs):

        for j in range(i,ncoeffs):

            at = (xx**exp[i])/yyunc
            a  = (xx**exp[j])/yyunc

            alpha[i,j] = np.sum(at*a)
            beta[i] = np.sum(at*b)

# Now transpose and add to get the other side

    alpha = alpha+np.transpose(alpha)

# Finally, divide the diagonal elements by 2 to make up for the addition
# in the transpose

    zdiag = np.arange(0,ncoeffs*ncoeffs,ncoeffs+1)
    zdiag = np.unravel_index(zdiag,(ncoeffs,ncoeffs))
    alpha[zdiag] = alpha[zdiag]/2.

# Solve things (need to remember what you are doing...)
    
    coeffs = np.linalg.solve(alpha,beta)
    covar = np.linalg.inv(alpha)

    var = np.diagonal(covar)

    yfit = np.polynomial.polynomial.polyval(x,coeffs)
    residual = y-yfit
    rms = np.std(residual[znonan])

    chi2 = np.sum( (residual[znonan]/yunc[znonan])**2)
    rchi2 = chi2/ndof
    
# Report results if requested

    if silent is not True:

        print(' ')
        print('             Number of points = ',len(x))
        print('               Number of Nans = ',nnan)
        print('         Number of parameters = ',order+1)
        print(' Number of degrees of freedom = ',ndof)
        print('                  Chi-Squared = ',chi2)
        print('          Reduced Chi-Squared = ',rchi2)
        print('         RMS deviation of fit = ',rms)
        print(' ')
        print('Coefficients:')
        print(' ')

        for i in range(0,order+1):
        
            print('Coeff #',str(i).zfill(2),': ',coeffs[i],'+-',\
                  np.sqrt(var[i]),sep='')

        print(' ')
        print('Covariance Matrix:')
        print(covar)
        print(' ')


    return({"coeffs":coeffs,"var":var,"covar":covar,"yfit":yfit,"nparm":order+1,"ndof":ndof,
            "chi2":chi2,"rchi2":rchi2,"rms":rms})
