import numpy as np

def polyfit1d(x,y,order,yunc=None,justfit=False,silent=True):

    '''
    Fits a polynomial of a given order to a set of 1-D data.

    Input Parameters
    ----------------
    x : numpy.ndarray
        an array of independent values

    y : numpy.ndarray
        an array of dependent values

    order : int
        the order of the polynomial

    yunc : numpy.ndarray, optional
        an array of uncertainties on the dependent values

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
    compute the alpha and beta matrices of the normal equations and then 
    solve np.linalg.solve.  See e.g., Numerical Recipes for details.  

    Examples
    --------
    > x = np.array([10. ,20. ,30. ,40. ,50., 60. ,70. ,80. ,90.,100.])
    > y = np.array([0.37,0.58,0.83,1.15,1.36,1.62,1.90,2.18,2.45,math.nan])
    > yerr = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
    > result = polyfit1d(x,y,1,yunc=yerr)

    Data from Chapter 6 of Data Reduction and Error Analysis for the 
    Physical Sciences by Bevington & Robinson


    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_polyfit1d.pro
        2022-06-07 - Added the justfit parameter
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

    if justfit is True:

        return({"coeffs":coeffs,"var":None,"covar":None,"yfit":None,\
            "nparm":None,"ndof":None,"chi2":None,"rchi2":None,"rms":None})   
        
    elif justfit is False:

# Keep going with other outputs and possible command line output
        
        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        yfit = np.polynomial.polynomial.polyval(x,coeffs)
        residual = y-yfit
        rms = np.std(residual[znonan])

        chi2 = np.sum( (residual[znonan]/yunc[znonan])**2)
        rchi2 = chi2/ndof
    
# Report results if requested

        if silent is False:
            
            print(' ')
            print('             Number of points = ',len(x))
            print('          Number of NaNs in y = ',nnan)
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
            
            
    return({"coeffs":coeffs,"var":var,"covar":covar,"yfit":yfit,\
            "nparm":order+1,"ndof":ndof,"chi2":chi2,"rchi2":rchi2,"rms":rms})

