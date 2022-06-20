import numpy as np

def polyfit1d(x,y,order,yunc=None,doalpha=False,justfit=False,silent=True):

    '''
    Fits a polynomial of a given order to a set of 1-D data.

    Input Parameters
    ----------------
    x : array_like
        An array of independent values

    y : array_like
        An array of dependent values

    order : int
        The order of the polynomial

    yunc : array_like, optional
        An array of uncertainties on the dependent values `y`

    doalpha: {False, True}, optional
        If true, then A**T * A is computed directly instead of A

    justfit : {False, True}, optional 
        Set to compute only the coefficents

    silent : {True, False}, optional
        If False, the result of the fit will be written to the command line

    Output Parameters
    -----------------
    dict
        A dict where with the following entries:

        coeffs : numpy.ndarray
            The polynomial coefficients

        var : numpy.ndarray
            The variances of the coefficients

        covar : numpy.ndarray
            The covariance matrix

        yfit : numpy.ndarray 
            The polynomial evaluated at `x`

        nparm : int
            The number of parameters of the fit

        ndof : int
            The number of degrees of freedom

        chi2 : float
            The chi^2 value of the fit

        rchi2 : float
            The reduced chi^2 value of the fit

        rms : flat
            The rms of the fit (1/N, not 1/(N-1))

    Notes
    -----
    Standard least squares fitting (e.g. Numerical Recipes).  
    The program fits the following function:

    z = f(x) = c_0 + c_1*x + c_2*x**2 + c_3*x**3 + ...

    By default, the design matrix A is constructed.  The alpha 
    (A^T ## A) and beta (A^T ## b) arrays of the normal equations 
    are then constructed and then passed to np.linalg.solve.  

    If the dataset to be fit is large and/or the number of coefficients 
    is large, the design matrix A can be large (ncoeff,ndat) and may 
    take too much memory and/or too much time to construct.  In this case,
    set `doalpha` to True and the alpha and beta arrays are constructed 
    directly.  There is a trade off since constructing the alpha and beta 
    arrays directly takes longer for smaller surfaces.


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
        2022-06-17 - Simplfied the creation of alpha is doalpha=True
    '''

# Construct the uncertainty array if need be
    
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
    
    exp = np.arange(0,order+1) # exoponents of the basis functions
    
    alpha = np.zeros((ncoeffs,ncoeffs))
    beta = np.zeros(ncoeffs)

    b = yy/yyunc

# Now start the linear algebra construction

    if doalpha is True:
    
        for i in range(0,ncoeffs):
            
                for j in range(i,ncoeffs):
        
                    at = (xx**exp[i])/yyunc
                    a  = (xx**exp[j])/yyunc
                    
                    alpha[i,j] = np.sum(at*a)
                    alpha[j,i] = np.sum(at*a)
                    beta[i] = np.sum(at*b)
                    
    elif doalpha is False:

        AT = np.empty((ncoeffs,ndat))
        for i in range(ncoeffs):

            AT[i,:] = xx**exp[i]/yyunc

        alpha = np.matmul(AT,np.transpose(AT))
        beta = np.matmul(AT,b)        
                
# Solve things
    
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
            print('         Number of parameters = ',ncoeffs)
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

