import numpy as np
from poly2d import poly2d

def polyfit2d(x,y,z,xorder,yorder,zunc=None,doalpha=False,silent=True,\
             justfit=False):

    '''
    Fits a 2D polynomial of a given order to a set of 2D data

    Standard least squares fitting (e.g. Numerical Recipes).  This 
    program is based on solving the equation A ## coeff = b.  Except 
    now, z = f(x,y) = (a_0 + a_1*x) * (b_0 + b_1*y + b_2*y^2)
    = b_0*a_0 + b_0*a_1*x + b_1*a_0*x + b_1*a_0 + b_1*a_1*x*y...

    The design matrix A is constructed.  The alpha (A^T ## A) and 
    beta (A^T ## b) arrays of the normal equations are constructed
    and then passed to np.linalg.solve.  

    If the surface to be fit is large and/or the number of coefficients 
    is large, set the DOALPHA keyword.  Since the design matrix A can be 
    large (ncoeff,ndat) it may take too much memory and or too much time.  
    If `doalpha` is True, the alpha and beta arrays are constructed directly 
    and passed to GAUSSJ.  There is a trade off since constructing the 
    alpha and beta arrays directly takes longer for smaller surfaces.

    Input Parameters
    ----------------
    x : array_like, int or float
        an array of independent values

    y : array_like, int of float
        an array of independent values

    z : array_like, int or float 
        an array of dependent values

    xorder : int
        the order of the polynomial for the x dimension

    yorder : int
        the order of the polynomial for the y dimension

    zunc : array_like, float, optional
        an array of uncertainties on the dependent values

    doalpha: {False, True}, optional
        If true, then A**T * A is computed directly instead of A

    silent : {True, False}, optional
        If False, the result of the fit will be written to the command line

    justfit : {False, True}, optional
        If True, only the coefficients are computed and returned

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

        zfit : numpy.ndarray 
            the polynomial evaluated at `x` and `y`

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

    Examples
    --------
    later?

    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_polyfit2d.pro
    '''    

    if zunc is None: zunc = np.full(len(x),1.0)

# Get rid of NaNs 

    znonan = ~np.isnan(z)
    nnan = np.size(znonan)-np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    zz = z[znonan]
    zzunc = zunc[znonan]

# Get basic parameters
    
    ndat = len(xx)      # number of idependent variables
    ncoeffs = (xorder+1)*(yorder+1)   # number of coefficients
    ndof = ndat-ncoeffs # number of degrees of freedom    

# Determine the exponents of the basis functions

    xexp = np.arange(xorder+1,dtype=int)
    yexp = np.arange(yorder+1,dtype=int)
    order = np.empty((ncoeffs,2),dtype=int)
    idx = 0

    for i in range(yorder+1):

        for j in range(xorder+1):

            order[idx,0] = xexp[j]
            order[idx,1] = yexp[i]
            idx +=1

# Create the b array

    b = zz/zzunc

# Now start the linear algebra construction
    
    if doalpha is True:

        alpha = np.empty((ncoeffs,ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i,ncoeffs):

                at = (xx**order[i,0] * yy**order[i,1])/zzunc
                a  = (xx**order[j,0] * yy**order[j,1])/zzunc
                    
                alpha[j,i] = np.sum(at*a)
                beta[i] = np.sum(at*b)                    
                
        alpha +=np.transpose(alpha)

        zdiag = np.arange(0,ncoeffs*ncoeffs,ncoeffs+1)
        zdiag = np.unravel_index(zdiag,(ncoeffs,ncoeffs))
        alpha[zdiag] = alpha[zdiag]/2.

    elif doalpha is False:

        AT = np.empty((ncoeffs,ndat))
        for i in range(ncoeffs):

            AT[i,:] = xx**order[i,0]*yy**order[i,1]/zzunc

        alpha = np.matmul(AT,np.transpose(AT))
        beta = np.matmul(AT,b)
        
# Solve things (need to remember what you are doing...)
    
    coeffs = np.linalg.solve(alpha,beta)

    if justfit is True:

# Just return the coefficients
        
        return({"coeffs":coeffs,"var":None,"covar":None,"yfit":None,\
            "nparm":None,"ndof":None,"chi2":None,"rchi2":None,"rms":None})   

    elif justfit is True:

# Keep going with other outputs and possible command line output
        
        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        zfit = poly2d(x,y,xorder,yorder,coeffs) 

        residual = zz-zzfit
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
                    npsqrt(var[i]),sep='')

            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')
            

        return({"coeffs":coeffs,"var":var,"covar":covar,"yfit":yfit,\
            "nparm":order+1,"ndof":ndof,"chi2":chi2,"rchi2":rchi2,"rms":rms})


        
        


