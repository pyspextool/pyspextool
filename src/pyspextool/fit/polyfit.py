import numpy as np


def poly2d(x, y, xorder, yorder, coeffs):

    """
    evaluates a polynomial function of two independent variables
    

    Parameters
    ----------------
    x : array_like, int or float
        an array of independent values

    y : array_like, int or float
        an array of independent values

    xorder : int
        the order of the polynomial for the x dimension

    yorder : int
        the order of the polynomial for the y dimension

    coeffs : np.ndarray
        an array of coefficients from polyfit2d

    Output Parameters
    -----------------
    numpy.ndarray
        the polynomial evaluated at `x` and `y`

    Examples
    --------
    later?

    Modification History
    --------------------
        2022-06-07 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_poly2d.pro
    """

    # Get set up with a bunch of necessary numbers and arrays
    
    ndat = x.size
    ncoeffs = (xorder+1)*(yorder+1)

    xexp = np.arange(xorder+1, dtype=int)
    yexp = np.arange(yorder+1, dtype=int)

    z = np.empty(ndat)

    # Now fill in the data
    
    for i in range(ndat):

        z[i] = np.sum(np.ravel(np.outer(y[i]**yexp, x[i]**xexp))*coeffs)

    return z


def polyfit1d(x, y, order, yunc=None, doalpha=False, justfit=False, silent=True):

    """
    Fits a polynomial of a given order to a set of 1-D data.

    Parameters
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

    Returns
    -------
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
    directly.  There is a trade-off since constructing the alpha and beta
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
    """

    # Construct the uncertainty array if need be
    
    if yunc is None:
        yunc = np.full(len(x), 1.0)

    # Get rid of NaNs 

    znonan = ~np.isnan(y)
    nnan = np.size(znonan)-np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    yyunc = yunc[znonan]

    # Get basic parameters
    
    ndat = len(xx)       # number of idependent variables
    ncoeffs = order+1    # number of coefficients
    ndof = ndat-ncoeffs  # number of degrees of freedom

    #  Construct the alpha and beta matrix of the normal equations.  
    
    exp = np.arange(0, order+1)  # exoponents of the basis functions
    
    alpha = np.zeros((ncoeffs, ncoeffs))
    beta = np.zeros(ncoeffs)

    b = yy/yyunc

    # Now start the linear algebra construction

    if doalpha is True:
    
        for i in range(0, ncoeffs):

            for j in range(i, ncoeffs):
        
                at = (xx**exp[i])/yyunc
                a = (xx**exp[j])/yyunc

                alpha[i, j] = np.sum(at*a)
                alpha[j, i] = np.sum(at*a)
                beta[i] = np.sum(at*b)

    elif doalpha is False:

        at = np.empty((ncoeffs, ndat))
        for i in range(ncoeffs):

            at[i, :] = xx**exp[i]/yyunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)
                
    # Solve things
    
    coeffs = np.linalg.solve(alpha, beta)

    if justfit is True:

        return({"coeffs": coeffs, "var": None, "covar": None, "yfit": None,
                "nparm": None, "ndof": None, "chi2": None, "rchi2": None, "rms": None})
        
    elif justfit is False:

        # Keep going with other outputs and possible command line output
        
        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        yfit = np.polynomial.polynomial.polyval(x, coeffs)
        residual = y-yfit
        rms = np.std(residual[znonan])

        chi2 = np.sum((residual[znonan]/yunc[znonan])**2)
        rchi2 = chi2/ndof
    
        # Report results if requested

        if silent is False:
            
            print(' ')
            print('             Number of points = ', len(x))
            print('          Number of NaNs in y = ', nnan)
            print('         Number of parameters = ', ncoeffs)
            print(' Number of degrees of freedom = ', ndof)
            print('                  Chi-Squared = ', chi2)
            print('          Reduced Chi-Squared = ', rchi2)
            print('         RMS deviation of fit = ', rms)
            print(' ')
            print('Coefficients:')
            print(' ')
            
            for i in range(0, order+1):
            
                print('Coeff #', str(i).zfill(2), ': ', coeffs[i], '+-',
                      np.sqrt(var[i]), sep='')
                    
            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')

    return({"coeffs": coeffs, "var": var, "covar": covar, "yfit": yfit,
            "nparm": order+1, "ndof": ndof, "chi2": chi2, "rchi2": rchi2,
            "rms": rms})


def polyfit2d(x, y, z, xorder, yorder, zunc=None, doalpha=False, silent=True,
              justfit=False):

    """
    Fits a 2D polynomial of a given order to a set of 2D data

    Standard least squares fitting (e.g. Numerical Recipes).  This 
    program is based on solving the equation A ## coeff = b.  Except 
    now, 

    
    Parameters
    ----------------
    x : array_like, int or float
        An array of independent values

    y : array_like, int of float
        An array of independent values

    z : array_like, int or float 
        An array of dependent values

    xorder : int
        The order of the polynomial for the x dimension

    yorder : int
        The order of the polynomial for the y dimension

    zunc : array_like, float, optional
        An array of uncertainties on the dependent values

    doalpha: {False, True}, optional
        If true, then A**T * A is computed directly instead of A

    silent : {True, False}, optional
        If False, the result of the fit will be written to the command line

    justfit : {False, True}, optional
        If True, only the coefficients are computed and returned

    Returns
    -------
    dict
        A dict where with the following entries:

        coeffs : numpy.ndarray
            The polynomial coefficients

        var : numpy.ndarray
            The variances of the coefficients

        covar : numpy.ndarray
            The covariance matrix

        zfit : numpy.ndarray 
            The polynomial evaluated at `x` and `y`

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

    z = f(x,y) = (a_0 + a_1*x) * (b_0 + b_1*y + b_2*y**2)
      = b_0*a_0 + b_0*a_1*x + b_1*a_0*x + b_1*a_0 + b_1*a_1*x*y + ...
      = c_0 + c_1*x + c_2*x*y + c_3*x*y**2 + ...

    By default, the design matrix A is constructed.  The alpha 
    (A^T ## A) and beta (A^T ## b) arrays of the normal equations 
    are then constructed and then passed to np.linalg.solve.  

    If the dataset to be fit is large and/or the number of coefficients 
    is large, the design matrix A can be large (ncoeff,ndat) and may 
    take too much memory and/or too much time to construct.  In this case,
    set `doalpha` to True and the alpha and beta arrays are constructed 
    directly.  There is a trade-off since constructing the alpha and beta
    arrays directly takes longer for smaller surfaces.

    Examples
    --------
    later?

    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_polyfit2d.pro
        2022-06-17 - Simplfied the creation of alpha is doalpha=True
    """

    # Construct the uncertainty array if necessary
    
    if zunc is None:
        zunc = np.full(len(x), 1.0)

    # Get rid of NaNs 

    znonan = ~np.isnan(z)
    nnan = np.size(znonan)-np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    zz = z[znonan]
    zzunc = zunc[znonan]

    # Get basic parameters
    
    ndat = len(xx)       # number of idependent variables
    ncoeffs = (xorder+1)*(yorder+1)   # number of coefficients
    ndof = ndat-ncoeffs  # number of degrees of freedom

    # Determine the exponents of the basis functions

    xexp = np.arange(xorder+1, dtype=int)
    yexp = np.arange(yorder+1, dtype=int)
    exp = np.empty((ncoeffs, 2), dtype=int)
    idx = 0

    for i in range(yorder+1):

        for j in range(xorder+1):

            exp[idx, 0] = xexp[j]
            exp[idx, 1] = yexp[i]
            idx += 1

    # Create the b array

    b = zz/zzunc

    # Now start the linear algebra construction
    
    if doalpha is True:

        alpha = np.empty((ncoeffs, ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i, ncoeffs):

                at = (xx**exp[i, 0] * yy**exp[i, 1])/zzunc
                a = (xx**exp[j, 0] * yy**exp[j, 1])/zzunc

                val = np.sum(at*a)
                    
                alpha[j, i] = val
                alpha[i, j] = val
                beta[i] = np.sum(at*b)                    
                
    elif doalpha is False:

        at = np.empty((ncoeffs, ndat))
        for i in range(ncoeffs):

            at[i, :] = xx**exp[i, 0]*yy**exp[i, 1]/zzunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)
        
    # Solve things (need to remember what you are doing...)
    
    coeffs = np.linalg.solve(alpha, beta)

    if justfit is True:

        # Just return the coefficients
        
        return({"coeffs": coeffs, "var": None, "covar": None, "zfit": None,
                "nparm": None, "ndof": None, "chi2": None, "rchi2": None,
                "rms": None})

    elif justfit is False:

        # Keep going with other outputs and possible command line output
        
        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        zfit = poly2d(x, y, xorder, yorder, coeffs)

        residual = z-zfit
        rms = np.std(residual[znonan])

        chi2 = np.sum((residual[znonan]/zunc[znonan])**2)
        rchi2 = chi2/ndof
    
        # Report results if requested

        if silent is False:

            print(' ')
            print('             Number of points = ', len(x))
            print('          Number of NaNs in y = ', nnan)
            print('         Number of parameters = ', ncoeffs)
            print(' Number of degrees of freedom = ', ndof)
            print('                  Chi-Squared = ', chi2)
            print('          Reduced Chi-Squared = ', rchi2)
            print('         RMS deviation of fit = ', rms)
            print(' ')
            print('Coefficients:')
            print(' ')
            
            for i in range(0, ncoeffs):
            
                print('Coeff #', str(i).zfill(2), ': ', coeffs[i], '+-',
                      np.sqrt(var[i]), sep='')

            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')

        return({"coeffs": coeffs, "var": var, "covar": covar, "zfit": zfit,
                "nparm": ncoeffs, "ndof": ndof, "chi2": chi2, "rchi2": rchi2,
                "rms": rms})


def robustpolyfit1d(x, y, order, thresh, eps, yunc=None, goodbad=None,
                    justfit=False, silent=True):

    """
    Fits a "robust" polynomial of a given order to a set of 1-D data.

    Parameters
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

    Returns
    -------
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

    Notes
    -----
 

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
    """

    # Check to see if the user passed an uncertainty array.
    
    if yunc is None:
        yunc = np.full(len(x), 1.0)

    # Set up the ogoodbad array

    ogoodbad = np.full(len(x), 2)
        
    # Get rid of NaNs

    znotnan = ~np.isnan(y)
    init_goodcnt = np.sum(znotnan)
    nnan = np.sum(~znotnan)
    xx = x[znotnan]
    yy = y[znotnan]
    yyunc = yunc[znotnan]

    ogoodbad[znotnan] = 1

    # Do a first pass of the data

    fit = polyfit1d(xx, yy, order, yunc=yyunc, justfit=justfit, silent=True)
    
    # Compute the residuals and the mean and sigma of the residuals
    
    residual = yy-np.polynomial.polynomial.polyval(xx, fit['coeffs'])

    mean = np.mean(residual)
    stddev = np.std(residual)

    # Now check for residual points that are outside the threshhold
    
    good = np.abs((residual-mean)/stddev) <= thresh
    goodcnt = np.sum(good)

    if goodcnt != init_goodcnt:

        itter = 0
        for i in range(0, 11):

            itter += 1

            oldstddev = stddev
            fit = polyfit1d(xx[good], yy[good], order, yunc=yyunc[good],
                            justfit=justfit, silent=True)
            residual = yy[good] - \
                       np.polynomial.polynomial.polyval(xx[good],
                                                        fit['coeffs'])
            mean = np.mean(residual)
            stddev = np.std(residual)

            # Check to see if the new stddev isn't much of a change from the
            # old one

            if (oldstddev-stddev)/oldstddev < eps:
                break

            # Now generate a new full residual array

            residual = yy-np.polynomial.polynomial.polyval(xx, fit['coeffs'])

            good = np.abs((residual-mean)/stddev) <= thresh
            goodcnt = np.sum(good)

    # Let's reconstuct the goodbad array

    tmp = np.full(init_goodcnt, 0)
    tmp[good] = 1
    ogoodbad[znotnan] = tmp

    # Create a full yfit

    yfit = np.polynomial.polynomial.polyval(x, fit['coeffs'])

    # Report results if requested

    if silent is not True:

        nbad = np.sum(good == 0)

        print(' ')
        print('    Number of original points = ', len(x))
        print('              Sigma threshold = ', thresh)
        print('              Number of loops = ', itter)
        print('           Number of outliers = ', nbad)
        print('          Number of NaNs in y = ', nnan)
        print('         Number of parameters = ', order+1)
        print(' Number of degrees of freedom = ', fit["ndof"])
        print('                  Chi-Squared = ', fit["chi2"])
        print('          Reduced Chi-Squared = ', fit["rchi2"])
        print('         RMS deviation of fit = ', fit["rms"])
        print(' ')
        print('Coefficients:')
        print(' ')

        for i in range(0, order+1):
        
            print('Coeff #', str(i).zfill(2), ': ', fit['coeffs'][i], '+-',
                  np.sqrt(fit['var'][i]), sep='')

        print(' ')
        print('Covariance Matrix:')
        print(fit['covar'])
        print(' ')

    return({"coeffs": fit['coeffs'], "var": fit['var'], "covar": fit['covar'],
            "ogoodbad": ogoodbad, "yfit": yfit, "nparm": order+1,
            "ndof": fit['ndof'], "chi2": fit['chi2'], "rchi2": fit['rchi2'],
            "rms": fit['rms']})


def image_poly(img, coeffs):

    """
    evaluates a polynomial function of a variable for images

    Parameters
    ----------
    img : numpy.ndarray
        an [nrows,ncols]

    coeffs : numpy.ndarray
        an [ncoeffs,nrows,ncols] array of polynomial coefficients
        [0,:,:] is the c_0 array
        [1,:,:] is the c_1 array
        etc.

    Returns
    -------
    numpy.ndarray
        an [nrows,ncols] array the 

    Procedure
    ---------
    Follows the IDL poly.pro technique to evaluate a polynomial but does 
    it for each value in the array at once

    Examples
    --------
    > import numpy as np
    > img = np.array([[1,1,1],[2,2,2],[3,3,3]])
    > coeffs0 = np.array([[2,2,2],[2,2,2],[2,2,2]])
    > coeffs1 = np.array([[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]])
    > coeffs = np.stack((coeffs0,coeffs1))

    [[2.5 2.5 2.5]
     [3.  3.  3. ]
     [3.5 3.5 3.5]]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.  
                 Based on the mc_polyimg.pro IDL program.
    2022-06-02 - changed indexing from [*,*,i] to [i,*,*]
    2022-06-03 - updated doc string and fixed bug with n and the indexing of y

    """

    n = coeffs.shape[0]-1
    
    y = coeffs[n, :, :]

    for i in range(n-1, -1, -1):

        y = y*img + coeffs[i, :, :]
    
    return y
