import numpy as np
from numpy.polynomial.polynomial import polyval

from pyspextool.io.check import check_parameter


def goodbad_init(*args, goodbad=None):

    """
    To set up the goodbad array for a robust fit.

    Parameters
    ----------
    *args : array_like, optional
        (ndat,) arrays to be used in a robust fit.  Typically, the
        dependent variable array and its associated uncertainty.

    goodbad : array_like, default=None, optional
         An (ndat,) array of int values where 0=bad, 1=good, 2=NaN.

    Returns
    -------
    numpy.ndarray
         A goodbad array with properties consistent with `goodbad` and
         `*args`.

    Notes
    -----
    The function identifies NaN values in *args and sets those pixels
    in either a user-passed goodbad array or a function-created goodbad 
    array to 2.

    Examples
    --------

    > x= [1,2,np.nan,5,6]
    > y = [0,1,2,3,np.nan]

    """

    if goodbad is None:

        # Create a good bad array

        goodbad = np.full(len(args[0]), 1, dtype=int)

    else:

        # Convert to a numpy array if need be

        goodbad = np.array(goodbad, dtype=int)

    # Find NaNs if need be

    for arg in args:

        if arg is not None:

            znan = np.isnan(arg)
            if np.sum(znan) > 0:
                goodbad[znan] = 2

    return goodbad

#
#==============================================================================
#
def image_poly(img, coeffs):

    """
    Evaluates a polynomial function of a variable for images

    Parameters
    ----------
    img : numpy.ndarray
        An (nrows, ncols) image of "independent" values.

    coeffs : numpy.ndarray
        An (ncoeffs, nrows, ncols) array of polynomial coefficients
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

    """

    n = coeffs.shape[0]-1
    
    y = coeffs[n, :, :]

    for i in range(n-1, -1, -1):

        y = y*img + coeffs[i, :, :]
    
    return y

#
#==============================================================================
#
def make_alphabeta_1d(x, y, yunc, order, doalpha=False):

    """
    Creates the alpha and beta arrays of the normal equations of
    the general linear least squares problem.

    Parameters
    ----------
    x : array_like
        (ndat,) array of independent values.

    y : array_like
        (ndat,) array of dependent values.

    yunc : array_like
        (ndat,) array of uncertainties.

    order : int
        The polynomial order.

    doalpha : {False, True}, optional
        Set to True to create the alpha directly.  See Notes.

    Returns
    -------
    alpha, beta : numpy.ndarray, numpy.ndarray

    Notes
    -----
    If `doalpha` is set, the alpha and beta arrays are constructed directly.
    If not, then the design matrix A and the vector b are constructed which
    are then used to compute the alpha and beta arrays.  There is a trade-off
    since constructing the alpha and beta arrays directly takes longer for
    smaller surfaces.

    Examples
    --------
    later

    """

    # How many coefficients
    
    ncoeffs = order+1

    # Compute the exponents for the basis functions
    
    exp = np.arange(0, order+1)  

    # Create the b vector

    b = y/yunc
    
    if doalpha is True:

        # Create empty alpha and beta arrays and compute the b vector
    
        alpha = np.zeros((ncoeffs, ncoeffs))
        beta = np.zeros(ncoeffs)

        # Loop to fill them in
        
        for i in range(0, ncoeffs):

            for j in range(i, ncoeffs):
        
                at = (x**exp[i])/yunc
                a = (x**exp[j])/yunc

                alpha[i, j] = np.sum(at*a)
                alpha[j, i] = np.sum(at*a)
                beta[i] = np.sum(at*b)

    elif doalpha is False:

        at = np.empty((ncoeffs, len(x)))
        for i in range(ncoeffs):

            at[i, :] = x**exp[i]/yunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)

    return alpha, beta

#
#==============================================================================
#
def make_alphabeta_2d(x, y, z, zunc, xorder, yorder, doalpha=False):

    """
    Creates the alpha and beta arrays of the normal equations of
    the general linear least squares problem.

    Parameters
    ----------
    x : array_like
        (ndat,) array of independent values.

    y : array_like
        (ndat,) array of independent values.

    z : array_like
        (ndat,) array of dependent values.

    zunc : array_like
        (ndat,) array of uncertainties.

    xorder : int
        The polynomial order of the first independent variable.

    yorder : int
        The polynomial order of the second indepenent variable.

    doalpha : {False, True}, optional
        Set to True to create the alpha directly.  See Notes.

    Returns
    -------
    alpha, beta : numpy.ndarray, numpy.ndarray

    Notes
    -----
    If `doalpha` is set, the alpha and beta arrays are constructed directly.
    If not, then the design matrix A and the vector b are constructed which
    are then used to compute the alpha and beta arrays.  There is a trade-off
    since constructing the alpha and beta arrays directly takes longer for
    smaller surfaces.

    Examples
    --------
    later

    """

    # How many coefficients
    
    ncoeffs = (xorder+1)*(yorder+1)

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

    # Create the b vector

    b = z/zunc
    
    if doalpha is True:

        # Create empty alpha and beta arrays and then fill them in.
        
        alpha = np.empty((ncoeffs, ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i, ncoeffs):

                at = (x**exp[i, 0] * y**exp[i, 1])/zunc
                a = (x**exp[j, 0] * y**exp[j, 1])/zunc

                val = np.sum(at*a)
                    
                alpha[j, i] = val
                alpha[i, j] = val
                beta[i] = np.sum(at*b)                    
                
    elif doalpha is False:

        at = np.empty((ncoeffs, len(x)))
        for i in range(ncoeffs):

            at[i, :] = x**exp[i, 0]*y**exp[i, 1]/zunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)

    return alpha, beta
#
#==============================================================================
#
def poly_1d(x, coeffs, covar=None):

    """
    Evaluates a polynomial function of an independent variables
    

    Parameters
    ----------
    x : array_like, int or float
        An (ndat,) array of independent values

    coeffs : np.ndarray
        An (ndat,) array of coefficients from poly_fit_1d

    covar : np.ndarray, default=None, optional
        An (ncoefs, ncoeffs) covariance array.

    Output Parameters
    -----------------
    numpy.ndarray, optional numpy.ndarray 
        The polynomial evaluated at `x`.
        The variances evaluated at `x`.

    Notes
    -----

    Examples
    --------
    later?

    """

    #
    # Check the parameters
    #
    check_parameter('poly_1d', 'x', x, ['list','ndarray'])
    
    check_parameter('poly_1d', 'coeffs', coeffs, 'ndarray')

    check_parameter('poly_1d', 'covar', covar, ['ndarray', 'NoneType'])

    #
    # Get set up
    #
    ncoeffs = len(coeffs)
    ndat = len(x)

    #
    # Compute the values at x
    #

    z = np.zeros(ndat, dtype=float)    
    for i in range(ncoeffs):

        z = z+coeffs[i]*x**i

    # Compute the variances if requested

    if covar is not None:

        var = np.zeros(ndat, dtype=float)
        for i in range(ncoeffs):

            for j in range(ncoeffs):

                var = var + covar[j,i]*x**i * x**j
                
        return (z, var)

    else:

        return z

#
#==============================================================================
#
def poly_2d(x, y, xorder, yorder, coeffs):

    """
    evaluates a polynomial function of two independent variables
    

    Parameters
    ----------
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

    """
    idx = 0
    z = 0
    for i in range(yorder+1):

        for j in range(xorder+1):

            z = z+x**j*y**i*coeffs[idx]
            idx +=1

    return z

#
#==============================================================================
#
def poly_fit_1d(x, y, order, yunc=None, goodbad=None, robust=None,
                doalpha=False, justfit=False, silent=True):

    """
    Fits a (robust) polynomial of a given order to a set of 1-D data.

    Parameters
    ----------------
    x : array_like
        (ndat,) array of independent values.

    y : array_like
        (ndat,) array of dependent values.

    order : int
        The order of the polynomial

    yunc : array_like, optional
        (ndat,) array of uncertainties on the dependent values `y`

    goodbad : array_like
        A goodbad array identifying each value.  0=bad, 1=good, 2=NaN.

    robust : dict, optional
        

    doalpha: {False, True}, optional
        If true, then A**T * A is computed directly instead of A

    justfit : {False, True}, optional 
        Set to compute only the coefficients.  Faster this way.

    silent : {True, False}, optional
        If False, the result of the fit will be written to the command line

    Returns
    -------
    dict
        A dict where with the following entries:

        coeffs : numpy.ndarray
            The polynomial coefficients

        coeff_var : numpy.ndarray
            The variances of the coefficients

        coeff_covar : numpy.ndarray
            The covariance matrix

        yfit : numpy.ndarray 
            The polynomial evaluated at `x`

        yfit_var : numpy.ndarray 
            The variances of `yfit` calculated using `coeff_covar`.

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
    > result = poly_fit1d(x,y,1,yunc=yerr)

    Data from Chapter 6 of Data Reduction and Error Analysis for the 
    Physical Sciences by Bevington & Robinson

    """

    # Construct the uncertainty array if need be
    
    if yunc is None:
        yunc = np.full(len(x), 1.0)

    # Create the goodbad array

    goodbad = goodbad_init(y, yunc, goodbad=goodbad)

    # Find the NaNs

    znan = goodbad == 2
    nnan = np.sum(znan)

    # Find the intial good points

    z_initgood = goodbad == 1
    n_initgood = np.sum(z_initgood)

    # Clip to the initial good data

    xx = x[z_initgood]
    yy = y[z_initgood]
    yyunc = yunc[z_initgood]

    # Create the alpha and beta arrays of the normal equatioon
    
    alpha, beta = make_alphabeta_1d(xx, yy, yyunc, order, doalpha=doalpha)

    # Solve things

    coeffs = np.linalg.solve(alpha, beta)

    # Are we doing a robust fit?  
    
    if robust is not None:

        # Compute the residuals and the mean and sigma of the residuals
    
        residual = yy-polyval(xx, coeffs)

        mean = np.mean(residual)
        stddev = np.std(residual)

        # Now check for residual points that are outside the threshhold
    
        z_good = np.abs((residual-mean)/stddev) <= robust['thresh']
        n_bad = np.sum(~z_good)

        if n_bad != 0:

            # Found a bad data point, so start looping

            itter = 0
            for i in range(11):

                itter += 1
                old_stddev = stddev
                alpha, beta = make_alphabeta_1d(xx[z_good], yy[z_good],
                                                yyunc[z_good], order,
                                                doalpha=doalpha)                

                coeffs = np.linalg.solve(alpha, beta)
                
                residual = yy[z_good]-polyval(xx[z_good], coeffs)

                mean = np.mean(residual)
                stddev = np.std(residual)

                # Check to see if the new stddev isn't much of a change
                # from the old one

                if (old_stddev-stddev)/old_stddev < robust['eps']:
                    break                

                # Create full residual array and search for bad ones
                
                residual = yy - polyval(xx, coeffs)                
                z_good = np.abs((residual-mean)/stddev) <= robust['thresh']

        # Let's reconstuct the goodbad array

        tmp = np.zeros(len(xx))
        tmp[z_good] = 1
        goodbad[z_initgood] = tmp

# Start reporting and returning results
    
    if justfit is True:

        return({"coeffs": coeffs, "coeff_var": None, "coeff_covar": None,
                "yfit": None, "yfit_var":None, "goodbad": goodbad,
                "nparm": None, "ndof": None, "chi2": None, "rchi2": None,
                "rms": None})
        
    else:

        # Keep going with other outputs and possible command line output

        # Compute the covariance array and get the variances of the parameters
        
        covar = np.linalg.inv(alpha)
        coeff_var = np.diagonal(covar)

        # Get the rms of the fit and the chi2 value

        yfit, yfit_var = poly_1d(x, coeffs, covar=covar)
        residual = y-yfit
        z_good = goodbad == 1
        n_good = np.sum(z_good)
        rms = np.std(residual[z_good])
        chi2 = np.sum((residual[z_good]/yunc[z_good])**2)

        # Compute the reduced chi2
        
        ncoeffs = order+1
        ndof = n_good-ncoeffs  # number of degrees of freedom
        rchi2 = chi2/ndof
    
        # Report results if requested

        if silent is False:
            
            print(' ')
            print('                  Number of points = ', len(x))
            print('               Number of NaNs in y = ', nnan)
            print('     Number of initial good points = ', n_initgood)       
            print('Total number of bad/outlier points = ', np.sum(goodbad == 0))
            print('              Number of parameters = ', ncoeffs)
            print('      Number of degrees of freedom = ', ndof)
            print('                       Chi-Squared = ', chi2)
            print('               Reduced Chi-Squared = ', rchi2)
            print('              RMS deviation of fit = ', rms)
            print(' ')
            print('Coefficients:')
            print(' ')
            
            for i in range(0, order+1):
            
                print('Coeff #', str(i).zfill(2), ': ', coeffs[i], '+-',
                      np.sqrt(coeff_var[i]), sep='')
                    
            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')

        # Return the results

    return({"coeffs": coeffs, "coeff_var": coeff_var, "coeff_covar": covar,
            "yfit": yfit, "yfit_var":yfit_var, "goodbad": goodbad,
            "nparm": order+1, "ndof": ndof, "chi2": chi2, "rchi2": rchi2,
            "rms": rms})

#
#==============================================================================
#
def poly_fit_2d(x, y, z, xorder, yorder, zunc=None, goodbad=None,
                robust=None, doalpha=False, silent=True, justfit=False):

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

    # Create the goodbad array

    goodbad = goodbad_init(z, zunc, goodbad=goodbad)

    # Find the NaNs

    znan = goodbad == 2
    nnan = np.sum(znan)

    # Find the intial good points

    z_initgood = goodbad == 1
    n_initgood = np.sum(z_initgood)        

    # Clip to the initial good data

    xx = x[z_initgood]
    yy = y[z_initgood]
    zz = z[z_initgood]    
    zzunc = zunc[z_initgood]

    # Create the alpha and beta arrays of the normal equatioon
    
    alpha, beta = make_alphabeta_2d(xx, yy, zz, zzunc, xorder, yorder,
                                    doalpha=doalpha)
    
    # Solve things
    
    coeffs = np.linalg.solve(alpha, beta)

    # Are we doing a robust fit?  
    
    if robust is not None:

        # Compute the residuals and the mean and sigma of the residuals

        residual = zz-poly_2d(xx, yy, xorder, yorder, coeffs)

        mean = np.mean(residual)
        stddev = np.std(residual)

        # Now check for residual points that are outside the threshhold
    
        z_good = np.abs((residual-mean)/stddev) <= robust['thresh']
        n_bad = np.sum(~z_good)

        if n_bad != 0:

            # Found a bad data point, so start looping

            itter = 0
            for i in range(11):

                itter += 1
                old_stddev = stddev
                alpha, beta = make_alphabeta_2d(xx[z_good], yy[z_good],
                                                zz[z_good], zzunc[z_good],
                                                xorder, yorder, 
                                                doalpha=doalpha)                

                coeffs = np.linalg.solve(alpha, beta)
                
                residual = zz[z_good] - poly_2d(xx[z_good], yy[z_good],
                                                xorder, yorder, coeffs)

                mean = np.mean(residual)
                stddev = np.std(residual)

                # Check to see if the new stddev isn't much of a change
                # from the old one

                if (old_stddev-stddev)/old_stddev < robust['eps']:
                    break                

                # Create full residual array and search for bad ones
                
                residual = zz - poly_2d(xx, yy, xorder, yorder, coeffs)    
                z_good = np.abs((residual-mean)/stddev) <= robust['thresh']

        # Let's reconstuct the goodbad array

        tmp = np.zeros(len(xx))
        tmp[z_good] = 1
        goodbad[z_initgood] = tmp

# Start reporting and returning results
    
    if justfit is True:

        return({"coeffs": coeffs, "var": None, "covar": None, "zfit": None,
                "goodbad": goodbad, "nparm": None, "ndof": None,
                "chi2": None, "rchi2": None, "rms": None})
        
    else:

        # Keep going with other outputs and possible command line output

        # Compute the covariance array and get the variances of the parameters
        
        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        # Get the rms of the fit and the chi2 value

        zfit = poly_2d(x, y, xorder, yorder, coeffs)
        residual = z-zfit
        z_good = goodbad == 1
        n_good = np.sum(z_good)
        rms = np.std(residual[z_good])
        chi2 = np.sum((residual[z_good]/zunc[z_good])**2)

        # Compute the reduced chi2
        
        ncoeffs = (xorder+1)*(yorder+1)
        ndof = n_good-ncoeffs  # number of degrees of freedom
        rchi2 = chi2/ndof
    
        # Report results if requested

        if silent is False:
            
            print(' ')
            print('                  Number of points = ', len(x))
            print('               Number of NaNs in y = ', nnan)
            print('     Number of initial good points = ', n_initgood)       
            print('Total number of bad/outlier points = ', np.sum(goodbad == 0))
            print('              Number of parameters = ', ncoeffs)
            print('      Number of degrees of freedom = ', ndof)
            print('                       Chi-Squared = ', chi2)
            print('               Reduced Chi-Squared = ', rchi2)
            print('              RMS deviation of fit = ', rms)
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

        # Return the results

    return({"coeffs": coeffs, "var": var, "covar": covar, "zfit": zfit,
            "goodbad": goodbad, "nparm": ncoeffs, "ndof": ndof, "chi2": chi2,
            "rchi2": rchi2, "rms": rms})        
        



    
