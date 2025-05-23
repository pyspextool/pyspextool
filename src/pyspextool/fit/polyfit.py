import numpy as np
import numpy.typing as npt
from numpy.polynomial.polynomial import polyval
import matplotlib.pyplot as pl

from pyspextool.io.check import check_parameter

def goodbad_init(*args,
                 goodbad:npt.ArrayLike=None):

    """
    To set up the goodbad array for a robust fit.

    Parameters
    ----------
    *args : ndarray
        A tuple of (ndat, ) arrays to be used in a robust fit.  Typically, the
        dependent variable array and its associated uncertainty.

    goodbad : ndarray, default None
         An (ndat, ) array of int values where 0=bad, 1=good, 2=NaN.

    Returns
    -------
    ndarray
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

    #
    # Check parameters
    #

    check_parameter('goodbad_init', 'args', args, 'tuple')

    check_parameter('goodbad_init', 'goodbad', goodbad, ['NoneType', 'ndarray'])

    #
    # Make sure the arrays are all the same size
    #

    


    
    #  Was a goodbad array passed?
    
    if goodbad is None:

        # Create a good bad array

        goodbad = np.full(len(args[0]), 1, dtype=int)

    else:

        # Convert to a numpy array if need be

        goodbad = np.array(goodbad, dtype=int)

    #
    # Find NaNs if need be
    #
    
    for arg in args:

        if arg is not None:

            znan = np.isnan(arg)
            if np.sum(znan) > 0:
                goodbad[znan] = 2

    return goodbad



def image_poly(img:npt.ArrayLike,
               coeffs:npt.ArrayLike):

    """
    Evaluates a polynomial function of a variable for images

    Parameters
    ----------
    img : ndarray
        An (nrows, ncols) image of "independent" values.

    coeffs : ndarray
        An (ncoeffs, nrows, ncols) array of polynomial coefficients
        [0,:,:] is the c_0 array
        [1,:,:] is the c_1 array
        etc.

    Returns
    -------
    ndarray
        An (rows,ncols) array where the each element is 

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

    n = coeffs.shape[0] - 1

    y = coeffs[n, :, :].copy()

    for i in range(n - 1, -1, -1):
        y = y * img + coeffs[i, :, :]

    return y



def make_alphabeta_1d(x:npt.ArrayLike,
                      y:npt.ArrayLike,
                      yunc:npt.ArrayLike,
                      order:int,
                      doalpha:bool=False):

    """
    Creates the alpha and beta arrays of the normal equations of
    the general linear least squares problem.

    Parameters
    ----------
    x : ndarray
        (ndat, ) array of independent values.

    y : ndarray
        (ndat, ) array of dependent values.

    yunc : ndarray
        (ndat,) array of uncertainties.

    order : int
        The polynomial order.

    doalpha : {False, True}
        Set to True to create the alpha directly.  See Notes.

    Returns
    -------
    alpha, beta : ndarray, ndarray

    Notes
    -----
    If `doalpha` is set, the alpha and beta arrays are constructed directly.
    If not, then the design matrix A and the vector b are constructed which
    are then used to compute the alpha and beta arrays.  There is a trade-off
    since constructing the alpha and beta arrays directly takes longer for
    smaller surfaces.

    """

    #
    # Check parameters
    #

    # How many coefficients

    ncoeffs = order + 1

    # Compute the exponents for the basis functions

    exp = np.arange(0, order + 1)

    # Create the b vector

    b = y / yunc

    if doalpha is True:

        # Create empty alpha and beta arrays and compute the b vector

        alpha = np.zeros((ncoeffs, ncoeffs))
        beta = np.zeros(ncoeffs)

        # Loop to fill them in

        for i in range(0, ncoeffs):

            for j in range(i, ncoeffs):
                at = (x ** exp[i]) / yunc
                a = (x ** exp[j]) / yunc

                alpha[i, j] = np.sum(at * a)
                alpha[j, i] = np.sum(at * a)
                beta[i] = np.sum(at * b)

    elif doalpha is False:

        at = np.empty((ncoeffs, len(x)))
        for i in range(ncoeffs):
            at[i, :] = x ** exp[i] / yunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)

    return alpha, beta



def make_alphabeta_2d(x:npt.ArrayLike,
                      y:npt.ArrayLike,
                      z:npt.ArrayLike,
                      zunc:npt.ArrayLike,
                      xorder:int,
                      yorder:int,
                      doalpha:bool=False):

    """
    Creates the alpha and beta arrays of the normal equations of
    the general linear least squares problem.

    Parameters
    ----------
    x : ndarray
        An (ndat, ) array of independent values.

    y : ndarray
        An (ndat, ) array of independent values.

    z : ndarray
        An (ndat ) array of dependent values.

    zunc : ndarray
        (ndat, ) array of uncertainties.

    xorder : int
        The polynomial order of the first independent variable.

    yorder : int
        The polynomial order of the second indepenent variable.

    doalpha : {False, True}, 
        Set to True to create the alpha directly.  See Notes.

    Returns
    -------
    ndarray, ndarray
    
    The alpha and beta arrays.  

    Notes
    -----
    If `doalpha` is set, the alpha and beta arrays are constructed directly.
    If not, then the design matrix A and the vector b are constructed which
    are then used to compute the alpha and beta arrays.  There is a trade-off
    since constructing the alpha and beta arrays directly takes longer for
    smaller surfaces.

    """

    #
    # Check parameters
    #

    
    check_parameter('make_alphabeta_2d', 'x', x, 'ndarray')

    check_parameter('make_alphabeta_2d', 'y', y, 'ndarray')

    check_parameter('make_alphabeta_2d', 'z', z, 'ndarray')

    check_parameter('make_alphabeta_2d', 'zunc', zunc, 'ndarray')

    check_parameter('make_alphabeta_2d', 'xorder', xorder, 'int')

    check_parameter('make_alphabeta_2d', 'yorder', yorder, 'int')    
    
    check_parameter('make_alphabeta_2d', 'doalpha', doalpha, 'bool')        
    
    # How many coefficients

    ncoeffs = (xorder + 1) * (yorder + 1)

    # Determine the exponents of the basis functions

    xexp = np.arange(xorder + 1, dtype=int)
    yexp = np.arange(yorder + 1, dtype=int)
    exp = np.empty((ncoeffs, 2), dtype=int)
    idx = 0

    for i in range(yorder + 1):

        for j in range(xorder + 1):
            exp[idx, 0] = xexp[j]
            exp[idx, 1] = yexp[i]
            idx += 1

    # Create the b vector

    b = z / zunc

    if doalpha is True:

        # Create empty alpha and beta arrays and then fill them in.

        alpha = np.empty((ncoeffs, ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i, ncoeffs):
                at = (x ** exp[i, 0] * y ** exp[i, 1]) / zunc
                a = (x ** exp[j, 0] * y ** exp[j, 1]) / zunc

                val = np.sum(at * a)

                alpha[j, i] = val
                alpha[i, j] = val
                beta[i] = np.sum(at * b)

    elif doalpha is False:

        at = np.empty((ncoeffs, len(x)))
        for i in range(ncoeffs):
            at[i, :] = x ** exp[i, 0] * y ** exp[i, 1] / zunc

        alpha = np.matmul(at, np.transpose(at))
        beta = np.matmul(at, b)

    return alpha, beta


def poly_1d(x:int | float | npt.ArrayLike,
            coeffs:npt.ArrayLike,
            covar:npt.ArrayLike=None,
            talk=False):
    
    """
    Evaluates a polynomial function of an independent variables
    

    Parameters
    ----------
    x : int, float, ndarray
        An (ndat, ) array of independent values

    coeffs : ndarray
        An (ncoeffs,) array of coefficients from polyfit_1d

    covar : ndarray, default None
        An (ncoeffs, ncoeffs) covariance array.

    Returns
    -------
    ndarray, optional ndarray 
        The polynomial evaluated at `x`.
        The variances evaluated at `x`.

    """

    #
    # Check the parameters
    #
    
    check_parameter('poly_1d', 'x', x,
                    ['int', 'float', 'float64', 'list', 'ndarray'])

    check_parameter('poly_1d', 'coeffs', coeffs, 'ndarray')

    check_parameter('poly_1d', 'covar', covar, ['ndarray', 'NoneType'])

    #
    # Get set up
    #

    x = x.astype(np.float64)
    ncoeffs = np.size(coeffs)
    ndat = np.size(x)
    if talk is True: print('first', ncoeffs)
    #
    # Compute the values at x
    #

    z = np.zeros(ndat, dtype=np.float64)
    for i in range(ncoeffs):
        z = z + coeffs[i] * x ** i

    # Compute the variances if requested

    if talk is True:  print('second',ncoeffs)
    if covar is not None:

        var = np.zeros(ndat, dtype=np.float64)
        for i in range(ncoeffs):

            for j in range(ncoeffs):
                var = var + covar[j, i] * x ** i * x ** j
                if talk is True:  print(i,j,covar[j, i])
#                print(covar[j, i] * x ** i * x ** j)
                
        return z, var

    else:

        return z


def poly_2d(x:int | float | npt.ArrayLike,
            y:int | float | npt.ArrayLike,
            xorder:int,
            yorder:int,
            coefficients:npt.ArrayLike):

    """
    Evaluates a polynomial function of two independent variables
    
    Parameters
    ----------
    x : ndarray
        An (ndat, ) array of independent values (int or float)

    y : ndarray
        An (ndat, ) array of independent values (int or float)

    xorder : int
        The order of the polynomial for the x dimension

    yorder : int
        The order of the polynomial for the y dimension

    coefficients : ndarray
        An (xorder+1 * yorder+1, ) array of coefficients from polyfit2d

    Returns
    -------
    ndarray
        The polynomial evaluated at `x` and `y`.


    """
    
    #
    # Check parameters
    #

    all = ['int', 'int8', 'float', 'float64', 'ndarray']
    check_parameter('poly_2d', 'x', x, all)

    check_parameter('poly_2d', 'y', y, all)

    check_parameter('poly_2d', 'xorder', xorder, 'int')

    check_parameter('poly_2d', 'yorder', yorder, 'int')            

    check_parameter('poly_2d', 'coefficients', coefficients, 'ndarray')    

    x = np.float64(x)
    y = np.float64(y)
    
    idx = 0
    z = 0
    
    for i in range(yorder + 1):

        for j in range(xorder + 1):

            z = z + x ** j * y ** i * coefficients[idx]
            idx += 1

    return z


def polyfit_1d(x:npt.ArrayLike,
               y:npt.ArrayLike,
               order:int,
               yunc:npt.ArrayLike=None,
               goodbad:npt.ArrayLike=None,
               robust:dict=None,
               doalpha:bool=False,
               justfit:bool=False,
               silent:bool=True):

    """
    Fits a (robust) polynomial to a set of 1D data.

    Parameters
    ----------
    x : ndarray
        An (ndat, ) array of independent values.

    y : ndarray
        An (ndat, ) array of dependent values.

    order : int
        The order of the polynomial fit.

    yunc : ndarray, default None
        An (ndat, ) array of uncertainties on the dependent values `y`

    goodbad : ndarray, deafult None
        An (ndat, ) array identifying each value as either, bad (0),
        good(1) or Nan (2). If passed, the the function will ignore "bad"
        values during the fit.

    robust : dict of {str:int or float, str:int or float}, optional

        `"thresh"` : int or float
            The standard deviation threshold over which to identify pixels as
            an outlier.

        `"eps"` : int or float
            The epsilon value to decide when to stop trying to identify
            outliers.  If (stddev_i-1 - stddev_i) / stddev_i < `"eps"` then
            the search is ended.

        If given, an attempt will be made to robustly determine the fit.  

    doalpha: {False, True}
        Set to False to compute the design matrix A, then alpha and beta.
        Set to True to compute the alpha and beta arrays directly.

    justfit : {False, True}
        Set to False to do everything.
        Set to True to only fit the coefficients.  Faster this way.
    
    silent : {True, False}
        Set to True to report the results to the command line.
        Set to False to not report the results to the command line.

    Returns
    -------
    dict
        `"coeffs"` : ndarray
            The polynomial coefficients

        `"coeff_var"`: ndarray
            The variances of the coefficients

        `"coeff_covar"` : ndarray
            The covariance matrix

        `"yfit"` : ndarray 
            The polynomial evaluated at `x`

        `"yfit_var"` : ndarray 
            The variances of `yfit` calculated using `coeff_covar`.

        `"nparm"` : int
            The number of parameters of the fit

        `"ndof"` : int
            The number of degrees of freedom

        `"chi2"` : float64
            The chi^2 value of the fit

        `"rchi2"`: float64
            The reduced chi^2 value of the fit

        `"rms"`: float64
            The rms of the fit (1/N, not 1/(N-1))

    Notes
    -----
    Standard least squares fitting (e.g. Numerical Recipes).  
    The program fits the following function:

    y = f(x) = c_0 + c_1*x + c_2*x**2 + c_3*x**3 + ...

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

    #
    # Check parameters
    #

    check_parameter('polyfit_1d','x', x, 'ndarray')

    check_parameter('polyfit_1d','y', y, 'ndarray')

    check_parameter('polyfit_1d','order', order, 'int')

    check_parameter('polyfit_1d','yunc', yunc, ['ndarray','NoneType'])

    check_parameter('polyfit_1d','goodbad', goodbad, ['ndarray','NoneType'])

    check_parameter('polyfit_1d','robust', robust, ['dict','NoneType'])

    check_parameter('polyfit_1d','doalpha', doalpha, 'bool')

    check_parameter('polyfit_1d','silent', silent, 'bool')    

    check_parameter('polyfit_1d','justfit', justfit, 'bool')                


    # Construct the uncertainty array if need be

    if yunc is None:
        yunc = np.full(len(x), 1.0, dtype=np.float64)

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

        residual = yy - polyval(xx, coeffs)

        mean = np.mean(residual)
        stddev = np.std(residual)

        # Now check for residual points that are outside the threshhold

#        if silent is False:
#            
#            pl.figure()
#            pl.plot(np.abs((residual - mean) / stddev),'or')
#            pl.axhline(y=robust['thresh'])
#            pl.title('test')
#            pl.show()
        
        z_good = np.abs((residual - mean) / stddev) <= robust['thresh']
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

                residual = yy[z_good] - polyval(xx[z_good], coeffs)

                mean = np.mean(residual)
                stddev = np.std(residual)
                
                if stddev == 0.0:
                    break

                # Check to see if the new stddev isn't much of a change
                # from the old one

                if (old_stddev - stddev) / old_stddev < robust['eps']:
                    break

                    # Create full residual array and search for bad ones

                residual = yy - polyval(xx, coeffs)
                z_good = np.abs((residual - mean) / stddev) <= robust['thresh']

        # Let's reconstuct the goodbad array

        tmp = np.zeros(len(xx))
        tmp[z_good] = 1
        goodbad[z_initgood] = tmp

    # Start reporting and returning results

    if justfit is True:

        return ({"coeffs": coeffs, "coeff_var": None, "coeff_covar": None,
                 "yfit": None, "yfit_var": None, "goodbad": goodbad,
                 "nparm": None, "ndof": None, "chi2": None, "rchi2": None,
                 "rms": None})

    else:

        # Keep going with other outputs and possible command line output

        # Compute the covariance array and get the variances of the parameters

        covar = np.linalg.inv(alpha)
        coeff_var = np.diagonal(covar)

        # Get the rms of the fit and the chi2 value

        yfit, yfit_var = poly_1d(x, coeffs, covar=covar)
        residual = y - yfit
        z_good = goodbad == 1
        n_good = np.sum(z_good)
        rms = np.std(residual[z_good])
        chi2 = np.sum((residual[z_good] / yunc[z_good]) ** 2)

        # Compute the reduced chi2

        ncoeffs = order + 1
        ndof = n_good - ncoeffs  # number of degrees of freedom
        rchi2 = chi2 / ndof

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

            for i in range(0, order + 1):
                print('Coeff #', str(i).zfill(2), ': ', coeffs[i], '+-',
                      np.sqrt(coeff_var[i]), sep='')

            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')

        # Return the results

    return ({"coeffs": coeffs, "coeffs_var": coeff_var, "coeffs_covar": covar,
             "yfit": yfit, "yfit_var": yfit_var, "goodbad": goodbad,
             "nparm": order + 1, "ndof": ndof, "chi2": chi2, "rchi2": rchi2,
             "rms": rms})


#
# ==============================================================================
#
def polyfit_2d(x:npt.ArrayLike,
               y:npt.ArrayLike,
               z:npt.ArrayLike,
               xorder:int,
               yorder:int,
               zunc:npt.ArrayLike=None,
               goodbad:npt.ArrayLike=None,
               robust:dict=None,
               doalpha:bool=False,
               silent:bool=True,
               justfit:bool=False):
    
    """
    Fits a 2D polynomial to a set of 2D data

    
    Parameters
    ----------
    x : ndarray
        An (ndat, ) array of independent values

    y : ndarray
        An (ndat, ) array of independent values

    z : ndarray
        An (ndat, ) array of dependent values

    xorder : int
        The order of the polynomial for the x dimension.  1=linear

    yorder : int
        The order of the polynomial for the y dimension.  1=linear

    zunc : nndarray, default None
        An (ndat, ) array of uncertainties on the dependent values

    goodbad : ndarray, deafult None
        An (ndat, ) array identifying each value as either, bad (0),
        good(1) or Nan (2). If passed, the the function will ignore "bad"
        values during the fit.

    robust : dict of {str:int or float, str:int or float}, optional

        `"thresh"` : int or float
            The standard deviation threshold over which to identify pixels as
            an outlier.

        `"eps"` : int or float
            The epsilon value to decide when to stop trying to identify
            outliers.  If (stddev_i-1 - stddev_i) / stddev_i < `"eps"` then
            the search is ended.

        If given, an attempt will be made to robustly determine the fit.  

    doalpha: {False, True}
        Set to False to compute the design matrix A, then alpha and beta.
        Set to True to compute the alpha and beta arrays directly.

    justfit : {False, True}
        Set to False to do everything.
        Set to True to only fit the coefficients.  Faster this way.
    
    silent : {True, False}
        Set to True to report the results to the command line.
        Set to False to not report the results to the command line.

    Returns
    -------
    dict

        `"coeffs"` : ndarray 
            The polynomial coefficients

        `"var"` : ndarray 
            The variances of the coefficients

        `"covar"` : ndarray
            The covariance matrix

        `"zfit"` : ndarray 
            The polynomial evaluated at `x` and `y`

        `"nparm"` : int
            The number of parameters of the fit

        `"ndof"` : int
            The number of degrees of freedom

        `"chi2"` : float64
            The chi^2 value of the fit

        `"rchi2"` : float64
            The reduced chi^2 value of the fit

        `"rms"` : float64
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

    """
    
    #
    # Check parameters
    #

    check_parameter('polyfit_2d','x', x, 'ndarray')

    check_parameter('polyfit_2d','y', y, 'ndarray')

    check_parameter('polyfit_2d','z', z, 'ndarray')

    check_parameter('polyfit_2d','xorder', xorder, 'int')

    check_parameter('polyfit_2d','yorder', yorder, 'int')

    check_parameter('polyfit_2d','zunc', zunc, ['ndarray','NoneType'])

    check_parameter('polyfit_2d','goodbad', goodbad, ['ndarray','NoneType'])

    check_parameter('polyfit_2d','robust', robust, 'dict')

    check_parameter('polyfit_2d','doalpha', doalpha, 'bool')

    check_parameter('polyfit_2d','silent', silent, 'bool')    

    check_parameter('polyfit_2d','justfit', justfit, 'bool')                

    #
    # Get set up
    # 
    
    # Convert to float64

    x = np.float64(x)
    y = np.float64(y)    
    z = np.float64(z)
    
    # Construct the uncertainty array if necessary

    if zunc is None:
        zunc = np.full(len(x), 1.0,dtype=np.float64)

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

    #
    # Now get the fit started
    #
    
    # Create the alpha and beta arrays of the normal equatioon

    alpha, beta = make_alphabeta_2d(xx,
                                    yy,
                                    zz,
                                    zzunc,
                                    xorder,
                                    yorder,
                                    doalpha=doalpha)

    # Solve things

    coeffs = np.linalg.solve(alpha, beta)

    # Are we doing a robust fit?  

    if robust is not None:

        # Compute the residuals and the mean and sigma of the residuals

        residual = zz - poly_2d(xx, yy, xorder, yorder, coeffs)

        mean = np.mean(residual)
        stddev = np.std(residual)

        # Now check for residual points that are outside the threshhold

        z_good = np.abs((residual - mean) / stddev) <= robust['thresh']
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

                if (old_stddev - stddev) / old_stddev < robust['eps']:
                    break

                    # Create full residual array and search for bad ones

                residual = zz - poly_2d(xx, yy, xorder, yorder, coeffs)
                z_good = np.abs((residual - mean) / stddev) <= robust['thresh']

        # Let's reconstuct the goodbad array

        tmp = np.zeros(len(xx))
        tmp[z_good] = 1
        goodbad[z_initgood] = tmp

    # Start reporting and returning results

    if justfit is True:

        return ({"coeffs": coeffs, "coeffs_var": None, "coeffs_covar": None,
                 "zfit": None, "zfit_var": None, "goodbad": goodbad,
                 "nparm": None, "ndof": None, "chi2": None, "rchi2":
                     None, "rms": None})

    else:

        # Keep going with other outputs and possible command line output

        # Compute the covariance array and get the variances of the parameters

        covar = np.linalg.inv(alpha)
        var = np.diagonal(covar)

        # Get the rms of the fit and the chi2 value

        zfit = poly_2d(x, y, xorder, yorder, coeffs)
        residual = z - zfit
        z_good = goodbad == 1
        n_good = np.sum(z_good)
        rms = np.std(residual[z_good])
        chi2 = np.sum((residual[z_good] / zunc[z_good]) ** 2)

        # Compute the reduced chi2

        ncoeffs = (xorder + 1) * (yorder + 1)
        ndof = n_good - ncoeffs  # number of degrees of freedom
        rchi2 = chi2 / ndof

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

            for i in range(0, ncoeffs):
                print('Coeff #', str(i).zfill(2), ': ', coeffs[i], '+-',
                      np.sqrt(var[i]), sep='')

            print(' ')
            print('Covariance Matrix:')
            print(covar)
            print(' ')

        # Return the results

    return ({"coeffs": coeffs, "coeffs_var": var, "coeffs_covar": covar,
             "zfit": zfit, "zfit_var": None, "goodbad": goodbad,
             "nparm": ncoeffs, "ndof": ndof, "chi2": chi2, "rchi2": rchi2,
             "rms": rms})
