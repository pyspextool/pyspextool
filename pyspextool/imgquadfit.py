import numpy as np

def imgquadfit(img,imgunc=None,doalpha=False):

    '''
    Fits a quadratic surface to an image

    Input Parameters
    ----------------
    img : array_like
        The image to to fit.

    imgunc : array_like, optional
        The uncertainty image.

    doalpha : {False, True}, optional
        If true, then A**T * A is computed directly instead of A.

    Returns
    --------
    numpy.ndarray
        The coefficients of the fit

    Notes
    -----
    Standard least squares fitting (e.g. Numerical Recipes).  
    The program fits the following function:

    z = c[0] + c[1]*x + c[2]*y + c[3]x^2 + c[4]y^2 + c[5]*x*y

    where x is the column number and y is the row number.

    By default, the design matrix A is constructed.  The alpha 
    (A^T ## A) and beta (A^T ## b) arrays of the normal equations 
    are then constructed and then passed to np.linalg.solve.  

    If the dataset to be fit is large and/or the number of coefficients 
    is large, the design matrix A can be large (ncoeff,ndat) and may 
    take too much memory and/or too much time to construct.  In this case,
    set `doalpha` to True and the alpha and beta arrays are constructed 
    directly.  There is a trade off since constructing the alpha and beta 
    arrays directly takes longer for smaller surfaces.

    Unlike polyfit1d and polyfit2d, this function ONLY computes the 
    coefficients.

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_quadfit.pro

    '''
    
# Do basic things and create basic things
    
    nrows, ncols = img.shape
    if imgunc is None: imgunc = np.full((nrows,ncols),1.0)

# Create coordinate images
    
    ximg = np.tile(np.arange(ncols,dtype=float),(nrows,1))
    yimg = np.tile(np.reshape(np.arange(nrows,dtype=float),\
                   (nrows,1)),(1,ncols))

# Ravel the images

    x = np.ravel(ximg)
    y = np.ravel(yimg)
    z = np.ravel(img)
    zunc = np.ravel(imgunc)
    
# Get rid of NaNs 

    znonan = ~np.isnan(z)
    nnan = np.size(znonan)-np.sum(znonan)
    xx = x[znonan]
    yy = y[znonan]
    zz = z[znonan]
    zzunc = zunc[znonan]    
    
# Set up to construct the equations

    ndat = len(xx)
    ncoeffs = 6
    exp = np.array([[0,0],[1,0],[0,1],[2,0],[0,2],[1,1]])

    b = zz/zzunc

# Now start the linear algebra construction
    
    if doalpha is True:

        alpha = np.empty((ncoeffs,ncoeffs))
        beta = np.empty(ncoeffs)

        for i in range(ncoeffs):

            for j in range(i,ncoeffs):

                at = (xx**exp[i,0] * yy**exp[i,1])/zzunc
                a  = (xx**exp[j,0] * yy**exp[j,1])/zzunc

                val = np.sum(at*a)
                    
                alpha[j,i] = val
                alpha[i,j] = val
                beta[i] = np.sum(at*b)                    
                
    elif doalpha is False:

        AT = np.empty((ncoeffs,ndat))
        for i in range(ncoeffs):

            AT[i,:] = xx**exp[i,0]*yy**exp[i,1]/zzunc

        alpha = np.matmul(AT,np.transpose(AT))
        beta = np.matmul(AT,b)
        
# Solve things 

    coeffs = np.linalg.solve(alpha,beta)

    return coeffs

    
    
