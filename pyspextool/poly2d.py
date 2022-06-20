import numpy as np

def poly2d(x,y,xorder,yorder,coeffs):

    '''
    evaluates a polynomial function of two independent variables
    

    Input Parameters
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
    '''    

# Get set up with a bunch of necessary numbers and arrays
    
    ndat = x.size
    ncoeffs = (xorder+1)*(yorder+1)

    xexp = np.arange(xorder+1,dtype=int)
    yexp = np.arange(yorder+1,dtype=int)    

    z = np.empty(ndat)

# Now fill in the data
    
    for i in range(ndat):

        z[i] = np.sum(np.ravel(np.outer(y[i]**yexp,x[i]**xexp))*coeffs)

    return z
