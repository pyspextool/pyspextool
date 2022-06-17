import numpy as np
from bicucof import bicucof

def bicuval(z,dz1,dz2,dz12,xl,xu,yl,yu,x,y):

    '''
    Evaluates a bicubic interpolation

    Input Parameters
    ----------------
    z : array_like
        The grid points values starting at the lower left and moving 
        counterclockwise.

    dz1 : float
        The gradient in dimension 1 evaluated at z.

    dz2 : float
        The gradient in dimension 2 evaluated at z.

    dz12 : float
        The cross derivative evaluated at z.

    xl : int
        The lower coordinate of the grid in the x direction.

    xu : int
        The upper coordinate of the grid in the x direction.

    yl : int
        The lower coordinate of the grid in the y direction.

    yu : int
        The upper coordinate of the grid in the y direction.

    x : array_like
        The x coordinates of the desired interpolation point(s)

    y : array_like
        The y values of the desired interpolation point(s)

    Returns
    --------
    numpy.ndarray
        The bicubic interpolation values at x and y.

    Notes
    ---------
    See Numerical Recipes bcucof

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_bicuval.pro.

    '''

    c = bicucof(z,dz1,dz2,dz12,xu-xl,yu-yl)
    
    t = (x-xl)/(xu-xl)
    u = (y-yl)/(yu-yl)

    nz = np.empty(x.shape)

    for i in range(3,-1,-1):

        nz = t*nz + c[i,0] + u*(c[i,0] + u*(c[i,0] + u*c[i,3]))
    
    return nz
    


    
    

