import numpy as np

def bicucof(z,dz1,dz2,dz12,nd1,nd2):

    '''
    Returns the coefficients used in bicuvals.py.

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

    nd1 : int
        Length of the grid cell in dimension 1.

    nd2 : int
        Length of the grid cell in dimension 2.

    Returns
    --------
    numpy.ndarray
        Returns a (4,4) set of coefficients used by bicuval.py.

    Notes
    ---------
    See bicucof in Numerical Recipes

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_bicucoeffs.pro.

    '''
    
    wt = np.array([[1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],\
        [0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.],\
        [-3.,0.,0.,3.,0.,0.,0.,0.,-2.,0.,0.,-1.,0.,0.,0.,0.],\
        [2.,0.,0.,-2.,0.,0.,0.,0.,1.,0.,0.,1.,0.,0.,0.,0.],\
        [0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],\
        [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.],\
        [0.,0.,0.,0.,-3.,0.,0.,3.,0.,0.,0.,0.,-2.,0.,0.,-1.],\
        [0.,0.,0.,0.,2.,0.,0.,-2.,0.,0.,0.,0.,1.,0.,0.,1.],\
        [-3.,3.,0.,0.,-2.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],\
        [0.,0.,0.,0.,0.,0.,0.,0.,-3.,3.,0.,0.,-2.,-1.,0.,0.],\
        [9.,-9.,9.,-9.,6.,3.,-3.,-6.,6.,-6.,-3.,3.,4.,2.,1.,2.],\
        [-6.,6.,-6.,6.,-4.,-2.,2.,4.,-3.,3.,3.,-3.,-2.,-1.,-1.,-2.],\
        [2.,-2.,0.,0.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],\
        [0.,0.,0.,0.,0.,0.,0.,0.,2.,-2.,0.,0.,1.,1.,0.,0.],
        [-6.,6.,-6.,6.,-3.,-3.,3.,3.,-4.,4.,2.,-2.,-2.,-2.,-1.,-1.],\
        [4.,-4.,4.,-4.,2.,2.,-2.,-2.,2.,-2.,-2.,2.,1.,1.,1.,1.]])

    x = np.concatenate((z,dz1*nd1,dz2*nd2,dz12*nd1*nd2))

    return np.reshape(np.matmul(wt,x),(4,4))
