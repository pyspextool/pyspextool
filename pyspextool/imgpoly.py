def imgpoly(img,coeffs):

    """
    evaluates a polynomial function of a variable for images

    Input Parameters
    ----------------
    img : numpy.ndarray
        an [nrows,ncols]

    coeffs : numpy.ndarray
        an [ncoeffs,nrows,ncols] array of polynomial coefficients
        [0,:,:] is the c_0 array
        [1,:,:] is the c_1 array
        etc.

    Returns
    --------
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
    
    y = coeffs[n,:,:]

    for i in range(n-1,-1,-1):

        y = y*img + coeffs[i,:,:]
    
    return(y)
