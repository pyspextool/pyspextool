import numpy as np

def findidx(xarr,x):

    '''
    Finds the effective index of a function value in an ordered array.

    Input Parameters
    ----------------
    xarr : array_like of float
        (N,) The array to be searched, must be monotonically increasing or 
        decreasing.

    x : array_like of float or float
        (M,) The value or values whose indices are required.

    Returns
    -------
    numpy.ndarray of float or float
        Effective indices in `xarr` for `x`.

    Notes
    -----
    Near line-by-line copy of IDL Astronomy User's Library routine tabinv.pro. 

    From its documentation:  A binary search is used to find the values XARR(I) 
    and XARR(I+1) where XARR(I) < X < XARR(I+1).  IEFF is then computed using 
    linear interpolation between I and I+1.  IEFF = I + (X-XARR(I)) / 
    (XARR(I+1)-XARR(I)).

    IEFF values outside of `xarr` are set to NaN.  

    Examples
    --------
    > x = [1,2.5,3,5.5]
    > findidx(x,[1.5,3.1,6])
      [0.33333333 2.04              nan]

    Modification History
    --------------------
    2022-06-30 - Written by M. Cushing, University of Toledo.

    '''

# Convert to numpy arrays and get basic things
    
    xarr = np.asarray(xarr)
    ndat = len(xarr)

# Deal with array_like versus int/float

    if hasattr(x,'__len__') is False:

        single = True
        x = np.asarray([x])        

    else:

        single = False
        x = np.asarray(x)                

# Initialize binary search area and compute number of divisions needed

    ileft = np.zeros(len(x),dtype=int)
    iright = np.zeros(len(x),dtype=int)

    ndivisions = int(np.log10(ndat)/np.log10(2) + 1)

# Test for monotonicity

    i = xarr - np.roll(xarr,1)
    i = i[1:]

# Is it increasing?
    
    a = i >= 0

    if np.sum(a) == ndat-1:

        iright +=ndat-1

    else:

        a = i <= 0

        if np.sum(a) == ndat-1:

            ileft +=ndat-1        

        else:

            raise ValueError('findidx:  arr is not monotonic.')

# Perform binary search by dividing search interval in half NDIVISIONS times

    for i in range(ndivisions):

        idiv = (ileft + iright)//2   # Split interval in half
        xval = xarr[idiv]            # Find function values at center
        greater = x > xval           # Determine which side of X is on
        less = x <= xval
        np.multiply(ileft, less, out=ileft)
        np.add(ileft, idiv*greater, out=ileft)
        np.multiply(iright, greater, out=iright)
        np.add(iright, idiv*less, out=iright)

# linearly interpolate 

    xleft = xarr[ileft]
    xright = xarr[iright]

    mask = xleft == xright
    ieff = (xright-x)*ileft
    np.add(ieff, (x-xleft)*iright, out=ieff)
    np.add(ieff, mask*ileft, out=ieff)
    np.divide(ieff, (xright - xleft + mask), out=ieff)

# Clip points beyond interpolation to NaN
    
    z = ieff < 0
    ieff[z] = np.nan

    z = ieff > ndat-1
    ieff[z] = np.nan

# Return numpy.ndarray of float or float    

    return ieff[0] if single is True else ieff
    
