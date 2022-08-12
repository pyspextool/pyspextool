import numpy as np


def find_index(x_array, x):

    """
    Finds the effective index of a function value in an ordered array.

    Parameters
    ----------------
    x_array : array_like
        (N,) The array of floats or integers to be searched,
        must be monotonically increasing or decreasing.
    x : array_like or float or int
        (M,) The value or values whose indices are required.

    Returns
    -------
    ndarray or float
        Effective indices in `x_array` for `x`.

    Raises
    ------
    ValueError
        If input arrays are not array_like or
        if input array is not monotonic

    Notes
    -----
    Near line-by-line copy of IDL Astronomy User's Library routine tabinv.pro. 

    From its documentation:  A binary search is used to find the values XARR(I) 
    and XARR(I+1) where XARR(I) < X < XARR(I+1).  IEFF is then computed using 
    linear interpolation between I and I+1.
    IEFF = I + (X-XARR(I)) / (XARR(I+1)-XARR(I)).

    IEFF values outside of `xarr` are set to NaN.  

    Examples
    --------
    >>> x_array = [1, 2.5, 3, 5.5]
    >>> x = [1.5, 3.1, 6]
    >>> find_index(x_array, x)
    array([0.33333333, 2.04      ,        nan])

    Modification History
    --------------------
    2022-06-30 - Written by M. Cushing, University of Toledo.
    """

    # Convert to numpy arrays and get basic things
    try:
        x_array = np.asarray(x_array, dtype='float')
    except ValueError:
        raise

    ndat = len(x_array)

    # Deal with array_like versus int/float and convert to array
    if isinstance(x, int) or isinstance(x, float) is True:
        try:
            x = np.asarray([x], dtype='float')
            single = True
        except ValueError:
            raise
    else:
        try:
            single = False
            x = np.asarray(x, dtype='float')
        except ValueError:
            raise

    # Initialize binary search area and compute number of divisions needed
    ileft = np.zeros(len(x), dtype=int)
    iright = np.zeros(len(x), dtype=int)

    n_divisions = int(np.log10(ndat)/np.log10(2) + 1)

    # Test array for monotonicity. If not, raise ValueError
    i = x_array - np.roll(x_array, 1)
    i = i[1:]
    # check if array is increasing or decreasing
    a = i >= 0
    if np.sum(a) == ndat-1:
        iright += ndat-1
    else:
        a = i <= 0
        if np.sum(a) == ndat-1:
            ileft += ndat-1
        else:
            raise ValueError('find_index:  array is not monotonic.')

    # Perform binary search by dividing search interval in half NDIVISIONS times
    for i in range(n_divisions):
        idiv = (ileft + iright)//2   # Split interval in half
        x_val = x_array[idiv]            # Find function values at center
        greater = x > x_val           # Determine which side of X is on
        less = x <= x_val
        np.multiply(ileft, less, out=ileft)
        np.add(ileft, idiv*greater, out=ileft)
        np.multiply(iright, greater, out=iright)
        np.add(iright, idiv*less, out=iright)

    # linearly interpolate
    x_left = x_array[ileft]
    x_right = x_array[iright]

    mask = x_left == x_right
    ieff = (x_right-x) * ileft
    np.add(ieff, (x-x_left) * iright, out=ieff)
    np.add(ieff, mask * ileft, out=ieff)
    np.divide(ieff, (x_right - x_left + mask), out=ieff)

    # Clip points beyond interpolation to NaN
    z = ieff < 0
    ieff[z] = np.nan

    z = ieff > ndat-1
    ieff[z] = np.nan

    # Return numpy.ndarray of float or float
    return ieff[0] if single is True else ieff


def make_image_indices(nrows, ncols, dtype=int):

    """

    To generate indice grids for an image

    Input Parameters
    ----------------
    nrows : int
        The number of rows in the image

    ncols : int
        The number of columns in the image

    dtype : dtype, default, `int`, optional
        The type of the output grids


    Returns
    --------
    ximg, yimg : numpy.ndnarray, numpy.ndnarray of `dtype`
         A list of integers giving the individual file numbers


    Notes
    -----
    https://stackoverflow.com/questions/1550130/\
    cloning-row-or-column-vectors

    Examples
    --------
    > ximg, yimg = mkimgidxs(3,5)
    > print(ximg)
      [[0 1 2 3 4]
       [0 1 2 3 4]
       [0 1 2 3 4]]

    > print(yimg)
      [[0 0 0 0 0]
      [1 1 1 1 1]
      [2 2 2 2 2]]

    Modification History
    --------------------
    2022-06-17 - Written by M. Cushing, University of Toledo.
    Just the old reform/rebin trick from IDL.

    """

    ximg = np.tile(np.arange(ncols, dtype=dtype), (nrows, 1))
    yimg = np.tile(np.reshape(np.arange(nrows, dtype=dtype),
                   (nrows, 1)), (1, ncols))

    return ximg, yimg


def nan_trim(arr, flag=0, trim=False):

    """
    To clip NaNs from an array.


    Input Parameters
    ----------------
    arr : array_like
        An (ndat,) array

    flag : {0,1,2,3}, optional
        0 - The trailing NaNs are removed (default).
        1 - The leading NaNs are removed.
        2 - Both leading and trailing NaNs are removed.
        3 -  All NaNs are removed.

    trim : {False, True}, optional
        If True, the program will return the trimmed array instead of a bool
        array.

    Returns
    -------
    numpy.ndarray
         A bool array where True indicates the number is not a NaN and False
         indicates the number is a NaN.

    Notes
    -----
    None

    Examples
    --------
    > x = np.array([np.nan,2,3,4,np.nan,7,89,90,np.nan])
    > nantrim(x,flag=0,trim=True)
      [nan  2.  3.  4. nan  7. 89. 90.]
    > nantrim(x,flag=1,trim=True)
      [ 2.  3.  4. nan  7. 89. 90. nan]
    > nantrim(x,flag=2,trim=True)
      [ 2.  3.  4. nan  7. 89. 90.]
    > nantrim(x,flag=3,trim=True)
      [ 2.  3.  4.  7. 89. 90.]

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_nantrim.pro.

    """

    # Test for correct flags

    if flag not in [0, 1, 2, 3]:

        raise ValueError('unknown flag.')

    # Convert to numpy array

    arr = np.array(arr)

    # Test for 1D

    if arr.ndim != 1:

        raise ValueError('must be 1D.')

    # test the easy cases of non NaNs

    test = np.isnan(arr)

    if flag == 3 or sum(test) == 0:

        z = np.invert(test)

    if flag == 0:

        cumsum = np.nancumsum(np.flip(arr))
        z = np.flip(cumsum) != 0

    elif flag == 1:

        cumsum = np.nancumsum(arr)
        z = cumsum != 0

    elif flag == 2:

        cumsum = np.nancumsum(arr)
        zleading = cumsum != 0

        cumsum = np.nancumsum(np.flip(arr))
        ztrailing = np.flip(cumsum) != 0

        z = zleading*ztrailing

    if trim is True:

        return arr[z]

    else:

        return z


def idl_rotate(img, direction):

    """
    Rotates and/or tranposes an image (rotate is in multiples of 90deg)

    Input Parameters
    ----------------
    img : numpy.ndarray
        an image

    direction : int

        Direction  Transpose?  Rotation Counterclockwise
        -------------------------------------------------

        0          No          None
        1          No          90 deg
        2          No          180 deg
        3          No          270 deg
        4          Yes         None
        5          Yes         90 deg
        6          Yes         180 deg
        7          Yes         270 deg

    The directions follow the IDL rotate function convention.


    Returns
    --------
    numpy.ndarray
        the rotated/transposed image

    Procedure
    ---------

    uses numpy.rot90, flipud, fliplr, ant transpose to rotate and/or 
    transpose the image as requested.

    Examples
    --------

    > import numpy as np
    > img = np.array([[1,2,3],[4,5,6],[7,8,9]])
    > idlrotate(img,6)

    [[9 6 3]
    [8 5 2]
    [7 4 1]]

    Modification History
    --------------------
    2022-06-02 - Written by M. Cushing, University of Toledo.
    Based on the IDL program rotate.pro

    """
    
    if direction == 0:

        return img

    elif direction == 1:

        return np.rot90(img, 3)

    elif direction == 2:

        return np.rot90(img, 2)

    elif direction == 3:

        return np.rot90(img, 1)
        
    elif direction == 4:

        return np.transpose(img)
        
    elif direction == 5:

        return np.fliplr(img)
        
    elif direction == 6:

        return np.fliplr(np.prot90(img, 1))
    
    elif direction == 7:               

        return np.flipud(img)
        
    else:

        print('idlrotate:  Unknown direction.')
        exit(1)
    




    
