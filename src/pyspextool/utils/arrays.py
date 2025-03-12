import numpy as np
import numpy.typing as npt
import re

from pyspextool.io.check import check_parameter


def find_index(x:npt.ArrayLike,
               x_want:npt.ArrayLike | int | float,
               ends_to_nan:bool=False):

    """
    Finds the effective index of a value(s) in an ordered array.

    Parameters
    ----------
    x : ndarray
        An (n,) array of floats or integers to be searched,
        must be monotonically increasing or decreasing.

    x_want : ndarray or float or int
        An int or float value or (m,) array whose indices are desired.  

    Returns
    -------
    ndarray or float
        Effective indices in `x` for `x_want`.

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

    IEFF values outside `x` are set to NaN.

    Examples
    --------
    >>> x_array = [1, 2.5, 3, 5.5]
    >>> x = [1.5, 3.1, 6]
    >>> find_index(x, x_want)
    array([0.33333333, 2.04      ,        nan])

    """

    # Convert to numpy arrays and get basic things
    try:
        x = np.asarray(x, dtype='float')
    except ValueError:
        raise

    ndat = len(x)

    # Deal with array_like versus int/float and convert to array
    if isinstance(x_want, int) or isinstance(x_want, float) is True:
        try:
            x_want = np.asarray([x_want], dtype='float')
            single = True
        except ValueError:
            raise
    else:
        try:
            single = False
            x_want = np.asarray(x_want, dtype='float')
        except ValueError:
            raise

    # Initialize binary search area and compute number of divisions needed
    ileft = np.zeros(len(x_want), dtype=int)
    iright = np.zeros(len(x_want), dtype=int)

    n_divisions = int(np.log10(ndat)/np.log10(2) + 1)

    # Test array for monotonicity. If not, raise ValueError
    i = x - np.roll(x, 1)
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
        x_val = x[idiv]            # Find function values at center
        greater = x_want > x_val           # Determine which side of X is on
        less = x_want <= x_val
        np.multiply(ileft, less, out=ileft)
        np.add(ileft, idiv*greater, out=ileft)
        np.multiply(iright, greater, out=iright)
        np.add(iright, idiv*less, out=iright)

    # linearly interpolate
    x_left = x[ileft]
    x_right = x[iright]

    mask = x_left == x_right
    ieff = (x_right - x_want) * ileft
    np.add(ieff, (x_want - x_left) * iright, out=ieff)
    np.add(ieff, mask * ileft, out=ieff)
    np.divide(ieff, (x_right - x_left + mask), out=ieff)

    # Clip points beyond interpolation to NaN

    z = x_want < float(x[0])
    ieff[z] = np.nan if ends_to_nan is True else 0

    z = x_want > float(x[-1])
    ieff[z] = np.nan if ends_to_nan is True else ndat-1

    # Return numpy.ndarray of float or float
    return ieff[0] if single is True else ieff


def make_image_indices(nrows:int,
                       ncols:int,
                       dtype=int):

    """

    To generate pixel indice arrays for an image

    Parameters
    ----------
    nrows : int
        The number of rows in the image

    ncols : int
        The number of columns in the image

    dtype : dtype, default, `int`, optional
        The data type of the output arrays

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


    """

    #
    # Check parameters
    #

    check_parameter('make_image_indices', 'nrows', nrows, 'int')

    check_parameter('make_image_indices', 'ncols', ncols, 'int')    

    #
    # The old IDL rebin/reform trick
    # 
        
    ximg = np.tile(np.arange(ncols, dtype=dtype), (nrows, 1))
    yimg = np.tile(np.reshape(np.arange(nrows, dtype=dtype),
                   (nrows, 1)), (1, ncols))

    return ximg, yimg



def trim_nan(arr:npt.ArrayLike,
             flag:int=0,
             trim:bool=False):

    """
    To clip NaNs from an array.


    Parameters
    ----------
    arr : ndarray
        An (ndat,) array

    flag : {0,1,2,3}
        0 - The trailing NaNs are removed (default).
        1 - The leading NaNs are removed.
        2 - Both leading and trailing NaNs are removed.
        3 - All NaNs are removed.
    
    trim : {False, True} 
        If True, the program will return the trimmed array instead of a bool
        array.

    Returns
    -------
    ndarray
         A bool array where True indicates the number is not a NaN and False
         indicates the number is a NaN.

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

        return np.fliplr(np.rot90(img, 1))
    
    elif direction == 7:               

        return np.flipud(img)
        
    else:

        raise ValueError('Unknown rotation direction.')


def idl_unrotate(img, direction):

    """
    Rotates and/or tranposes an image (rotate is in multiples of 90deg).

    The program effectively undoes what the idl_rotate program does.

    Input Parameters
    ----------------
    img : numpy.ndarray
        an image

    direction : int
        This is the direction passed to the idl_rotate function that is 
        to be undone.  

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

        return np.rot90(img, 1)

    elif direction == 2:

        return np.rot90(img, 2)

    elif direction == 3:

        return np.rot90(img, 3)
        
    elif direction == 4:

        return np.transpose(img)
        
    elif direction == 5:

        return np.fliplr(img)
        
    elif direction == 6:

        return np.fliplr(np.rot90(img, 1))
    
    elif direction == 7:               

        return np.flipud(img)
        
    else:

        raise ValueError('Unknown rotation direction.')
        
        
def numberList(inp,sort=True):
    '''
    Purpose
    -------
    Converts between a string listing of numbers and an array of integers
    
    Parameters
    ----------

    inp : string or list/numpy.array
        input as a number list ('45,50-67,69,72-90') or array of integers ([1,2,3,7,8,9])

    sort : bool, default = False
        set to True to sort the output list (string -> array)

    Outputs
    -------
        if input is a string, returns the ist of integers specified by string
        if input is an array of integers, returns a compressed list of the numbers

    Example
    -------
        >>> from pyspextool.utils.arrays import numberList
        >>> print(numberList([1,2,6,8,9,12,15]))
            1-2,6,8-9,12,15
        >>> print(numberList('1-4,9-12'))
            [1, 2, 3, 4, 9, 10, 11, 12]

    Dependencies
    ------------
    numpy
    re

    '''
# string --> number list
    if isinstance(inp,str):
        numlist = []
        tmp1 = inp.replace(' ','')
        tmp2 = re.split('[,;]',tmp1)
        for a in tmp2:
            tmp4 = a.split('-')
            if len(tmp4) > 1: numlist.extend(list(range(int(tmp4[0]),int(tmp4[1])+1)))
            else: numlist.append(int(tmp4[0]))
        
        if sort==True: numlist = sorted(numlist)

# number list --> string
    elif isinstance(inp,list) or isinstance(inp,np.ndarray) or isinstance(inp,set):
        tmp = [int(x) for x in inp]
        tmp = np.array(sorted(list(set(tmp))))
        diff = np.abs(tmp-np.roll(tmp,1))-1
        w = np.where(diff > 0)[0]
        lnum = tmp[0]
        numlist = str(lnum)
        if len(w)>1: 
            for i,j in enumerate(w[1:]):
                if tmp[j-1]>lnum: numlist+='-{},{}'.format(str(tmp[j-1]),str(tmp[j]))
                else: numlist+=',{}'.format(str(tmp[j]))
                lnum = tmp[j]
        if tmp[-1]!=lnum: numlist+='-{}'.format(str(tmp[-1]))

# don't recognize
    else: raise ValueError('Could not interpret input {} which is type {}'.format(inp,type(inp)))

# return output
    return numlist







    
