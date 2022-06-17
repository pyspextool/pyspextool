import numpy as np

def mkimgidxs(nrows,ncols,dtype=int):

    '''

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

    '''
        
    ximg = np.tile(np.arange(ncols,dtype=dtype),(nrows,1))
    yimg = np.tile(np.reshape(np.arange(nrows,dtype=dtype),\
                   (nrows,1)),(1,ncols))

    return (ximg, yimg)
