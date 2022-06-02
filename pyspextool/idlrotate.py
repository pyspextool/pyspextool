from sys import exit
from numpy import rot90 as nprot90
from numpy import flipud as npflipud
from numpy import fliplr as npfliplr
from numpy import transpose as nptranspose

def idlrotate(img,direction):

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

        return(img)        

    elif direction == 1:

        return(nprot90(img,3))

    elif direction == 2:

        return(nprot90(img,2))

    elif direction == 3:

        return(nprot90(img,1))
        
    elif direction == 4:

        return(nptranspose(img))
        
    elif direction == 5:

        return(npfliplr(img))        
        
    elif direction == 6:

        return(npfliplr(nprot90(img,1)))
    
    elif direction == 7:               

        return(npflipud(img))        
        
    else:

        print('idlrotate:  Unknown direction.')
        exit(1)
