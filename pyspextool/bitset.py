from numpy import zeros_like as zeros_like
from numpy import int8 as int8

def bitset(array, bits):
    """
    To determine if the given bits are set in an array.

    Input Parameters
    ----------------
    array : numpy.ndarray
        numpy array to search

    bits : int, list, numpy.ndarray
        bit values to search in `array`

    Returns
    --------
    numpy.ndarray
        Returns an byte array of the same size as `array`.  An element 
        is set if any of the bits requested are set in the same element
        of `array`.

    Procedure
    ---------
    Uses the Gumley IDL ishft technique.  Note that the "first" bit 
    is denoted as zero, while the "second" bit is denoted as 1.


    Example
    --------

    > import numpy as np
    > bitset(np.array([3,4,1]),0)

    [1, 0, 1]

    > import numpy as np
    > bitset(np.array([3,4,1]),[0,3])

    [1, 0, 1]

    > import numpy as np
    > bitset(np.array([3,4,1]),[2,3])

    [0, 1, 0]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.
                 Based on the Spextool IDL mc_bitset.pro program.
    """
    
#  Define empty mask

    mask = zeros_like(array, dtype=int8)

#  test to see if bits is iterable

    try:

        iter(bits)

    except TypeError:

        bits = [bits]    

#  Loop over every bit requested and identify those pixels for
#  which that bit is set.

    for val in bits:
        tmp = (array >> val) & 1
        mask = mask | tmp    

    return mask
