import numpy as np


def bitset(array, bits):
    """
    To determine if the given bits are set in an array.

    Input Parameters
    ----------------
        array : array
            A numpy array to search.
        bits : list or array
            A list or numpy array of bits to search.
            Note that the "first" bit is denoted as zero,
            while the "second" bit is denoted as 1.

    Optional Parameters:
        None

    Returns
    --------
    array
        Returns a byte array of the same size as array.  A pixel is
        set if any of the bits requested are set in the same pixel
        of array.

    Procedure
    ---------
        Uses the Gumley IDL ishft technique.

    Example
    --------
        print(mc_bitset(np.array([3,4,1]),[0]))
        [1 0 1]

    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.
                     Based on the mc_bitset.pro IDL program.
    """
    
    #  Define empty mask
    mask = np.zeros_like(array, dtype=np.int8)

    #  Loop over every bit requested and identify those pixels for which that bit is set.
    for val in bits:
        tmp = (array >> val) & 1
        mask = mask | tmp    

    return mask
