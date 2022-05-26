import numpy as np
from bitset import bitset

def combflagstack(stack,nbits=8):

    '''
    To combine bit-set flag arrays together.

    Input Parameters:
        stack - The stack of bit-set flag arrays to combine.  The stack
                can either be a stack of spectra [ndat,nspec] or a stack
                of images [nx,ny,nimg].   

    Optional Parameters:
        nbits - The number of bits that can potentially be set.  This
                routine assumes the bits are set sequentially, starting
                with the zeroth bit.  So if nbits is 2, then it will
                check the 0th and 1st bit.  The default is to check all
                eight bits

    Output Parameters:
        A bit-set flag array that reflects the bit-set flags from all of
        the spectra or images.

    Procedure:
        Just some basic math.

    Example:
        TBD

    Modification History:
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the mc_combflagstack.pro IDL program.
    '''

# Determine whether this is a spectral stack or image stack
    
    ndim = stack.ndim
    shape = stack.shape

# Set up the output array
# Replace with match case statement when you upgrade to 3.10.

    if ndim == 2: comb = np.zeros(shape[1],dtype=np.int8)
    if ndim == 3: comb = np.zeros(shape[0:2],dtype=np.int8)        
        
# Now just loop over each bit requested.

    for i in range(0,nbits):

#  Identify the pixels with the particular bit set

        set = mc_bitset(stack,[i])
        
#  Collapse everything down one dimension

        sum = np.sum(set,axis=ndim-1)

#  Identify which pixels are set

        mask = sum > 0

#  Set them to the proper bit value and add to the comb

        comb = comb+mask*2**i 

    
    return(comb)
    
