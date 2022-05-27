from numpy import zeros as npzeros
from numpy import sum as npsum
from numpy import int8 as npint8
from bitset import bitset

def combflagstack(stack,nbits=8):

    '''
    To combine bit-set flag arrays together.

    Input Parameters
    ----------------
    stack : numpy.ndarray 
        The stack of bit-set flag arrays to combine.  The stack
        can either be a stack of spectra [nspec,ndat] or a stack
        of images [nimg,nx,ny].   

    nbits : int, optional
        The number of bits that can potentially be set.  This
        routine assumes the bits are set sequentially, starting
        with the zeroth bit.  So if nbits is 2, then it will
        check the 0th and 1st bit.  The default is to check all
        eight bits

    Output Parameters
    ------------------
    numpy.ndarray
        A bit-set flag array that reflects the bit-set flags from all of
        the spectra or images.

    Procedure
    ---------
    Just some basic math.

    Example
    -------
    Consider a two spectra masks
    
    > spec1 = np.array([0,4,4,2])
    > spec2 = np.array([0,0,3,1])
    > stack = np.stack((spec1,spec2))
    > combflagstack(stack)

    [0 4 7 3]

    Consider two image masks

    > img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
    > img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
    > stack = np.stack((img1,img2))     
    > combflagstack(stack)

    [[1 2 0]
     [3 0 4]
     [0 0 0]]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.  
    Based on the mc_combflagstack.pro IDL program.

    '''

    # Determine whether this is a spectral stack or image stack
    
    ndim = stack.ndim
    shape = stack.shape

# Set up the output array
# Replace with match case statement when you upgrade to 3.10.

    if ndim == 2: comb = npzeros(shape[1],dtype=npint8)
    if ndim == 3: comb = npzeros(shape[1:2],dtype=npint8)        
        
# Now just loop over each bit requested.

    for i in range(0,nbits):

#  Identify the pixels with the particular bit set

        set = bitset(stack,i)

#  Collapse everything down one dimension

        sum = npsum(set,axis=0)

#  Identify which pixels are set

        mask = sum > 0

#  Set them to the proper bit value and add to the comb

        comb = comb+mask*2**i 

    return(comb)
    
