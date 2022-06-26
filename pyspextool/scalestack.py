import numpy as np
import sys

def scalestack(stack,*var,mask=None,idx=None):

    '''
    Scales a stack of spectra or images to a common intensity level.             

    Input Parameters
    ----------------
    stack : array_like
        (nspec,nwave) - a stack of spectra.

        (nimgs,nrows,cols) - a stack of images.


    *var : array_like, optional
        The variances associated with the spectra or images.


    mask : array_like, optional
        A mask array with the same shape as `stack`.  
        0 = bad, 1=good


    idx : int, optional 
        A mask array with the same shape as `data`.  
        0 = bad, 1=good


    Returns
    --------
    sstack : numpy.ndarray
        The scaled stack.


    svar : numpy.ndarray or None
        The scaled variance, if *var is passed.


    scales : numpy.ndarray
        The scale factors.


    Notes
    -----
    Computes the median of each spectrum or image and then determines 
    scale factors to scale each spectrum or image to either the median 
    value of the median values or the median value associated 
    with a particular spectrum or image.


    Examples
    --------
    Spectral stack, scale to median:

        > spec1 = np.array([1,1,1,1])
        > spec2 = np.array([2,2,2,2])
        > spec3 = np.array([4,4,4,4])
        > specstack = np.stack((spec1,spec2,spec3))
        > scaledstack, var, scales = scalestack(specstack,None)
        > print(scaledstack)
        > print(scales)
          [[2. 2. 2. 2.]
           [2. 2. 2. 2.]
           [2. 2. 2. 2.]]
          [2.  1.  0.5]

    Spectral stack, scale to first spectrum:

        > spec1 = np.array([1,1,1,1])
        > spec2 = np.array([2,2,2,2])
        > spec3 = np.array([4,4,4,4])
        > specstack = np.stack((spec1,spec2,spec3))
        > scaledstack, var, scales = scalestack(specstack,None,idx=0)
        > print(scaledstack)
        > print(scales)
          [[1. 1. 1. 1.]
           [1. 1. 1. 1.]
           [1. 1. 1. 1.]]
          [1.   0.5  0.25]

    Image stack, scale to median:

        > img1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
        > img2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
        > img3 = np.array([[4,4,4],[4,4,4],[4,4,4]])
        > imgstack = np.stack((img1,img2,img3))
        > scaledstack, var, scales = scalestack(imgstack, None)
        > print(scaledstack)
        > print(scales)
          [[[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]

           [[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]

           [[2. 2. 2.]
            [2. 2. 2.]
            [2. 2. 2.]]]
          [2.  1.  0.5]

    Image stack, scale to first image:

        > img1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
        > img2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
        > img3 = np.array([[4,4,4],[4,4,4],[4,4,4]])
        > imgstack = np.stack((img1,img2,img3))
        > scaledstack, var, scales = scalestack(imgstack, None)
        > print(scaledstack)
        > print(scales)

          [[[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]

           [[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]

           [[1. 1. 1.]
            [1. 1. 1.]
            [1. 1. 1.]]]
          [1.   0.5  0.25]

    Modification History
    --------------------
    2022-06-18 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_getspecscales.pro

    '''

# Get array dimensions

    ndimen = np.ndim(stack)
    shape = np.shape(stack)

# If no mask passed, create one.

    if mask is None:

        mask = np.ones(shape,dtype=int)

# Now search and replace any masked pixels with NaNs
        
    stack = np.where(mask != 0, stack, np.nan)

# Set up median, reshape, and axis variables
    
    if ndimen == 2:

# Spectral stack
        
        axis_info = 1
        reshape_info = (shape[0],1)
        tile_info = (1,shape[1])

    elif ndimen == 3:

# Image stack

        axis_info = (1,2)
        reshape_info = (shape[0],1,1)
        tile_info = (1,shape[1],shape[2])        
        
    else:

        print('scalestack:  Unknown data shape.')
        sys.exit(1)

# Do the calculation

    medvals = np.nanmedian(stack,axis=axis_info)

# Check whether you are scaling to the median or spectrum/image 

    if idx is None:

# Scale to the median 
            
        scales = np.median(medvals)/medvals
                        
    else:

# Scale to a spectrum or image
            
        scales = medvals[idx]/medvals
                    
# Build the final scale array and do the math

    sclarr = np.tile(np.reshape(scales,reshape_info),tile_info)
    
    sstack = stack*sclarr

# Now return the scaled stack and potentially the scaled variance

    if len(var) == 1:

        return (sstack,var*sclarr**2,scales)

    else: return (sstack,None,scales)
