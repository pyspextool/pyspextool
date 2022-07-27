"""Functions for image utlitiy management."""


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


def ten(val):

    """
    Converts a sexigesmal number to a decimal


    Input Parameters
    ----------------
    val : str, list, numpy.ndarray 
          A sexigesimal number that is a colon-delimited string or 
          3-element list of numpy.npdarray

    Returns
    --------
    float
        The decimal number

    Procedure
    ---------
    Based on the IDL Astronomy User's Library sixty program.  
    Basic manipulation and formating

    Examples
    --------
    > x = '-00:00:40.04424'
    > ten(x)

    -0.0111234

    > x = [-0.0,0.0,40.04424]
    > ten(x)

    -0.0111234

    > x = [0.0,0.0,-40.04424]
    > ten(x)

    -0.0111234

    > import numpy as np
    > x = np.array([-0.0,0.0,40.04424])
    > ten(x)

    -0.0111234

    > import numpy as np
    > x = np.array([0.0,0.0,-40.04424])
    > ten(x)

    -0.0111234


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """
    
    # Figure out what the name is 
    
    typ = type(val).__name__

    # String input
    
    if typ == 'str':

        hms=(val.split(':'))

        decimal = abs(float(hms[0]))+float(hms[1])/60.+float(hms[2])/3600.

        # Grab the first element of the string to test for positivity 
        
        posneg = hms[0][0]

        if posneg == '+':
            
            return(decimal)
        
        elif posneg == '-':
            
            return(-1*decimal)        
    
        else:

            return(decimal)        

    # A list of numpy array
        
    elif typ == 'list' or typ == 'ndarray':

        # Convert to positive (to deal with things like [-0.0,0.0.40]
        
        decimal = abs(float(val[0])+float(val[1])/60.+float(val[2])/3600.)

        # Check for negative
        
        prod = val[0]*val[1]*val[2]
        if prod == -0.0:

            decimal *=-1

        return(decimal)
        

def sixty(val,colons=None,trailsign=False):

    """
    Converts a decimal number to sexigesmal


    Input Parameters
    ----------------
    val : float 
        A decimal number.

    colons: dict of {'dec':int, 'plus':int}, optional
        If given, the output is a colon-separated string.  'dec' gives 
        the number of decimal places on the third value set 'plus' 
        for plus sign if the first value is positive.  

    trailsign: {True, False}
        If `colons` is given, this has no affect.  By default (False), 
        the first non-zero value has the negative sign.  Setting 
        `trailsign` forces the first element to have the negative sign.  

    Returns
    --------
    list, str
        list:  a 3-element list giving the sexigesimal values
        str: a string of the form hh:mm:ss 

    Procedure
    ---------
    Based on the IDL Astronomy User's Library sixty program.  
    Basic manipulation and formating

    Examples
    --------
    > x = -0.0111234
    > sixty(x)
    
    [0.0, 0.0, -40.04424]

    > x = -0.0111234
    > sixty(x,trainsign=True)
    
    [-0.0, 0.0, 40.04424]

    > x = -0.0111234
    > sixty(x,trainsign=True)
    
    [-0.0, 0.0, 40.04424]

    > x = -0.0111234
    > sixty(x,colons={'dec':2,'plus':1},trailsign=True)

    -00:00:40.04

    > x = 0.0111234
    > sixty(x)

    [0.0, 0.0, 40.04424]

    > x = 0.0111234
    > sixty(x,trailsign=True)

    [0.0, 0.0, 40.04424]

    > x = 0.0111234
    > sixty(x,colons={'dec':2,'plus':1},trailsign=True)

    +00:00:40.04

    > x = 0.0111234
    > sixty(x,colons={'dec':3,'plus':0},trailsign=True)

    00:00:40.044


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """    

    # First check to see if the value is negative
    
    neg = [0,1][val < 0]

    # Convert the value to positive degrees, minutes, and seconds
    
    ss = abs(3600.*val)
    mm = abs(60.*val)
    dd = abs(val)

    # Now determine the positive sexigesimal values
    
    sexg = [int(dd)]
    sexg.append(int(mm-60.*sexg[0]))
    sexg.append(ss-3600.*sexg[0]-60*sexg[1])

    # Check to see whether a colon-separated string is requested.
    
    if colons is not None:

        # Create the positive colon-separated string.
        
        fmt = '{:.'+str(colons['dec'])+'f}'

        result = str(sexg[0]).zfill(2)+':'+\
                 str(sexg[1]).zfill(2)+':'+\
                 fmt.format(sexg[2]).zfill(3+colons['dec'])

        # Now deal with the positive/negative issue

        if neg == 1:

            result = '-'+result

        else:

            result = [result,'+'+result][colons['plus'] == 1]
            
    if colons is None:

        sexg = [float(x) for x in sexg]
        if neg == 1:

            if trailsign is True:

                sexg[0]*=-1

            else:

                if sexg[0] != 0:

                    sexg[0]*=-1

                elif sexg[1] != 0:

                    sexg[1]*=-1

                elif sexg[2] != 0:

                    sexg[2]*=-1                    
                    
        result = sexg

    return(result)
    