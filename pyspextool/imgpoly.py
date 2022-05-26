import numpy

def imgpoly(img,coeffs):
    '''
    Evaluate a polynomial where the independent variable is an image.

    Input Parameters:
       img    - a 2D numpy array of size (NAXIS1,NAXIS2)
       coeffs - a 3D numpy array of coefficients (NAXIS1, NAXIS2, ncoeffs)

    Optional Parameters:
       None

    Output Parameters:
       Returns a 2D numpy array evaluted at each pixel in img.

    Procedure:
       Follows the IDL poly.pro technique to evaluate a polynomial

    Example:
       NA

    Modification History:
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the mc_polyimg.pro IDL program.
    '''

    n = coeffs.shape[2]-1

    y = coeffs[:,:,n]

    for i in range(n-1,-1,-1):

        y = y*img + coeffs[:,:,i]
    
    return(y)
