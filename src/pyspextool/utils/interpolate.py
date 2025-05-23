import numpy as np
import numpy.typing as npt

from pyspextool.io.check import check_parameter
from pyspextool.utils.arrays import find_index
from pyspextool.utils.math import bit_set
from pyspextool.utils.for_print import for_print
from pyspextool.pyspextoolerror import pySpextoolError

def linear_interp1d(input_x:npt.ArrayLike,
                    input_y:npt.ArrayLike,
                    output_x:npt.ArrayLike,
                    input_u:npt.ArrayLike=None,
                    leave_nans:bool=False):

    """
    To preform 1D linear interpolation with optional error propagation.

    Parameters
    ----------
    input_x : ndarray
        An (ndat,) array of independent values.  Can include NaNs, but should 
        be monotonically increasing.

    input_y : ndarray
        An (ndat,) array of dependent values.  Can include NaNs.

    output_x : ndarray, int, float
        An (ndat,) array or scalar of requested independent values.  Can 
        include NaNs, but should be monotonically increasing.

    input_u : ndarray, optional
        An (ndat,) array of uncertainty values.  Can include NaNs.

    leave_nans : {False, True}
        Set to True allow the NaN to propagate through the interpolation.  


    Returns
    -------
    

    """

    #
    # Check parameters
    #

    check_parameter('linear_interp1d', 'input_x', input_x, ['ndarray', 'list'])

    check_parameter('linear_interp1d', 'input_y', input_y, ['ndarray', 'list'])

    check_parameter('linear_interp1d', 'output_x', output_x,
                    ['ndarray', 'list', 'int', 'float'])

    check_parameter('linear_interp1d', 'input_u', input_u,
                    ['NoneType', 'ndarray', 'list'])    

    check_parameter('linear_interp1d', 'leave_nans', leave_nans, 'bool')    

    #
    #  Convert possible lists to numpy arrays
    #
    
    if isinstance(input_x, list) is True:

        input_x = np.array(input_x)

    if isinstance(input_y, list) is True:

        input_y = np.array(input_y)

    if isinstance(input_u, list) is True:

        input_u = np.array(input_u)                                

    # Now deal with the output_x array.  First check to see if it is a scalar,
    # and then check if it is a list.

    if np.isscalar(output_x):

        scalar = True
        output_x = np.array([output_x])        

    else:

        scalar = False
            
    if isinstance(output_x, list) is True:

        output_x = np.array(output_x)

    #
    #  Remove NaNs
    #

    # Deal with the inputs first

    if leave_nans is True:

        # Only worry about the input_x array
        
        z_input_nonan = ~np.isnan(input_x)        

    else:

        # Deal with the input_y and potentially the input_u arrays
        
        z_input_nonan = ~np.isnan(input_y)*~np.isnan(input_x)

        if input_u is not None:

            z_input_u_nonan = ~np.isnan(input_u)            
            np.multiply(z_input_nonan,z_input_u_nonan, out=z_input_nonan)
    
    input_x = input_x[z_input_nonan]
    input_y = input_y[z_input_nonan]

    if input_u is not None:

        input_u = input_u[z_input_nonan]        
    
    # Now deal with the output x array.

    z_output_nonan = ~np.isnan(output_x)
    x = output_x[z_output_nonan]    

    #
    # Do the interpolation
    #
        
    result = nonan_interp1d(input_x, input_y, x)

    # Now create an output_y array the same size as output_x and fill with
    # result.
    
    output_y = np.full_like(output_x, np.nan, dtype=np.float64)
    output_y[z_output_nonan] = result
    
    #
    # Do the error propagation
    #

    if input_u is not None:

        result = nonan_interp1d(input_x, input_u**2, x, variance=True)    

        output_v = np.full_like(output_x, np.nan, dtype=np.float64)
        output_v[z_output_nonan] = result

    #
    # Return the values accordingly
    #

    if input_u is not None:

        if scalar is True:

            return output_y[0], np.sqrt(output_v[0])

        else:

            return output_y, np.sqrt(output_v)            
        
    else:

        if scalar is True:
        
            return output_y[0]

        else:

            return output_y


        
def linear_bitmask_interp1d(input_x, input_y, output_x, nbits=8):

    """
    To preform 1D linear interpolation on a bit mask.

    The function let's bitset pixels to bleed to adjacent pixels if they
    they are used in the linear interpolation.  

    Parameters
    ----------

    input_x : ndarray 
        An (ndat,) array of independent values.  Can include NaNs, but should 
        be monotonically increasing.

    input_y : ndarray of int
        An (ndat,) array of dependent values. 

    output_x : 
        An (ndat,) array or scalar of requested independent values.  Can 
        include NaNs, but should be monotonically increasing.

    nbits : int, default=8
        The number of bits to loop over.

    Notes
    -----
    Out of range interpolates are not marked with NaNs but rather 0.  

    Returns
    -------
    ndarray of int
        
    """

    #
    # Check parameters
    #

    check_parameter('linear_bitset_interp1d', 'input_x', input_x,
                    ['ndarray', 'list'])

    check_parameter('linear_bitset_interp1d', 'input_y', input_y,
                    ['ndarray', 'list'])

    check_parameter('linear_bitset_interp1d', 'output_x', output_x,
                    ['ndarray', 'list', 'int', 'float'])

    check_parameter('linear_bitset_interp1d', 'nbits', nbits, 'int')

    #
    #  Convert any list inputs to numpy arrays
    #

    if isinstance(input_x, list) is True:

        input_x = np.array(input_x)

    if isinstance(input_y, list) is True:

        input_y = np.array(input_y)        
    
    # Now deal with the output_x array.  First check to see if it is a scalar,
    # then deal with the list

    if np.isscalar(output_x):

        scalar = True
        output_x = np.array([output_x])        

    else:

        scalar = False
            
    if isinstance(output_x, list) is True:

        output_x = np.array(output_x)

    #
    # Now squeeze them all in case you get passed a 1D slice
    #

    input_x = np.squeeze(input_x)
    input_y = np.squeeze(input_y)
    output_x = np.squeeze(output_x)

    #
    #  Remove NaNs
    #

    z_input_nonan = ~np.isnan(input_x)        
    input_x = input_x[z_input_nonan]

    z_output_nonan = ~np.isnan(output_x)        
    x = output_x[z_output_nonan]

    #
    # Loop over each bit requested
    #

    # Create an output_y array the same size as output_x and fill with zeros.
    
    output_y = np.full_like(output_x,0)

    for i in range(nbits):

        # Create a mask for each requested bit.

        is_set = bit_set(input_y, i)

        # Do the interpolation
        
        result = nonan_interp1d(input_x, is_set, x)
    
        # Convert NaNs to zero as they are out of range.

        z_nan = np.isnan(result)
        if np.sum(z_nan) != 0:

            result[z_nan] = 0.0

        # Create the interpolated mask and store
            
        mask = np.ceil(result)*2**i
        
        output_y[z_output_nonan] = output_y[z_output_nonan] + mask
                                   
    #
    # Return the results
    #
    
    # Convert to integers
        
    tmp = np.empty_like(output_y, dtype=np.int16)
    np.ceil(output_y, out=tmp, casting='unsafe')    

    # Check whether the input was a scalar and act accordingly
    
    if scalar is True:
        
        return tmp[0]

    else:

        return tmp



def nonan_interp1d(input_x, input_y, x, variance=False):

    """
    To perform a linear interpolation assuming non NaNs present.

    Parameters
    ----------
    input_x : ndarray 
        An (ndat1,) array of independent values.  Can include NaNs, but should 
        be monotonically increasing.

    input_y : ndarray
        An (ndat1,) array of dependent values. 

    x : ndarray or float or int
        An (ndat2,) array or scalar of requested independent values.  Should
        be monotonically increasing.

    variance : {False, True}
        Set to True to interpolate assuming the values are variances.

    Returns
    -------
    ndarray
          
    """

    # Determine the indices of x on input_x.  Things that don't overlap will
    # be set to Nan (on purpose for ease of identification).

    idx = find_index(input_x, x, ends_to_nan=True)
    
    # Check to see if there are at least some non NaN values. If not, then
    # the requested 'x' values are not in range of input_x.

    z_idx_nonan = ~np.isnan(idx)
    if np.sum(z_idx_nonan) == 0:

        message = 'nonan_interp1d from linear_inter1d:  `x` is not in range of `input_x`'
        raise pySpextoolError(message)
        
    # Trim the idx array and x array to ensure only values within
    # input_x.  Then create an array for the outputs.  
        
    idx = idx[z_idx_nonan]
    inrange_x = x[z_idx_nonan]
    inrange_y = np.zeros(np.sum(z_idx_nonan))

    # Determine the floor and ceil of each value for later use.  have to do it
    # this way to avoid using astype()
    
    floor = np.empty_like(idx, dtype=np.int64)
    np.floor(idx, out=floor, casting='unsafe')

    ceil = np.empty_like(idx, dtype=np.int64)
    np.ceil(idx, out=ceil, casting='unsafe')    
    
    # See which points in inrange_output_x_nonan land directly on points in
    # input_x_nonan.
    
    z_idx_same = ceil == floor

    if np.sum(z_idx_same) !=0:

        inrange_y[z_idx_same] = input_y[floor[z_idx_same]]

    # Now deal with the points that don't.

    if np.sum(z_idx_same) != np.sum(z_idx_nonan):

        # Find the pixels that aren't identical
        
        z_idx_differ = floor-ceil != 0

        # Compute the slope

        Deltay = input_y[ceil[z_idx_differ]] - input_y[floor[z_idx_differ]]
        Deltax = input_x[ceil[z_idx_differ]] - input_x[floor[z_idx_differ]]
        m = Deltay/Deltax

        # Get the dx

        dx = inrange_x[z_idx_differ] - input_x[floor[z_idx_differ]]
        
        if variance is False:

            # y = y_1 + m * dx where m = (y2-y1)/(x2-x1), dx=x-x1        

            inrange_y[z_idx_differ] = m*dx + input_y[floor[z_idx_differ]]

        else:

            # var_y = [1-dx/Deltax]**2 * var_y1 + [dx/Deltax]**2 * var_y2

            input_y = input_y**2
            
            term1 = (1-dx/Deltax)**2 * input_y[floor[z_idx_differ]]**2
            term2 = (dx/Deltax)**2 * input_y[ceil[z_idx_differ]]**2

            inrange_y[z_idx_differ] = term1 + term2
                
    # Make an output array that is the same size as x.

    y = np.full_like(x, np.nan,dtype=np.float64)

    # Now fill the inrange values back into the entire y array

    y[z_idx_nonan] = inrange_y
    
    return y 


def sinc_interpolation_fft(x: np.ndarray, s: np.ndarray, u: np.ndarray) -> np.ndarray:
    """
    Fast Fourier Transform (FFT) based sinc or bandlimited interpolation.
    
    Args:
        x (np.ndarray): signal to be interpolated, can be 1D or 2D
        s (np.ndarray): time points of x (*s* for *samples*) 
        u (np.ndarray): time points of y (*u* for *upsampled*)
        
    Returns:
        np.ndarray: interpolated signal at time points *u*
    """
    num_output = len(u)

    # Compute the FFT of the input signal
    X = np.fft.rfft(x)

    # Create a new array for the zero-padded frequency spectrum
    X_padded = np.zeros(num_output // 2 + 1, dtype=complex)

    # Copy the original frequency spectrum into the zero-padded array
    X_padded[:X.shape[0]] = X

    # Compute the inverse FFT of the zero-padded frequency spectrum
    x_interpolated = np.fft.irfft(X_padded, n=num_output)

    return x_interpolated * (num_output / len(s))

def sinc_interpolation(x, s, u):
            """Whittaker–Shannon or sinc or bandlimited interpolation.
            Args:
                x (NDArray): signal to be interpolated, can be 1D or 2D
                s (NDArray): time points of x (*s* for *samples*) 
                u (NDArray): time points of y (*u* for *upsampled*)
            Returns:
                NDArray: interpolated signal at time points *u*
            Reference:
                This code is based on https://gist.github.com/endolith/1297227
                and the comments therein.
            TODO:
                * implement FFT based interpolation for speed up
            """
            sinc_ = np.sinc((u - s[:, None])/(s[1]-s[0]))

            return np.dot(x, sinc_)
