import numpy as np

from pyspextool.io.check import check_parameter

def apertures_to_array(list:list,
                       norders:int,
                       napertures:int):

    """
    To convert a list of extraction apertures to an array suitable for FITS

    Parameters
    ----------
    list : list
        A (norders*napertures) list of ndarrays arrays.  Each array is a 
        (norders*napertures, 4, nwavelengths).

    norders : int
        The number of orders.

    napertures : int
        The number of apertures

    Returns
    -------
    ndarray

    
    """
    
    #
    # Check parameters
    #

    check_parameter('apertures_to_array', 'list', list, 'list')

    check_parameter('apertures_to_array', 'norders', norders, 'int')

    check_parameter('apertures_to_array', 'napertures', napertures, 'int')

    #
    # Start the process
    #
    
    # Determine the maximum size of each aperture array

    npixels = []
    for slice in list:
        npixels.append(np.shape(slice)[1])

    max_npixels = np.max(np.asarray(npixels))

    # Now create arrays into which the slices will be placed
    # Fill with NaN and zeros.

    array = np.full((norders * napertures, 4, max_npixels), np.nan)
    array[:,3,:] = 0.0

    # Now fill in the arrays

    l = 0
    for slice in list:

        array[l, :, 0:npixels[l]] = slice
        l += 1


    return array
