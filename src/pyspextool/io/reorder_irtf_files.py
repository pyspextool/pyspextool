import numpy as np

from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.check import check_parameter

def reorder_irtf_files(fullpaths:str | list):

    """
    reorder_irtf_files


    Parameters
    ----------
    fullpaths : str or list
        The fullpaths of IRTF FITS files

    Returns
    -------
    tuple
         tuple[0] : A list fullpaths names ordered such that the beams are 
                    ababab...

         tuple[1] : A npdarray index array ordered such that the beams are 
                    ababab...
        
    """

    #
    # Check parameters
    #
    
    check_parameter('reorder_irtf_files', 'fullpaths', fullpaths,
                    ['str','list'])

    #
    # Check a few things
    #
    
    if isinstance(fullpaths, str):
        fullpaths = [fullpaths]
    
    nfiles = len(fullpaths)

    if nfiles % 2 !=0:

        message = '`fullpaths` must contain an even number of images.'
        raise pySpextoolError(message)

    #
    # Create index array
    #

    indices = np.arange(nfiles)

    #
    # Now do the loop
    #
    
    for i in range(nfiles//2):

        # Pull the beam label from the file name
        
        position = fullpaths[i*2].rfind('.fits')
        beam = fullpaths[i*2][position-1]

        
        # If it is 'b', switch the filesnames
        
        if beam == 'b':

            tmp1 = fullpaths[i*2] 
            fullpaths[i*2] = fullpaths[i*2+1]
            fullpaths[i*2+1] = tmp1

            tmp2 = indices[i*2]
            indices[i*2] = indices[i*2+1]
            indices[i*2+1] = tmp2
            
    return fullpaths, indices
