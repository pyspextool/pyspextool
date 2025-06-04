import numpy as np
import logging

from pyspextool.extract import config
from pyspextool.io.check import check_parameter
from pyspextool.io.check import check_qakeywords

from pyspextool.pyspextoolerror import pySpextoolError


def override_aperturesigns(aperture_signs:str,
                           verbose=None):

    """
    To force the apertures signs to the users choice.

    This is useful of the automated functions do not properly identivy the 
    aperture signs.
    
    Parameters
    ----------
    apsigns : str
        A string of comma-separated + and - signs giving the aperture signs.
        E.g. '+,-'

    Returns
    -------
    ndarray
        A (naps,) integer array of 1, -1 values.
    

    Updates the config.state['aperture_signs']
    
    """

    #
    # Check parameters and keywords
    #

    check_parameter('override_aperturesigns', 'aperture_signs', aperture_signs,
                    'str')

    qa = check_qakeywords(verbose=verbose)

    #
    # Determine whether the requested number of signs matches the number being
    # extracted.
    #

    naps = len(config.state['average_aperturesigns'])
    user_naps = len(aperture_signs.split(','))

    # Does the request have the same number?

    if naps != user_naps:

        message = 'You have requested to set the signs of '+str(user_naps)+\
        ' apertures but the number of apertures being extracted is '+\
        str(naps)+'.'
        raise pySpextoolError(message)

    label = aperture_signs

    # 
    # Set the aperture signs
    #

    aperture_signs = aperture_signs.replace('+', '1')    
    aperture_signs = aperture_signs.replace('-', '-1')    


    message = ' Aperture signs are (' + label + ').'
    logging.info(message)

    aperture_signs = np.array(aperture_signs.split(','),dtype=int)

    # Update the state variable

    config.state['average_aperturesigns'] = aperture_signs

    # Return the value for possible testing

    return aperture_signs

    
