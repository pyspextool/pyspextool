import numpy as np

from pyspextool.extract import config
from pyspextool.io.check import check_parameter


def override_aperture_signs(apsigns):
    """
    To force the apertures signs to the users choice.


    Parameters
    ----------
    apsigns : numpy.dnarray
        a (naps,) array of desired aperture signs (1, or -1).


    Returns
    -------
    None
    

    Notes
    -----
    Simply updates the state variable if the number of apertures requested 
    matches the number of apertures stored in the state variable.

    Examples
    --------
    later

    """

    #
    # Check parameter
    #

    check_parameter('override_aperture_sign', 'apsigns', apsigns,
                    ['list', 'ndarray'])

    # Determine the number of apertures already identified

    naps = len(config.state['apsigns'])

    # Does the request have the same number?

    if naps != len(apsigns):
        message = 'The number of elements in `apsigns` does not match the ' + \
                  'number of apertures.'
        raise ValueError(message)

    # Swap the apertures

    apsigns = np.array(apsigns).astype(int)

    # Update the state variable

    config.state['apsigns'] = apsigns
