from numpy import exp
import numpy.typing as npt

from pyspextool.io.check import check_parameter

def irplanck(wavelength:int | float | npt.ArrayLike,
             temperature:int | float):

    """
    Gives the Planck function in useful infrared units, W m-2 um-1 sr-1.

    Parameters
    ----------
    wavelength : ndarray
        An array of wavelengths in microns.

    temperature : int or float
        The temperature in Kelvin

    Returns
    -------
    int, float, or ndarray
        The Planck function in W m-2 um-1 sr-1 at `wavelength` for a temperature
        of `temperature`.

    Notes
    -----
    Uses the formula in the chapter "Infrared Astronomy Fundamentals" from
    Planets, Stars, and Stellar Systems Volume 2: Astronomical Techniques,
    Software, Data
    
    """

    #
    # Check parameters
    #

    check_parameter('irplanck', 'wavelength', wavelength,
                    ['int', 'float', 'ndarray'])

    check_parameter('irplanck', 'temperature', temperature, ['int', 'float'])

    #
    # Compute the Planck function
    #

    planck = 1.1910e8 / wavelength**5 / (exp(14387.7/wavelength/temperature)-1)

    return planck
