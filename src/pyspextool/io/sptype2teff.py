import numpy as np
from os.path import join as osjoin

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter

def sptype2teff(sptype:str):

    """
    To obtain an effective temperature from a spectral type.

    The function uses a lookup table in pyspextool to convert a spectral type
    (really spectralclass+np.floor(subtype)) to an effective temperature.

    Parameters
    ----------
    sptype : str
        A spectral type.  The first character must be one of 'OBFGK' an the
        second character must be one of '123456789'.

    Returns
    -------
    int

        An effective temperature.
    
    """

    #
    # Check parameters
    #

    check_parameter('sptype2teff', 'sptype', sptype, 'str')

    #
    # Read the look up table into memory
    #

    fullpath = osjoin(setup.state["package_path"], "data",
                      "sptype_effectivetemperature.dat")
    
    sptypes, teffs = np.loadtxt(fullpath, comments='#', unpack=True,
                                dtype='str')

    z = np.where(sptypes == sptype[0:2])[0].item()

    return float(teffs[z])

    
