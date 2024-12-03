import numpy as np
from os.path import join as osjoin, basename as osbasename, splitext
import glob
from astropy.io import fits

from pyspextool.io.check import check_parameter
from pyspextool import config as setup

def load_atmosphere(resolving_power:int):

    """
    To load an atmospheric transmission curve.

    Parmameters
    -----------
    resolving_power : int
        The resolving power requested.

    Returns
    -------
    ndarray, ndarray
    
    """

    #
    # Check parameters
    #

    check_parameter('load_atmosphere','resolving_power', resolving_power, 'int')

    #
    # Start the process
    #
    
    # First we have to get the possible file names
    
    fullpath = glob.glob(osjoin(setup.state['package_path'],
                                'data', 'atran*.fits'))

    # Then strip the paths off
    
    basenames = [osbasename(x) for x in fullpath]
    
    # Now get the resolving powers
    
    rps = np.array([int(x[5:x.find('.')]) for x in basenames])

    # Find the closest one

    deltas = rps - resolving_power
    z = deltas == np.min(np.abs(deltas))
        
    # Load that file and return the results
    
    array = fits.getdata(np.array(fullpath)[z][0])
    
    return np.squeeze(array[0, :]), np.squeeze(array[1, :]) 
    
