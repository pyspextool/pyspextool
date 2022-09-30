import numpy as np

from pyspextool.cl import config
from pyspextool.io.check_parameter import check_parameter
from pyspextool.spectroscopy.find_peaks import find_peaks

def define_aperture_positions(method, apertures, fwhm=0.8):

    """
    To determine the locations of spectral extraction apertures


    Parameters
    ----------
    method : {'auto', 'guess', fixed'}
        The method by which the apertures are identified (see Notes).

    apertures : int, list, or numpy.ndarray
        `method` = 'auto'
            the number of apertures to search
        `method` = 'guess'
            a list of numpy.darray of apertures positions to seaerch around.
        `method` = 'fixed'
            a list of numpy.darray of apertures positions.

    fwhm: float, default 0.8 (arcseconds).
        The approximate FWHM of the peak to be identified.  Only used 
        if `method` is 'auto' or 'guess'.


    Returns
    -------
    Bone


    Notes
    -----
    Just does some organizing and calls find_peaks.  

    Examples
    --------
    later
    
    """

    #
    # Check parameters
    #
    
    check_parameter('define_aperture_positions', 'method', method, 'str',
                    possible_values=['auto', 'fixed', 'guess'])

    check_parameter('define_aperture_positions', 'apertures', apertures,
                    ['int', 'list', 'ndarray'])        

    #
    # Go based off of the method
    

    if method == 'auto':

        apertures, apsigns = find_peaks(config.state['profiles'],
                                        {'method':'auto', 'peaks':apertures},
                                        fwhm=fwhm)

    if method == 'guess':

        norders = len(config.state['profiles'])
        naps = len(apertures)
        apertures = np.tile(apertures, (norders, 1))
        
        apertures, apsigns = find_peaks(config.state['profiles'],
                                        {'method':'guess', 'peaks':apertures},
                                        fwhm=fwhm)

    if method == 'fixed':

        norders = len(config.state['profiles'])
        naps = len(apertures)
        apertures = np.tile(apertures, (norders, 1))        

        apertures, apsigns = find_peaks(config.state['profiles'],
                                        {'method':'fixed', 'peaks':apertures},
                                        fwhm=fwhm)                

    #
    # Now check to see if we are doing an extended source extraction.  If so,
    # then set all aperture signs to 1.
    #
    
    if config.state['exttype'] == 'xs':

        # Set all aperture signs to positive

        apsigns = np.abs(apsigns)

    #
    # Determine the average apsign
    #

    if norders > 1:

        average_apsign = np.sum(apsigns, axis=0)/np.sum(np.abs(apsigns), axis=0)
        apsigns = np.empty(naps, dtype=int)

        for i in range(naps):

            apsigns[i] = 1 if average_apsign[i] > 0 else -1

    else:

        apsigns = np.squeeze(apsigns)
        
    #
    # Store the results into the config variable
    #
    
    config.state['apertures'] = apertures
    config.state['apsigns'] = apsigns
