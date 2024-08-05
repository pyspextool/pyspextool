import copy
import logging
import numpy as np

from pyspextool import config as setup
from pyspextool.telluric import config as telluric
from pyspextool.telluric.core import correct_spectrum
from pyspextool.pyspextoolerror import pySpextoolError

def correct_spectra():

    """
    Corrects the object spectra with the telluric correction spectra

    Parameters
    ----------
    None

    Returns
    -------
    None
    Loads data into memory.

        telluric.state['corrected_spectra']
        telluric.state['correct_done']
    
    """

    #
    # Check the make_done variable
    #
    
    if telluric.state['make_done'] is False:

        message = "The telluric correction spectra have not been created.  " \
            "Please run make_telluric_spectra.py."
        raise pySpextoolError(message)
    
    logging.info(f" Correcting spectra... ")    

    corrected_spectra = copy.deepcopy(telluric.state['object_spectra'])

    for i in range(telluric.state['object_norders']):

        # Find the order

        z_order = np.where(telluric.state['object_orders'][i] == \
                           telluric.state['standard_orders'])

       # Now loop over the apertures
        
        for j in range(telluric.state['object_napertures']):

            k =i*telluric.state['object_napertures']+j
            
            # Interpolate the correction spectrum 

            tc_w= np.squeeze(telluric.state['telluric_spectra'][z_order,0,:])
            tc_f = np.squeeze(telluric.state['telluric_spectra'][z_order,1,:])
            tc_u = np.squeeze(telluric.state['telluric_spectra'][z_order,2,:]) 
            tc_m = np.squeeze(telluric.state['telluric_spectra'][z_order,3,:])
            
            obj_w = telluric.state['object_spectra'][k,0,:]
            obj_f = telluric.state['object_spectra'][k,1,:]
            obj_u = telluric.state['object_spectra'][k,2,:]
            obj_m = telluric.state['object_spectra'][k,3,:]

            fd, fd_u, fd_m = correct_spectrum(obj_w, obj_f, tc_w, tc_f,
                                              uncertainties=[obj_u, tc_u],
                                              masks=[obj_m, tc_m])
            
            # Store the results

            corrected_spectra[k,1,:] = fd
            corrected_spectra[k,2,:] = fd_u
            corrected_spectra[k,3,:] = fd_m

    # Store all the corrected_spectra
            
    telluric.state['corrected_spectra'] = corrected_spectra
    
    telluric.state['correct_done'] = True
