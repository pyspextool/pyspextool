import copy
import logging
import numpy as np

from pyspextool.telluric import config as tc
from pyspextool.utils.math import combine_flag_stack
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

        tc.state['corrected_spectra']
        tc.state['correct_done']
    
    """
    #
    # Check the make_done variable
    #
    
    if tc.state['make_done'] is False:

        message = "The telluric correction spectra have not been created.  " \
            "Please run make_telluric_spectra.py."
        raise pySpextoolError(message)
    
    logging.info(" Correcting spectra.")    

    #
    # Start the process
    #
    
    corrected_spectra = copy.deepcopy(tc.state['object_spectra'])

    # Start the loop over the object orders and apertures

    for i in range(tc.state['object_norders']):

        # Find the order in the standard spectrum

        z_order = np.where(tc.state['object_orders'][i] == \
                           tc.state['standard_orders'])

       # Now loop over the apertures
        
        for j in range(tc.state['object_napertures']):

            k =i*tc.state['object_napertures']+j
            
            tc_f = np.squeeze(tc.state['shiftedtc_spectra'][z_order,1,:])
            tc_u = np.squeeze(tc.state['shiftedtc_spectra'][z_order,2,:]) 
            tc_m = np.squeeze(tc.state['shiftedtc_spectra'][z_order,3,:])
            
            obj_f = tc.state['object_spectra'][k,1,:]
            obj_u = tc.state['object_spectra'][k,2,:]
            obj_m = tc.state['object_spectra'][k,3,:]

            # Do the correction and propagate uncertainties
            
            fd = obj_f*tc_f
            fd_u = np.sqrt(obj_f**2 * tc_u**2 + tc_f**2 * obj_u**2)

            # Combine the masks

            stack = np.stack((obj_m.astype(np.uint8),tc_m.astype(np.uint8)))
            fd_m = combine_flag_stack(stack)
                        
            # Store the results

            corrected_spectra[k,1,:] = fd
            corrected_spectra[k,2,:] = fd_u
            corrected_spectra[k,3,:] = fd_m

    # Store all the corrected_spectra
            
    tc.state['corrected_spectra'] = corrected_spectra
    
    tc.state['correct_done'] = True
