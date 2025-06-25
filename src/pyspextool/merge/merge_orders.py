import numpy as np
import logging

from pyspextool import config as setup
from pyspextool.merge import config as config
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.merge.core import merge_spectra

def merge_orders(verbose:bool=None,
                 qa_show:bool=None,
                 qa_showscale:float=None,
                 qa_showblock:bool=None,
                 qa_write:bool=None):

    """
    To merge the orders.

    Parameters
    ----------
    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].


    Returns
    -------
    None

    """

    #
    # Check to ensure previous steps completed.
    #

    if config.state['load_done'] is False:

        message = "Spectra have not been loaded.  Please run load_spectra.py."
        raise pySpextoolError(message)

    #
    # Check the parameters and keywords
    #

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Log the action
    #

    message = ' Merging the orders.'
    logging.info(message)    

    #
    # Start the merging
    #

    idx = config.state['merge_aperture']-1
    merged_wavelength = config.state['spectra'][idx,0,:]
    merged_intensity = config.state['spectra'][idx,1,:]
    merged_uncertainty = config.state['spectra'][idx,2,:]
    merged_bitmask = config.state['spectra'][idx,3,:].astype(np.uint8)
    
    for i in range(1,config.state['norders']):


        add_idx = (i*config.state['napertures']+ 
                   (config.state['merge_aperture']-1))

        add_wavelength = config.state['spectra'][add_idx,0,:]
        add_intensity = config.state['spectra'][add_idx,1,:]
        add_uncertainty = config.state['spectra'][add_idx,2,:]
        add_bitmask = config.state['spectra'][add_idx,3,:].astype(np.uint8)

        result = merge_spectra(merged_wavelength,
                               merged_intensity,
                               add_wavelength,
                               add_intensity,
                               anchor_uncertainty=merged_uncertainty,
                               anchor_bitmask=merged_bitmask,
                               add_uncertainty=add_uncertainty,
                               add_bitmask=add_bitmask)

        merged_wavelength = result['wavelength']
        merged_intensity = result['intensity']
        merged_uncertainty = result['uncertainty']
        merged_bitmask = result['bitmask']


    #
    # Create spectrum to store
    #

    spectrum = np.vstack((merged_wavelength,merged_intensity, 
                          merged_uncertainty, merged_bitmask))

    #
    # Store the results and set the done variable
    # 

    config.state['merged_spectrum'] = spectrum

    config.state["merge_done"] = False
        
