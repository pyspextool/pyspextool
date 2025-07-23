import numpy as np
import logging
import matplotlib.pyplot as pl
import os

from pyspextool import config as setup
from pyspextool.merge import config
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.merge.core import merge_spectra
from pyspextool.merge.qaplots import plot_merges

def merge_orders(
    merge_apertures:int=None,
    scale_orders:bool=True,
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    To merge the orders.

    Parameters
    ----------
    merge_apertures : int, list, default None
        The aperture number or a list of aperture numbers to merge.  
        The apetures are index starting with 1.  

    scale_orders : {True, False}
        Set to True to scale the flux density of each additional order to match
        the flux density level of the merged spectrum.
        Set to False to merge orders with no scaling.

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
        Set variable config.state['merged_spectrum']

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

    check_parameter("merge_orders", "merge_apertures", merge_apertures, 
                    ["NoneType", "int", "list"])

    check_parameter("merge_orders", "scale_orders", scale_orders, 'bool')


    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Determine which apertures we are merging
    #
    
    possible_apertures = np.arange(1, config.state['napertures']+1)

    if merge_apertures is None:

        merge_apertures = possible_apertures

    else:

        # Convert the user input into a numpy array

        merge_apertures = np.array(merge_apertures)

        # Check to make sure they are requesting apertures that actually exist.

        result = np.all(np.isin(merge_apertures, possible_apertures))
        
        if result.item() is False:

            string = ', '.join(str(x) for x in list(possible_apertures))
            message = 'Only apertures '+string+' are avaiable for merging.'
            raise pySpextoolError(message)

    # Store the result

    config.state['merge_apertures'] = merge_apertures

    message = ' Merging apertures '+\
        ', '.join(str(x) for x in list(merge_apertures))+'.'
    logging.info(message)

                                                  
    #
    # Start the loop over requested apertures
    #
    
    merged_spectra = []
    for aperture in merge_apertures:

        # Get the first order

        idx = aperture-1
        merged_wavelength = config.state['rawspectra'][idx,0,:]
        merged_intensity = config.state['rawspectra'][idx,1,:]
        merged_uncertainty = config.state['rawspectra'][idx,2,:]
        merged_bitmask = config.state['rawspectra'][idx,3,:].astype(np.uint8)

        # Now add in other orders.

        for i in range(1,config.state['norders']):

            add_idx = (aperture-1) + i*config.state['napertures']

            add_wavelength = config.state['rawspectra'][add_idx,0,:]
            add_intensity = config.state['rawspectra'][add_idx,1,:]
            add_uncertainty = config.state['rawspectra'][add_idx,2,:]
            add_bitmask = config.state['rawspectra'][add_idx,3,:].astype(np.uint8)

            # Scale the add order?


            # Merge the results

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
                        
        # Create spectrum to store

        spectrum = np.vstack((merged_wavelength,merged_intensity, 
                              merged_uncertainty, merged_bitmask))

        merged_spectra.append(spectrum)

    #
    # Make QA plot
    #

    if qa['show'] is True:

        plot_merges(setup.plots['shifts'],
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    qa['showscale'],
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    config.state['xlabel'],
                    config.state['rawspectra'],
                    config.state['orders'],
                    1,
                    merged_spectra)
    
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False:

            pl.pause(1)
        
    if qa['write'] is True:

        plot_merges(None,
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    qa['showscale'],
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    config.state['xlabel'],
                    config.state['rawspectra'],
                    config.state['orders'],
                    1,
                    merged_spectra)
        
        pl.savefig(os.path.join(setup.state['qa_path'],
                                config.state['outputfile_root']+ \
                                '_merges' + \
                                setup.state['qa_extension']))
        pl.close()

    #
    # Store the results and set the done variable
    # 

    config.state['merged_spectrum'] = merged_spectra

    config.state["merge_done"] = False
        
