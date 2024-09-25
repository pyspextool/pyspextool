import logging
import numpy as np
from math import floor

from pyspextool import config as setup
from pyspextool.telluric import config as telluric
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils.arrays import find_index
from pyspextool.telluric.core import measure_linerv

def get_radialvelocity(fwhm_scale:float,
                       resolving_power:bool | float |int =None, 
                       verbose:bool=None,
                       qa_show:bool=None,
                       qa_scale:float | int=None,
                       qa_block:bool=None,
                       qa_write:bool=None):

    """
    To determine the radial velocity of the A0 V star using a model.

    Parameters
    ----------
    fwhm_scale : float
        The number of FWHMs of the line around which to do the cross correlation

    resolving_power : int, float, NoneType
        If given, an attempt is made to smooth the model to the resolving
        power of the data before the cross correlation is run.  
        The default is None and no smoothing is attempted.
    
    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_block : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_scale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].
    
    Returns
    -------
    None
    Loads data into the config state variable.

        telluric.load['rv_nfwhm']
        telluric.plotwindows['radial_velocity']
        telluric.state['ew_scale'] 
        telluric.state['standard_rv']
        telluric.state['standard_redshift']
        telluric.state['1+z']
        telluric.state['rv_done']
    
    """

    #
    # Check the load_done variable
    #
    
    if telluric.state['prepare_done'] is False:

        message = "Line has not been prepared.  Please run prepare_line.py."
        raise pySpextoolError(message)
    
    #
    # Check parameters and keywords
    #

    check_parameter('get_radialvelocity', 'fwhm_scale', fwhm_scale,
                    ['int','float'])

    check_parameter('get_radialvelocity', 'resolving_power', resolving_power,
                    ['float','bool','NoneType'])
    
    check_parameter('get_radialvelocity', 'verbose', verbose,
                    ['NoneType','bool'])

    check_parameter('get_radialvelocity', 'qa_show', qa_show,
                    ['NoneType','bool'])

    check_parameter('get_radialvelocity', 'qa_scale', qa_scale,
                    ['NoneType','float','int'])
    
    check_parameter('get_radialvelocity', 'qa_block', qa_block,
                    ['NoneType','bool'])    
    
    check_parameter('get_radialvelocity', 'qa_write', qa_write,
                    ['NoneType','bool'])

    keywords = check_keywords(verbose=verbose, qa_show=qa_show,
                              qa_scale=qa_scale, qa_block=qa_block,
                              qa_write=qa_write)

    #
    # Log the action
    #
    
    logging.info(f" Computing radial velocity... ")

    telluric.load['rv_nfwhm'] = fwhm_scale
    
    #
    # Determine the wavelength window over which to measure the radial velocity
    #
    
    range = telluric.state['line_fwhm']*fwhm_scale
    window = np.array([telluric.state['line_center']-range/2.,
                      telluric.state['line_center']+range/2.])
    
    # Now ensure the range doesn't go beyond what is available.

    wavelength = telluric.state['normalized_line_wavelength']

    window[0] = max(window[0],np.nanmin(wavelength))
    window[1] = min(window[1],np.nanmax(wavelength))                   
    
    telluric.state['rv_window'] = window
        
    #
    # Smooth model to approximate resolving power of line
    #
    
    # Figure out resolving power at the wavelength in question

    if resolving_power is None:
    
        delta = wavelength-np.roll(wavelength,
                                   floor(telluric.state['slitw_pix']))
        resolving_powers = wavelength/delta
    
        line = np.sum(np.array(telluric.state['normalization_window']))/2
        
        idx = floor(find_index(wavelength, telluric.state['line_center']))
        resolving_power = floor(resolving_powers[idx])

    if resolving_power is not False:

        # 

        # Find the dispersion of the Vega model at this wavelength
        
        idx = int(find_index(telluric.state['vega_wavelength'],
                             telluric.state['line_center']))

        model_dispersion = telluric.state['vega_wavelength'][idx]- \
            telluric.state['vega_wavelength'][idx-1] 
        
        # Determine the number of pixels for the resolving power.

        fwhm_kernel = int(telluric.state['line_center']/resolving_power/ \
                          model_dispersion)

        # Create a Gaussian kernel

        npixels_kernel = fwhm_kernel*5
        if npixels_kernel % 2 == 0:

            npixels_kernel +=1

        x = np.arange(npixels_kernel)-npixels_kernel//2
        
        sigma = fwhm_kernel/2.354

        gaussian = np.exp(-(x/sigma)**2 / 2)
        gaussian /= np.sum(gaussian)

        model = np.convolve(telluric.state['vega_normalized_fluxdensity'],
                            gaussian, mode='same')

    else:

        model = telluric.state['vega_normalized_fluxdensity']
        
        
    #
    # Get QA set up
    #

    xlabel = telluric.state['standard_hdrinfo']['LXLABEL'][0]
    title = telluric.state['standard_name']+', '+\
        telluric.state['mode']+' Order '+\
        str(telluric.state['normalization_order'])

    if keywords['qa_show'] is True:

        # Build the qashow_info dictionary.

        qashow_info = {'number':telluric.plotwindows['radial_velocity'],
                       'scale':keywords['qa_scale'],
                       'block':keywords['qa_block'],
                       'xlabel':xlabel,
                       'title':title}

    else:

        qashow_info = None

    if keywords['qa_write'] is True:

        # Build the qafile_info dictionary.        

        qafile_info = {'filepath':setup.state['qa_path'],
                       'filename':telluric.load['output_filename']+'_rv',
                       'extension':setup.state['qa_extension'],
                       'xlabel':xlabel,
                       'title':title}

    else:

        qafile_info = None

    #
    # Compute the RV
    #

    # Locate the region for the object and Vega

    zd = np.where((telluric.state['normalized_line_wavelength'] >= window[0]) &
                  (telluric.state['normalized_line_wavelength'] <= \
                   window[1]))[0]
    
    zm = np.where((telluric.state['vega_wavelength'] >= window[0]) & 
                  (telluric.state['vega_wavelength'] <= window[1]))[0]
    
    result = measure_linerv(telluric.state['normalized_line_wavelength'][zd],
                            telluric.state['normalized_line_flux'][zd],
                            telluric.state['vega_wavelength'][zm],
                            model[zm],
                            qashow_info=qashow_info,
                            qafile_info=qafile_info)

    # Did we show it?

    if keywords['qa_show'] is True:

        telluric.plotwindows['radial_velocity'] = result['plotnum']

    logging.info(f" Radial velocity = "+\
                 '{:g}'.format(float('{:.4g}'.format(result['rv'])))+" km s-1")

    #
    # Store the results
    #

    telluric.state['ew_scale'] = 1.0
    telluric.state['standard_rv'] = result['rv']
    telluric.state['standard_redshift'] = result['z']
    telluric.state['1+z'] = (1+result['z'])
        
    #
    # Set the done variable
    #
        
    telluric.state['rv_done'] = True
