import logging
import numpy as np
from math import floor
from os.path import join as osjoin

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils.arrays import find_index
from pyspextool.telluric.core import measure_linerv
from pyspextool.pyspextoolerror import pySpextoolError

def measure_radialvelocity(
    *args:int | float,
    resolving_power:bool | float |int =None, 
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float | int=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    To determine the radial velocity of the A0 V star using a model.

    Parameters
    ----------
    args : int, float
        If given, the number of FWHMs of the line to do the cross correlation
        If not given, then the default value is used.

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
    Later

    """

    #
    # Setup output dictionary
    #

    output = {}
    
    #
    # Check the load_done variable
    #
    
    if tc.state['normalize_done'] is False:

        message = "Line has not been normalized.  Please run normalize_line.py."
        raise pySpextoolError(message)
    
    #
    # Did the user pass a value for number of FWHMs around line center?
    #

    if len(args) == 0:

        # No.  We are using the default value

        output['deconvolution_rvfit_nfwhm'] = [tc.state['deconvolution_rvfit_nfwhm'][0],
                                               tc.state['deconvolution_rvfit_nfwhm'][1]]

    else:

        # Yes.  Check and store it.

        check_parameter('measure_radialvelocity', args[0], 
                        args[0], 'int')

        output['deconvolution_rvfit_nfwhm'] = [tc.state['deconvolution_rvfit_nfwhm'][0],
                                               args[0]]

    #
    # Check the remaining parameters and keywords
    #

    check_parameter('get_radialvelocity', 'resolving_power', 
                    resolving_power, ['float','bool','NoneType'])
    
    check_parameter('get_radialvelocity', 'verbose', 
                    verbose, ['NoneType','bool'])

    check_parameter('get_radialvelocity', 'qa_show', 
                    qa_show, ['NoneType','bool'])

    check_parameter('get_radialvelocity', 'qa_showscale', 
                    qa_showscale, ['NoneType','float','int'])
    
    check_parameter('get_radialvelocity', 'qa_showblock', 
                    qa_showblock, ['NoneType','bool'])    
    
    check_parameter('get_radialvelocity', 'qa_write', 
                    qa_write, ['NoneType','bool'])

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)
    
    #
    # Log the action and setup outout
    #
    
    logging.info(" Computing radial velocity.")



#    tc.state['rv_nfwhm'] = fwhm_scale
    
    #
    # Determine the wavelength window over which to measure the radial velocity
    #
    
    wrange = tc.state['normalizedline_fwhm']*\
        output['deconvolution_rvfit_nfwhm'][1]
    window = np.array([tc.state['normalizedline_linecenter']-wrange/2.,
                      tc.state['normalizedline_linecenter']+wrange/2.])
    
    # Now ensure the range doesn't go beyond what is available.

    wavelength = tc.state['normalizedline_wavelengths']

    window[0] = max(window[0],np.nanmin(wavelength))
    window[1] = min(window[1],np.nanmax(wavelength))                   
    
    output['deconvolution_rvwindow'] = window
        
    #
    # Smooth model to approximate resolving power of line
    #
    
    # Figure out resolving power at the wavelength in question

    if resolving_power is None:
    
        delta = wavelength-np.roll(wavelength,
                                   floor(tc.state['slitw_pix']))
        resolving_powers = wavelength/delta
    
        idx = floor(find_index(wavelength, tc.state['normalizedline_linecenter']))
        resolving_power = floor(resolving_powers[idx])

    if resolving_power is not False:

        # 

        # Find the dispersion of the Vega model at this wavelength
        
        idx = int(find_index(tc.state['vega_wavelength'],
                             tc.state['normalizedline_linecenter']))

        model_dispersion = tc.state['vega_wavelength'][idx]- \
            tc.state['vega_wavelength'][idx-1] 
        
        # Determine the number of pixels for the resolving power.

        fwhm_kernel = int(tc.state['normalizedline_linecenter']/resolving_power/ \
                          model_dispersion)

        # Create a Gaussian kernel

        npixels_kernel = fwhm_kernel*5
        if npixels_kernel % 2 == 0:

            npixels_kernel +=1

        x = np.arange(npixels_kernel)-npixels_kernel//2
        
        sigma = fwhm_kernel/2.354

        gaussian = np.exp(-(x/sigma)**2 / 2)
        gaussian /= np.sum(gaussian)

        model = np.convolve(tc.state['vega_normalized_fluxdensity'],
                            gaussian, mode='same')

    else:

        model = tc.state['vega_normalized_fluxdensity']
        
        
    #
    # Get QA set up
    #

    xlabel = tc.state['latex_xlabel']
    title = tc.state['standard_name']+', '+\
        tc.state['instrument_mode']+' Order '+\
        str(tc.state['deconvolution_line_order'])

    if qa['show'] is True:

        # Build the qashow_info dictionary.

        figure_size = (setup.plots['portrait_size'][0]*qa['showscale'],
                       setup.plots['portrait_size'][1]*qa['showscale'])
        
        font_size = setup.plots['font_size']*qa['showscale']

        
        qashow_info = {
            'plot_number':setup.plots['radial_velocity'],
            'figure_size':figure_size,
            'font_size':font_size,
            'spectrum_linewidth':setup.plots['zoomspectrum_linewidth'],
            'spine_linewidth':setup.plots['spine_linewidth'],
            'block':qa['showblock'],
            'xlabel':xlabel,
            'title':title}

    else:

        qashow_info = None

    if qa['write'] is True:

        fullpath = osjoin(setup.state['qa_path'],
                          tc.state['telluric_output_filename']+\
                          '_rv'+setup.state['qa_extension'])    
        
        # Build the qafile_info dictionary.        

        qafile_info = {
            'figure_size':setup.plots['portrait_size'],
            'font_size':setup.plots['font_size'],
            'spectrum_linewidth':setup.plots['zoomspectrum_linewidth'],
            'spine_linewidth':setup.plots['spine_linewidth'],
            'file_fullpath':fullpath,
            'xlabel':xlabel,
            'title':title}

    else:

        qafile_info = None
        
    #
    # Compute the RV
    #

    # Locate the region for the object and Vega

    zd = np.where((tc.state['normalizedline_wavelengths'] >= window[0]) &
                  (tc.state['normalizedline_wavelengths'] <= \
                   window[1]))[0]
    
    zm = np.where((tc.state['vega_wavelength'] >= window[0]) & 
                  (tc.state['vega_wavelength'] <= window[1]))[0]
    
    # Do the fit

    result = measure_linerv(
        tc.state['normalizedline_wavelengths'][zd],
        tc.state['normalizedline_fluxes'][zd],
        tc.state['vega_wavelength'][zm],
        model[zm],
        qashow_info=qashow_info,
        qafile_info=qafile_info)

    # Report the value
    
    logging.info(" Radial velocity = "+\
                 '{:g}'.format(float('{:.4g}'.format(result['rv'])))+" km s-1")

    #
    # Store the results
    #

    output['standard_rv'] = result['rv']
    output['standard_redshift'] = result['z']
    output['1+z'] = (1+result['z'])

    
    # Update the telluric config file

    keys = list(output.keys())
    vals = list(output.values())

    for i in range(len(keys)):

        tc.state[keys[i]] = vals[i]
        
    #
    # Set the done variable and return output
    #
        
    tc.state['rv_done'] = True

    return output
