import logging
import numpy as np
from os.path import join as osjoin

from pyspextool.utils.interpolate import linear_interp1d
from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric.core import deconvolve_line
from pyspextool.telluric.core import make_instrument_profile

def get_kernels(*args:int | float,
                verbose:bool=None,
                qa_show:bool=None,
                qa_showscale:float | int=None,
                qa_showblock:bool=None,
                qa_write:bool=None):

    """
    Determines the convolution kernel.
    
    Parameters
    ----------
    args : int, float
        If given, it is tbe number of FWHMs around the line to do the
        deconvolution.
        If not given, then the instrument profile is loaded instead.

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
    
    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].
    
    Returns
    -------
    None
    Loads data into memory and writes QA plots to disk.  

        tc.state['method']
        tc.state['deconvolution_window']
        tc.state['rms_deviation']
        tc.state['max_deviation']
        tc.state['ew_scale']
        tc.state['kernels']
        tc.state['kernel_done']
    
    """
    
    #
    # Check the _done variables
    #

    if len(args) == 0:

        # We are doing the IP method

        tc.state['method'] = 'ip'

        if tc.state['load_done'] is False:

            message = "The spectra have not been loaded.  Please run "\
            "load_spectra.py."
            raise pySpextoolError(message)

        
    else:

        tc.state['method'] = 'deconvolution'        
              
        if tc.state['rv_done'] is False:

            message = "Radial velocity has not been measured.  Please run "\
            "get_radialvelocity.py."
            raise pySpextoolError(message)
        
    #
    # Check parameters and keywords
    #

    check_parameter('get_kernels', 'verbose', verbose, ['NoneType','bool'])
    
    check_parameter('get_kernels', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('get_kernels', 'qa_scale', qa_showscale,
                    ['NoneType','float','int'])
    
    check_parameter('get_kernels', 'qa_block', qa_showblock,
                    ['NoneType','bool'])
        
    check_parameter('get_kernels', 'qa_write', qa_write, ['NoneType','bool'])

    if tc.state['method'] == 'deconvolution':

        check_parameter('get_kernels', 'args', args[0], ['float', 'int'])
            
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    logging.info(" Telluric method = "+tc.state['method'])

    #
    # Create and load the kernels
    #
    
    if tc.state['method'] == 'deconvolution':

        logging.info(" Deconvolving line.")        

        # Get the deconvolution range

        line_range = tc.state['line_fwhm']*args[0]
        tc.state['deconvolution_window'] = \
            [tc.state['line_center']-line_range/2.,
             tc.state['line_center']+line_range/2.]

        tc.state['deconvolution_window'] = \
            np.array(tc.state['deconvolution_window'])
        
        #
        # Get QA set up
        #

        xlabel = tc.state['standard_hdrinfo']['LXLABEL'][0]
        title = tc.state['mode']+' Order '+\
            str(tc.state['normalization_order'])

        if qa['show'] is True:

        # Build the qashow_info dictionary.

            figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                           setup.plots['landscape_size'][1]*qa['showscale'])
        
            font_size = setup.plots['font_size']*qa['showscale']

        
            qashow_info = {'plot_number':setup.plots['deconvolution'],
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

            # Build the qafile_info dictionary.        

            fullpath = osjoin(setup.state['qa_path'],
                              tc.state['output_filename']+\
                              '_decon'+setup.state['qa_extension'])    
        
        # Build the qafile_info dictionary.        
        
            qafile_info = {'figure_size':setup.plots['landscape_size'],
                           'font_size':setup.plots['font_size'],
                    'spectrum_linewidth':setup.plots['zoomspectrum_linewidth'],
                           'spine_linewidth':setup.plots['spine_linewidth'],
                           'file_fullpath':fullpath,
                           'xlabel':xlabel,
                           'title':title}

            
        else:

            qafile_info = None
            
        #    
        # Do the deconvolution
        #
    
        result = deconvolve_line(tc.state['normalized_line_wavelength'],
                                 tc.state['normalized_line_flux'],
                                 tc.state['vega_wavelength']*\
                                 tc.state['1+z'],
                                 tc.state['vega_normalized_fluxdensity'],
                                 tc.state['deconvolution_window'],
                                 qashow_info=qashow_info,
                                 qafile_info=qafile_info,
                                 verbose=verbose)
        
        #
        # Store the result
        #

        tc.state['rms_deviation'] = result['rms_deviation']
        tc.state['max_deviation'] = result['max_deviation']
        tc.state['ew_scale'] = result['ew_scale']            

        #
        # Now generate the kernels for each order
        #

        logging.info(" Generating the kernels.")        
        
        kernels = []        
        for i in range(tc.state['standard_norders']):

        # Generate data wavelengths based on the kernel pixels
            
            data_wavelengths = result['data_pixels']*\
                tc.state['standard_dispersions'][i]
            
            # Now generate model wavelengths
            
            min = np.min(data_wavelengths)
            max = np.max(data_wavelengths)
            
            delta = max-min
            
            npixels = int(delta/tc.state['vega_dispersions'][i])
            
            # enforce oddity
            
            if npixels % 2 == 0:

                npixels += 1
            
            vega_wavelengths = np.arange(-npixels//2+1,npixels//2+1)*\
                tc.state['vega_dispersions'][i]
            
            # Interpolate the data kernel onto the vega wavelegths
            
            rkernel= linear_interp1d(data_wavelengths, result['kernel'],
                                     vega_wavelengths)

            # Set NaNs to zero
            
            znan = np.isnan(rkernel)
            rkernel[znan] = 0.0

            # Renormalize the kernel
            
            rkernel /= np.sum(rkernel)
            
            kernels.append(rkernel)

    if tc.state['method'] == 'ip':
            
        logging.info(" Generating the kernels.")        

        tc.state['rms_deviation'] = np.nan
        tc.state['max_deviation'] = np.nan
        tc.state['ew_scale'] = 1.0
                
        kernels = []        
        for i in range(tc.state['standard_norders']):

            # Get the min/max wavelengths of the data
        
            min = np.nanmin(tc.state['standard_spectra'][i,0,:])
            max = np.nanmax(tc.state['standard_spectra'][i,0,:])
            
            # Determine the Vega model dispersion over these wavelengths
            
            z = np.where((tc.state['vega_wavelength'] >= min) &
                         (tc.state['vega_wavelength'] <= max))[0]
            
            # Compute the dispersion
            
            dispersion = tc.state['vega_wavelength'][z] - \
                np.roll(tc.state['vega_wavelength'][z],1)
            
            # Determine the median dispersion, ignoring the first pixel
            
            vega_dispersion = np.median(dispersion[1:-1])
            
            # Figure out the number of pixels required.
            
            nkernel = np.round(10*tc.state['standard_fwhm'][i]/\
                               vega_dispersion).astype(int)
            
            # enforce oddity
            
            if nkernel % 2 == 0:

                nkernel += 1
            
            # Create x values
            
            x = np.arange(-1*(nkernel//2),nkernel//2+1)*vega_dispersion/\
                tc.state['standard_dispersions'][i]
            
            # Create the profile
            
            p = make_instrument_profile(x,tc.state['ip_coefficients'])
            
            kernels.append(p)
            
    # Store the results
        
    tc.state['kernels'] = kernels

    #
    # Generate the control points for the EW scale factors
    #

    control_points = []        
    for i in range(tc.state['standard_norders']):

        minw = tc.state['standard_wavelengthranges'][i][0]
        maxw = tc.state['standard_wavelengthranges'][i][1]

        zlines = np.where((tc.state['H_wavelengths'] >= minw) &
                          (tc.state['H_wavelengths'] <= maxw))[0]

        
        points = tc.state['H_wavelengths'][zlines]
        points = np.insert(points,0, minw)
        points = np.append(points, maxw)        

        scales = np.full_like(points, tc.state['ew_scale'])

        stack = np.stack((points, scales))
        control_points.append(stack)
        
    tc.state['control_points'] = control_points
        
    # Set the done variable
    
    tc.state['kernel_done'] = True
