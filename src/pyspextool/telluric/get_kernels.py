import logging
import numpy as np

from pyspextool.utils.interpolate import linear_interp1d
from pyspextool import config as setup
from pyspextool.telluric import config as telluric
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric.core import deconvolve_line
from pyspextool.telluric.core import make_instrument_profile

def get_kernels(*args:int | float,
                verbose:bool=None,
                qa_show:bool=None,
                qa_scale:float | int=None,
                qa_block:bool=None,
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
    Loads data into memory and writes QA plots to disk.  

        telluric.state['method']
        telluric.state['deconvolution_window']
        telluric.state['rms_deviation']
        telluric.state['max_deviation']
        telluric.state['ew_scale']
        telluric.state['kernels']
        telluric.state['kernel_done']
    
    """
    
    #
    # Check the _done variables
    #

    if len(args) == 0:

        # We are doing the IP method

        telluric.state['method'] = 'ip'

        if telluric.state['load_done'] is False:

            message = "The spectra have not been loaded.  Please run "\
            "load_spectra.py."
            raise pySpextoolError(message)

        
    else:

        telluric.state['method'] = 'deconvolution'        
              
        if telluric.state['rv_done'] is False:

            message = "Radial velocity has not been measured.  Please run "\
            "get_radialvelocity.py."
            raise pySpextoolError(message)
        
    #
    # Check parameters and keywords
    #

    check_parameter('get_kernels', 'verbose', verbose, ['NoneType','bool'])
    
    check_parameter('get_kernels', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('get_kernels', 'qa_scale', qa_scale,
                    ['NoneType','float','int'])
    
    check_parameter('get_kernels', 'qa_block', qa_block, ['NoneType','bool'])
        
    check_parameter('get_kernels', 'qa_write', qa_write, ['NoneType','bool'])

    if telluric.state['method'] == 'deconvolution':

        check_parameter('get_kernels', 'args', args[0], ['float', 'int'])
            
    keywords = check_keywords(verbose=verbose, qa_show=qa_show,
                              qa_scale=qa_scale, qa_block=qa_block,
                              qa_write=qa_write)


    logging.info(f" Telluric method = "+telluric.state['method'])

    #
    # Create and load the kernels
    #
    
    if telluric.state['method'] == 'deconvolution':

        logging.info(f" Deconvolving line... ")        

        # Get the deconvolution range

        line_range = telluric.state['line_fwhm']*args[0]
        telluric.state['deconvolution_window'] = \
            [telluric.state['line_center']-line_range/2.,
             telluric.state['line_center']+line_range/2.]

        telluric.state['deconvolution_window'] = \
            np.array(telluric.state['deconvolution_window'])
        
        #
        # Get QA set up
        #

        xlabel = telluric.state['standard_hdrinfo']['LXLABEL'][0]
        title = telluric.state['mode']+' Order '+\
            str(telluric.state['normalization_order'])

        if keywords['qa_show'] is True:

        # Build the qashow_info dictionary.

            qashow_info = {'number':telluric.plotwindows['deconvolution'],
                           'scale':keywords['qa_scale'],
                           'xlabel':xlabel,
                           'title':title,
                           'block':keywords['qa_block']}
            
        else:
        
            qashow_info = None

        if keywords['qa_write'] is True:

            # Build the qafile_info dictionary.        

            qafile_info = {'filepath':setup.state['qa_path'],
                           'filename':telluric.load['output_filename']+\
                           '_deconvolution',
                           'extension':setup.state['qa_extension'],
                           'xlabel':xlabel,
                           'title':title}
            
        else:

            qafile_info = None

        #    
        # Do the deconvolution
        #
    
        result = deconvolve_line(telluric.state['normalized_line_wavelength'],
                                 telluric.state['normalized_line_flux'],
                                 telluric.state['vega_wavelength']*\
                                 telluric.state['1+z'],
                                 telluric.state['vega_normalized_fluxdensity'],
                                 telluric.state['deconvolution_window'],
                                 qashow_info=qashow_info,
                                 qafile_info=qafile_info,
                                 verbose=verbose)
        
        # Did we show it?

        if keywords['qa_show'] is True:
            telluric.plotwindows['deconvolution'] = result['plot_number']

        #
        # Store the result
        #

        telluric.state['rms_deviation'] = result['rms_deviation']
        telluric.state['max_deviation'] = result['max_deviation']
        telluric.state['ew_scale'] = result['ew_scale']            

        #
        # Now generate the kernels for each order
        #

        logging.info(f" Generating the kernels... ")        
        
        kernels = []        
        for i in range(telluric.state['standard_norders']):

        # Generate data wavelengths based on the kernel pixels
            
            data_wavelengths = result['data_pixels']*\
                telluric.state['standard_dispersions'][i]
            
            # Now generate model wavelengths
            
            min = np.min(data_wavelengths)
            max = np.max(data_wavelengths)
            
            delta = max-min
            
            npixels = int(delta/telluric.state['vega_dispersions'][i])
            
            # enforce oddity
            
            if npixels % 2 == 0: npixels += 1
            
            vega_wavelengths = np.arange(-npixels//2+1,npixels//2+1)*\
                telluric.state['vega_dispersions'][i]
            
            # Interpolate the data kernel onto the vega wavelegths
            
            rkernel= linear_interp1d(data_wavelengths, result['kernel'],
                                     vega_wavelengths)

            # Set NaNs to zero
            
            znan = np.isnan(rkernel)
            rkernel[znan] = 0.0

            # Renormalize the kernel
            
            rkernel /= np.sum(rkernel)
            
            kernels.append(rkernel)

    if telluric.state['method'] == 'ip':
            
        logging.info(f" Generating the kernels... ")        

        telluric.state['rms_deviation'] = np.nan
        telluric.state['max_deviation'] = np.nan
        telluric.state['ew_scale'] = 1.0
                
        kernels = []        
        for i in range(telluric.state['standard_norders']):

            # Get the min/max wavelengths of the data
        
            min = np.nanmin(telluric.state['standard_spectra'][i,0,:])
            max = np.nanmax(telluric.state['standard_spectra'][i,0,:])
            
            # Determine the Vega model dispersion over these wavelengths
            
            z = np.where((telluric.state['vega_wavelength'] >= min) &
                         (telluric.state['vega_wavelength'] <= max))[0]
            
            # Compute the dispersion
            
            dispersion = telluric.state['vega_wavelength'][z] - \
                np.roll(telluric.state['vega_wavelength'][z],1)
            
            # Determine the median dispersion, ignoring the first pixel
            
            vega_dispersion = np.median(dispersion[1:-1])
            
            # Figure out the number of pixels required.
            
            nkernel = np.round(10*telluric.state['standard_fwhm'][i]/\
                               vega_dispersion).astype(int)
            
            # enforce oddity
            
            if nkernel % 2 == 0: nkernel += 1
            
            # Create x values
            
            x = np.arange(-1*(nkernel//2),nkernel//2+1)*vega_dispersion/\
                telluric.state['standard_dispersions'][i]
            
            # Create the profile
            
            p = make_instrument_profile(x,telluric.state['ip_coefficients'])
            
            kernels.append(p)
            
    # Store the results
        
    telluric.state['kernels'] = kernels

    # Set the done variable
    
    telluric.state['kernel_done'] = True
