import logging
import numpy as np
from os.path import join as osjoin

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric.core import deconvolve_line
from pyspextool.telluric.core import make_instrument_profile

def make_kernel(
    *args:int | float,
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float | int=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    Build  the convolution kernel.
    
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
    Later

    """

    #
    # Setup output dictionary
    #

    output = {}

    #
    # Check parameters and keywords sans args
    #

    check_parameter('build_kernel', 'verbose', 
                    verbose, ['NoneType','bool'])
    
    check_parameter('build_kernel', 'qa_show', 
                    qa_show, ['NoneType','bool'])

    check_parameter('build_kernel', 'qa_scale', 
                    qa_showscale, ['NoneType','float','int'])
    
    check_parameter('build_kernel', 'qa_block', 
                    qa_showblock, ['NoneType','bool'])
        
    check_parameter('build_kernel', 'qa_write', 
                    qa_write, ['NoneType','bool'])

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)
    
    logging.info(" Telluric method = "+tc.state['telluric_method'])

    #
    # Start the process by asking what kind of IP
    #

    if tc.state['telluric_method'] == 'deconvolution':

        # Check previous steps

        if tc.state['rv_done'] is False:

            message = "Radial velocity has not been measured.  Please run "\
            "measure_radialvelocity.py."
            raise pySpextoolError(message)

        # Did the user pass for number of FWHMs around line center?

        if len(args) == 0:

            # No. Use the default value

            output['deconvolution_nfwhm'] = [tc.state['deconvolution_nfwhm'][0],
                                             tc.state['deconvolution_nfwhm'][1]]
            
        else:

            # Yes.  Check and store.

            check_parameter('build_kernel', 'args', 
                            args[0], ['int'])

            output['deconvolution_nfwhm'] = [tc.state['deconvolution_nfwhm'][0],
                                             args[0]]

        # Get on with it.

        logging.info(" Deconvolving line.")        

        # Get the deconvolution range
        
        line_range = tc.state['normalizedline_fwhm']*\
            output['deconvolution_nfwhm'][1]
        
        window = np.array([tc.state['normalizedline_linecenter']-line_range/2.,
                           tc.state['normalizedline_linecenter']+line_range/2.])

        output['deconvolution_window'] = window
        
        #
        # Get QA set up
        #

        xlabel = tc.state['latex_xlabel']
        title = tc.state['instrument_mode']+' Order '+\
            str(tc.state['deconvolution_line_order'])

        if qa['show'] is True:

        # Build the qashow_info dictionary.

            figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                           setup.plots['landscape_size'][1]*qa['showscale'])
        
            font_size = setup.plots['font_size']*qa['showscale']
        
            qashow_info = {
                'plot_number':setup.plots['deconvolution'],
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

            fullpath = osjoin(
                setup.state['qa_path'],
                tc.state['telluric_output_filename']+\
                '_decon'+setup.state['qa_extension'])    
            
        # Build the qafile_info dictionary.        
        
            qafile_info = {
                'figure_size':setup.plots['landscape_size'],
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
    
        result = deconvolve_line(
            tc.state['normalizedline_wavelengths'],
            tc.state['normalizedline_fluxes'],
            tc.state['vega_wavelength']*\
            tc.state['1+z'],
            tc.state['vega_normalized_fluxdensity'],
            output['deconvolution_window'],
            qashow_info=qashow_info,
            qafile_info=qafile_info,
            verbose=verbose)
        
        #
        # Store the result
        #

        output['kernel'] = result['kernel']
        output['rms_deviation'] = result['rms_deviation']
        output['max_deviation'] = result['max_deviation']
        output['ew_scale'] = result['ew_scale']            

    if tc.state['telluric_method'] == 'IP':    
        
        logging.info(" Generating the kernel.")        

        output['rms_deviation'] = np.nan
        output['max_deviation'] = np.nan
        output['ew_scale'] = 1.0

        FWHM = np.ceil(tc.state["slitw_pix"])

        # enforce oddity
            
        if FWHM % 2 == 0:
            
            FWHM += 1
                
        x = np.arange(-1*(FWHM*4//2),FWHM*4//2+1)
            
        # Create the profile
            
        kernel = make_instrument_profile(x,tc.state['ip_coefficients'])

        output['kernel'] = kernel

    #
    # Generate the default control points for the EW scale factors
    #

    control_points = []
    info = tc.state['ew_scale_info']

    # Are there EWs to adjust?

    if info is None:

        # No.  Just create control points at the ends of the orders.

        for i in range(tc.state['standard_norders']):
        
            minw = tc.state['standard_wavelengthranges'][i][0]
            maxw = tc.state['standard_wavelengthranges'][i][1]

            points = np.array([minw,maxw])
            labels = np.array(['start','end'])
            scales = np.full_like(points, output['ew_scale'])

            package = [list(points), list(labels), list(scales)]
            control_points.append(package)

    else:

        # Yes.  Parse the user input to create control points. 
        
        orders = np.array(info['Order Numbers'])

        # Loop over each order

        for i in range(tc.state['standard_norders']):
        
            points = []
            labels = []
            
            minw = tc.state['standard_wavelengthranges'][i][0]
            maxw = tc.state['standard_wavelengthranges'][i][1]

            z = np.where(orders == tc.state['standard_orders'][i])[0]

            # Did the user pass something?

            if len(z) != 0:

                # Yes

                for fit in info['Fits'][z[0]]:

                    points.extend([fit['Fit Range'][0],
                                  *fit['Lines'],
                                  fit['Fit Range'][1]])

                    labels.extend(['edge',
                                   *fit['Line IDs'],
                                   'edge'])

            points = np.array([minw]+points+[maxw])
            labels = np.array(['start']+labels+['end'])
            scales = np.full_like(points, output['ew_scale'])
            
            s = np.argsort(points)
            points = points[s]
            labels = labels[s]

            package = [list(points), list(labels), list(scales)]
            control_points.append(package)

    output['control_points'] = control_points

    # Update the telluric config file

    keys = list(output.keys())
    vals = list(output.values())

    for i in range(len(keys)):

        tc.state[keys[i]] = vals[i]

    tc.state['kernel_done'] = True

    return output

