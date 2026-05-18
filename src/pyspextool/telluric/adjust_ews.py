import numpy as np
import logging
import os

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.telluric.core import estimate_ewscales
from pyspextool.telluric.core import find_modelshift
from pyspextool.utils.loop_progress import loop_progress

def adjust_ews(
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    To adjust the EWs of absorption lines in the standard star spectrum.

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
    later

    """

    if tc.state['ew_scale_info'] is None:

        return


    #
    # Check the parameters and keywords
    #

    check_parameter('adjust_ews', 'verbose', 
                    verbose, ['NoneType','bool'])

    check_parameter('adjust_ews', 'qa_show', 
                    qa_show, ['NoneType','bool'])

    check_parameter('adjust_ews', 'qa_showscale', 
                    qa_showscale, ['NoneType','float','int'])

    check_parameter('adjust_ews', 'qa_showblock', 
                    qa_showblock, ['NoneType','bool'])
    
    check_parameter('adjust_ews', 'qa_write', 
                    qa_write, ['NoneType','bool'])

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)

    #
    # Check to see whether this is a SpeX prism spectrum. 
    #

    tc.state['model_pixelshift'] = 0
    if tc.state['pixel_shift_info'] is not None:

        # Yes.  We need to determine the shift in pixels of Vega relative to 
        # the data.  

        # Build the qafile_info dictionary.

        qafile_info = None

        if qa['write'] is True:

            qafile_info = setup.plots['qafile_info']
            qafile_info['file_fullpath'] = os.path.join(
                setup.state['qa_path'],
                tc.state['telluric_output_filename']+\
                '_modelshift'+setup.state['qa_extension'])

            qafile_info['xlabel'] = tc.state['latex_xlabel']

        # Get necessary information

        order =  tc.state['pixel_shift_info']['Order Number']
        wrange =  tc.state['pixel_shift_info']['Wavelength Range']
    
        z = np.where(tc.state['standard_orders'] == order)[0]
        standard_wavelength = np.squeeze(tc.state['standard_spectra'][z,0,:])
        standard_fluxdensity = np.squeeze(tc.state['standard_spectra'][z,1,:])
        
        # Find the shift

        result = find_modelshift(
            standard_wavelength,
            standard_fluxdensity,
            tc.state['vega_wavelength'],
            tc.state['vega_fluxdensity'],
            wrange,
            qafile_info=qafile_info)

        # Store the results and log

        tc.state['vega_pixelshift'] = float(result)

        logging.info(" Telluric standard model shift="+\
                     '%+.2f' % tc.state['vega_pixelshift']+" pixels.")

    # 
    # Now get set up for the adjustments
    #

    # Rename things for ease of calling.

    standard_bmag = tc.state['standard_bmag']
    standard_vmag = tc.state['standard_vmag']
    standard_rv = float(tc.state['standard_rv'])
    
    vega_wavelength = tc.state['vega_wavelength']
    vega_fluxdensity = tc.state['vega_fluxdensity']
    vega_continuum = tc.state['vega_continuum']
    vega_fitted_continuum = tc.state['vega_fitted_continuum']

    kernel = tc.state['kernel']
    ew_scale = float(tc.state['ew_scale'])
    
    # Log the action

    message = ' Adjusting the EWs of absorption lines in the Vega model.'
    logging.info(message)

    # Get the information needed for the adjustments

    adjustments_orders = tc.state['ew_scale_info']['Order Numbers']
    adjustments_info = tc.state['ew_scale_info']['Fits']

    #
    # Get QA set up.  
    #

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
            'xlabel':tc.state['latex_xlabel'],
            'title':''}
    
    else:
                
        qashow_info = None
                
    if qa['write'] is True:
                
        # Build the qafile_info dictionary.  The fullpath and title will
        # be filled in in the loop so each correction has a unique file name.
        
        qafile_info = {
            'figure_size':setup.plots['landscape_size'],
            'font_size':setup.plots['font_size'],
            'spectrum_linewidth':setup.plots['zoomspectrum_linewidth'],
            'spine_linewidth':setup.plots['spine_linewidth'],
            'file_fullpath':'',
            'xlabel':tc.state['latex_xlabel'],
            'title':''}
        
    else:
        
        qafile_info = None

    #
    # Start the loop over each standard order.  
    #

    loop_progress(-1,0,tc.state['standard_norders'])

    for i in range(tc.state['standard_norders']):


        # Does this order have any adjustments?

        z_adjustments = np.where(adjustments_orders==tc.state['standard_orders'][i])[0]

        if z_adjustments.size == 0:

            # No it does not.

            loop_progress(i,0,tc.state['standard_norders'])
            continue

        #
        # Yes it does.  Adjust the EWs.
        #

        # Get the standard spectrum and adjustment_info for that order

        standard_wavelength = np.squeeze(tc.state['standard_spectra'][i,0,:])
        standard_fluxdensity = np.squeeze(tc.state['standard_spectra'][i,1,:])
        atmosphere = np.squeeze(tc.state['atmospheric_transmission'][i,1,:])
        
        info = adjustments_info[z_adjustments.item()]

        # Start the loop over the number of adjustments for this order

        points = []
        labels = []
        scales = []
        
        for j in range(len(info)):

            # Get information to update qaplot_info and qafile_info

            title = tc.state['instrument_mode']+' Order '+\
                str(tc.state['standard_orders'][i])+', '+str(info[j]['Fit Range'][0])+\
                '-'+str(info[j]['Fit Range'][1])+' '+setup.state['lxunits']+\
                ', poly degree='+str(info[j]['Poly Degree'])+\
                ', tolerance='+str(info[j]['Fit Tolerance'])

            if qa['show'] is True: 

                qashow_info['title'] = title
            
            if qa['write'] is True: 

                qafile_info['title'] = title
                
                filelabel = tc.state['instrument_mode']+'Order'+\
                    str(tc.state['standard_orders'][i])+'_'+\
                    str(info[j]['Fit Range'][0])+\
                    '-'+str(info[j]['Fit Range'][1])+setup.state['xunits']
                
                fullpath = os.path.join(setup.state['qa_path'],
                                  tc.state['telluric_output_filename']+\
                                  '_ewadjustments_'+\
                                  filelabel+\
                                  setup.state['qa_extension'])    

                qafile_info['file_fullpath'] = fullpath

            # Pull out the particular range for this order and estimate

            z = np.where((standard_wavelength > info[j]['Fit Range'][0]) & 
                         (standard_wavelength < info[j]['Fit Range'][1]))[0]

            
            result = estimate_ewscales(
                standard_wavelength[z],
                standard_fluxdensity[z],
                standard_rv,
                standard_vmag,
                standard_bmag,
                vega_wavelength,
                vega_fluxdensity,
                vega_continuum,
                vega_fitted_continuum,
                kernel,
                tc.state['vega_pixelshift'],
                ew_scale,
                np.array(info[j]['Lines']),
                atmosphere[z],
                info[j]['Poly Degree'],
                tolerance=info[j]['Fit Tolerance'],
                qashow_info=qashow_info,
                qafile_info=qafile_info)

            # Store the results

            new = ['start',*info[j]['Line IDs'],'end']
            labels.extend(new)

            new = [info[j]['Fit Range'][0],*info[j]['Lines'],info[j]['Fit Range'][1]]
            points.extend(new)

            new = [ew_scale,*list(result[2]),ew_scale]
            scales.extend(new)

            loop_progress(i,0,tc.state['standard_norders'])

        # Now add the edges of the orders and then finally sort

        labels.insert(0,'edge')
        labels.insert(-1,'edge')
        
        points.insert(0,tc.state["standard_wavelengthranges"][i][0])
        points.insert(-1,tc.state["standard_wavelengthranges"][i][1])
        
        scales.insert(0,ew_scale)
        scales.insert(-1,ew_scale)
        
        labels = np.array(labels)
        points = np.array(points)
        scales = np.array(scales)
        
        s = np.argsort(points)
        
        labels = labels[s]
        points = points[s]
        scales = scales[s]
        
        package = [list(points), list(labels), list(scales)]
        tc.state['control_points'][i] = package

    return tc.state['control_points']

