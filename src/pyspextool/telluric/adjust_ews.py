import numpy as np
import logging
from os.path import join as osjoin

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords, check_file
from pyspextool.telluric.core import estimate_ewscales
from pyspextool.utils.loop_progress import loop_progress

def adjust_ews(verbose:bool=None,
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
    None
    Updates data into memory:
    
        tc.state['control_points']
    
    """

    #
    # Check the parameters and keywords
    #

    check_parameter('adjust_ews', 'verbose', verbose, ['NoneType','bool'])

    check_parameter('adjust_ews', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('adjust_ews', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])

    check_parameter('adjust_ews', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
    
    check_parameter('adjust_ews', 'qa_write', qa_write, ['NoneType','bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    #
    # Log the action
    #

    message = ' Adjusting the EWs of standard star absorption lines.'
    logging.info(message)

    # 
    # Load the telluric_ewadjustments.dat file
    #

    # Load the file

    file = osjoin(setup.state['instrument_path'],
                  'telluric_ewadjustments.dat')
    
    result = check_file(file, raise_error=False)

    if result is None:

        message = ' File telluric_ewadjustments.dat not found.  '
        logging.info( message)
        return

    result = np.loadtxt(file, 
                        delimiter='|', 
                        comments='#',
                        dtype='str',
                        ndmin=1)
    
    # Parse the file to deal with single-line files

    if result.ndim == 1:
        
        modes = np.char.strip(np.array([result[0]]).astype(str))
        orders = np.array([result[1]]).astype(int)
        ranges = np.array([np.array([result[2],result[3]]).astype(float)])
        degrees = np.array([result[4]]).astype(int)
        tolerances = np.array([result[5]]).astype(float)
        
    else:
        
        modes = np.char.strip(np.array(result[:,0]).astype(str))
        orders = np.array(result[:,1]).astype(int)
        ranges = np.zeros((modes.size,2))
        ranges[:,0] = result[:,2]
        ranges[:,1] = result[:,3]
        degrees = np.array(result[:,4]).astype(int)
        tolerances = np.array(result[:,5]).astype(float)
    
    #
    # Determine if the mode being corrected even requires adjustments
    #

    z = np.where(modes == tc.state['mode'])[0]

    if z.size == 0:

        # No modifications required for this mode.

        return 

    ntotal = z.size
    # 
    # It does.  So get set up
    #

    # Clip the values for this mode.

    orders = orders[z]
    ranges = ranges[z]
    degrees = degrees[z]
    tolerances = tolerances[z]

    # Rename things for ease of calling.

    standard_bmag = tc.state['standard_bmag']
    standard_vmag = tc.state['standard_vmag']
    standard_rv = float(tc.state['standard_rv'])
    
    vega_wavelength = tc.state['vega_wavelength']
    vega_fluxdensity = tc.state['vega_fluxdensity']
    vega_continuum = tc.state['vega_continuum']
    vega_fitted_continuum = tc.state['vega_fitted_continuum']

    ew_scale = float(tc.state['ew_scale'])
     
    # 
    # Get QA set up
    #
    
    xlabel = tc.state['standard_hdrinfo']['LXLABEL'][0]

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
                       'title':''}
            
    else:
        
        qashow_info = None

    if qa['write'] is True:
            
        # Build the qafile_info dictionary.  The fullpath and title will
        # be filled in in the loop so each correction has a unique file name.
        
        qafile_info = {'figure_size':setup.plots['landscape_size'],
                       'font_size':setup.plots['font_size'],
                     'spectrum_linewidth':setup.plots['zoomspectrum_linewidth'],
                       'spine_linewidth':setup.plots['spine_linewidth'],
                       'file_fullpath':'',
                       'xlabel':xlabel,
                       'title':''}
            
    else:

        qafile_info = None
               
    #
    # Loop over each order checking for adjustments
    #

    idx = 0
    for i in range(tc.state['standard_norders']):

        # Does this order have any adjustments?

        z_adjustments = np.where(orders == tc.state['standard_orders'][i])[0]

        if z_adjustments.size == 0:

            # No it does not.

            continue

        # Yes it does.

        subranges = ranges[z_adjustments]
        subdegrees = degrees[z_adjustments]
        subtolerances = tolerances[z_adjustments]

        n_adjustments=z_adjustments.size
        
        standard_wavelength = tc.state['standard_spectra'][i,0,:]
        standard_fluxdensity = tc.state['standard_spectra'][i,1,:]
        atmosphere = tc.state['atmospheric_transmission'][i,1,:]
        kernel = tc.state['kernels'][i]        

        control_points = tc.state['control_points'][i]
        
        for j in range(n_adjustments):

            # Get the H lines in the the requested range
        
            z_lines = np.where((control_points[0,:] > subranges[j][0]) & 
                               (control_points[0,:] < subranges[j][1]))[0]
            
            lines = control_points[0,z_lines]


            # Grab the correct portion of the standard spectrum

            z_standard = np.where((standard_wavelength > subranges[j][0]) & 
                         (standard_wavelength < subranges[j][1]))[0]

            # Get the title for the QA plot

            title = tc.state['mode']+' Order '+\
                str(tc.state['standard_orders'][i])+', '+str(subranges[j][0])+\
                '-'+str(subranges[j][1])+' '+setup.state['lxunits']+\
                ', poly degree='+str(subdegrees[j].item())+\
                ', tolerance='+str(subtolerances[j].item())
            
            if qa['show'] is True: 

                qashow_info['title'] = title
            
            if qa['write'] is True: 

                qafile_info['title'] = title
                
                label = tc.state['mode']+'Order'+\
                    str(tc.state['standard_orders'][i])+'_'+\
                    str(subranges[j][0])+\
                '-'+str(subranges[j][1])+setup.state['xunits']

                fullpath = osjoin(setup.state['qa_path'],
                                  tc.state['output_filename']+\
                                  '_ewadjustments_'+\
                                  label+\
                                  setup.state['qa_extension'])    

                qafile_info['file_fullpath'] = fullpath

            # Do the estimation

            result = estimate_ewscales(standard_wavelength[z_standard],
                                       standard_fluxdensity[z_standard],
                                       standard_rv,
                                       standard_vmag,
                                       standard_bmag,
                                       vega_wavelength,
                                       vega_fluxdensity,
                                       vega_continuum,
                                       vega_fitted_continuum,
                                       kernel,
                                       ew_scale,
                                       lines,
                                       atmosphere[z_standard],
                                       subdegrees[j].item(),
                                       include_edges=False,
                                       tolerance=subtolerances[j].item(),
#                                       tolerance=1.0,
                                       qashow_info=qashow_info,
                                       qafile_info=qafile_info)

            if qa['verbose'] is True:
                
                loop_progress(idx, 0, ntotal)

            idx += 1
            control_points[1,z_lines] = result[2]
           
        tc.state['control_points'][i] = control_points
