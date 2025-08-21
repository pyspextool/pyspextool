import matplotlib.pyplot as pl
from astropy.io import fits
import copy
import logging
import numpy as np
import os
from astropy.coordinates import angular_separation

from pyspextool import config as setup
from pyspextool.telluric import config

from pyspextool.io.check import check_parameter, check_qakeywords, check_file
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io.files import inoutfiles_to_fullpaths
from pyspextool.io.files import make_full_path
from pyspextool.telluric.core import find_shift
from pyspextool.telluric.qaplots import plot_shifts
from pyspextool.plot.plot_spectra import plot_spectra
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.coords import ten
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.units import get_latex_fluxdensity

def correct_spectra(
    object_filenames:str | list,
    telluric_fullfilename:str,
    output_filenames:str,
    default_shiftranges:bool=True,
    user_shiftranges:tuple | list=None,
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float | int=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    Performs telluric correction and "flux calibration" on raw spectra.

    The function uses a previously created telluric correction/flux calibration 
    spectrum and applies it to raw data.  In addition, it shifts the correction 
    until noise in default/user defined wavelength ranges is minimized.

    Parameters
    ----------
    object_filenames : str or list
        If type is str, then a comma-separated string of full pyspextool file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits' or a single file 'spc00001.a.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.

        e.g. ['spectra', '1-2']

    telluric_fullfilename : str
        The full name (include the suffix '.fits') of the telluric correction file.

    output_filenames : str
        If type `object_filenames` is str, then a comma-separated string of the full 
        output file names, e.g. 'spectra00001, spectra00002'.  The number of input 
        files must equal the the number of output files.

        If type `object_filenames` is list, then the prefix, e.g. 'spectra'.  

    default_shiftranges : {True, False}
        Set to True to use the default ranges.
        Set to False to not use the default ranges.

    user_shiftranges : tuple, list, deafult None
        If a tuple, then a (3,) tuple as
        (order number, lower wavelength, upper wavelength).

        If a list, then a (norders,) list of (3,) tuples as
        (order number, lower wavelength, upper wavelength).

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
        the plot window which is (9,6).  This does NOT affect plots written
        to disk.  Set to None to default to setup.state['qa_scale'].    
      
    Returns
    -------
    None
        Writes QA plots to disk and a pyspextool FITS file.
    
    """

    #
    # Check the parameters and QA keywords
    #
    
    check_parameter('correct_spectra', 'object_filenames', object_filenames, 
                    ['str', 'list'], list_types=['str','str'])

    check_parameter('correct_spectra', 'telluric_fullfilename', telluric_fullfilename, 
                    'str')

    check_parameter('correct_spectra', 'output_filenames', output_filenames, 
                    ['str', 'list'])

    check_parameter('correct_spectra', 'default_shiftranges',
                    default_shiftranges, 'bool')

    check_parameter('correct_spectra', 'user_shiftranges', user_shiftranges,
                    ['tuple', 'list', 'NoneType'])    

    check_parameter('correct_spectra', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('correct_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)
    
    #
    # Log the process
    #
    
    logging.info(" Telluric Correction\n--------------------------\n")

    #
    # Clear the state variable
    #

    config.state.clear()

    #
    # Create the full input and output file names
    #

    result = inoutfiles_to_fullpaths(
        setup.state['proc_path'],
        object_filenames,
        setup.state['nint'],
        '',
        '.fits',
        setup.state['proc_path'],
        output_filenames)

    input_fullpaths = result['input_fullpaths']
    input_filenames = result['input_filenames']

    output_fullpaths = result['output_fullpaths']
    output_filenames = result['output_filenames']

    #
    # Now get the telluric correction file
    #

    logging.info(" Loading the telluric correction spectrum.\n")

    fullpath = make_full_path(setup.state["proc_path"], 
                              telluric_fullfilename, 
                              exist=True)

    telluric_spectra, telluric_dict = read_spectra_fits(fullpath)

    # Ensure the standard has only one aperture.

    if telluric_dict["napertures"] != 1:

        message = "The telluric file must have only one aperture."
        raise pySpextoolError(message)

    # Store results and useful information

    config.state['telluric_rawspectra'] = telluric_spectra
    config.state['telluric_dictionary'] = telluric_dict
    config.state['telluric_astropyheader'] = telluric_dict['astropyheader']
    config.state["telluric_orders"] = telluric_dict["orders"]

    #
    # Now start the loop over each `input_file`.
    #

    for i in range(len(input_fullpaths)):

        logging.info(' Loading file '+input_filenames[i]+".")
        
        object_spectra, object_dict = read_spectra_fits(input_fullpaths[i])

        # 
        # Check to ensure the telluric spectrum is compatible with the object.
        #

        airmass, angle = _check_telluric(object_dict['astropyheader'], 
                                         telluric_dict['astropyheader'])

        config.state['airmass_difference'] = airmass
        config.state['angular_separation'] = angle

        #
        # Store object information
        #

        config.state['object_spectra'] = object_spectra
        config.state['object_dictionary'] = object_dict
        config.state['object_astropyheader'] = object_dict['astropyheader']
        config.state['object_norders'] = object_dict['norders']
        config.state['object_orders'] = object_dict['orders']
        config.state['object_napertures'] = object_dict['napertures']
        config.state['instrument_mode'] = object_dict['obsmode']
        config.state['latex_xlabel'] = object_dict['lxlabel']

        #
        # Interpolate the telluric correction spectra onto the object wavelengths
        #

        telluric_spectra = np.copy(object_spectra)
        for j in range(config.state['object_norders']):

            tidx = np.where(config.state['telluric_orders'] == \
                            config.state['object_orders'][j])[0]

            for k in range(config.state['object_napertures']):
                
                idx = j*config.state['object_napertures'] + k        

                rtc_f, rtc_u = linear_interp1d(
                    config.state['telluric_rawspectra'][tidx,0,:],
                    config.state['telluric_rawspectra'][tidx,1,:],
                    config.state['object_spectra'][idx,0,:],
                    input_u=config.state['telluric_rawspectra'][tidx,2,:],)

                telluric_spectra[idx,1,:] = rtc_f
                telluric_spectra[idx,2,:] = rtc_u


        config.state['telluric_spectra'] = telluric_spectra
        config.state['shifted_telluric_spectra'] = np.copy(telluric_spectra)


        # 
        # Shift the telluric spectra to minimize residual telluric noise
        #

        if default_shiftranges is True or user_shiftranges is not None:

            _shift_spectra(default_shiftranges,
                           user_shiftranges,
                           output_filenames[i],
                           qa['verbose'],
                           qa['show'],
                           qa['showscale'],
                           qa['showblock'],
                           qa['write'])

        #
        # Now do the correction
        #

        message = " Correcting spectra for telluric absorption and flux calibrating."
        logging.info(message)

        #
        # Start the process
        #
    
        corrected_spectra = copy.deepcopy(config.state['object_spectra'])

        # Start the loop over the object orders and apertures

        for j in range(config.state['object_norders']):

            for k in range(config.state['object_napertures']):

                idx =j*config.state['object_napertures']+k
                
                tc_f = np.squeeze(config.state['shifted_telluric_spectra'][idx,1,:])
                tc_u = np.squeeze(config.state['shifted_telluric_spectra'][idx,2,:]) 
                tc_m = np.squeeze(config.state['shifted_telluric_spectra'][idx,3,:])
                
                obj_f = config.state['object_spectra'][idx,1,:]
                obj_u = config.state['object_spectra'][idx,2,:]
                obj_m = config.state['object_spectra'][idx,3,:]
                
                # Do the correction and propagate uncertainties
                
                fd = obj_f*tc_f
                fd_u = np.sqrt(obj_f**2 * tc_u**2 + tc_f**2 * obj_u**2)
                
                # Combine the masks
                
                stack = np.stack((obj_m.astype(np.uint8),tc_m.astype(np.uint8)))
                fd_m = combine_flag_stack(stack)
                
                # Store the results
                
                corrected_spectra[idx,1,:] = fd
                corrected_spectra[idx,2,:] = fd_u
                corrected_spectra[idx,3,:] = fd_m

        #
        #  Write the results to disk
        #

        output_astropyheader = copy.deepcopy(config.state['object_astropyheader'])

        # Add new keywords to the object header

        # Get telluric basics

        keywords = ['TC_SFILE', 'TC_STDID', 'TC_STDST', 'TC_TYPE', 'TC_METH']

        for keyword in keywords:

            new = (keyword, config.state['telluric_astropyheader'][keyword], \
                   config.state['telluric_astropyheader'].comments[keyword])
            output_astropyheader.append(new)
            
        # Update important keywords

        output_astropyheader['FILENAME'] = output_filenames[i]+'.fits'
        output_astropyheader['MODULE'] = 'telluric'
            

        # Deal with the units

        yunits = config.state['telluric_astropyheader']['YUNITS'].split('/')[0].strip()

        result = get_latex_fluxdensity(yunits)
        
        output_astropyheader['LYUNITS'] = result[0]
        output_astropyheader['LYLABEL'] = result[1]
        output_astropyheader['LuLABEL'] = result[2]

        fits.writeto(output_fullpaths[i]+'.fits',
                     corrected_spectra,
                     output_astropyheader,
                     overwrite=True)
        
        logging.info(' Wrote file '+output_filenames[i]+'.fits' + \
                     " to the proc directory.\n")

        #
        # Do the QA plotting
        #

        if qa['show'] is True:

            figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                           setup.plots['landscape_size'][1]*qa['showscale'])

            font_size = setup.plots['font_size']*qa['showscale']
        
            plot_spectra(
                output_fullpaths[i]+'.fits',
                ytype='flux and uncertainty',
                spectrum_linewidth=setup.plots['spectrum_linewidth'],
                spine_linewidth=setup.plots['spine_linewidth'],            
                title=output_filenames[i]+'.fits',
                showblock=qa['showblock'],
                plot_number=setup.plots['abeam_spectra'],
                figure_size=figure_size,
                font_size=font_size,
                colors=['green','black'])
            
            plot_spectra(
                output_fullpaths[i]+'.fits',
                ytype='snr',
                spectrum_linewidth=setup.plots['spectrum_linewidth'],
                spine_linewidth=setup.plots['spine_linewidth'],            
                title=output_filenames[i]+'.fits',
                showblock=qa['showblock'],
                plot_number=setup.plots['abeam_snr'],
                figure_size=figure_size,
                font_size=font_size,
                colors=['green','black'])
            

        if qa['write'] is True:
                
            file_fullpath = os.path.join(
                setup.state['qa_path'],
                output_filenames[i]+\
                setup.state['qa_extension'])
            
            plot_spectra(
                output_fullpaths[i]+'.fits',
                ytype='flux and uncertainty',
                spectrum_linewidth=setup.plots['spectrum_linewidth'],
                spine_linewidth=setup.plots['spine_linewidth'],            
                title=output_filenames[i]+'.fits',
                showblock=qa['showblock'],
                output_fullpath=file_fullpath,
                figure_size=setup.plots['landscape_size'],
                font_size=setup.plots['font_size'],
                colors=['green','black'])
            
            file_fullpath = os.path.join(
                setup.state['qa_path'],
                output_filenames[i]+\
                '_snr'+\
                setup.state['qa_extension'])
            
            plot_spectra(
                output_fullpaths[i]+'.fits',
                ytype='snr',
                spectrum_linewidth=setup.plots['spectrum_linewidth'],
                spine_linewidth=setup.plots['spine_linewidth'],            
                title=output_filenames[i]+'.fits',
                showblock=qa['showblock'],
                output_fullpath=file_fullpath,
                figure_size=setup.plots['landscape_size'],
                font_size=setup.plots['font_size'],
                colors=['green','black'])

            
def _shift_spectra(
    default_shiftranges:bool,
    user_shiftranges:tuple | list | None,
    output_filename:str,
    verbose:bool,
    qa_show:bool | None,
    qa_showscale:float | None,
    qa_showblock:bool | None,
    qa_write:bool | None):

    """
    To shift the telluric spectra to minimize residual telluric noise.
    
    The telluric spectrum is shifted relative to the object spectrum by sub pixel 
    shifts in order to minimize the noise over a particular wavelength range.

    Parameters
    ----------
    default_shiftranges : {True, False}
        Set to True to use the default ranges.
        Set to False to not use the default ranges.

    user_shiftranges : tuple, list, deafult None
        If a tuple, then a (3,) tuple as
        (order number, lower wavelength, upper wavelength).

        If a list, then a (norders,) list of (3,) tuples as
        (order number, lower wavelength, upper wavelength).
            
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
    Loads data into memory:

        config.state['shift_ranges']
        config.state['shifts']
        config.state['shiftedtc_spectra']
    
    """

    #
    # Check the parameters and keywords
    #

    check_parameter('shift_spectra', 'default_shiftranges',
                    default_shiftranges, 'bool')

    check_parameter('shift_spectra', 'user_shiftranges', user_shiftranges,
                    ['tuple', 'list', 'NoneType'])    
                   
    check_parameter('shift_spectra', 'verbose', verbose, 'bool')

    check_parameter('shift_spectra', 'qa_show', qa_show, 'bool')

    check_parameter('shift_spectra', 'qa_showscale', qa_showscale,
                    ['float','int'])

    check_parameter('shift_spectra', 'qa_showblock', qa_showblock, 'bool')
    
    check_parameter('shift_spectra', 'qa_write', qa_write, 'bool')

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    logging.info(" Shifting the telluric spectra to minimize residual telluric noise.")

    #
    # Create an empty shift_ranges array
    #

    shift_orders = config.state['object_orders']
    shift_ranges = np.full((config.state['object_norders'],2), np.nan)

    #
    # Now update with the default ranges if requested.
    #
    
    if default_shiftranges is True:

        # Read the telluric_shiftinfo.dat file.
        
        file = os.path.join(setup.state['instrument_path'],
                            'telluric_shiftinfo.dat')

        result = check_file(file, raise_error=False)

        if result is None:

            message = ' File telluric_shiftinfo.dat not found.  '+\
                'Skipping default shifts.'
            logging.info( message)

        else:
            
            mode, order, minwave, maxwave = np.loadtxt(file,
                                                       comments='#',
                                                       dtype='str',
                                                       unpack=True)
            
            order = order.astype(int)
            minwave = minwave.astype(float)
            maxwave = maxwave.astype(float)    
            
            #
            # yank out the pertinent information for the mode
            #
            
            z = mode == config.state['instrument_mode']
            if not np.sum(z):

                message = ' Default shift ranges for mode '+\
                    config.state['instrument_mode']+\
                    ' not found.  Skipping default shifts.'
                logging.info(message)
            
            default_orders = order[z]
            minwave = minwave[z]
            maxwave = maxwave[z]
        
            for i in range(len(shift_orders)):

                z = np.where(default_orders == shift_orders[i])[0]
                if len(z) == 1:
                    
                    shift_ranges[i,0] = minwave[z][0]
                    shift_ranges[i,1] = maxwave[z][0]                    
                                       
    #
    # Update the shifts_ranges for the user values
    #

    if isinstance(user_shiftranges,tuple):

        z = np.where(shift_orders == user_shiftranges[0])[0]
        if len(z) == 1:

            user_range = _check_inrange(user_shiftranges,
                                        verbose=qa['verbose'])
            if user_range is not None:

                shift_ranges[z,:] = user_range

    config.state['shift_ranges'] = shift_ranges
                
    #
    # Now do the shifts
    #

    shifts = np.zeros((config.state['object_norders'],
                       config.state['object_napertures']))

    for i in range(config.state['object_norders']):

        if np.isnan(shift_ranges[i][1]):

            continue
        
        for j in range(config.state['object_napertures']):

            idx = i*config.state['object_napertures'] + j

            shifts[i,j] = find_shift(config.state['object_spectra'][idx,0,:],
                                     config.state['object_spectra'][idx,1,:],
                                     config.state['telluric_spectra'][idx,1,:],
                                    list(shift_ranges[i][0:2]))

    config.state['shifts'] = np.round(shifts, decimals=2)

    #
    # Perform the shifts
    #

    for i in range(config.state['object_norders']):
        
        for j in range(config.state['object_napertures']):

            idx = i*config.state['object_napertures'] + j

            x = np.arange(len(config.state['object_spectra'][idx,0,:]))
            tc_f = config.state['telluric_spectra'][idx,1,:]
            tc_u = config.state['telluric_spectra'][idx,2,:]            

            rtc_f, rtc_u = linear_interp1d(
                x+config.state['shifts'][i,j].item(),
                tc_f,
                x,
                input_u=tc_u)
            
            config.state['shifted_telluric_spectra'][idx,1,:] = rtc_f
            config.state['shifted_telluric_spectra'][idx,2,:] = rtc_u

    #
    # Do the QA plots
    #

    if qa['show'] is True and np.sum(config.state['shifts']) != 0:

        plot_shifts(setup.plots['shifts'],
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    qa['showscale'],
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    config.state['latex_xlabel'],
                    config.state['object_orders'],
                    config.state['object_spectra'],
                    config.state['telluric_spectra'],
                    config.state['shifted_telluric_spectra'],
                    config.state['shift_ranges'],
                    config.state['shifts'])
        
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False:

            pl.pause(1)
        
    if qa['write'] is True  and np.sum(config.state['shifts']) != 0:

        plot_shifts(None,
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    1,
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    config.state['latex_xlabel'],
                    config.state['object_orders'],
                    config.state['object_spectra'],
                    config.state['telluric_spectra'],
                    config.state['shifted_telluric_spectra'],
                    config.state['shift_ranges'],
                    config.state['shifts'],
                    reverse_order=True)
        
        pl.savefig(os.path.join(setup.state['qa_path'],
                                output_filename+ \
                                '_shifts' + \
                                setup.state['qa_extension']))
        pl.close()

            
def _check_inrange(
    test:tuple,
    verbose:bool=True):

    """
    To determine if a requested user shift range falls within the order

    Parameters
    ----------
    test : tuple
        A (3,) tuple where test[0] is the order number, test[1] is the
        lower wavelength limit and test[2] is the upper wavelength limit.

    verbose : {True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
    
    Returns
    -------
    None, tuple
        If user wavelength range falls in the order, a numpy array with the
        wavelength range.

        If user wavelength range does not fall in the order, None.
    
    """

    #
    # Check parameters
    #

    check_parameter('_check_inrange', 'test', test, 'tuple')

    check_qakeywords(verbose=verbose)
    
    #
    # Do the check
    #
        
    # Rename the object wavelength ranges for ease.
    
    object_ranges = config.state['standard_wavelengthranges']
    
    z = np.where(config.state['object_orders'] == test[0])[0]

    # Does the requested order exist?
    
    if len(z) == 0:

        message = ' Requested telluric shift '+str(config.state['instrument_mode'])+\
            ' order '+str(test[0])+' not extracted.  Skipping this order.'
        logging.info(message)

        return None
           
    # Is the range in order?

    z = int(z)
    lower = np.logical_and(test[1] >= object_ranges[z][0], 
                           test[0] <= object_ranges[z][1])
    
    upper = np.logical_and(test[2] >= object_ranges[z][0], 
                           test[1] <= object_ranges[z][1])
    
    if not lower*upper:
        
        message = ' Requested telluric shift range for '+\
            str(config.state['instrument_mode'])+\
            ' Order '+str(test[0])+' of '+\
            ' and '.join([str(test[1]),str(test[2])])+\
            ' is out of range.  Skipping this order.'

        logging.info(message)

        return None
        
    return np.array(test[1:3])


def _check_telluric(object_astropyheader,
                    telluric_astropyheader):

    """
    To ensure the object and telluric files are compatible.

    The function checks that the observing mode and slit width are the same for the 
    telluric and object files and confirm that the telluric has the orders needed.  
    In addition, the angular and airmass separation are computed and returned.

    Parameters
    ----------
    object_astropyheader : Header
        The astropy header of the object spectra.

    telluric_astropyheader : Header
        The astropy header of the telluric correction spectra.

    Returns
    -------
    float : The angular separation between the object and standard.

    float : The difference in air mass (obj-std).
      
    """

    #
    # Check input parameters
    #

    check_parameter('_check_telluric', 'object_astropyheader', object_astropyheader,
                    'Header')

    check_parameter('_check_telluric', 'telluric_astropyheader', telluric_astropyheader,
                    'Header')

    #
    # Check to ensure object and telluric files are compatible
    #

    # Has the telluric been processed through the combine module?

    if telluric_astropyheader['MODULE'] != 'telluric':

        message = "The telluric file must have been processed through the telluric"+\
            " module."
        raise pySpextoolError(message)

    # Are the modes the same?

    if telluric_astropyheader['MODE'] != object_astropyheader['MODE']:

        message = "The telluric observing mode of "+telluric_astropyheader['MODE']+\
            " does not match the object observing mode of "+\
            object_astropyheader['MODE']+'.'
        raise pySpextoolError(message)
    
    # is the slitw_pix the same?
    
    if telluric_astropyheader['SLTW_ARC'] != object_astropyheader['SLTW_ARC']:
        
        message = "The telluric slit width of "+\
            str(telluric_astropyheader['SLTW_ARC'])+\
            " does not match the object slit width of "+\
            str(object_astropyheader['SLTW_ARC'])+'.'
        raise pySpextoolError(message)
    
    # Are all the orders in the object file present in the telluric file?
    
    object_orders = np.array([int(x) for x in 
                              object_astropyheader['ORDERS'].split(',')])
    telluric_orders = np.array([int(x) for x in 
                                telluric_astropyheader['ORDERS'].split(',')])


    intersection = np.intersect1d(object_orders, telluric_orders)
    
    if np.size(intersection) != np.size(object_orders):
        
        message = "The telluric file lacks an order the object file has."
        raise pySpextoolError(message)
    
    #
    # Now compute angular and airmass separationb
    #
    
    obj_long = 15 * ten(object_astropyheader["RA"], toradians=True)
    obj_lat = ten(object_astropyheader["DEC"], toradians=True)
    std_long = 15 * ten(telluric_astropyheader["RA"], toradians=True)
    std_lat = ten(telluric_astropyheader["DEC"], toradians=True)
    
    angle = (angular_separation(obj_long, obj_lat, std_long, std_lat) 
             * 180 / np.pi)

    # Do the airmass.  You must account for the fact that the object spectra
    # can be either a single spectrum with a single airmass or a combined spectrum 
    # with an average airmass.  

    if object_astropyheader['MODULE'] == 'extract':

        object_airmass = object_astropyheader['AM']

    elif object_astropyheader['MODULE'] == 'combine':

        object_airmass = object_astropyheader['AVE_AM']

    else:

        message = 'Module '+object_astropyheader['MODULE']+' unrecognized.'
        raise pySpextoolError(message)

    delta_airmass = object_airmass - telluric_astropyheader["AVE_AM"]


    return delta_airmass, angle
