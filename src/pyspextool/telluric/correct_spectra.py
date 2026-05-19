import matplotlib.pyplot as pl
from astropy.io import fits
import copy
import logging
import numpy as np
import os
from astropy.coordinates import angular_separation

from pyspextool import config as setup
from pyspextool.telluric import config

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io.files import inoutfiles_to_fullpaths
from pyspextool.io.files import make_full_path
from pyspextool.telluric.core import find_shift
from pyspextool.telluric.qaplots import plot_shifts
from pyspextool.plot.plot_spectra import plot_spectra
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io import read
from pyspextool.utils.coords import ten
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.units import get_latex_fluxdensity

def correct_spectra(
    object_filenames:str | list,
    telluric_fullfilename:str,
    output_filenames:str,
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

    user_shiftranges : tuple, list, default None
        If a (3,) tuple, then 
        (order number, lower wavelength, upper wavelength).

        If a (0,) tuple, then no shifts are applied.
  
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
    
    check_parameter('correct_spectra', 'object_filenames', 
                    object_filenames, ['str', 'list'], list_types=['str','str'])

    check_parameter('correct_spectra', 'telluric_fullfilename', 
                    telluric_fullfilename, 'str')

    check_parameter('correct_spectra', 'output_filenames', 
                    output_filenames, ['str', 'list'])

    check_parameter('correct_spectra', 'user_shiftranges', 
                    user_shiftranges, ['tuple', 'list', 'NoneType'])    

    check_parameter('correct_spectra', 'verbose', 
                    verbose, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_write', 
                    qa_write, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_show', 
                    qa_show, ['NoneType', 'bool'])

    check_parameter('correct_spectra', 'qa_showscale', 
                    qa_showscale, ['int', 'float', 'NoneType'])

    check_parameter('correct_spectra', 'qa_showblock', 
                    qa_showblock, ['NoneType', 'bool'])
    
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

    logging.info(" Loading the telluric correction file "+telluric_fullfilename+".")

    fullpath = make_full_path(
        setup.state["proc_path"], 
        telluric_fullfilename, 
        exist=True)
    
    telluric_spectra, telluric_dict = read_spectra_fits(fullpath)

    # Ensure the telluric file has only one aperture.

    if telluric_dict["napertures"] != 1:

        message = "The telluric file must have only one aperture."
        raise pySpextoolError(message)

    # Store results and useful information

    config.state['telluric_rawspectra'] = telluric_spectra
    config.state['telluric_dictionary'] = telluric_dict
    config.state['telluric_astropyheader'] = telluric_dict['astropyheader']
    config.state["telluric_orders"] = telluric_dict["orders"]
    config.state["telluric_norders"] = telluric_dict["norders"]
    config.state["instrument_mode"] = telluric_dict["obsmode"]

    #
    # Get the wavelength ranges for the telluric correction spectra
    #

    wavelength_ranges = []
    for i in range(telluric_dict["norders"]):

        # Do min max first

        min = np.nanmin(telluric_spectra[i, 0, :])
        max = np.nanmax(telluric_spectra[i, 0, :])

        wavelength_ranges.append(np.array([min, max]))

    # Store the results

    config.state["telluric_wavelengthranges"] = wavelength_ranges

    #
    # Get set up for minimizing residual telluric noise
    #

    logging.info(' Loading residuals telluric noise minimization ranges. ')
    shift_requests = _get_shiftinfo(        
        user_shiftranges)
    
    #
    # Now start the loop over each `input_file`.
    #

    for i in range(len(input_fullpaths)):

        logging.info(' Loading object file '+input_filenames[i]+".")
        
        object_spectra, object_dict = read_spectra_fits(input_fullpaths[i])

        # 
        # Check to ensure the telluric spectrum is compatible with the object.
        #

        airmass, angle = _check_telluric(
            object_dict['astropyheader'], 
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

            tidx = np.where(
                config.state['telluric_orders'] == \
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

        shifts = _shift_spectra(
            shift_requests,
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

        #
        #  Add the shift ranges and values
        #

        shifts = config.state['shifts']
        np.nan_to_num(shifts,copy=False,nan=0.0)
        for j in range(config.state['object_norders']):

            order = str(config.state['object_orders'][j]).zfill(3)
            keyword = 'SHFTO'+order
            string = ", ".join(map(str,shifts[j,:]))

            output_astropyheader[keyword] = (string, 
                                             ' shift of order '+order+\
                                             ' telluric spectrum (pixels).')

        fits.writeto(
            output_fullpaths[i]+'.fits',
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

                    
def _check_inrange(
    test_order,
    test_range:list,
    verbose:bool=True):

    """
    To determine if a requested user shift range falls within the order

    Parameters
    ----------
    test_order : int
        The order number associated with `test_range`.

    test_range : list
        A (2,) list where test_range[0] is the lower wavelength limit and 
        test[1] is the upper wavelength limit.

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

    check_parameter('_check_inrange', 'test_order', test_order, 'int')

    check_parameter('_check_inrange', 'test_range', test_range, 'list')

    check_qakeywords(verbose=verbose)
    
    #
    # Do the check
    #
        
    # Rename the object wavelength ranges for ease.
    
    telluric_ranges = config.state['telluric_wavelengthranges']
    
    z = np.where(config.state['telluric_orders'] == test_order)[0]
           
    # Is the range in order?

    z = int(z)
    lower = np.logical_and(test_range[1] >= telluric_ranges[z][0], 
                           test_range[0] <= telluric_ranges[z][1])
    
    upper = np.logical_and(test_range[1] >= telluric_ranges[z][0], 
                           test_range[0] <= telluric_ranges[z][1])
    
    if not lower*upper:
        
        message = ' Requested telluric shift range for '+\
            str(config.state['instrument_mode'])+\
            ' Order '+str(test_order)+' of '+\
            ' and '.join([str(test_range[0]),str(test_range[1])])+\
            ' is out of range.  Skipping this order.'

        logging.info(message)

        return None
        
    return list(test_range[0:2])


def _check_telluric(
    object_astropyheader,
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

    check_parameter('_check_telluric', 'object_astropyheader', 
                    object_astropyheader, 'Header')

    check_parameter('_check_telluric', 'telluric_astropyheader', 
                    telluric_astropyheader, 'Header')

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

def _get_shiftinfo(
    user_shifts:tuple | list | None,
    verbose:bool=None):

    """
    To merge the default shift information and user shift information

    Parameters
    ----------
    user_shifts : tuple | list | None
        If a (3,) tuple, then
        (order number, lower wavelength, upper wavelength).

        If a (0,) tuple, then no shifts

        If a list, then a (norders,) list of (3,) tuples as
        (order number, lower wavelength, upper wavelength).
   
    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].

    Returns
    -------
    list or None
        A (nshiftorders,) list of dictionaries:

        'Order Number' : int
            The order number.

        'Wavelength Range' : list
            A (2,) list giving the wavelength range over which to determine shifts.

    """

    #
    # Check parameters and QA
    #

    check_parameter('_merge_shiftinfo', 'user_shifts',
                    user_shifts, ['NoneType','list','tuple'])

    check_parameter('merge_shiftinfo', 'verbose', 
                    verbose, ['NoneType', 'bool'])

    qa = check_qakeywords(
        verbose=verbose)

    shifts = None
    
    if user_shifts is not None:

        # The user did pass shift information. 
        
        if len(user_shifts) == 0:

            # The user is requesting no shifts

            return shifts
    
        if isinstance(user_shifts,tuple) is True:

            # embed in a list for ease of looping.

            user_shifts = [user_shifts]

        # Loop over user shifts and create dictionary

        shifts = []
        for element in user_shifts:
                
            range = _check_inrange(
                element[0],
                [element[1],element[2]],
                verbose=qa['verbose'])
            
            if range is not None:
                    
                dict = {'Order Number':element[0],'Wavelength Range':range}
                shifts.append(dict)

    else:

        # Get and store the mode telluric information
        
        fullpath = os.path.join(
            setup.state["instrument_path"], 
            config.state['telluric_dictionary']['obsmode']+"_telluric.dat")
        
        result = read.read_telluric_file(fullpath)
        shifts = result['telluric_noise_info']


    return shifts


def _shift_spectra(
    shift_requests:list | None,
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
    shift_requests : dict, None
            
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

    check_parameter('shift_spectra', 'shift_requests', 
                    shift_requests, ['list', 'NoneType'])    
                   
    check_parameter('shift_spectra', 'verbose', 
                    verbose, 'bool')

    check_parameter('shift_spectra', 'qa_show', 
                    qa_show, 'bool')

    check_parameter('shift_spectra', 'qa_showscale', 
                    qa_showscale, ['float','int'])

    check_parameter('shift_spectra', 'qa_showblock', 
                    qa_showblock, 'bool')
    
    check_parameter('shift_spectra', 'qa_write', 
                    qa_write, 'bool')

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)

    #
    # Create an empty shift and ranges array
    #

    shifts = np.full(
        (config.state['object_norders'],config.state['object_napertures']),np.nan)

    ranges = np.full((config.state['object_norders'],2),np.nan)

    
    if shift_requests is None:

        config.state['shifted_telluric_spectra'] = config.state['telluric_spectra']
        config.state['shifts'] = shifts
        config.state['shift_ranges'] = ranges

        return

    #
    # Determine the shifts
    #

    message = " Determining subpixel shifts to minimize residual telluric noise."
    logging.info(message)

    # Start the loop over the object number

    order_requests = np.array([x['Order Number'] for x in shift_requests])

    for i in range(config.state['object_norders']):

        z = np.where(order_requests == config.state['object_orders'][i])[0]

        if z.size == 0: # no match

            continue

        else:

            z = z[0]

        ranges[i,:] = shift_requests[z]['Wavelength Range']

        for j in range(config.state['object_napertures']):
            
            idx = i*config.state['object_napertures'] + j
            
            shifts[i,j] = find_shift(
                config.state['object_spectra'][idx,0,:],
                config.state['object_spectra'][idx,1,:],
                config.state['telluric_spectra'][idx,1,:],
                shift_requests[z]['Wavelength Range'])
            
    config.state['shifts'] = np.round(shifts, decimals=2)
    config.state['shift_ranges'] = ranges
                    
    #
    # Perform the shifts
    #
        
    for i in range(config.state['object_norders']):

        for j in range(config.state['object_napertures']):

            if np.isnan(config.state['shifts'][i,j]):

                continue

            else:

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

    if qa['show'] is True:

        plot_shifts(
            setup.plots['shifts'],
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

        plot_shifts(
            None,
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
        
        pl.savefig(
            os.path.join(
                setup.state['qa_path'],
                output_filename+ \
                '_shifts' + \
                setup.state['qa_extension']))

        pl.close()

