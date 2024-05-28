import os
import numpy as np
import numpy.typing as npt
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table.table import Table
import logging
from astropy.io import fits
import typing
import matplotlib.pyplot as pl
import copy
from astropy.coordinates import SkyCoord
import astropy.units as u

from pyspextool import config as setup
from pyspextool.setup_utils import pySpextoolError
from pyspextool.telluric import config as telluric
from pyspextool.io.check import *
from pyspextool.io.files import *
from pyspextool.io.fitsheader import get_header_info
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.plot.plot_spectra import plot_spectra
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.interpolate import linear_bitmask_interp1d
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.split_text import split_text
from pyspextool.utils import units
from pyspextool.utils.spectra import normalize_continuum
from pyspextool.utils.spectra import model_xcorrelate
from pyspextool.telluric.telluric_utils import make_instrument_profile
from pyspextool.telluric.telluric_utils import deconvolve_line
from pyspextool.telluric.telluric_utils import make_telluric_spectrum
from pyspextool.utils.for_print import for_print
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.utils.units import convert_fluxdensity
from pyspextool.utils.units import get_latex_fluxdensity

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def telluric_correction(object_file:str,
                        standard_file:str,
                        standard_info:list | str | dict,
                        output_name:str,
                        fluxdensity_units:str='W m-2 um-1',
                        use_instrument_profile:bool=False,
                        reflectance:bool=False,
                        input_path:typing.Optional[str]=None,
                        output_path:typing.Optional[str]=None,
                        write_vegamodel:bool=False,
                        write_telluric:bool=True,                        
                        qa_path:typing.Optional[str]=None,
                        qa_show:typing.Optional[bool]=None,
                        qashow_scale:float=1.0,
                        qa_write:typing.Optional[bool]=None,
                        qawrite_extension:typing.Optional[str]=None,
                        qa_block:typing.Optional[bool]=None,
                        verbose:typing.Optional[bool]=None,
                        overwrite:bool=True):

    """
    To correct spectra for telluric absorption and flux calibrate

    Parameters
    ----------
    object_file : str 
        The name of the FITS file containing the object star spectra.

    standard_file : str 
        The name of the FITS file containing the standard star spectra.

    standard_info : str, list, dict
        If a string is passed, it is assumed to be the name of the standard.
        SIMBAD is queried for the spectral type and B- and V-band magnitudes.

        If a list is passed, it is assumed to contain the coordinates of the
        standard.  standard_info[0] = RA, standard_info[1] = Dec.  Each element
        is a string which gives the sexigesimal coordidates, separated by a
        blank space, :, or hms, dms.  SIMBAD is then queried for the id,
        spectral type, and B- and V-band magnitudes.

        If a dict is passed, it is assumed to contain the standard information.

        `"id"` : str
            The name of the standard.
    
        `"sptype"` : str
            The standard star spectral type.

        `"bmag"` : int or float
            The standard star B-band magnitude

        `"vmag"` : int or float
            The standard star V-band magnitude.
      
    output_name : str
        The output file name sans the suffix.

    fluxdensity_units : {'W m-2 um-1', 'erg s-1 cm-2 A-1', 'W m-2 Hz-1', 
                         'ergs s-1 cm-2 Hz-1', 'Jy', 'mJy', 'uJy'}
        The flux density units of the output.  

    use_instrument_profile : {False, True}
        Set to True to force use of the insrument profile.
    
    refelectence : {False, True} 
        Set to simply divide by the standard and return a relative 
        reflectance spectrum.

    standard_data: dict or None
        An optional dictionary giving the standard star information.  If given, 
        the program will not query SIMBAD for the information.

        `"sptype"` : str
            The standard star spectral type.

        `"bmag"` : int or float
            The standard star B-band magnitude

        `"vmag"` : int or float
            The standard star V-band magnitude.

    input_path : str or None
        An optional path with the data.  Otherwise the default is the proc/
        directory.

    output_path : str or None
        An optional output path.  Otherwise the default is the proc/
        directory.

    write_vegamodel : {False, True}
        Set to True to write the modified Vega model to disk.

    write_telluric : {False, True}
        Set to True to write the telluric correction spectrum to disk.
    
    qa_path : str or None
        An optional qa path.  Otherwise the default is the qa/
        directory.

    qa_show : {None, True, False}
        Set to True/False to override config.state['qa_show'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qashow_scale : tuple, default=(10, 6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}
        Set to True/False to override config.state['qa_write'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

    qawrite_extension : {None, '.pdf', '.png'}
        If qa_write is True, set this to override config.state['qa_extension']
        in the pyspextool config file.
    
    qa_block : {False, True}, optional
        Set to make the plot block access to the command line, e.g. pl.ioff().

    verbose : {None, True, False}
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file.  

    overwrite: {True, False}, optional
        Set to False to to not overwrite a file already on disk.

    Returns
    -------
    None
        Writes a pyspextool FITS file to disk.

    """

    #
    # Check the parameters
    #

    check_parameter('telluric_correction', 'object_file', object_file, 'str')

    check_parameter('telluric_correction', 'standard_file', standard_file,
                    'str')

    check_parameter('telluric_correction', 'standard_info', standard_info,
                    ['str','list','dict'])
        
    check_parameter('telluric_correction', 'output_name', output_name, 'str')

    check_parameter('telluric_correction', 'fluxdensity_units',
                    fluxdensity_units, 'str',
                    possible_values=['W m-2 um-1', 'erg s-1 cm-2 A-1',
                                     'W m-2 Hz-1', 'ergs s-1 cm-2 Hz-1',
                                     'Jy', 'mJy', 'uJy'])    

    check_parameter('telluric_correction', 'use_instrument_profile',
                    use_instrument_profile, 'bool')    
    
    check_parameter('telluric_correction', 'reflectance', reflectance, 'bool')
    
    check_parameter('telluric_correction', 'input_path', input_path,
                    ['NoneType', 'str'])

    check_parameter('telluric_correction', 'output_path', output_path,
                    ['NoneType', 'str'])

    check_parameter('telluric_correction', 'write_vegamodel', write_vegamodel,
                    'bool')

    check_parameter('telluric_correction', 'write_telluric', write_telluric,
                    'bool')
        
    check_parameter('telluric_correction', 'qa_path', qa_path,
                    ['NoneType', 'str'])            

    check_parameter('telluric_correction', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'qashow_scale', qashow_scale,
                    ['float','int'])
    
    check_parameter('telluric_correction', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'qawrite_extension',
                    qawrite_extension, ['NoneType', 'str'])
    
    check_parameter('telluric_correction', 'qa_block', qa_block,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'overwrite', overwrite, 'bool')    

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_write is None:

        qa_write = setup.state['qa_write']

    if qa_show is None:

        qa_show = setup.state['qa_show']

    if qa_block is None:

        qa_block = setup.state['qa_block']

    if verbose is None:
        verbose = setup.state['verbose']

    if overwrite is None:
        overwrite = setup.state['overwrite']

    if qawrite_extension is None:
        qawrite_extension = setup.state['qa_extension']
    
    if verbose is True:
        logging.getLogger().setLevel(logging.INFO)
        setup.state["verbose"] = True

    elif verbose is False:
        logging.getLogger().setLevel(logging.ERROR)
        setup.state["verbose"] = False

    # Get user paths if need be.

    if input_path is None:

        input_path = setup.state['proc_path']

    if output_path is None:

        output_path = setup.state['proc_path']

    if qa_path is None:

        qa_path = setup.state['qa_path']                

    check_path(input_path)
    check_path(output_path)
    check_path(qa_path)        
       
    #
    # Store user inputs
    #

    telluric.load['object_file'] = object_file
    telluric.load['standard_file'] = standard_file
    telluric.load['standard_info'] = standard_info    
    telluric.load['output_name'] = output_name
    telluric.load['fluxdensity_units'] = fluxdensity_units
    telluric.load['reflectance'] = reflectance
    telluric.load['input_path'] = input_path
    telluric.load['output_path'] = output_path
    telluric.load['write_vegamodel'] = write_vegamodel
    telluric.load['write_telluric'] = write_telluric
    telluric.load['qa_path'] = qa_path      
    telluric.load['qa_show'] = qa_show
    telluric.load['qashow_scale'] = qashow_scale
    telluric.load['qa_write'] = qa_write
    telluric.load['qawrite_extension'] = qawrite_extension
    telluric.load['qa_block'] = qa_block
    telluric.load['verbose'] = verbose
    telluric.load['overwrite'] = overwrite
    telluric.load['use_instrument_profile'] = use_instrument_profile

    #
    # Load the data
    #

    load_data()

    #
    # Start the process, first checking whether we are correcting a
    # solar system object or not.
    #

    if reflectance is False:

        #
        # Load the Vega model
        #

        load_vega()


        #
        # Normalize the order with the RV/deconvolution line

        normalize_order()

        #
        # Determine the radial velocity
        #

        get_radialvelocity()
    
        
        #
        # Get the convolution kernel
        #

        get_kernel()
        
        #
        # Load the kernels for each order
        #

        load_kernels()

    else:

        # Set things to that need to be set.
        
        telluric.state['standard_rv'] = 0
        telluric.state['max_deviation'] = 0
        telluric.state['rms_deviation'] = 0
            
    #
    # Make the telluric correction spectra
    #
        
    make_telluric_spectra()

    #
    # Correct the data for telluric absorption and flux calibrate
    #
    
    correct_spectra()

    #
    # Write files to disk
    #
    
    write_files()


    
def correct_spectra():

    """
    Corrects the object spectra with the telluric correction spectra

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    logging.info(f" Correcting spectra... ")    

    corrected_spectra = copy.deepcopy(telluric.state['object_spectra'])
    
    for i in range(telluric.state['object_norders']):

        # Find the order

        z_order = np.where(telluric.state['object_orders'][i] == \
                           telluric.state['standard_orders'])

       # Now loop over the apertures
        
        for j in range(telluric.state['object_napertures']):

            k =i*telluric.state['object_napertures']+j
            
            # Interpolate the correction spectrum 

            tc_w= np.squeeze(telluric.state['telluric_spectra'][z_order,0,:])
            tc_f = np.squeeze(telluric.state['telluric_spectra'][z_order,1,:])
            tc_u = np.squeeze(telluric.state['telluric_spectra'][z_order,2,:]) 
            tc_m = np.squeeze(telluric.state['telluric_spectra'][z_order,3,:])
            
            obj_w = telluric.state['object_spectra'][k,0,:]
            obj_f = telluric.state['object_spectra'][k,1,:]
            obj_u = telluric.state['object_spectra'][k,2,:]
            obj_m = telluric.state['object_spectra'][k,3,:]
            
            rtc_f, rtc_u = linear_interp1d(tc_w, tc_f, obj_w, input_u=tc_u)
            
            # Do the correction
            
            corrected_flux = obj_f*rtc_f

            corrected_var = obj_f**2 * rtc_u**2 + rtc_f**2 * obj_u**2

            # Interpolate the masks and combine

            mask = linear_bitmask_interp1d(tc_w, tc_m.astype(np.uint8), obj_w)

            stack = np.stack((obj_m.astype(np.uint8),mask))
            corrected_mask = combine_flag_stack(stack)

            # Store the results

            corrected_spectra[k,1,:] = corrected_flux
            corrected_spectra[k,2,:] = np.sqrt(corrected_var)
            corrected_spectra[k,3,:] = corrected_mask

    # Store all the corrected_spectra
            
    telluric.state['corrected_spectra'] = corrected_spectra


    
def get_kernel():

    """
    Determines the convolution kernel.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    #
    # Only run if the method is deconvolution
    #

    if telluric.state['method'] == 'ip':

        return

    logging.info(f" Deconvolving line... ")
            
    #
    # Get QA set up
    #

    xlabel = telluric.state['standard_hdrinfo']['LXLABEL'][0]
    title = telluric.state['mode']+' Order '+\
        str(telluric.state['normalized_order'])
    normalization_range = [np.min(telluric.state['normalization_windows']),
                           np.max(telluric.state['normalization_windows'])]




    if telluric.load['qa_show'] is True:

        # Build the qashow_info dictionary.

        qashow_info = {'plot_number':setup.plotwindows['telluric_deconvolution'],
                       'figure_scale':telluric.load['qashow_scale'],
                       'normalization_wavelength_range':normalization_range,
                       'plot_xlabel':xlabel,
                       'plot_title':title,
                       'block':telluric.load['qa_block']}

    else:

        qashow_info = None

    if telluric.load['qa_write'] is True:

        # Build the qafile_info dictionary.        

        qafile_info = {'filepath':telluric.load['qa_path'],
                       'filename':telluric.load['output_name']+'_deconv',
                       'extension':telluric.load['qawrite_extension'],
                       'normalization_wavelength_range':normalization_range,
                       'plot_title':title,
                       'plot_xlabel':xlabel}
    else:

        qafile_info = None

    #    
    # Do the deconvolution
    #
    
    result = deconvolve_line(telluric.state['normalized_order_wavelength'],
                             telluric.state['normalized_order_flux'],
                             telluric.state['vega_wavelength_rvshifted'],
                             telluric.state['vega_normalized_fluxdensity'],
                             telluric.state['deconvolution_window'],
                             qafile_info=qafile_info,
                             qashow_info=qashow_info,
                             verbose=setup.state['verbose'])
    
    # Did we show it?

    if telluric.load['qa_show'] is True:
        setup.plotwindows['telluric_deconvolution'] = result['plot_number']

    #
    # Store the result
    #

    telluric.state['kernel'] = result['kernel']
    telluric.state['pixels_kernel'] = result['object_pixels']
    telluric.state['rms_deviation'] = result['rms_deviation']
    telluric.state['max_deviation'] = result['max_deviation']    

    

def get_modeinfo(mode:str):

    """
    Reads the telluric_modeinfo file and returns values matching `mode`.

    Parameters
    ----------
    mode : str
        The instrument mode

    Returns
    -------
    

    """

    #
    # Check parameters
    #
    
    check_parameter('get_modeinfo', 'mode', mode, 'str')

    #
    # Get the file name and read it
    #
    
    file = os.path.join(setup.state['instrument_path'],'telluric_modeinfo.dat')

    values = np.loadtxt(file, comments='#', delimiter='|', dtype='str')

    # Parse the results
    
    modes = list(values[:,0])
    methods = list(values[:,1])
    models = list(values[:,2])
    orders = list(values[:,3])    
    normalization_windows = list(values[:,4])
    normalization_degrees = list(values[:,5])
    radialvelocity_windows = list(values[:,6])
    deconvolution_windows = list(values[:,7])


    #
    # Modify each vector accordingly.
    #

    modes = np.array([x.strip() for x in modes])
    methods = np.array([x.strip().lower() for x in methods])
    models = np.array([x.strip() for x in models])        
    orders = np.array([int(x) for x in orders])        

    normalization_windows = [x.split() for x in normalization_windows]
    for i in range(len(normalization_windows)):

        normalization_windows[i] = \
            np.array([float(x) for x in normalization_windows[i]])

    normalization_windows = np.array(normalization_windows)

    normalization_degrees = np.array([int(x) for x in normalization_degrees])
    
    radialvelocity_windows = [x.split() for x in radialvelocity_windows]
    for i in range(len(radialvelocity_windows)):

        radialvelocity_windows[i] = \
            np.array([float(x) for x in radialvelocity_windows[i]])

    radialvelocity_windows = np.array(radialvelocity_windows)

    deconvolution_windows = [x.split() for x in deconvolution_windows]
    for i in range(len(deconvolution_windows)):

        deconvolution_windows[i] = \
            np.array([float(x) for x in deconvolution_windows[i]])

    deconvolution_windows = np.array(deconvolution_windows)
      
    #
    # Find the matching mode and return results
    #

    z = np.where(modes == mode)

    return {'method':methods[z][0],
            'model':models[z][0], 
            'order':orders[z][0],
            'normalization_windows':normalization_windows[z][0],
            'normalization_degree':int(normalization_degrees[z][0]),
            'radialvelocity_window':radialvelocity_windows[z][0],
            'deconvolution_window':deconvolution_windows[z][0]}




def get_radialvelocity():

    """
    To determine the radial velocity of the A0 V star

    Parameters
    ----------
    None

    Returns
    -------
    None - Loads data into the config state variable.

    """

    logging.info(f" Computing radial velocity... ")

    #
    # Get QA set up
    #

    xlabel = telluric.state['standard_hdrinfo']['LXLABEL'][0]
    title = telluric.state['standard_name']+', '+\
        telluric.state['mode']+' Order '+\
        str(telluric.state['normalized_order'])


    if telluric.load['qa_show'] is True:

        # Build the qashow_info dictionary.

        qashow_info = {'plot_number':setup.plotwindows['telluric_rv'],
                       'plot_scale':telluric.load['qashow_scale'],
                       'block':telluric.load['qa_block'],
                       'plot_xlabel':xlabel,
                       'plot_title':title}

    else:

        qashow_info = None

    if telluric.load['qa_write'] is True:

        # Build the qafile_info dictionary.        

        qafile_info = {'filepath':telluric.load['qa_path'],
                       'filename':telluric.load['output_name']+'_rv',
                       'extension':telluric.load['qawrite_extension'],
                       'plot_xlabel':xlabel,
                       'plot_title':title}

    else:

        qafile_info = None

    # Compute the RV
        
    rv, z, num= model_xcorrelate(telluric.state['normalized_order_wavelength'],
                                telluric.state['normalized_order_flux'],
                                telluric.state['vega_wavelength'],
                                telluric.state['vega_normalized_fluxdensity'],
                             telluric.state['radialvelocity_window'][0].item(),
                             telluric.state['radialvelocity_window'][1].item(),
                            resolving_power=telluric.state['resolving_power'],
                                qashow_info=qashow_info,
                                qafile_info=qafile_info)


    logging.info(f" Radial velocity = "+\
                 '{:g}'.format(float('{:.4g}'.format(rv)))+" km s-1")

    #
    # Store the results
    #
    
    telluric.state['standard_rv'] = rv
    telluric.state['standard_redshift'] = z

    telluric.state['vega_wavelength_rvshifted'] = \
        telluric.state['vega_wavelength']*(1+z)

    # Did we show it?

    if telluric.load['qa_show'] is True:

        setup.plotwindows['telluric_rv'] = num

        
        
def load_data():

    """
    Loads the users requested data into memory

    Parameters
    ----------
    None

    Returns
    -------
    None - Loads data into the config state variable.

    """
    
    #
    # Load the files into memory
    #


    
    logging.info(f" Telluric Correction\n-------------------------\n")
    logging.info(f" Loading the data...")
            
    # Object first
    
    fullpath = make_full_path(telluric.load['input_path'],
                              telluric.load['object_file'], exist=True)

    object_spectra, object_info = read_spectra_fits(fullpath)

    header = fits.getheader(fullpath)

    object_hdrinfo = get_header_info(header,
                                     keywords=setup.state['telluric_keywords'])

    # Now the standard

    fullpath = make_full_path(telluric.load['input_path'],
                              telluric.load['standard_file'], exist=True)

    standard_spectra, standard_info = read_spectra_fits(fullpath)

    header = fits.getheader(fullpath)

    standard_hdrinfo = \
        get_header_info(header, keywords=setup.state['telluric_keywords'])

    
    #
    # Get standard star information from the user or SIMBAD
    #

    if isinstance(telluric.load['standard_info'], str):

        # user has passed the name of the standard

        telluric.load['standard_name'] = telluric.load['standard_info']
    
        # Get SIMBAD information of the standard

        if telluric.load['verbose'] is True:
            logging.info(f' Querying SIMBAD for standard star information...')
        
        Simbad.add_votable_fields('sptype', 'flux(B)', 'flux(V)')
        table = Simbad.query_object(telluric.load['standard_name'])

        if isinstance(table,Table):
            standard_name = table['MAIN_ID'][0]
            standard_sptype = table['SP_TYPE'][0]
            standard_vmag = float(table['FLUX_V'][0])
            standard_bmag = float(table['FLUX_B'][0])
            
        else:
            
            message = 'Standard name "{}"" was not found in SIMBAD; provide '\
                'the correct name, correct coordinates, or a dictionary '\
                'containing keys "id", "sptype", "bmag", '\
                'and "vmag".'.format(telluric.load['standard_name'])
            
            raise pySpextoolError(message)
        
    if isinstance(telluric.load['standard_info'], list):
    
        # user has passed the coordinates of the standard

        if telluric.load['verbose'] is True:
            logging.info(f' Querying SIMBAD for standard star information...')
        
        Simbad.add_votable_fields('id','sptype', 'flux(B)', 'flux(V)')

        c = SkyCoord(telluric.load['standard_info'][0],
                     telluric.load['standard_info'][1],
                     unit=(u.hourangle, u.deg))
        

        table = Simbad.query_region(c,radius='0d1m0s')

        if isinstance(table,Table):
            standard_name = table['MAIN_ID'][0]
            standard_sptype = table['SP_TYPE'][0]
            standard_vmag = float(table['FLUX_V'][0])
            standard_bmag = float(table['FLUX_B'][0])
            
        else:
            
            message = 'Standard coordiantes "{}"" were not found in SIMBAD; '\
                'provide the correct name, correct coordinates, or a '\
                'dictionary containing keys "id", "sptype", "bmag", '\
                'and "vmag".'.format(telluric.load['standard_info'])
            
            raise pySpextoolError(message)

    if isinstance(telluric.load['standard_info'], dict):

        # user has passed the standard information        

        standard_name = telluric.load['standard_info']['id']
        standard_sptype = telluric.load['standard_info']['sptype']
        standard_vmag = float(telluric.load['standard_info']['bmag'])
        standard_bmag = float(telluric.load['standard_info']['vmag'])
         
    #
    # Let's do some checks to ensure 1) the standard has only one aperture and
    # 2) the standard has the orders required.
    #

    if standard_info['napertures'] != 1:

        message = 'The standard has more than one aperture.'
        raise pySpextoolError(message)

    intersection = np.intersect1d(object_info['orders'],
                                  standard_info['orders'])

    if np.size(intersection) != object_info['norders']:

        message = 'The standard lacks an order the object has.'
        raise pySpextoolError(message)

    #
    # Store the results
    #

    telluric.state['standard_spectra'] = standard_spectra
    telluric.state['standard_info'] = standard_info
    telluric.state['standard_hdrinfo'] = standard_hdrinfo
    telluric.state['standard_name'] = standard_name
    telluric.state['standard_vmag'] = standard_vmag
    telluric.state['standard_bmag'] = standard_bmag
    telluric.state['delta_airmass'] = object_hdrinfo['AVE_AM'][0]-\
        standard_hdrinfo['AVE_AM'][0] 
    telluric.state['object_spectra'] = object_spectra
    telluric.state['object_info'] = object_info    
    telluric.state['object_hdrinfo'] = object_hdrinfo
    telluric.state['mode'] = standard_hdrinfo['MODE'][0]
    telluric.state['resolving_power'] = standard_hdrinfo['RP'][0]

    #
    # Get other information from the headers
    #

    telluric.state['object_norders'] = object_info['norders']
    telluric.state['object_napertures'] = object_info['napertures']    
    telluric.state['standard_norders'] = standard_info['norders']

    telluric.state['object_orders'] = object_info['orders']
    telluric.state['standard_orders'] = standard_info['orders']

    #
    # Compute the minimum and maximum wavelengths and dispersions for each
    # order (This may go away eventually once you implement mikeline_convolve).
    #
       
    pixels = np.arange(len(telluric.state['standard_spectra'][0,0,:]))
    dispersions = np.empty(telluric.state['standard_norders'])            
    wavelength_ranges = []
    
    for i in range(telluric.state['standard_norders']):

        # Do min max first

        min = np.nanmin(telluric.state['standard_spectra'][i,0,:])
        max = np.nanmax(telluric.state['standard_spectra'][i,0,:])

        wavelength_ranges.append(np.array([min,max]))

        # Now do the dispersion
                
        fit = poly_fit_1d(pixels,telluric.state['standard_spectra'][i,0,:],1)
        dispersions[i] = fit['coeffs'][1]

    # Store the results

    telluric.state['standard_wavelengthranges'] = wavelength_ranges
    telluric.state['standard_dispersions'] = dispersions
    telluric.state['standard_fwhm'] = dispersions*standard_info['slitw_pix']

    #
    # Finally, get the mode info and store pertinent information
    #

    mode_info = get_modeinfo(telluric.state['mode'])

    telluric.state['method'] = mode_info['method']
    telluric.state['model'] = mode_info['model']
    telluric.state['normalized_order'] = mode_info['order']
    telluric.state['normalization_windows'] = mode_info['normalization_windows']
    telluric.state['normalization_degree'] = mode_info['normalization_degree']
    telluric.state['radialvelocity_window'] = mode_info['radialvelocity_window']
    telluric.state['deconvolution_window'] = mode_info['deconvolution_window']

    # Override type if requested

    if telluric.load['use_instrument_profile'] is True:

        telluric.state['method'] = 'ip'    
    
    #
    # Load the instrument profile parameters
    #

    # Get the file name
        
    file = os.path.join(setup.state['instrument_path'],'IP_coefficients.dat')

    # Read the file
        
    slitw_arc, c0, c1, c2 = np.loadtxt(file,comments='#', unpack=True)

    # Find the right slit width and store the coefficients
    
    z = np.where(slitw_arc == \
                 telluric.state['standard_hdrinfo']['SLTW_ARC'][0])[0]

    ip_coefficients = np.squeeze(np.array((c0[z],c1[z],c2[z])))

    telluric.state['ip_coefficients'] = ip_coefficients

    #
    # Clear out variables
    #

    telluric.state['rms_deviation'] = np.nan
    telluric.state['max_deviation'] = np.nan    


def load_vega():

    """
    Loads the proper Vega model 

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    logging.info(f' Loading Vega model...')

    #
    # Determine which Vega model to use and load
    #

    root = 'Vega'+telluric.state['model']+'.fits'
    file = os.path.join(setup.state['package_path'],'data',root)
    
    hdul = fits.open(file) 
    data  = hdul[1].data
    
    vega_wavelength = data['wavelength']
    vega_flux = data['flux density']
    vega_continuum = data['continuum flux density']
    vega_fitted_continuum = data['fitted continuum flux density']
    
    hdul.close()
            
    # Store the results

    telluric.state['vega_wavelength'] = vega_wavelength
    telluric.state['vega_fluxdensity'] = vega_flux
    telluric.state['vega_continuum'] = vega_continuum
    telluric.state['vega_fitted_continuum'] = vega_fitted_continuum
    telluric.state['vega_normalized_fluxdensity'] = vega_flux/\
        vega_fitted_continuum

    #
    # Compute the dispersions over the order wavelengths
    #

    dispersions = np.empty(telluric.state['standard_norders'])            
    for i in range(telluric.state['standard_norders']):

        zleft = (vega_wavelength > \
                 telluric.state['standard_wavelengthranges'][i][0])
        zright = (vega_wavelength < \
                  telluric.state['standard_wavelengthranges'][i][1])
        
        zselection = np.logical_and(zleft,zright)

        pixels = np.arange(np.sum(zselection))
        
        fit = poly_fit_1d(pixels,vega_wavelength[zselection],1)
        dispersions[i] = fit['coeffs'][1]

    telluric.state['vega_dispersions'] = dispersions
    

def load_kernels():


    """
    To calculate the convolution kernel for each order 

    Parameters
    ----------
    None

    Returns
    -------
    None
    
    """

    logging.info(f' Loading the convolution kernels...')
    
    kernels = []
    

    if telluric.state['method'] == 'deconvolution':

        for i in range(telluric.state['standard_norders']):

            # Generate data wavelengths based on the kernel pixels
            
            data_wavelengths = telluric.state['pixels_kernel']*\
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

            rkernel= linear_interp1d(data_wavelengths, telluric.state['kernel'],
                                     vega_wavelengths)
            
            znan = np.isnan(rkernel)
            rkernel[znan] = 0.0

            kernels.append(rkernel)

    if telluric.state['method'] == 'ip':
    
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
    


def make_telluric_spectra():

    """
    To make the telluric correction spectra 

    Parameters
    ----------
    None

    Returns
    -------
    None - Loads data into the config state variable.

    """

    logging.info(f' Making telluric correction spectra...')


    telluric_spectra = copy.deepcopy(telluric.state['standard_spectra'])

    
    if telluric.load['reflectance'] is False:

        vega_spectra = copy.deepcopy(telluric.state['standard_spectra'])
        vega_spectra[:,2,:] = np.nan
        vega_spectra[:,3,:] = 0
        
        for i in range(telluric.state['standard_norders']):    
            
            standard_wavelength = telluric.state['standard_spectra'][i,0,:]
            standard_fluxdensity = telluric.state['standard_spectra'][i,1,:]
            standard_uncertainty = telluric.state['standard_spectra'][i,2,:]
            
            standard_bmag = telluric.state['standard_bmag']
            standard_vmag = telluric.state['standard_vmag']
            standard_rv = telluric.state['standard_rv']                
            
            vega_wavelength = telluric.state['vega_wavelength']
            vega_fluxdensity = telluric.state['vega_fluxdensity']
            vega_continuum = telluric.state['vega_continuum']
            vega_fitted_continuum = telluric.state['vega_fitted_continuum']
            kernel = telluric.state['kernels'][i]
            
            result = make_telluric_spectrum(standard_wavelength,
                                            standard_fluxdensity,
                                            standard_uncertainty,
                                            standard_rv,
                                            standard_bmag,
                                            standard_vmag,
                                            vega_wavelength,
                                            vega_fluxdensity,
                                            vega_continuum,
                                            vega_fitted_continuum,
                                            kernel)
            
            #
            # Change units to those requested by the user
            #

            flux = result[0]
            unc = result[1]
            vega = result[2]
            
            flux = convert_fluxdensity(standard_wavelength,
                                       result[0],'um','erg s-1 cm-2 A-1',
                                       telluric.load['fluxdensity_units'])
            
            unc = convert_fluxdensity(standard_wavelength,
                                      result[1],'um','erg s-1 cm-2 A-1',
                                      telluric.load['fluxdensity_units'])
            
            vega = convert_fluxdensity(standard_wavelength,
                                       result[2],'um','erg s-1 cm-2 A-1',
                                       telluric.load['fluxdensity_units'])
            
            # Updates labels
            
            telluric_spectra[i,1,:] = flux
            telluric_spectra[i,2,:] = unc
            vega_spectra[i,1,:] = vega
            
        #
        # Store the results
        #

        telluric.state['telluric_spectra'] = telluric_spectra
        telluric.state['vega_spectra'] = vega_spectra    


    else:

        for i in range(telluric.state['object_norders']):    
            
            standard_fluxdensity = telluric.state['standard_spectra'][i,1,:]
            standard_uncertainty = telluric.state['standard_spectra'][i,2,:]

        
            telluric_spectra[i,1,:] = 1/standard_fluxdensity
            telluric_spectra[i,2,:] = 1/standard_fluxdensity**2 * \
                standard_uncertainty
            
        #
        # Store the results
        #

        telluric.state['telluric_spectra'] = telluric_spectra

        

     
def normalize_order():

    """
    To normalize an order for radial velocity/deconvolution operations

    Parameters
    ----------
    None

    Returns
    -------
    None - Loads data into the config state variable.

    """

    #
    # log the operation
    #
    
    logging.info(f" Normalizing continuum in order "+\
                 str(telluric.state['normalized_order'])+"...")

    #
    # Get set up for the normalization
    #
    
    # Find the order given the modeinfo file

    z_order = np.where(telluric.state['standard_orders'] == \
                       telluric.state['normalized_order'])

    # Store values in shorter variable names for ease
    
    wavelength = np.squeeze(telluric.state['standard_spectra'][z_order,0,:])
    flux = np.squeeze(telluric.state['standard_spectra'][z_order,1,:])
    windows = telluric.state['normalization_windows']
    degree = telluric.state['normalization_degree']    
    robust = {'threshold':5, 'epsilon':0.1}
    xlabel = telluric.state['standard_hdrinfo']['LXLABEL'][0]
    title = telluric.state['standard_name']+', '+\
        telluric.state['mode']+' Order '+\
    str(telluric.state['normalized_order'])+', degree='+str(degree)
    
    #
    # Get QA set up
    #
    
    if telluric.load['qa_show'] is True:

        # Build the qashow_info dictionary.

        qashow_info = {'plot_number':setup.plotwindows['telluric_normalize'],
                       'plot_scale':telluric.load['qashow_scale'],
                       'block':telluric.load['qa_block'],
                       'plot_xlabel':xlabel,
                       'plot_title':title}

    else:

        qashow_info = None

    if telluric.load['qa_write'] is True:

        # Build the qafile_info dictionary.        

        qafile_info = {'filepath':telluric.load['qa_path'],
                       'filename':telluric.load['output_name']+'_normalize',
                       'extension':telluric.load['qawrite_extension'],
                       'plot_xlabel':xlabel,
                       'plot_title':title}

    else:

        qafile_info = None

        
    # Normalize the order
        
    nspectrum, plotnum = normalize_continuum(wavelength, flux, windows,
                                             degree, robust=robust,
                                             qashow_info=qashow_info,
                                             qafile_info=qafile_info)

    #
    # Store the results
    #

    telluric.state['normalized_order_wavelength'] = wavelength
    telluric.state['normalized_order_flux'] = nspectrum

    # Did we show it?

    if telluric.load['qa_show'] is True:

        setup.plotwindows['telluric_normalize'] = plotnum



        
def write_files():

    """
    To write FITS files to disk.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

#    if telluric.load['write_vegamodel'] is True:
#
#        x = 1
#
    if telluric.load['write_telluric'] is True:

        hdr = telluric.state['standard_hdrinfo']

        # Store the history

        old_history = hdr['HISTORY']

        # remove it from the avehdr

        hdr.pop('HISTORY')

        # Add things to it
        
        hdr['MODULE'][0] = 'telluric'
        
        hdr['FILENAME'][0] = telluric.load['output_name']+'_telluric.fits'
        
        hdr['STDFILE'] = [telluric.load['standard_file'],
                          'Telluric standard input file']
        
        hdr['STD_ID'] = [telluric.load['standard_name'],'Telluric standard']

        f = '{:.3f}'
        hdr['STD_BMAG'] = [float(f.format(telluric.state['standard_bmag'])),
                           'Telluric standard B mag']

        hdr['STD_VMAG'] = [float(f.format(telluric.state['standard_vmag'])),
                       'Telluric standard V mag']

        # Deal with the units and plot labels
    
        units = telluric.load['fluxdensity_units']+' / DN s-1'
        hdr['YUNITS'][0] = units

        lylabel = get_latex_fluxdensity(telluric.load['fluxdensity_units'])[0]

        hdr['LYLABEL'][0] = 'Ratio ('+lylabel+' / DN s$^{-1}$)'
        
        # Create the basic headers

        phdu = fits.PrimaryHDU()
        newhdr = phdu.header

        # Add our keywords
    
        keys = list(hdr.keys())
    
        for i in range(len(keys)):
            
            if keys[i] == 'COMMENT':
                
                junk = 1
                
            else:
                
                newhdr[keys[i]] = (hdr[keys[i]][0], hdr[keys[i]][1])
                
        # Do the history
        
        for hist in old_history:

            newhdr['HISTORY'] = hist

        #
        # Write the file out
        #

        full_path = os.path.join(telluric.load['output_path'],
                                 telluric.load['output_name']+'_telluric.fits')
    

        fits.writeto(full_path, telluric.state['telluric_spectra'], newhdr,
                 overwrite=telluric.load['overwrite'])
        
        logging.info(' Wrote file '+os.path.basename(full_path) + ' to disk.')
        
    
    #
    # Write the corrected spectra to disk
    #
    
    # Rename header for ease of reading
            
    hdr = telluric.state['object_hdrinfo']
                
    # Store the history

    old_history = hdr['HISTORY']

    # remove it from the avehdr

    hdr.pop('HISTORY')

    #
    # Add things to the header
    #

    hdr['MODULE'][0] = 'telluric'

    hdr['FILENAME'][0] = telluric.load['output_name']+'.fits'

    hdr['TC_OFILE'] = [telluric.load['object_file'],
                      'Telluric object input file']

    hdr['TC_SFILE'] = [telluric.load['standard_file'],
                      'Telluric standard input file']

    hdr['TC_STDID'] = [telluric.load['standard_name'],'Telluric standard']


    hdr['TC_STDB'] = [float('{:.3f}'.format(telluric.state['standard_bmag'])),
                       'Telluric standard B mag']

    hdr['TC_STDV'] = [float('{:.3f}'.format(telluric.state['standard_vmag'])),
                       'Telluric standard V mag']

    hdr['TC_STDRV'] =  [float('{:.2f}'.format(telluric.state['standard_rv'])),
                      'Telluric standard radial velocity (km s-1)']

    hdr['TC_dAM'] = [float('{:.2f}'.format(telluric.state['delta_airmass'])),\
                       'Telluric Average airmass difference']    

    hdr['TC_METH'] = [telluric.state['method'],'Telluric method']

    hdr['TC_MXDEV'] = [float('{:.5f}'.format(telluric.state['max_deviation'])),
                       'Telluric maxmimum % deviation of Vega-data']

    hdr['TC_RMS'] = [float('{:.5f}'.format(telluric.state['rms_deviation'])),
                       'Telluric RMS deviation of Vega-data']

    # Deal with the units and plot labels
    
    units = telluric.load['fluxdensity_units']
    hdr['YUNITS'][0] = units

    latex = get_latex_fluxdensity(units) 

    hdr['LYUNITS'][0] = latex[0]    
    hdr['LYLABEL'][0] = latex[1]
    
    #
    # Create the header
    #

    # Create the basic headers

    phdu = fits.PrimaryHDU()
    newhdr = phdu.header

    # Add our keywords
    
    keys = list(hdr.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            newhdr[keys[i]] = (hdr[keys[i]][0], hdr[keys[i]][1])
    
    # Do the history

    for hist in old_history:

        newhdr['HISTORY'] = hist

    #
    # Write the file out
    #

    full_path = os.path.join(telluric.load['output_path'],
                             telluric.load['output_name']+'.fits')
    
    fits.writeto(full_path, telluric.state['corrected_spectra'], newhdr,
                 overwrite=telluric.load['overwrite'])

    logging.info(' Wrote file '+os.path.basename(full_path) + ' to disk.')
        
    #
    # Do the QA plotting
    #

    if telluric.load['qa_show'] is True:

        plot_spectra(full_path, ytype='flux and uncertainty',
                     line_width=0.5,
                     title=os.path.basename(full_path))

        
    if telluric.load['qa_write'] is True:

        qafileinfo = {'figsize': (6*telluric.load['qashow_scale'],
                                  4*telluric.load['qashow_scale']),
                      'filepath': setup.state['qa_path'],
                      'filename': telluric.load['output_name'],
                      'extension': setup.state['qa_extension']}

        plot_spectra(full_path, ytype='uncertainty',
                     order_numbers=False,
                     line_width=0.5,colors=['green','black'],
                     title=os.path.basename(full_path), file_info=qafileinfo)
        









    

