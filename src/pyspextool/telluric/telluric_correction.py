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

from pyspextool import config as setup
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

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def telluric_correction(object_file:str,
                        standard_name:str,
                        standard_file:str,
                        output_name:str,
                        fluxdensity_units:str='W m-2 um-1',
                        use_instrument_profile:bool=False,
                        reflectance:bool=False,
                        standard_data:typing.Optional[dict]=None,
                        input_path:typing.Optional[str]=None,
                        output_path:typing.Optional[str]=None,
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

    standard_name : str 
        The name of the standard star

    standard_file : str 
        The name of the FITS file containing the standard star spectra.

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

    check_parameter('telluric_correction', 'standard_name', standard_name,
                    'str')
    
    check_parameter('telluric_correction', 'standard_file', standard_file,
                    'str')

    check_parameter('telluric_correction', 'output_name', output_name, 'str')

    check_parameter('telluric_correction', 'fluxdensity_units',
                    fluxdensity_units, 'str')    

    check_parameter('telluric_correction', 'use_instrument_profile',
                    use_instrument_profile, 'bool')    
    
    check_parameter('telluric_correction', 'reflectance', reflectance, 'bool')
    
    check_parameter('telluric_correction', 'standard_data', standard_data,
                    ['dict', 'NoneType'])    

    check_parameter('telluric_correction', 'input_path', input_path,
                    ['NoneType', 'str'])

    check_parameter('telluric_correction', 'output_path', output_path,
                    ['NoneType', 'str'])

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
    telluric.load['standard_name'] = standard_name
    telluric.load['standard_file'] = standard_file
    telluric.load['standard_data'] = standard_data    
    telluric.load['output_name'] = output_name
    telluric.load['fluxdensity_units'] = fluxdensity_units
    telluric.load['reflectance'] = reflectance
    telluric.load['input_path'] = input_path
    telluric.load['output_path'] = output_path
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

        #
        # Make the telluric correction spectra
        #

        make_telluric_spectra()

        #
        # Correct the data for telluric absorption and flux calibrate
        #

        correct_spectra()

        

        

            

    
#        kernel_info = get_kernel(np.squeeze(standard_spectra[z,0,:]),
#                                 np.squeeze(standard_spectra[z,1,:]),
#                                 standard_info['resolvingpower'],
#                                 vega_wavelength, vega_flux, vega_continuum,
#                                 mode_info['normalization_windows'],
#                                 mode_info['polynomial_degree'],
#                                 mode_info['rv_range'])
#                    
#    else:
#
#        yunits = ''
#        lyunits = ''
#        lylabel = 'ratio'
#          
#
#    return
#        
#    #
#    # Start the loop over object order number
#    #
#
#    for i in range(object_info['norders']):
#
#        # Find the order
#        
#        z = np.where(object_info['orders'][i] == standard_info['orders'])
#
#        # Now loop over the apertures
#        
#        for j in range(object_info['napertures']):
#
#            k =z[0][0]*object_info['napertures']+j
#            
#            # Interpolate the standard and bit mask onto the object
#                
#            rspectrum, runc = linear_interp1d(standard_spectra[z,0,:],
#                                              standard_spectra[z,1,:],
#                                              object_spectra[k,0,:],
#                                              input_u=standard_spectra[z,2,:])
#
#            if reflectance is False:
#
#                rvega = linear_interp1d(vega_wavelength, vega_flux,
#                                        object_spectra[k,0,:])
#
#                correction = rvega/rspectrum
#                correction_unc =  rvega/rspectrum**2 * runc               
#                
#
#            else:
#
#                correction = 1/rspectrum
#                correction_unc = 1/rspectrum**2 * runc
#                
#            
#            # Do the correction and error propagation
#
#            value = object_spectra[k,1,:] * correction
#
#            var = object_spectra[k,1,:]**2 * correction_unc**2 + \
#                  correction**2 * object_spectra[k,2,:]**2
#            
#            object_spectra[k,1,:] = value
#            object_spectra[k,2,:] = np.sqrt(var)
#
#            # Combine the masks
#
#            mask = linear_bitmask_interp1d(standard_spectra[z,0,:],
#                                 standard_spectra[z,3,:].astype(np.uint8),
#                                           object_spectra[k,0,:])
#
#
#            stack = np.stack((object_spectra[k,3,:].astype(np.uint8),mask))
#            mask = combine_flag_stack(stack)
#
#            object_spectra[k,3,:] = mask
#
#    #
#    # Write the file to disk
#    #
#
#    # Update the object hdrinfo 
#
#    object_hdrinfo['MODULE'][0] = 'telluric'
#    object_hdrinfo['FILENAME'][0] = output_name+'.fits'
#    object_hdrinfo['YUNITS'][0] = yunits
#    object_hdrinfo['LYUNITS'][0] = lyunits
#    object_hdrinfo['LYLABEL'][0] = lylabel
#
#    object_hdrinfo['OBJFILE'] = (object_file, ' Object file name')
#    object_hdrinfo['STANDARD'] = (standard_name, ' Standard')
#    object_hdrinfo['STDSPTYP'] = (standard_sptype, ' Standard spectral type')
#    object_hdrinfo['STDBMAG'] = (standard_bmag, ' Standard spectral B mag')
#    object_hdrinfo['STDVMAG'] = (standard_vmag, ' Standard spectral V mag')    
#    object_hdrinfo['STDFILE'] = (standard_file, ' Standard file name')        
#
#    # Store the history
#
#    old_history = object_hdrinfo['HISTORY']
#    
#    # remove it from the avehdr
#
#    object_hdrinfo.pop('HISTORY')
#    
#    # Create a new header
#
#    phdu = fits.PrimaryHDU()
#    hdr = phdu.header
#
#    # Add our keywords
#    
#    keys = list(object_hdrinfo.keys())
#    
#    for i in range(len(keys)):
#
#        if keys[i] == 'COMMENT':
#
#            junk = 1
#
#        else:
#
#            hdr[keys[i]] = (object_hdrinfo[keys[i]][0],
#                            object_hdrinfo[keys[i]][1])
#
#    history = split_text(old_history, length=65)            
#
#    for hist in history:
#
#        hdr['HISTORY'] = hist
#
#    output_fullpath = os.path.join(output_path, output_name+'.fits')
#    
#    fits.writeto(output_fullpath, object_spectra, hdr, overwrite=overwrite)
#
#    #
#    # Plot the results
#    #
#
##    if telluric.load['qa_show'] is True:
##
##        number = plot_spectra(output_fullpath, title=output_name+'.fits',
##                              plot_size=qa_showsize,
##                              plot_number=telluric.state['spectra_plotnum'])
##        telluric.state['spectra_plotnum'] = number
##                    
##    if telluric.load['qa_write'] is True:
##
##        qafileinfo = {'figsize': telluric.load['qa_showsize'],
##                      'filepath': telluric.load['qa_path'],
##                      'filename': telluric.load['output_name'],
##                      'extension': setup.state['qa_extension']}
##
##        plot_spectra(output_fullpath, file_info=qafileinfo)
#
#    #
#    # Update the user
#    #
#
#    message = ' Wrote '+os.path.basename(output_fullpath)+' to disk.'
#    logging.info(message)
#
#

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


    #
    # Start the loop over object order number
    #
        
    for i in range(telluric.state['object_norders']):

        # Find the order

        z_order = np.where(telluric.state['object_orders'][i] == \
                           telluric.state['standard_orders'])

       # Now loop over the apertures
        
        for j in range(telluric.state['object_napertures']):

            k =z_order[0][0]*telluric.state['object_napertures']+j
#            
#            # Interpolate the standard and bit mask onto the object
#                
#            rspectrum, runc = linear_interp1d(standard_spectra[z,0,:],
#                                              standard_spectra[z,1,:],
#                                              object_spectra[k,0,:],
#                                              input_u=standard_spectra[z,2,:])
#

    
    

    
    x = 1

    # Convert to the users requested units
    
#    vega_flux = units.convert_fluxdensity(vega_wavelength, vega_flux,
#                                          'um', 'erg s-1 cm-2 A-1',
#                                          telluric.load['fluxdensity_units']) 
#    
#    # Set the units
#        
#    yunits = telluric.load['fluxdensity_units']
#    latex_yunits, latex_ylabel = \
#        units.get_latex_fluxdensity(telluric.load['fluxdensity_units'])
#
#    
#    telluric.state['yunits'] = yunits
#    telluric.state['latex_yunits'] = latex_yunits
#    telluric.state['latex_ylabel'] = latex_ylabel
    


    




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
                       'title':title,
                       'xlabel':xlabel}
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

    if telluric.load['standard_data'] is None:
    
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

        standard_sptype = standard_data['sptype']
        standard_vmag = float(standard_data['bmag'])
        standard_bmag = float(standard_data['vmag'])


    #
    # error check standard information - what do we do if there is no
    # standard data?
    #
    
    try:
        type(standard_name)

    except:
        
        message = 'Standard name "{}"" was not found in SIMBAD; provided '\
            'correct name or provide the optional standard_data keyword with '\
            'a dictionary containing keys '\
            '"sptype", "bmag", and "vmag"'.format(standard)

        
        raise pySpextoolError(message)

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
       
    # Scale it by the vband magnitude
    
    vega_vmagnitude = 0.03
    vega_bminv = 0.00
    scale = 10.0**(-0.4*(telluric.state['standard_vmag']-vega_vmagnitude))
    vega_flux /= scale
    vega_continuum /= scale
    vega_fitted_continuum /= scale        
    
        
    # Store the results

    telluric.state['vega_wavelength'] = vega_wavelength
    telluric.state['vega_fluxdensity'] = vega_flux
    telluric.state['vega_continuum'] = vega_continuum
    telluric.state['vega_fitted_continuum'] = vega_fitted_continuum
    telluric.state['vega_normalized_fluxdensity'] = vega_flux/vega_fitted_continuum

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
    
    for i in range(telluric.state['object_norders']):    

        telluric_spectra = telluric.state['standard_spectra']
        
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
        # Change units to those requested by users
        #

        

        
        telluric_spectra[i,1,:] = result[0]
        telluric_spectra[i,2,:] = result[1]
        


    #
    # Store the results
    #

    telluric.state['standard_spectra'] = telluric_spectra
        

    
    
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
    
                




                













    

