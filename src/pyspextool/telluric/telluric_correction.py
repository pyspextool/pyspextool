import os
import numpy as np
import numpy.typing as npt
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table.table import Table
import logging
from astropy.io import fits
import typing

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
from pyspextool.telluric.kernel import make_instrument_profile
from pyspextool.utils.for_print import for_print

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def telluric_correction(object_file:str, standard:str, standard_file:str,
                        output_name:str, fluxdensity_units:str='W m-2 um-1',
                        reflectance:bool=False,
                        standard_data:typing.Optional[dict]=None,
                        input_path:typing.Optional[str]=None,
                        output_path:typing.Optional[str]=None,
                        qa_path:typing.Optional[str]=None,
                        qa_show:typing.Optional[bool]=None,
                        qa_showscale:float=1.0,
                        qa_write:typing.Optional[bool]=None,
                        verbose:typing.Optional[bool]=None,
                        overwrite:bool=True):

    """
    To correct spectra for telluric absorption and flux calibrate

    Parameters
    ----------
    object_file : str 
        The name of the FITS file containing the object star spectra.

    standard : str 
        The name of the standard star

    standard_file : str 
        The name of the FITS file containing the standard star spectra.

    output_name : str
        The output file name sans the suffix.

    fluxdensity_units : {'W m-2 um-1', 'erg s-1 cm-2 A-1', 'W m-2 Hz-1', 
                         'ergs s-1 cm-2 Hz-1', 'Jy', 'mJy', 'uJy'}
        The flux density units of the output.  
        
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

    qa_showscale : tuple, default=(10, 6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}
        Set to True/False to override config.state['qa_write'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

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

    check_parameter('telluric_correction', 'standard', standard, 'str')
    
    check_parameter('telluric_correction', 'standard_file', standard_file,
                    'str')

    check_parameter('telluric_correction', 'output_name', output_name, 'str')

    check_parameter('telluric_correction', 'fluxdensity_units',
                    fluxdensity_units, 'str')    

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

    check_parameter('telluric_correction', 'qa_showscale', qa_showscale,
                    'float')
    
    check_parameter('telluric_correction', 'qa_write', qa_write,
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

    if verbose is None:
        verbose = setup.state['verbose']

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
    telluric.load['standard'] = standard
    telluric.load['standard_file'] = standard_file
    telluric.load['standard_data'] = standard_data    
    telluric.load['output_name'] = output_name
    telluric.load['fluxdensity_units'] = fluxdensity_units
    telluric.load['reflectance'] = reflectance
    telluric.load['input_path'] = input_path
    telluric.load['output_path'] = output_path
    telluric.load['qa_path'] = qa_path      
    telluric.load['qa_show'] = qa_show
    telluric.load['qa_showscale'] = qa_showscale
    telluric.load['qa_write'] = qa_write
    telluric.load['verbose'] = verbose
    telluric.load['overwrite'] = overwrite

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
        # Get the kernel
        #

        telluric.state['mode_info']['method'] = 'IP'

        if telluric.state['mode_info']['method'] == 'IP':

            logging.info(f" Mike must deal with RV shift.")            

        
        if telluric.state['mode_info']['method'] == 'deconvolution':

            deconvolve_line()

            

    
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
    vega_models = list(values[:,1])
    method = list(values[:,2])
    order = list(values[:,3])    
    norm_ranges = list(values[:,4])
    poly_degree = list(values[:,5])
    rv_range = list(values[:,6])    

    #
    # Modify each vector accordingly.
    #
    
    modes = np.array([x.strip() for x in modes])
    
    vega_models = np.array([x.strip() for x in vega_models])
    method = np.array([x.strip() for x in method])
    order = np.array([int(x) for x in order])    

    norm_ranges = [x.split() for x in norm_ranges]
    for i in range(len(norm_ranges)):

        norm_ranges[i] = np.array([float(x) for x in norm_ranges[i]])

    norm_ranges = np.array(norm_ranges)
    poly_degree = np.array([int(x) for x in poly_degree])        

    rv_range = [x.split() for x in rv_range]
    for i in range(len(rv_range)):

        rv_range[i] = np.array([float(x) for x in rv_range[i]])

    rv_range = np.array(rv_range)
      
    #
    # Find the matching mode and return results
    #
    
    z = modes == mode

    return {'vega_model':vega_models[z][0], 'method':method[z][0],
            'order':order[z][0], 'normalization_windows':norm_ranges[z][0],
            'polynomial_degree':int(poly_degree[z][0]),
            'rv_range':rv_range[z][0]}



    

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
    # Get standard star information
    #

    if telluric.load['standard_data'] is None:
    
        # Get SIMBAD information of the standard

        logging.info(f' Querying SIMBAD for standard star information...')
        
        Simbad.add_votable_fields('sptype', 'flux(B)', 'flux(V)')
        table = Simbad.query_object(telluric.load['standard'])

        if isinstance(table,Table):
            standard_name = table['MAIN_ID'][0]
            standard_sptype = table['SP_TYPE'][0]
            standard_vmag = table['FLUX_V'][0]
            standard_bmag = table['FLUX_B'][0]         

    else:

        standard_name = standard
        standard_sptype = standard_data['sptype']
        standard_vmag = standard_data['bmag']
        standard_bmag = standard_data['vmag']

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

        
        raise ValueError(message)

    #
    # Let's do some checks to ensure 1) the standard has only one aperture and
    # 2) the standard has the orders required.
    #

    if standard_info['napertures'] != 1:

        message = 'The standard has more than one aperture.'
        raise ValueError(message)
    
    intersection = np.intersect1d(object_info['orders'],
                                  standard_info['orders'])

    if np.size(intersection) != object_info['norders']:

        message = 'The standard lacks an order the object has.'
        raise ValueError(message)

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

    #
    # Get other information from the headers
    #
    
    telluric.state['object_norders'] = object_hdrinfo['NORDERS'][0]
    telluric.state['standard_norders'] = standard_hdrinfo['NORDERS'][0]

    telluric.state['object_orders'] = object_hdrinfo['ORDERS'][0]
    telluric.state['standard_orders'] = standard_hdrinfo['ORDERS'][0]

    telluric.state['standard_dispersions'] = standard_info['dispersions']

    telluric.state['standard_fwhm'] = standard_info['dispersions']*\
        standard_info['slitw_pix']
            
    #
    # Finally, get the mode info 
    #
    
    mode_info = get_modeinfo(telluric.state['standard_hdrinfo']['MODE'][0])
    telluric.state['mode_info'] = mode_info

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

    root = 'Vega'+telluric.state['mode_info']['vega_model']+'.fits'
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
    
    # Convert to the users requested units
    
    vega_flux = units.convert_fluxdensity(vega_wavelength, vega_flux,
                                          'um', 'erg s-1 cm-2 A-1',
                                          telluric.load['fluxdensity_units']) 
    
    # Set the units
        
    yunits = telluric.load['fluxdensity_units']
    latex_yunits, latex_ylabel = \
        units.get_latex_fluxdensity(telluric.load['fluxdensity_units'])
        
    # Store the results

    telluric.state['vega_wavelength'] = vega_wavelength
    telluric.state['vega_flux'] = vega_flux
    telluric.state['vega_continuum'] = vega_continuum
    telluric.state['vega_fitted_continuum'] = vega_fitted_continuum
    telluric.state['yunits'] = yunits
    telluric.state['latex_yunits'] = latex_yunits
    telluric.state['latex_ylabel'] = latex_ylabel
    
    #
    # Generate the instrument profiles so they are only created once
    #

    instrument_profiles = []
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

        instrument_profiles.append(p)

    # Store the results
        
    telluric.state['convolution_kernels'] = instrument_profiles
    
        



    
def get_kernel(object_wavelength:npt.ArrayLike, object_flux:npt.ArrayLike,
               resolving_power:float, vega_wavelength:npt.ArrayLike,
               vega_flux:npt.ArrayLike, vega_continuum:npt.ArrayLike,
               normalization_windows:npt.ArrayLike, polynomial_degree:float,
               rv_range:npt.ArrayLike):


    """
    To calculate the convolution kernel

    Parameters
    ----------

    Returns
    -------


    """


    #
    # Normalized the spectrum
    #


#        z = np.where(telluric.state['standard_info']['orders'] == \
#                     telluric.state['mode_info']['order'])[0]

    
    print('hi')
    object_nflux = normalize_continuum(object_wavelength, object_flux,
                                       normalization_windows,
                                       polynomial_degree)

    #
    # Find the velocity shift
    #

#    z1 = object_wavelength > rv_range[0]
#    z2 = object_wavelength < rv_range[1]
#    z = np.logical_and(z1,z2)
#    
#    vshift = model_xcorrelate(object_wavelength[z], object_nflux[z],
#                              vega_wavelength, vega_flux/vega_continuum,
#                              resolving_power=resolving_power,qa_show=True)
#
#    print(vshift)
    

    
    

    

    x = 1

    return 1
    


    


    
