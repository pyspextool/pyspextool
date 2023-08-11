import os
import numpy as np
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table.table import Table
import logging
from astropy.io import fits

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

from pyspextool.utils.for_print import for_print


def telluric_correction(object_file, standard, standard_file, output_name,
                        fluxdensity_units='W m-2 um-1', reflectance=False,
                        standard_data=None, input_path=None, output_path=None,
                        qa_path=None, qa_plot=None, qa_plotsize=(10, 6),
                        qa_file=None, verbose=None, overwrite=True):

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

    qa_plot : {None, True, False}
        Set to True/False to override config.state['qa_plot'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qa_plotsize : tuple, default=(10, 6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_file : {None, True, False}
        Set to True/False to override config.state['qa_file'] in the 
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

    check_parameter('basic_tellcor', 'object_file', object_file, 'str')

    check_parameter('basic_tellcor', 'standard', standard, 'str')
    
    check_parameter('basic_tellcor', 'standard_file', standard_file, 'str')

    check_parameter('basic_tellcor', 'output_name', output_name, 'str')

    check_parameter('basic_tellcor', 'fluxdensity_units', fluxdensity_units,
                    'str')    

    check_parameter('basic_tellcor', 'reflectance', reflectance, 'bool')    

    check_parameter('basic_tellcor', 'standard_data', standard_data,
                    ['dict', 'NoneType'])    

    check_parameter('basic_tellcor', 'input_path', input_path,
                    ['NoneType', 'str'])

    check_parameter('basic_tellcor', 'output_path', output_path,
                    ['NoneType', 'str'])

    check_parameter('basic_tellcor', 'qa_path', qa_path, ['NoneType', 'str'])            
    check_parameter('basic_tellcor', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('basic_tellcor', 'qa_plotsize', qa_plotsize, 'tuple')    

    check_parameter('basic_tellcor', 'qa_file', qa_file, ['NoneType', 'bool'])

    check_parameter('basic_tellcor', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('basic_tellcor', 'overwrite', overwrite, 'bool')    

    #
    # Check the qa and verbose variables and set to system default if need be.
    #
        
    if qa_file is None:

        qa_file = setup.state['qa_file']

    if qa_plot is None:

        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

    if verbose is True:
        logging.getLogger().setLevel(logging.INFO)
        setup.state["verbose"] = True
        
    elif verbose is False:
        logging.getLogger().setLevel(logging.ERROR)
        setup.state["verbose"] = False

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)            
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
    telluric.load['output_name'] = output_name
    telluric.load['fluxdensity_units'] = fluxdensity_units
    telluric.load['reflectance'] = reflectance
    telluric.load['input_path'] = input_path
    telluric.load['output_path'] = output_path
    telluric.load['qa_path'] = qa_path      
    telluric.load['qa_plot'] = qa_plot
    telluric.load['qa_plotsize'] = qa_plotsize
    telluric.load['qa_file'] = qa_file
    telluric.load['verbose'] = verbose
    telluric.load['overwrite'] = overwrite

    #
    # Load the files into memory
    #

    logging.info(f" Telluric Correction\n-------------------------\n")

    logging.info(f" Loading the data...")
            
    # Object first
    
    fullpath = make_full_path(input_path, object_file, exist=True)

    object_spectra, object_info = read_spectra_fits(fullpath)

    header = fits.getheader(fullpath)

    object_hdrinfo = get_header_info(header,
                                    keywords=setup.state['telluric_keywords'])

    # Now the standard

    fullpath = make_full_path(input_path, standard_file, exist=True)

    standard_spectra, standard_info = read_spectra_fits(fullpath)

    header = fits.getheader(fullpath)

    standard_hdrinfo = get_header_info(header,
                                    keywords=setup.state['telluric_keywords'])

    #
    # Get standard star information
    #

    if standard_data is None:
    
        # Get SIMBAD information of the standard

        logging.info(f'Querying SIMBAD for standard star information...')
        
        Simbad.add_votable_fields('sptype', 'flux(B)', 'flux(V)')
        table = Simbad.query_object(standard)

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

    # error check standard information - what do we do if there is no standard data?
    try:
        type(standard_name)
    except:
        raise ValueError('Standard name "{}"" was not found in SIMBAD; provided correct name or provide the optional standard_data keyword with a dictionary containing keys "sptype", "bmag", and "vmag"'.format(standard))

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
    # Now figure out which Vega model to use and load if need be
    #

    if reflectance is False:

        model = get_modeinfo(standard_hdrinfo['MODE'][0])

        file = os.path.join(setup.state['package_path'],'Vega'+model+'.fits')

        hdul = fits.open(setup.state['package_path']+'/data/Vega5000.fits') 
        data  = hdul[1].data

        vega_wavelength = data['wavelength']
        vega_flux = data['flux density']
        vega_continuum = data['continuum flux density']
        vega_fitted_continuum = data['fitted continuum flux density']
        
        hdul.close()

        # Scale it by the vband magnitude

        vega_vmagnitude = 0.03
        vega_bminv = 0.00
        scale = 10.0**(-0.4*(standard_vmag-vega_vmagnitude))
        vega_flux /= scale

        # Convert to the users requested units
        
        vega_flux = units.convert_fluxdensity(vega_wavelength, vega_flux,
                                              'um', 'erg s-1 cm-2 A-1',
                                              fluxdensity_units) 
        
        # Set the units

        yunits = fluxdensity_units
        lyunits, lylabel = units.get_latex_fluxdensity(fluxdensity_units)

        
    else:

        yunits = ''
        lyunits = ''
        lylabel = 'ratio'
        
    #
    # Start the loop over object order number
    #

    for i in range(object_info['norders']):

        # Find the order
        
        z = np.where(object_info['orders'][i] == standard_info['orders'])

        # Now loop over the apertures
        
        for j in range(object_info['napertures']):

            k =z[0][0]*object_info['napertures']+j
            
            # Interpolate the standard and bit mask onto the object
                
            rspectrum, runc = linear_interp1d(standard_spectra[z,0,:],
                                              standard_spectra[z,1,:],
                                              object_spectra[k,0,:],
                                              input_u=standard_spectra[z,2,:])

            if reflectance is False:

                rvega = linear_interp1d(vega_wavelength, vega_flux,
                                        object_spectra[k,0,:])

                correction = rvega/rspectrum
                correction_unc =  rvega/rspectrum**2 * runc               
                

            else:

                correction = 1/rspectrum
                correction_unc = 1/rspectrum**2 * runc
                
            
            # Do the correction and error propagation

            value = object_spectra[k,1,:] * correction

            var = object_spectra[k,1,:]**2 * correction_unc**2 + \
                  correction**2 * object_spectra[k,2,:]**2
            
            object_spectra[k,1,:] = value
            object_spectra[k,2,:] = np.sqrt(var)

            # Combine the masks

            mask = linear_bitmask_interp1d(standard_spectra[z,0,:],
                                 standard_spectra[z,3,:].astype(np.uint8),
                                           object_spectra[k,0,:])


            stack = np.stack((object_spectra[k,3,:].astype(np.uint8),mask))
            mask = combine_flag_stack(stack)

            object_spectra[k,3,:] = mask

    #
    # Write the file to disk
    #

    # Update the object hdrinfo 

    object_hdrinfo['MODULE'][0] = 'telluric'
    object_hdrinfo['FILENAME'][0] = output_name+'.fits'
    object_hdrinfo['YUNITS'][0] = yunits
    object_hdrinfo['LYUNITS'][0] = lyunits
    object_hdrinfo['LYLABEL'][0] = lylabel

    object_hdrinfo['OBJFILE'] = (object_file, ' Object file name')
    object_hdrinfo['STANDARD'] = (standard_name, ' Standard')
    object_hdrinfo['STDSPTYP'] = (standard_sptype, ' Standard spectral type')
    object_hdrinfo['STDBMAG'] = (standard_bmag, ' Standard spectral B mag')
    object_hdrinfo['STDVMAG'] = (standard_vmag, ' Standard spectral V mag')    
    object_hdrinfo['STDFILE'] = (standard_file, ' Standard file name')        

    # Store the history

    old_history = object_hdrinfo['HISTORY']
    
    # remove it from the avehdr

    object_hdrinfo.pop('HISTORY')
    
    # Create a new header

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add our keywords
    
    keys = list(object_hdrinfo.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (object_hdrinfo[keys[i]][0],
                            object_hdrinfo[keys[i]][1])

    history = split_text(old_history, length=65)            

    for hist in history:

        hdr['HISTORY'] = hist

    output_fullpath = os.path.join(output_path, output_name+'.fits')
    
    fits.writeto(output_fullpath, object_spectra, hdr, overwrite=overwrite)

    #
    # Plot the results
    #

    if telluric.load['qa_plot'] is True:

        number = plot_spectra(output_fullpath, title=output_name+'.fits',
                              plot_size=qa_plotsize,
                              plot_number=telluric.state['spectra_plotnum'])
        telluric.state['spectra_plotnum'] = number
                    
    if telluric.load['qa_file'] is True:

        qafileinfo = {'figsize': telluric.load['qa_plotsize'],
                      'filepath': telluric.load['qa_path'],
                      'filename': telluric.load['output_name'],
                      'extension': setup.state['qa_extension']}

        plot_spectra(output_fullpath, file_info=qafileinfo)

    #
    # Update the user
    #

    message = ' Wrote '+os.path.basename(output_fullpath)+' to disk.'
    logging.info(message)


def get_modeinfo(mode):

    """
    Reads the telluric_modeinfo file and returns values for a given mode.

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

    # Get the file name and read it
    
    file = os.path.join(setup.state['instrument_path'],'telluric_modeinfo.dat')
    
    modes, vega_models = np.loadtxt(file, comments='#', unpack=True,
                                    delimiter='|', dtype='str')

    # Convert to list and compress strings
    
    modes = list(modes)
    vega_models = list(vega_models)

    modes = np.array([m.strip() for m in modes])
    vega_models = np.array([m.strip() for m in vega_models]    )

    # Find the matching mode and return results
    
    z = modes == mode

    return vega_models[z][0]

    

    
