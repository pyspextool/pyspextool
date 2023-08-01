import os
import numpy as np
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table.table import Table

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

from pyspextool.utils.for_print import for_print


def telluric_correction(object_file, standard, standard_file, output_name,
                        standard_data=None, input_path=None, output_path=None,
                        qa_path=None, qa_plot=None, qa_plotsize=(10, 6),
                        qa_file=None, verbose=None, overwrite=True):

    """
    To combine raw extracted spectra

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

    if verbose is True:

        print('Reflectance Telluric Correction')
        print('-------------------------------')
        print('Loading the data...')        
        
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

        if verbose is True:

            print('Querying SIMBAD for standard star information...')
    
    
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
    # Start the loop over object order number
    #

    for i in range(object_info['norders']):

        # Find the order
        
        z = np.where(object_info['orders'][i] == standard_info['orders'])

        # Now loop over the apertures
        
        for j in range(object_info['napertures']):

            k =z[0][0]*object_info['napertures']+j
            
            # Interpolate the data and the bit mask
                
            spectrum, unc = linear_interp1d(standard_spectra[z,0,:],
                                            standard_spectra[z,1,:],
                                            object_spectra[k,0,:],
                                            input_u=standard_spectra[z,2,:])

            # Do the correction and error propagation

            value = object_spectra[k,1,:] / spectrum

            var = value**2 * ((object_spectra[k,2,:]/object_spectra[k,1,:])**2+
                              (unc/spectrum)**2)     

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
    object_hdrinfo['YUNITS'][0] = ''
    object_hdrinfo['LYUNITS'][0] = ''
    object_hdrinfo['LYLABEL'][0] = 'ratio'

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

    if verbose is True:

            print('Wrote', os.path.basename(output_fullpath), 'to disk.')
