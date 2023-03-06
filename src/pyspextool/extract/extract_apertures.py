import numpy as np
import os

from pyspextool.extract import config
from pyspextool.extract.check_continue import check_continue
from pyspextool.extract.extract_pointsource_1dxd import extract_pointsource_1dxd
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.io.check import check_parameter
from pyspextool.io.write_apertures_fits import write_apertures_fits



def extract_apertures(output_filenames=None, verbose=True):

    """
    User function to extract spectra.

    Parameters
    ----------
    output_filenames : list of str, optional
        A (nfiles,) list of output file names sans the suffix.  Only required 
        if the data are read in in the file read mode.

    verbose : {True, False}, optional
        Set to True for command line updates during execution.

    Returns 
    -------
    None

    """

    #
    # Check parameters.  Just do it now to not waste extracting.
    #

    check_parameter('extract_apertures', 'output_filenames',
                    output_filenames, ['str', 'list', 'NoneType'])

    check_parameter('extract_apertures', 'verbose',
                    verbose, 'bool')


    #
    # Store user inputs
    #

    config.user['extract']['verbose'] = verbose
    
    #
    # The procedure depends on the extraction type, point or extended source
    #

    xsbginfo = None
    psbginfo = None
    
    if config.state['exttype'] == 'ps':

        #
        #========================= Point Source ============================
        #

        # Grab which orders are being extracted
            
        z = config.state['psdoorders'] == 1

        # Create the background information dictionary
        
        if config.state['bgradius'] is not None:

            psbginfo = {'radius': config.state['bgradius'],
                        'width': config.state['bgwidth'],
                        'degree': config.state['bgfitdeg']}
        
        # Create the aperture radii list

        apradii = np.full(config.state['naps'], config.state['apradii'])
        
        # Do the extraction
            
        spectra = extract_pointsource_1dxd(config.state['workimage'],
                                           config.state['varimage'],
                                           config.state['ordermask'],
                                           config.state['orders'][z],
                                           config.state['wavecal'],
                                           config.state['spatcal'],
                                           config.state['tracecoeffs'],
                                           config.state['apradii'],
                                           config.state['apsigns'],                                                        bginfo=psbginfo)

    else:

        #
        #======================= Extended Source ===========================
        #

        #
        # Create the full path
        #
        
        # Make a list if it isn't one.
        
        full_filenames = output_filenames
        if isinstance(output_filenames, str):

            output_filenames = os.path.join(config.state['procpath'],
                                            output_filenames)

        else:
            
            output_filenames = [os.path.join(config.state['procpath'],
                                output_filenames) for e in output_filenames]

            
        #
        # Now do the extracton
        #
        
        # Grab which orders are being extracted
            
        z = config.state['xsdoorders'] == 1
        
        # Create background information dictionary

        if config.state['bgregions'] is not None:

            xsbginfo = {'regions': config.state['bgregions'],
                        'degree': config.state['bgfitdeg']}

        # Do the extraction
        
        spectra = extract_extendedsource_1dxd(config.state['workimage'],
                                              config.state['varimage'],
                                              config.state['ordermask'],
                                              config.state['orders'][z],
                                              config.state['wavecal'],
                                              config.state['spatcal'],
                                              config.state['apertures'][z],
                                              config.state['apradii'],
                                              bginfo=xsbginfo)

    #
    # Write the results to disk
    #
        
    write_apertures(spectra, psbginfo=psbginfo, xsbginfo=xsbginfo,
                    output_filenames=output_filenames, verbose=verbose)

        
def write_apertures(spectra, psbginfo=None, xsbginfo=None,
                    output_filenames=None, verbose=True):

    """
    To write extracted spectra to disk.

    Parameters
    ----------
    spectra : dict
        The dictionary returned from any of the extract_ functions.
        
        `"spectra"` : list
            (norders*naps,) list of spectral planes of the extracted target.

        `"background"` : list
            (norders*naps,) list of spectral planes of the background of the 
            target.

    psbginfo : dict, optional
        `"radius"` : float
            The background start radius in arcseconds

        `"width"` : float
            The background annulus width in arcseconds

        `"degree"` : int
            The background polynomial fit degree.

    xsbginfo : dict, optional

        `"regions"` : str
            A comma-separated list of background regions, e.g. '1-2,14-15'.

        `"degree"` : int
            The background polynomial fit degree.

    output_filenames : list of str, optional
        A (nfiles,) list of output file names sans the suffix.  Only required 
        if the data are loaded using the file read mode.

    verbose : {True, False}, optional
        Set to True for command line updates during execution.

    Returns 
    -------
    None

       Writes FITS files to disk.

    """
    
    # Unpack the spectra dictionary
    
    if spectra['background'] is not None:
        background = spectra

    else:

        background = None

    spectra = spectra['spectra']
    
    # Get the wavelengh calibration info if used.

    if config.user['load']['wavecalfile'] != None:

        wavecalinfo = {'file':config.user['load']['wavecalfile'],
                       'wavecaltype':config.state['wavecaltype'],
                       'wavetype':'vacuum'}                           

    else:

        wavecalinfo = None
    
    # Things proceed depending on the extraction mode

    if config.state['exttype'] == 'ps':

        #
        #========================= Point Source ============================
        #

        # Grab things all modes need

        xranges = config.state['xranges'][config.state['psdoorders']]
        orders = config.state['orders'][config.state['psdoorders']]    
        apertures = config.state['apertures'][config.state['psdoorders']]
        norders = np.sum(config.state['psdoorders'])

        # Now go mode by mode
        
        if config.state['reductionmode'] == 'A':

            # Get file information
             
            hdrinfo = config.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'
            output_fullpath = config.state['output_files'][0]

            # Write ther results
            
            write_apertures_fits(spectra, xranges, aimage, sky,
                                 config.user['load']['flatfile'],
                                 config.state['naps'],
                                 orders, hdrinfo, apertures,
                                 config.state['apradii'],
                                 config.state['plate_scale'],
                                 config.state['slith_pix'],
                                 config.state['slith_arc'],
                                 config.state['slitw_pix'],
                                 config.state['slitw_arc'],
                                 config.state['resolvingpower'],
                                 'um', 'DN s-1', 'microns', 'flux',
                                 config.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 psbginfo=psbginfo,
                                 verbose=verbose)


            
        elif config.state['reductionmode'] == 'A-B':       

            # Are there positive apertures?

            z_pos = config.state['apsigns'] == 1
            if np.sum(z_pos) != 0:

                # Get file information
             
                hdrinfo = config.state['hdrinfo'][0]
                aimage = hdrinfo['FILENAME'][0]
                sky = config.state['hdrinfo'][1]['FILENAME'][0]
                output_fullpath = config.state['output_files'][0]

                # Now get the indices for the spectra array

                full_apsign = np.tile(config.state['apsigns'],norders)
                z = (np.where(full_apsign == 1))[0]
                pos_spectra = [spectra[i] for i in z]

                write_apertures_fits(pos_spectra, xranges, aimage, sky,
                                     config.user['load']['flatfile'],
                                     config.state['naps'], orders,
                                     hdrinfo, apertures,
                                     config.state['apradii'],
                                     config.state['plate_scale'],
                                     config.state['slith_pix'],
                                     config.state['slith_arc'],
                                     config.state['slitw_pix'],
                                     config.state['slitw_arc'],
                                     config.state['resolvingpower'],
                                     'um', 'DN s-1', 'microns', 'flux',
                                     config.state['version'],
                                     output_fullpath, wavecalinfo=wavecalinfo,
#                                     background_spectra=background[z],
                                     psbginfo=psbginfo,
                                     verbose=verbose)

            # Are there negative apertures?

            z_pos = config.state['apsigns'] == -1
            if np.sum(z_pos) != 0:

                # Get file information
             
                hdrinfo = config.state['hdrinfo'][1]
                aimage = hdrinfo['FILENAME'][1]
                sky = config.state['hdrinfo'][0]['FILENAME'][0]
                output_fullpath = config.state['output_files'][1]

                # Now get the indices for the spectra array

                full_apsign = np.tile(config.state['apsigns'],norders)
                z = (np.where(full_apsign == -1))[0]
                neg_spectra = [spectra[i] for i in z]

                write_apertures_fits(neg_spectra, xranges, aimage, sky,
                                     config.user['load']['flatfile'],
                                     config.state['naps'], orders,
                                     hdrinfo, apertures,
                                     config.state['apradii'],
                                     config.state['plate_scale'],
                                     config.state['slith_pix'],
                                     config.state['slith_arc'],
                                     config.state['slitw_pix'],
                                     config.state['slitw_arc'],
                                     config.state['resolvingpower'],
                                     'um', 'DN s-1', 'microns', 'flux',
                                     config.state['version'],
                                     output_fullpath, wavecalinfo=wavecalinfo,
#                                     background_spectra=background[z],
                                     psbginfo=psbginfo,
                                     verbose=verbose)                
                
            
        elif config.state['reductionmode'] == 'A-Sky/Dark':

            x = 1
            
        else:

            print('Reduction Mode Unknown.')     
            
    else:

        #
        #======================= Extended Source ===========================
        #

        # Grab things all modes need
        
        xranges = config.state['xranges'][config.state['psdoorders']]
        orders = config.state['orders'][config.state['psdoorders']]    
        apertures = config.state['apertures'][config.state['psdoorders']]
        norders = np.sum(config.state['psdoorders'])

        if config.state['reductionmode'] == 'A':

            hdrinfo = config.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 config.user['load']['flatfile'],
                                 config.state['naps'],
                                 orders, hdrinfo, apertures,
                                 config.state['apradii'],
                                 config.state['plate_scale'],
                                 config.state['slith_pix'],
                                 config.state['slith_arc'],
                                 config.state['slitw_pix'],
                                 config.state['slitw_arc'],
                                 config.state['resolvingpower'],
                                 'um', 'DN s-1', 'microns', 'flux',
                                 config.state['version'],
                                 output_filenames, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

        
        elif config.state['reductionmode'] == 'A-B':

            hdrinfo = config.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = config.state['hdrinfo'][1]['FILENAME'][0]

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 config.user['load']['flatfile'],
                                 config.state['naps'],
                                 orders, hdrinfo, apertures,
                                 config.state['apradii'],
                                 config.state['plate_scale'],
                                 config.state['slith_pix'],
                                 config.state['slith_arc'],
                                 config.state['slitw_pix'],
                                 config.state['slitw_arc'],
                                 config.state['resolvingpower'],
                                 'um', 'DN s-1', 'microns', 'flux',
                                 config.state['version'],
                                 output_filenames, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

        elif config.state['reductionmode'] == 'A-Sky/Dark':

            x = 1
            
        else:

            print('Reduction Mode Unknown.')

    print(' ')

