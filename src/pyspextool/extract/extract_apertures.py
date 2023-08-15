import numpy as np
import os

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.extract_pointsource_1dxd import extract_pointsource_1dxd
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.io.check import check_parameter
from pyspextool.io.write_apertures_fits import write_apertures_fits
from pyspextool.plot.plot_spectra import plot_spectra


def extract_apertures(qa_show=None, qa_showsize=(10, 6), qa_write=None,
                      fix_bad_pixels=True, use_mean_profile=False,
                      bad_pixel_thresh=extract.state['bad_pixel_thresh'],
                      verbose=None):
    
    """
    User function to extract spectra.

    Parameters
    ----------
    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, default=(6,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    fix_bad_pixels : {True, False}, optional
        Set to True to fix bad pixels using the 2D model profiles.  

    use_mean_profile : {False, True}, optional
        Set to True to use the mean profile at all wavelengths.  

    Returns 
    -------
    list
    A list of str giving the names of the files successfully written to disk.

    """

    #
    # Check if we can proceed.
    #

    if extract.state['parameters_done'] is False:

        message = "extract.state['parameters_done']=False.  "+\
          "Previous steps not complete."        
        raise ValueError(message)

    #
    # Check parameters.  Just do it now to not waste extracting.
    #

    check_parameter('extract_apertures', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_showsize', qa_showsize,
                    'tuple')

    check_parameter('extract_apertures', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'fix_bad_pixels',
                    fix_bad_pixels, 'bool')

    check_parameter('extract_apertures', 'use_mean_profile',
                    use_mean_profile, 'bool')        
    
    check_parameter('extract_apertures', 'verbose',
                    verbose, ['NoneType', 'bool'])

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_write is None:
        qa_write = setup.state['qa_write']

    if qa_show is None:
        qa_show = setup.state['qa_show']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Store user inputs
    #

    extract.extract['qafile'] = qa_write
    extract.extract['qaplot'] = qa_show
    extract.extract['verbose'] = verbose
    extract.extract['fix_bad_pixels'] = fix_bad_pixels
    extract.extract['use_mean_profile'] = use_mean_profile

    #
    # The procedure depends on the extraction type, point or extended source
    #

    xsbginfo = None
    psbginfo = None

    optimal_extraction = True if extract.state['psfradius'] is not None else \
      False
        
    if extract.state['type'] == 'ps':

        #
        # ========================= Point Source ============================
        #

        # Grab which orders are being extracted

        z = (extract.state['psdoorders'] == 1)

        #
        # Now build the various information dictionaries
        #

        # Background first
        
        if extract.state['bgradius'] is not None:

            psbginfo = {'radius': extract.state['bgradius'],
                        'width': extract.state['bgwidth'],
                        'degree': extract.state['bgfitdeg']}

        else:

            psbginfo = None

        # Now bad pixels
        
        if fix_bad_pixels is True:

            badpixelinfo = {'images':np.asarray(extract.state['rectorders'])[z],
                            'usemeanprofile':use_mean_profile,
                            'mask':extract.state['bad_pixel_mask'],
                            'thresh':bad_pixel_thresh}

            # Add the atmospheric transmission if a wavecal file was used.
                
            if extract.load['wavecalfile'] is not None:

                badpixelinfo['atmosphere'] = extract.state['atmosphere']

        else:

            badpixelinfo = None

        # Now optimal extraction
            
        if optimal_extraction is True:

            # Set badpixelinfo to None since it isn't redundant
            
            badpixelinfo = None

            optimalinfo = {'images':np.asarray(extract.state['rectorders'])[z],
                           'psfradius':extract.parameters['psfradius'],
                           'usemeanprofile':use_mean_profile,
                           'mask':extract.state['bad_pixel_mask'],
                           'thresh':bad_pixel_thresh}

            # Add the atmospheric transmission if a wavecal file was used.
                
            if extract.load['wavecalfile'] is not None:

                optimalinfo['atmosphere'] = extract.state['atmosphere']
                    
        else:

            optimalinfo = None
        
        # Do the extraction
        
        spectra = extract_pointsource_1dxd(extract.state['workimage'],
                                           extract.state['varimage'],
                                           extract.state['ordermask'],
                                           extract.state['orders'][z],
                                           extract.state['wavecal'],
                                           extract.state['spatcal'],
                                           extract.state['tracecoeffs'],
                                           extract.state['apradii'],
                                           extract.state['apsigns'],
                                     linmax_bitmask=extract.state['maskimage'],
                                           badpixel_info=badpixelinfo,
                                           optimal_info=optimalinfo,
                                           background_info=psbginfo,
                                           verbose=verbose)

    else:

        #
        # ======================= Extended Source ===========================
        #

        # Grab which orders are being extracted

        z = extract.state['xsdoorders'] == 1

        # Create background information dictionary

        if extract.state['bgregions'] is not None:
            xsbginfo = {'regions': extract.state['bgregions'],
                        'degree': extract.state['bgfitdeg']}

        # Do the extraction

        spectra = extract_extendedsource_1dxd(extract.state['workimage'],
                                              extract.state['varimage'],
                                              extract.state['ordermask'],
                                              extract.state['orders'][z],
                                              extract.state['wavecal'],
                                              extract.state['spatcal'],
                                              extract.state['apertures'][z],
                                              extract.state['apradii'],
                                              bginfo=xsbginfo, verbose=verbose)

    #
    # Write the results to disk
    #

    filenames = write_apertures(spectra, psbginfo=psbginfo, xsbginfo=xsbginfo,
                                optimal_info=optimalinfo,
                                badpixel_info=badpixelinfo, qa_write=qa_write,
                                qa_show=qa_show, qa_showsize=qa_showsize,
                                verbose=verbose)
    
    #
    # Set the done variable
    #

    extract.state['extract_done'] = True

    #
    # Return the filenames written to disk
    #

    return filenames
    

def write_apertures(spectra, psbginfo=None, xsbginfo=None, optimal_info=None,
                    badpixel_info=None, qa_show=None, qa_write=None,
                    qa_showsize=(10, 6), verbose=None):

    """
    To write extracted spectra to disk.

    Parameters
    ----------
    spectra : list
        An (norders*naps,) list of spectral planes of the extracted target.

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

    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, default=(10,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    Returns 
    -------
    list
    A list of str giving the names of the files successfully written to disk.

    Also writes FITS files to disk.

    """

    #
    # Get the plotting ready
    #

    qafileinfo = {'figsize': (11, 8.5), 'filepath': setup.state['qa_path'],
                  'filename': extract.state['qafilename'],
                  'extension': setup.state['qa_extension']}

    # Get the wavelengh calibration info if used.

    if extract.load['wavecalfile'] is not None:

        wavecalinfo = {'file': extract.load['wavecalfile'],
                       'wavecaltype': extract.state['wavecaltype'],
                       'wavetype': 'vacuum'}

    else:

        wavecalinfo = None

    # Things proceed depending on the extraction mode

    if extract.state['type'] == 'ps':

        #
        # ========================= Point Source ============================
        #

        # Grab things all modes need

        xranges = extract.state['xranges'][extract.state['psdoorders']]
        orders = extract.state['orders'][extract.state['psdoorders']]
        apertures = extract.state['apertures'][extract.state['psdoorders']]
        norders = np.sum(extract.state['psdoorders'])

        # Now go mode by mode

        if extract.state['reductionmode'] == 'A':

            # Get file information

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'
            output_fullpath = extract.state['output_files'][0]

            # Write ther results

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 setup.state['version'],
                                 output_fullpath,
                                 wavecalinfo=wavecalinfo,
                                 psbginfo=psbginfo,
                                 optimal_info=optimal_info,
                                 badpixel_info=badpixel_info,
                                 verbose=verbose)

            # Plot the spectra

            filename = os.path.basename(output_fullpath)            
            if qa_show is True:
                
                number = plot_spectra(output_fullpath + '.fits',
                                      title=filename+'.fits',                                                         plot_size=qa_showsize,
                                      flag_linearity=True,
                            plot_number=extract.state['spectra_a_plotnum'])

                extract.state['spectra_a_plotnum'] = number

            if qa_write is True:

                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits',
                             flat_linearity=True,
                             file_info=qafileinfo)

            # Store name of successfully written files 
                
            filenames = filename

        elif extract.state['reductionmode'] == 'A-B':

            filenames = []

            # Are there positive apertures?

            z_pos = extract.state['apsigns'] == 1
            pos_naps = int(np.sum(z_pos))
            if pos_naps != 0:

                # Get file information

                hdrinfo = extract.state['hdrinfo'][0]
                aimage = hdrinfo['FILENAME'][0]
                sky = extract.state['hdrinfo'][1]['FILENAME'][0]
                output_fullpath = extract.state['output_files'][0]

                # Now get the indices for the spectra array

                full_apsign = np.tile(extract.state['apsigns'], norders)
                z = (np.where(full_apsign == 1))[0]
                pos_spectra = [spectra[i] for i in z]

                write_apertures_fits(pos_spectra, xranges, aimage, sky,
                                     extract.load['flatfile'], pos_naps,
                                     orders, hdrinfo, apertures,
                                     extract.state['apradii'],
                                     extract.state['plate_scale'],
                                     extract.state['slith_pix'],
                                     extract.state['slith_arc'],
                                     extract.state['slitw_pix'],
                                     extract.state['slitw_arc'],
                                     extract.state['resolvingpower'],
                                     'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                     'Wavelength ($\mu$m)',
                                     'Count Rate (DN s$^{-1}$)',
                                     setup.state['version'],
                                     output_fullpath,
                                     wavecalinfo=wavecalinfo,
                                     psbginfo=psbginfo,
                                     optimal_info=optimal_info,
                                     badpixel_info=badpixel_info,
                                     verbose=verbose)

                # Plot the spectra

                a_filename = os.path.basename(output_fullpath)

                if qa_show is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                 title=a_filename+'.fits',
                                 plot_size=qa_showsize,
                                 flag_linearity=True,
                                 plot_number=extract.state['spectra_a_plotnum'])

                    extract.state['spectra_a_plotnum'] = number

                    
                if qa_write is True:

                    qafileinfo['filename'] = a_filename
                    plot_spectra(output_fullpath + '.fits',
                                 flag_linearity=True,                 
                                 file_info=qafileinfo)

                filenames.append(a_filename+'.fits')
                
            # Are there negative apertures?

            z_neg = extract.state['apsigns'] == -1
            neg_naps = int(np.sum(z_neg))
            if neg_naps != 0:

                # Get file information

                hdrinfo = extract.state['hdrinfo'][1]
                aimage = hdrinfo['FILENAME'][1]
                sky = extract.state['hdrinfo'][0]['FILENAME'][0]
                output_fullpath = extract.state['output_files'][1]

                # Now get the indices for the spectra array

                full_apsign = np.tile(extract.state['apsigns'], norders)
                z = (np.where(full_apsign == -1))[0]
                neg_spectra = [spectra[i] for i in z]

                write_apertures_fits(neg_spectra, xranges, aimage, sky,
                                     extract.load['flatfile'],
                                     neg_naps, orders, hdrinfo, apertures,
                                     extract.state['apradii'],
                                     extract.state['plate_scale'],
                                     extract.state['slith_pix'],
                                     extract.state['slith_arc'],
                                     extract.state['slitw_pix'],
                                     extract.state['slitw_arc'],
                                     extract.state['resolvingpower'],
                                     'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                     'Wavelength ($\mu$m)',
                                     'Count Rate (DN s$^{-1}$)',
                                     setup.state['version'],
                                     output_fullpath,
                                     wavecalinfo=wavecalinfo,
                                     psbginfo=psbginfo,
                                     verbose=verbose)

                # Plot the spectra

                b_filename = os.path.basename(output_fullpath)
                if qa_show is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                          title=b_filename+'.fits',
                                          flag_linearity=True,                
                                          plot_size=qa_showsize,
                                plot_number=extract.state['spectra_b_plotnum'])
                    extract.state['spectra_b_plotnum'] = number
                    
                if qa_write is True:

                    qafileinfo['filename'] = b_filename
                    plot_spectra(output_fullpath + '.fits',flag_linearity=True,
                                 file_info=qafileinfo)

                filenames.append(b_filename+'.fits')
                
            # Store name of successfully written files 
            
        elif extract.state['reductionmode'] == 'A-Sky/Dark':

            x = 1

        else:

            print('Reduction Mode Unknown.')

    else:

        #
        # ======================= Extended Source ===========================
        #

        # Grab things all modes need

        xranges = extract.state['xranges'][extract.state['psdoorders']]
        orders = extract.state['orders'][extract.state['psdoorders']]
        apertures = extract.state['apertures'][extract.state['psdoorders']]

        if extract.state['reductionmode'] == 'A':

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'
            output_fullpath = extract.state['output_files'][0]
            
            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 setup.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

            # Plot the spectra

            if qa_show is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_showsize)

            if qa_write is True:
                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-B':

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = extract.state['hdrinfo'][1]['FILENAME'][0]
            output_fullpath = extract.state['output_files'][0]            

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 extract.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

            # Plot the spectra

            if qa_show is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_showsize)

            if qa_write is True:
                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-Sky/Dark':

            x = 1

        else:

            print('Reduction Mode Unknown.')

    print(' ')

    return filenames
    
