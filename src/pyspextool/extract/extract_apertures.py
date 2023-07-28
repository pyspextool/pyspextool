import numpy as np
import os
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.extract_pointsource_1dxd import extract_pointsource_1dxd
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.io.check import check_parameter
from pyspextool.io.write_apertures_fits import write_apertures_fits
from pyspextool.plot.plot_spectra import plot_spectra
from pyspextool.io.check import check_parameter
from pyspextool.utils.split_text import split_text


def extract_apertures(qa_plot:bool=None, qa_plotsize:tuple=(10, 6), qa_file:bool=None,
                      fix_bad_pixels:bool=True, use_mean_profile:bool=False,
                      bad_pixel_thresh=extract.state['bad_pixel_thresh'],
                      verbose:bool=None):
    
    """
    User function to extract spectra.

    Parameters
    ----------
    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(6,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.state['qa_file'] in the
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
    None

    """

    #
    # Check if we can proceed.
    #

    if extract.state['parameters_done'] is False:
        message = 'Previous steps not completed.'
        print(message)
        return

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_file is None:
        qa_file = setup.state['qa_file']

    if qa_plot is None:
        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Check parameters.  Just do it now to not waste extracting.
    #

    check_extract_aperture_parameters(qa_plot=qa_plot, qa_plotsize=qa_plotsize, qa_file=qa_file, 
                                      fix_bad_pixels=fix_bad_pixels, use_mean_profile=use_mean_profile,
                                      bad_pixel_thresh=bad_pixel_thresh, verbose=verbose)

    #
    # Store user inputs
    #

    extract.extract['qafile'] = qa_file
    extract.extract['qaplot'] = qa_plot
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
        spectra = extract_pointsource(fix_bad_pixels=fix_bad_pixels, use_mean_profile=use_mean_profile,
                                      bad_pixel_thresh=bad_pixel_thresh, optimal_extraction=optimal_extraction,
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

    write_apertures(spectra, psbginfo=psbginfo, xsbginfo=xsbginfo,
                    optimal_info=optimalinfo, badpixel_info=badpixelinfo,
                    qa_file=qa_file, qa_plot=qa_plot, qa_plotsize=qa_plotsize,
                    verbose=verbose)

    #
    # Set the done variable
    #

    extract.state['extract_done'] = True


def check_extract_aperture_parameters(qa_plot, qa_plotsize, qa_file, fix_bad_pixels,
                              use_mean_profile, bad_pixel_thresh, verbose):
    
    check_parameter('extract_apertures', 'qa_plot', qa_plot,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_plotsize', qa_plotsize,
                    'tuple')

    check_parameter('extract_apertures', 'qa_file', qa_file,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'fix_bad_pixels',
                    fix_bad_pixels, 'bool')

    check_parameter('extract_apertures', 'use_mean_profile',
                    use_mean_profile, 'bool')        
    
    check_parameter('extract_apertures', 'verbose',
                    verbose, ['NoneType', 'bool'])
    
    return


def extract_pointsource(fix_bad_pixels=None, use_mean_profile=None, bad_pixel_thresh=None, 
                        optimal_extraction=None, verbose=None):
    '''
     ========================= Point Source ============================
    
    Returns
    -------
    list
        A (norders,) list. Each entry is a (4, nwave) numpy.ndarray where:
          wave = (0,:)
          intensity = (1,:)
          uncertainty = (2,:)
          flags = (3,:)
    '''

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


def write_apertures(spectra, psbginfo=None, xsbginfo=None, optimal_info=None,
                    badpixel_info=None, qa_plot=None, qa_file=None,
                    qa_plotsize=(10, 6), verbose=None):

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

    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(10,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.state['qa_file'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    Returns 
    -------
    None

       Writes FITS files to disk.

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
            if qa_plot is True:
                
                number = plot_spectra(output_fullpath + '.fits',
                                      title=filename+'.fits',                                                         plot_size=qa_plotsize,
                                      flag_linearity=True,
                            plot_number=extract.state['spectra_a_plotnum'])

                extract.state['spectra_a_plotnum'] = number

            if qa_file is True:

                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits',
                             flat_linearity=True,
                             file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-B':

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

                filename = os.path.basename(output_fullpath)
                if qa_plot is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                 title=filename+'.fits',
                                 plot_size=qa_plotsize,
                                 flag_linearity=True,
                                 plot_number=extract.state['spectra_a_plotnum'])

                    extract.state['spectra_a_plotnum'] = number

                    
                if qa_file is True:

                    qafileinfo['filename'] = filename
                    plot_spectra(output_fullpath + '.fits',
                                 flag_linearity=True,                 
                                 file_info=qafileinfo)

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

                filename = os.path.basename(output_fullpath)
                if qa_plot is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                          title=filename+'.fits',
                                          flag_linearity=True,                
                                          plot_size=qa_plotsize,
                                plot_number=extract.state['spectra_b_plotnum'])
                    extract.state['spectra_b_plotnum'] = number
                    
                if qa_file is True:

                    qafileinfo['filename'] = filename
                    plot_spectra(output_fullpath + '.fits',flag_linearity=True,
                                 file_info=qafileinfo)

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

            if qa_plot is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_plotsize)

            if qa_file is True:
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

            if qa_plot is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_plotsize)

            if qa_file is True:
                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-Sky/Dark':

            x = 1

        else:

            print('Reduction Mode Unknown.')

    print(' ')



def write_apertures_fits(spectra: list, xranges: np.ndarray, aimage: str, sky: str, 
                         flat: str, naps: int, 
                         orders: list[int | list[int] | np.ndarray],
                         header_info: dict, aperture_positions: np.ndarray, 
                         aperture_radii: list[int | float | np.ndarray],
                         plate_scale: float, slith_pix: float, slith_arc: float, slitw_pix: float,
                         slitw_arc: float, resolving_power: float, xunits: str, yunits: str,
                         latex_xunits: str, latex_yunits: str, latex_xlabel: str,
                         latex_ylabel: str, version, output_fullpath: str,
                         wavecalinfo: list[dict | None] = None, 
                         psbginfo=None, xsbginfo=None,
                         optimal_info=None, 
                         badpixel_info: list[dict | None]=None,
                         lincormax: None=None, overwrite: bool=True, verbose:bool=True):
                         

    check_apertures_parameters(aimage, sky, flat, naps, orders,
                           header_info, aperture_positions, aperture_radii,
                         plate_scale, slith_pix, slith_arc, slitw_pix,
                         slitw_arc, resolving_power, xunits, yunits,
                         latex_xunits, latex_yunits, latex_xlabel,
                         latex_ylabel, version, output_fullpath,
                         wavecalinfo=wavecalinfo, psbginfo=psbginfo, xsbginfo=xsbginfo,
                         optimal_info=optimal_info, badpixel_info=badpixel_info,
                         lincormax=lincormax, overwrite=overwrite, verbose=verbose)

    aperature_array = create_aperture_array(spectra, naps, orders)

    aperature_header = create_aperature_header(aimage, sky, flat, naps, orders,
                                               header_info, aperture_positions, aperture_radii,
                         plate_scale, slith_pix, slith_arc, slitw_pix,
                         slitw_arc, xunits, yunits,
                         latex_xunits, latex_yunits, latex_xlabel,
                         latex_ylabel, version, output_fullpath,
                         wavecalinfo, xsbginfo,
                         optimal_info, badpixel_info)

    write_fits(aperature_array, aperature_header, output_fullpath, xsbginfo=xsbginfo, 
               overwrite=overwrite, verbose=verbose)

    return


def check_apertures_parameters(spectra, xranges, aimage, flat, naps, orders,
                         header_info, aperture_positions, aperture_radii,
                         plate_scale, slith_pix, slith_arc, slitw_pix,
                         slitw_arc, resolving_power, xunits, yunits,
                         latex_xunits, latex_yunits, latex_xlabel,
                         latex_ylabel, version, output_fullpath,
                         wavecalinfo=None, psbginfo=None, xsbginfo=None,
                         optimal_info=None, badpixel_info=None,
                         lincormax=None, overwrite=True, verbose=True):
    
    check_parameter('write_apertures_fits', 'spectra', spectra, 'list')

    check_parameter('write_apertures_fits', 'xranges', xranges, 'ndarray')

    check_parameter('write_apertures_fits', 'aimage', aimage, 'str')

    check_parameter('write_apertures_fits', 'flat', flat, 'str')

    check_parameter('write_apertures_fits', 'naps', naps, 'int')

    check_parameter('write_apertures_fits', 'orders', orders,
                    ['int', 'list', 'ndarray'])

    check_parameter('write_apertures_fits', 'header_info', header_info,
                    'dict')

    check_parameter('write_apertures_fits', 'aperture_positions',
                    aperture_positions, 'ndarray')

    check_parameter('write_apertures_fits', 'aperture_radii', aperture_radii,
                    ['int', 'float', 'ndarray'])

    check_parameter('write_apertures_fits', 'plate_scale', plate_scale, 'float')

    check_parameter('write_apertures_fits', 'slith_pix', slith_pix, 'float')

    check_parameter('write_apertures_fits', 'slith_arc', slith_arc, 'float')

    check_parameter('write_apertures_fits', 'slitw_pix', slitw_pix, 'float')

    check_parameter('write_apertures_fits', 'slitw_arc', slitw_arc, 'float')

    check_parameter('write_apertures_fits', 'resolving_power',
                    resolving_power, 'float')

    check_parameter('write_apertures_fits', 'xunits', xunits, 'str')

    check_parameter('write_apertures_fits', 'yunits', yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xunits', latex_xunits, 'str')

    check_parameter('write_apertures_fits', 'latex_yunits', latex_yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xlabel', latex_xlabel, 'str')

    check_parameter('write_apertures_fits', 'latex_ylabel', latex_ylabel, 'str')

    check_parameter('write_apertures_fits', 'version', version, 'str')

    check_parameter('write_apertures_fits', 'output_fullpath',
                    output_fullpath, 'str')

    check_parameter('write_apertures_fits', 'psbginfo', psbginfo,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'xsbginfo', xsbginfo,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'optimal_info', optimal_info,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'badpixel_info', badpixel_info,
                    ['NoneType', 'dict'])    

    check_parameter('write_apertures_fits', 'wavecalinfo', wavecalinfo,
                    ['dict', 'NoneType'])

    check_parameter('write_apertures_fits', 'lincormax', lincormax,
                    'NoneType')

    check_parameter('write_apertures_fits', 'overwrite', overwrite,
                    'bool')

    check_parameter('write_apertures_fits', 'verbose', verbose,
                    'bool')

    return


def create_aperture_array(spectra,  naps, orders) -> np.ndarray:

    """
    To write a spextool spectral FITS file to disk

    Parameters
    ----------

    Returns
    -------
    None.  Writes FITS files to disk.


    """

    orders = np.asarray(orders)
    norders = len(orders)


    #
    # Get the list of spectra into an array format
    #

    # Determine the maximum size of each aperture array

    npixels = []
    for slice in spectra:
        npixels.append(np.shape(slice)[1])

    max_npixels = np.max(np.asarray(npixels))

    # Now create arrays into which the slices will be placed

    array = np.full((norders * naps, 4, max_npixels), np.nan)

    # Now fill in the arrays

    l = 0
    for slice in spectra:
        array[l, :, 0:npixels[l]] = slice
        l += 1

    return array


def create_aperature_header(aimage, sky, flat, naps, orders,
                         header_info, aperture_positions, aperture_radii,
                         plate_scale, slith_pix, slith_arc, slitw_pix,
                         slitw_arc, xunits, yunits,
                         latex_xunits, latex_yunits, latex_xlabel,
                         latex_ylabel, version, output_fullpath,
                         wavecalinfo, xsbginfo,
                         optimal_info, badpixel_info):


    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the original file header keywords.  Cut out any history first
    # if need be.

    is_history_present = 'HISTORY' in header_info

    if is_history_present is True:

        old_history = header_info['HISTORY']

        # remove it from the header_info dictionary

        header_info.pop('HISTORY')

    keys = list(header_info.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (header_info[keys[i]][0], header_info[keys[i]][1])

    # Add spextool keywords

    hdr['MODULE'] = ('extract', ' Creation module')
    hdr['VERSION'] = (version, ' pySpextool version')
    hdr['AIMAGE'] = (aimage, ' A image')
    hdr['SKYORDRK'] = (sky, ' Sky or dark image')
    hdr['FLAT'] = (flat, ' Flat field image')

    if wavecalinfo is not None:
        hdr['WAVECAL'] = (wavecalinfo['file'], ' Wavecal file')
        hdr['WCTYPE'] = (wavecalinfo['wavecaltype'],
                         ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (wavecalinfo['wavetype'], ' Wavelength type')

    else:

        hdr['WAVECAL'] = (None, ' Wavecal file')
        hdr['WCTYPE'] = (None, ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (None, ' Wavelength type')
        
    orders = np.asarray(orders)
    norders = len(orders)

    hdr['NORDERS'] = (norders, ' Number of orders')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), ' Orders')
    hdr['NAPS'] = (naps, ' Number of apertures')

    hdr['PLTSCALE'] = (plate_scale, ' Plate scale (arcseconds pixel-1)')
    hdr['SLTH_ARC'] = (slith_arc, ' Nominal slit height (arcseconds)')
    hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
    hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit height (pixels)')
    hdr['SLTW_PIX'] = (slitw_pix, ' Slit width (pixels)')

    # Add the aperture positions

    for i in range(norders):
        name = 'APOSO' + str(orders[i]).zfill(3)
        comment = ' Aperture positions (arcseconds) for order ' + \
                  str(orders[i]).zfill(3)
        val = ','.join([str(round(elem, 2)) for elem in aperture_positions[i]])
        hdr[name] = (val, comment)

    if optimal_info is not None:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')
        hdr['THRESH'] = (optimal_info['thresh'], ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (True, ' Optimal extraction?')
        hdr['SUMEXT'] = (False, 'Sum extraction?')
        hdr['PSFRAD'] = (optimal_info['psfradius'], ' PSF radius (arcseconds)')

                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        
        hdr['THRESH'] = (None, ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (False, ' Optimal extraction?')
        hdr['SUMEXT'] = (True, ' Sum extraction?') 
        hdr['PSFRAD'] = (None, ' PSF radius (arcseconds)')


    if badpixel_info is not None:

        hdr['BDPXFIX'] = (True, ' Bad pixels fixed?')        
        hdr['THRESH'] = (badpixel_info['thresh'], ' Bad pixel sigma threshold')
                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        
                
    # Add the aperture radii

    if isinstance(aperture_radii, np.ndarray):

        val = ','.join([str(elem) for elem in aperture_radii])
        hdr['APRADII'] = (val, ' Aperture radii (arcseconds)')

    else:

        hdr['APRADII'] = (aperture_radii, ' Aperture radii (arcseconds)')

        # Add the background info

    if xsbginfo is not None:

        tmplist = []
        for val in xsbginfo['regions']:
            region = str(val[0]) + '-' + str(val[1])
            tmplist.append(region)

        hdr['BGREGS'] = (','.join(tmplist),
                         ' Background regions (arcseconds)')
        hdr['BGDEGREE'] = (xsbginfo['degree'],
                           ' Background polynomial fit degree')

    hdr['XUNITS'] = (xunits, ' Units of the x axis')
    hdr['YUNITS'] = (yunits, ' Units of the y axis')

    hdr['LXUNITS'] = (latex_xunits, ' LateX units of the x axis')
    hdr['LYUNITS'] = (latex_yunits, ' LateX units of the y axis')

    hdr['LXLABEL'] = (latex_xlabel, ' LateX x axis label')
    hdr['LYLABEL'] = (latex_ylabel, ' LateX Y axis label')

    # Now add the HISTORY

    history = 'Spextool FITS files contain an array of size ' \
              '[nwaves,4,norders*naps]. The ith image (array[*,*,i]) ' \
              'contains the data for a single extraction aperture within ' \
              'an order such that, lambda=array[*,0,i], flux=array[*,1,i], ' \
              'uncertainty=array[*,2,i],flag=array[*,3,i].  The zeroth ' \
              'image (array[*,*,0]) contains the data for the aperture in ' \
              'the order closest to the bottom of the detector that is ' \
              'closest to the bottom of the slit (i.e. also closest to the ' \
              'bottom of the detector).  Moving up the detector, the FITS ' \
              'array is filled in with subsequent extraction apertures.  ' \
              'If no orders have been deselected in the extraction process, ' \
              'the contents of the ith aperture in order j can be found as ' \
              'follows: lambda=array[*,0,{j-min(orders)}*naps + (i-1)], ' \
              'flux=array[*,1,{j-min(orders)}*naps + (i-1)], ' \
              'uncertainty=array[*,2,{j-min(orders)}*naps + (i-1)], ' \
              'flag=array[*,3,{j-min(orders)}*naps + (i-1)].'

    history = split_text(history, length=65)

    if is_history_present is True:

        # Join the two histories

        history = old_history + [' '] + history
    
    for hist in history:
        hdr['HISTORY'] = hist


    # Set the file name for the spectra file

    hdr['FILENAME'] = (os.path.basename(output_fullpath)+'.fits', ' File name')

    return hdr


def write_fits(spectrum, hdr, output_fullpath, overwrite=False, verbose=True):

    # Put the header data in the primary HDU
    hdu0 = fits.PrimaryHDU(header = hdr)
    
    # Put the spectrum in the second HDU
    hdu1 = fits.BinTableHDU(data=spectrum)
    hdu1.header['EXTNAME'] = 'SPECTRUM'

    spectrum_mef =  fits.HDUList([hdu0, hdu1])

    fits_filename = output_fullpath + '.fits'

    try:
        spectrum_mef.writeto(fits_filename, overwrite=overwrite)
        if verbose:
            print(f'Wrote {os.path.basename(fits_filename)} to disk.')
    except:
        raise Exception(f'Could not write {os.path.basename(fits_filename)} to disk.')

    return
