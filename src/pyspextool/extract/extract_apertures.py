import numpy as np
from os.path import join, basename as osbasename
import logging

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.extraction import extract_1dxd, write_apertures_fits
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.plot.plot_spectra import plot_spectra
from pyspextool.pyspextoolerror import pySpextoolError


def extract_apertures(fix_badpixels:bool=True,
                      use_meanprofile:bool=False,
                      badpixel_thresh:int | float=7,
                      verbose:bool=None,
                      qa_show:bool=None,
                      qa_showscale:float | int=None,
                      qa_showblock:bool=None,
                      qa_write:bool=None):
    
    """
    Command line function to extract spectra.

    Parameters
    ----------
    fix_bad_pixels : {True, False}, optional
        Set to True to fix bad pixels using the 2D model profiles.
        Set to False to not fix bad pixels using the 2D model profile.

        If optimal extraction is done, bad pixels are ignored and so this
        keyword has no effect.

    use_meanprofile : {False, True}, optional
        Set to True to use the mean profile for bad pixel identification,
        fixing, and optimal extraction instead of the 2D model profile.

        Set to False to use the 2D model profile for bad pixel identification,
        fixing, and optimal extraction.

    badpixel_thresh : int or float, default 6
        The sigma threshold over which to pixels are identified as bad.
    
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
    list
    A list of str giving the names of the files successfully written to disk.

    """

    #
    # Check if we can proceed.
    #

    if extract.state['parameters_done'] is False:

        message = "Previous steps not complete.  Please run "+ \
            'define_aperture_parameters.'
        raise pySpextoolError(message)

    #
    # Check parameters and QA keywords
    #

    check_parameter('extract_apertures', 'fix_badpixels', fix_badpixels, 'bool')

    check_parameter('extract_apertures', 'use_meanprofile', use_meanprofile,
                    'bool')    

    check_parameter('extract_apertures', 'badpixel_thresh', badpixel_thresh,
                   ['int', 'float'])    
    
    check_parameter('extract_apertures', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('extract_apertures', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)    

    #
    # Get set up
    #

    optimal_extraction = True if extract.state['psf_radius'] else False 

    z = (extract.state['doorders'] == 1)

    # Set up the bad pixel information
    
    if fix_badpixels:

        badpixelinfo = {'images':np.asarray(extract.state['rectorders'])[z],
                        'usemeanprofile':use_meanprofile,
                        'mask':extract.state['bad_pixel_mask'],
                        'thresh':badpixel_thresh}
        
        # Add the atmospheric transmission if a wavecal file was used.
        
        if extract.state['wavecalfile'] is not None:
            
            badpixelinfo['atmosphere'] = extract.state['atmosphere']
            
        else:
            
            badpixelinfo['atmosphere'] = None

    else:

        badpixelinfo = None
            
            
    # Set up the optimal extraction information
        
    if optimal_extraction:

        # Set badpixelinfo to None since it is now redundant
            
        badpixelinfo = None
        
        optimalinfo = {'images':np.asarray(extract.state['rectorders'])[z],
                       'psfradius':extract.state['psf_radius'],
                       'usemeanprofile':use_meanprofile,
                       'mask':extract.state['bad_pixel_mask'],
                       'thresh':badpixel_thresh}
        
        # Add the atmospheric transmission if a wavecal file was used.
        
        if extract.state['wavecalfile'] is not None:
            
            optimalinfo['atmosphere'] = extract.state['atmosphere']
            
        else:
            
            optimalinfo['atmosphere'] = None
            
    else:
        
        optimalinfo = None

    #
    # Do the extraction
    #

    naps = len(extract.state['average_aperturesigns'])
    norders = len(extract.state['tracecoeffs']) // naps

    message = ' Optimal ' if optimal_extraction is True else ' Sum '
    
    message = message+('extracting '+str(naps)+' apertures in '+ \
                       str(norders)+' orders')

    if extract.state['bg_annulus'] is not None or \
       extract.state['bg_regions'] is not None:
        
        logging.info(message + ' (with background subtraction).')
            
    else:

        logging.info(message + ' (without background subtraction).')


    z = extract.state['doorders'] == True
    extract_orders = extract.state['orders'][z]
              
    spectra = extract_1dxd(extract.state['workimage'],
                           extract.state['varimage'],
                           extract.state['ordermask'],
                           extract.state['wavecal'],
                           extract.state['spatcal'],
                           extract_orders,
                           extract.state['tracecoeffs'],
                           extract.state['aperture_radii'],
                           extract.state['average_aperturesigns'],
                           linmax_bitmask=extract.state['maskimage'],
                           bg_annulus=extract.state['bg_annulus'],
                           bg_regions=extract.state['bg_regions'],
                           bg_fitdegree=extract.state['bg_fitdegree'],
                           badpixel_info=badpixelinfo,
                           optimal_info=optimalinfo,
                           progressbar=qa['verbose'])
        
    #
    # Write the results to disk
    #

    if extract.state['wavecalfile'] is not None:

        wavecalinfo = {'file': extract.state['wavecalfile'],
                       'wavecaltype': extract.state['wavecaltype'],
                       'wavetype': 'vacuum'}

    else:

        wavecalinfo = None

    # Grab things all modes need

    z = extract.state['doorders'] == 1    
    xranges = extract.state['xranges'][z]
    apertures = extract.state['aperture_positions'][z,:]
    norders = np.sum(z)
    
    # Now go mode by mode

    if extract.state['reductionmode'] == 'A':
        
        # Get file information

        hdrinfo = extract.state['hdrinfo'][0]
        aimage = hdrinfo['FILENAME'][0]
        sky = 'None'
        abeam_fullpath = join(setup.state['proc_path'],
                              extract.state['output_files'][0]+'.fits')
        bbeam_fullpath = None
        
        # Write ther results
        
        write_apertures_fits(spectra,
                             xranges,
                             aimage,
                             sky,
                             extract.state['flat_filename'],
                             extract.state['flat_fielded'],                    
                             extract.state['naps'],
                             extract_orders,
                             hdrinfo,
                             apertures,
                             extract.state['aperture_radii'],
                             extract.state['plate_scale'],
                             extract.state['slith_pix'],
                             extract.state['slith_arc'],
                             extract.state['slitw_pix'],
                             extract.state['slitw_arc'],
                             extract.state['resolvingpower'],
                             extract.state['xunits'],
                             extract.state['yunits'],
                             extract.state['latex_xunits'],
                             extract.state['latex_yunits'],
                             extract.state['latex_xlabel'],
                             extract.state['latex_ylabel'],
                             extract.state['latex_ulabel'],                  
                             setup.state['version'],
                             abeam_fullpath,
                             bg_annulus=extract.state['bg_annulus'],
                             bg_regions=extract.state['bg_regions'],
                             bg_fitdegree=extract.state['bg_fitdegree'],
                             wavecalinfo=wavecalinfo,
                             optimal_info=optimalinfo,
                             badpixel_info=badpixelinfo,
                             verbose=qa['verbose'])

    # Now go mode by mode

    if extract.state['reductionmode'] == 'A-B':

        abeam_fullpath = None
        bbeam_fullpath = None        
        
        # Are there positive apertures?

        z_pos = extract.state['average_aperturesigns'] == 1
        pos_naps = int(np.sum(z_pos))
        if pos_naps != 0:
            
            # Get file information
            
            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = extract.state['hdrinfo'][1]['FILENAME'][0]
            abeam_fullpath = join(setup.state['proc_path'],
                                  extract.state['output_files'][0]+'.fits')

            # Now get the indices for the spectra array
            
            full_apsign = np.tile(extract.state['average_aperturesigns'], 
                                  norders)
            z = (np.where(full_apsign == 1))[0]
            pos_spectra = [spectra[i] for i in z]
            
            # Write ther results
        
            write_apertures_fits(pos_spectra,
                                 xranges,
                                 aimage,
                                 sky,
                                 extract.state['flat_filename'],
                                 extract.state['flat_fielded'],                
                                 pos_naps,
                                 extract_orders,
                                 hdrinfo,
                                 apertures,
                                 extract.state['aperture_radii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 extract.state['xunits'],
                                 extract.state['yunits'],
                                 extract.state['latex_xunits'],
                                 extract.state['latex_yunits'],
                                 extract.state['latex_xlabel'],
                                 extract.state['latex_ylabel'],
                                 extract.state['latex_ulabel'],
                                 setup.state['version'],
                                 abeam_fullpath,
                                 bg_annulus=extract.state['bg_annulus'],
                                 bg_regions=extract.state['bg_regions'],
                                 bg_fitdegree=extract.state['bg_fitdegree'],
                                 wavecalinfo=wavecalinfo,
                                 optimal_info=optimalinfo,
                                 badpixel_info=badpixelinfo,
                                 verbose=qa['verbose'])

        # Are there negative apertures?

        z_neg = extract.state['average_aperturesigns'] == -1
        neg_naps = int(np.sum(z_neg))
        if neg_naps != 0:

            # Get file information
            
            hdrinfo = extract.state['hdrinfo'][1]
            aimage = hdrinfo['FILENAME'][1]
            sky = extract.state['hdrinfo'][0]['FILENAME'][0]
            bbeam_fullpath = join(setup.state['proc_path'],
                                  extract.state['output_files'][1]+'.fits')

            # Now get the indices for the spectra array
            
            full_apsign = np.tile(extract.state['average_aperturesigns'], 
                                  norders)
            z = (np.where(full_apsign == -1))[0]
            neg_spectra = [spectra[i] for i in z]
            
            # Write ther results
        
            write_apertures_fits(neg_spectra,
                                 xranges,
                                 aimage,
                                 sky,
                                 extract.state['flat_filename'],
                                 extract.state['flat_fielded'],
                                 neg_naps,
                                 extract_orders,
                                 hdrinfo,
                                 apertures,
                                 extract.state['aperture_radii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 extract.state['xunits'],
                                 extract.state['yunits'],
                                 extract.state['latex_xunits'],
                                 extract.state['latex_yunits'],
                                 extract.state['latex_xlabel'],
                                 extract.state['latex_ylabel'],
                                 extract.state['latex_ulabel'],                 
                                 setup.state['version'],
                                 bbeam_fullpath,
                                 bg_annulus=extract.state['bg_annulus'],
                                 bg_regions=extract.state['bg_regions'],
                                 bg_fitdegree=extract.state['bg_fitdegree'],
                                 wavecalinfo=wavecalinfo,
                                 optimal_info=optimalinfo,
                                 badpixel_info=badpixelinfo,
                                 verbose=qa['verbose'])

    #
    # Report results
    #

    if qa['verbose'] is True:

        message = ' Wrote file(s) '

        files = []
        if abeam_fullpath is not None:

            files.append(osbasename(abeam_fullpath))

        if bbeam_fullpath is not None:

            files.append(osbasename(bbeam_fullpath))

        message += ', '.join(files)+' to proc/.\n'

        logging.info(message)
            
            
    #
    # Do the QA plotting
    #
    
    if qa['show'] is True:

        figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                       setup.plots['landscape_size'][1]*qa['showscale'])

        font_size = setup.plots['font_size']*qa['showscale']

    # Get the blocking right if you are showing both beams
        
        if abeam_fullpath is not None and bbeam_fullpath is not None:

            a_showblock = False
            b_showblock = qa['showblock']

        else:  

            a_showblock = qa['showblock']
            b_showblock = qa['showblock']
        
        if abeam_fullpath is not None:
            
            plot_spectra(abeam_fullpath,
                         plot_number=setup.plots['abeam_spectra'],
                         figure_size=figure_size,
                         font_size=font_size,
                         title=osbasename(abeam_fullpath),
                         showblock=a_showblock)
            
        if bbeam_fullpath is not None:

            plot_spectra(bbeam_fullpath,
                         plot_number=setup.plots['bbeam_spectra'],
                         figure_size=figure_size,
                         font_size=font_size,
                         title=osbasename(abeam_fullpath),
                         showblock=b_showblock)                         
                    
    if qa['write'] is True:

        if abeam_fullpath is not None:

            fullpath = join(setup.state['qa_path'],
                            extract.state['output_files'][0]+\
                            setup.state['qa_extension'])

            plot_spectra(abeam_fullpath,
                         figure_size=setup.plots['landscape_size'],
                         font_size=setup.plots['font_size'],
                         spectrum_linewidth=setup.plots['spectrum_linewidth'],
                         spine_linewidth=setup.plots['spectrum_linewidth'],
                         title=osbasename(abeam_fullpath),
                         output_fullpath=fullpath)
            
        if bbeam_fullpath is not None:

            fullpath = join(setup.state['qa_path'],
                            extract.state['output_files'][1]+\
                            setup.state['qa_extension'])

            plot_spectra(bbeam_fullpath,
                         figure_size=setup.plots['landscape_size'],
                         font_size=setup.plots['font_size'],
                         spectrum_linewidth=setup.plots['spectrum_linewidth'],
                         spine_linewidth=setup.plots['spectrum_linewidth'],
                         title=osbasename(abeam_fullpath),
                         output_fullpath=fullpath)    

    #
    # Set the done variable
    #

    extract.state['extract_done'] = True
