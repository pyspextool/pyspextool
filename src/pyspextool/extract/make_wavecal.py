import importlib
import os
import numpy as np
import logging

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.extract.flat import read_flat_fits
from pyspextool.io.files import files_to_fullpath
from pyspextool.extract import wavecal
from pyspextool.extract.extraction import extract_1dxd
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils.math import scale_data_stack
from pyspextool.utils.math import median_data_stack
from pyspextool.utils.math import combine_flag_stack
from pyspextool.pyspextoolerror import pySpextoolError

def make_wavecal(arc_files:str | list,
                 flat_file:str,
                 output_filename:str,
                 sky_files:str=None,
                 linearity_correction=True,
                 detector_info:dict=None,
                 use_stored_solution:bool=False,
                 verbose:bool=None,
                 saturation_fraction:float=0.3,
                 ignore_saturation_error:bool=False,
                 qa_show:bool=None,
                 qa_showscale:float | int=None,
                 qa_showblock:bool=None,
                 qa_write_findlines:bool=False,
                 qa_write:bool=None):


    """
    To create a spextool wavecal file.

    Parameters
    ----------
    arc_files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc-00001.a.fits, spc-00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
        
    flat_file : str
        The full path of a pySpextool flat field file.

    output_filename : str
        The filename of the wavecal file to written to disk.      

    sky_files : str or list or None
        If type NoneType, then only arc files are required.
    
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc-00001.a.fits, spc-00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].

    linearity_correction : {True, False}
        Set to True to correct for non-linearity.
        Set to False to not correct for non-linearity.

    detector_info : dict, deafult None
        A dictionary with any information that needs to be passed to the
        instrument-specific readfits program.  
    
    use_stored_solution : {False, True}
        Set to True to use  the solution stored in the pyspextool cal file.  

    saturation_fraction : float, default 0.3
        The fraction of pixels identified as saturated within in order beyond 
        which a pySpextoolError will be thrown.  The wavelength calibration 
        is most likely bad.  
   
    ignore_saturation_error : {False, True}
        Set to True to ignore issues with saturation of lines.
        Set to False to throw a pySpextoolError error.  

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
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].


    Returns
    -------
     None
        Writes a pySpextool wavecal FITS file to disk.
        
    """
    
    #
    # Check parameters and keywords
    #

    check_parameter('make_wavecal', 'arc_files', arc_files, ['str', 'list'],
                    list_types=['str','str'])

    check_parameter('make_wavecal', 'flat_file', flat_file, 'str')

    check_parameter('make_wavecal', 'output_filename', output_filename, 'str')

    check_parameter('make_wavecal', 'sky_files', sky_files,
                    ['str', 'list','NoneType'], list_types=['str','str'])

    check_parameter('make_wavecal', 'linearity_correction',
                    linearity_correction, 'bool')

    check_parameter('make_wavecal', 'detector_info', detector_info,
                    ['dict', 'NoneType'])

    check_parameter('make_wavecal', 'use_stored_solution',
                    use_stored_solution, 'bool')

    check_parameter('get_wavecal', 'verbose', verbose, ['NoneType','bool'])

    check_parameter('get_wavecal', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('get_wavecal', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])
    
    check_parameter('get_wavecal', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
        
    check_parameter('get_wavecal', 'qa_write', qa_write, ['NoneType','bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Let the user know what you are doing.
    #

    message = ' Generating Wavelength Solution'
    logging.info(message+'\n'+'-'*(len(message)+6)+'\n')

    extract.load['wavecal']['output_filename'] = output_filename
    extract.load['wavecal']['use_stored_solution'] = use_stored_solution

    #
    # Create the image
    #

    # Create the file names

    flat_fullpath = os.path.join(setup.state['cal_path'], flat_file)
    
    result = files_to_fullpath(setup.state['raw_path'],
                               arc_files,
                               setup.state['nint'],
                               setup.state['suffix'],
                               setup.state['search_extension'])

    arcs_fullpath = result[0]

    if sky_files is not None:

        result = files_to_fullpath(setup.state['raw_path'],
                                   sky_files,
                                   setup.state['nint'],
                                   setup.state['suffix'],
                                   setup.state['search_extension'])

        skys_fullpath = result[0]

        skys_filenames = ','.join(result[2])
                      
    else:

        skys_fullpath = None
        skys_filenames = None

    #
    # Create the image
    #

    result = _make_wavecal_image(arcs_fullpath,
                                 skys_fullpath,
                                 linearity_correction,
                                 detector_info,
                                 flat_fullpath,
                                 qa['verbose'],
                                 setup.plots['wavecal_image'],
                                 setup.plots['square_size'],
                                 setup.plots['font_size'],
                                 qa['show'],
                                 qa['showscale'],
                                 qa['showblock'],    
                                 qa['write'])

    flatinfo = result[0]
    wavecalinfo = result[1]
    wavecal_image = result[2]
    wavecal_mask = result[3]
    
    #
    # Check for saturation of lines
    #

    fractions = []
    for i in range(flatinfo['norders']):

        zorder = np.where(flatinfo['ordermask'] == flatinfo['orders'][i])[0]
        zsat   = np.where((wavecal_mask == 1) & 
                          (flatinfo['ordermask'] == flatinfo['orders'][i]))[0]
       
        if len(zsat)/len(zorder) > saturation_fraction:

            fractions.append(True)

        else:

            fractions.append(False)

#    fractions[1] = True
    saturated_orders = flatinfo['orders'][fractions]

    if len(saturated_orders) > 0:

        filename1 = extract.load['wavecal']['output_filename'] + \
            '_image'+setup.state['qa_extension']

        filename2 = extract.load['wavecal']['output_filename'] + \
            '_findlines.pdf'


        message = '>'+str(saturation_fraction*100)+'% of pixels in Order(s) '+\
            ", ".join([str(x) for x in saturated_orders])+' are saturated.  '+\
            'The wavelength calibration will most likely be useless so '+\
            'please check the QA plot '+filename1+'.  Also '+\
            "consider setting `qa_write_findlines` to True and inspecting "+\
            'the QA plot '+filename2+' to ensure all lines are identified '+\
            'properly.  Please choose unsaturated files for the wavelength '+\
            'calibration To proceed without doing so, '+\
            "set the `ignore_saturation_error` keyword to True."

        if ignore_saturation_error is False:
            
            raise pySpextoolError(message)

    #
    # Let's do the extraction of the wavecal image
    #

    # Create wavecal and spatcal images

    result = wavecal.simulate_wavecal_1dxd(flatinfo['ncols'],
                                           flatinfo['nrows'],
                                           flatinfo['edgecoeffs'],
                                           flatinfo['xranges'],
                                           flatinfo['slith_arc'])

    wavecal_pixels = result[0]
    spatcal = result[1]

    # Extract the spectra     

    
    tracecoeffs = np.full((flatinfo['norders'], 2), 0)
    tracecoeffs[:, 0] = flatinfo['slith_arc'] / 2
    aperture_radii = np.full((flatinfo['norders'], 1), wavecalinfo['apradius'])

    
    message = ' Sum extracting 1 aperture in '+str(flatinfo['norders'])+\
        ' orders (without background subtraction).'
    logging.info(message)

    spectra = extract_1dxd(wavecal_image,
                           wavecal_mask,
                           flatinfo['ordermask'],
                           wavecal_pixels,
                           spatcal,
                           flatinfo['orders'],
                           tracecoeffs,
                           aperture_radii,
                           np.full((1),1),
                           progressbar=qa['verbose'])

    #
    # Find the pixel offset between these spectra and the disk spectra
    #

    # Get the anchor order and spectra.

    z = np.sum(np.where(flatinfo['orders'] == wavecalinfo['xcororder']))

    xanchor = np.arange(int(wavecalinfo['xranges'][z, 0]),
                        int(wavecalinfo['xranges'][z, 1] + 1), dtype=int)
    fanchor = np.squeeze(wavecalinfo['spectra'][z, 1, :])

    # Get the source order

    xsource = np.squeeze(spectra[z][0, :])
    fsource = np.squeeze(spectra[z][1, :])

    if qa['write'] is True:

            filename = extract.load['wavecal']['output_filename'] + \
                '_pixelshift'+setup.state['qa_extension']
            fullpath = os.path.join(setup.state['qa_path'],filename)

    else:

        fullpath = None

    offset = wavecal.get_spectral_pixelshift(xanchor,
                                     fanchor,
                                     xsource,
                                     fsource,
                                     qa_plotnumber=setup.plots['pixel_shift'],
                                     qa_figuresize=setup.plots['portrait_size'],
                                     qa_fontsize=setup.plots['font_size'],
                                     qa_show=qa['show'],
                                     qa_showscale=qa['showscale'],
                                     qa_showblock=qa['showblock'],
                                     qa_fullpath=fullpath,
                                     anchor_label='Stored Spectrum', 
                                     source_label='Observed Spectrum')

    #
    # Are we using the stored solution?
    #

    # Let's check to see what the wavecal file says as a function of the slit
    # width
    
    z = wavecalinfo['slits'] == flatinfo['slitw_arc'] 

    use_stored_solution_wavecal = bool(wavecalinfo['usestored'][z])


    # Now compare to the user request
       
    use_stored_solution = use_stored_solution_wavecal or use_stored_solution

    if use_stored_solution is False:

        #
        # Locate the line positions
        #

        # Get the line list to search for lines

        filename = os.path.join(setup.state['instrument_path'], 
                                wavecalinfo['linelist'])

        lineinfo = wavecal.read_line_list(filename,
                                          delta_to_microns=True)
        
        #  Determine the guess position and search range for each
    
        lineinfo = wavecal.get_line_guess_position(wavecalinfo['spectra'],
                                                   wavecalinfo['orders'],
                                                   wavecalinfo['xranges'],
                                                   lineinfo)
        
        # Add the shift offset to the results

        lineinfo['xguess'] = lineinfo['xguess'] + offset
        lineinfo['range_min_xguess'] = lineinfo['range_min_xguess'] + offset
        lineinfo['range_max_xguess'] = lineinfo['range_max_xguess'] + offset

        # Now find the lines

        logging.info(' Finding the lines.')

        if qa_write_findlines is True:

            filename = extract.load['wavecal']['output_filename'] + \
                '_findlines.pdf'
            fullpath = os.path.join(setup.state['qa_path'],filename)
            
        else:

            fullpath = None

        # Find the lines
    
        lineinfo = wavecal.find_lines_1dxd(spectra,
                                   wavecalinfo['orders'],
                                   lineinfo,
                                   flatinfo['slitw_pix'],
                                   qa_figuresize=setup.plots['portrait_size'],
                                   qa_fontsize=setup.plots['font_size'],
                                   qa_fullpath=fullpath,
                                   verbose=qa['verbose'])

        #
        # Let's do the actual calibration
        #

        logging.info(' Determining the wavelength solution.')

        # Get set up for either 1d of 1dxd

        if wavecalinfo['wcaltype'] == '1d':

            figure_size = setup.plots['landscape_size']
            xdinfo = None

        else:

            figure_size = setup.plots['portrait_size']
            xdinfo = {'homeorder': wavecalinfo['homeorder'],
                      'orderdeg': wavecalinfo['ordrdeg']}

        
        if qa['write'] is True:

            filename = extract.load['wavecal']['output_filename'] + \
                '_residuals'+setup.state['qa_extension']
            fullpath = os.path.join(setup.state['qa_path'],filename)
            
        else:
            
            fullpath = None

        # Find the solution

        solution = wavecal.wavecal_solution_1d(wavecalinfo['orders'],
                                       lineinfo,
                                       wavecalinfo['dispdeg'],
                                       xd_info=xdinfo,
                                       verbose=qa['verbose'],
                            qa_plotnumber=setup.plots['wavecal_residuals'],
                                       qa_figuresize=figure_size,
                                       qa_fontsize=setup.plots['font_size'],
                                       qa_show=qa['show'],
                                       qa_showblock=qa['showblock'],
                                       qa_fullpath=fullpath)
        stored_solution_offset = 0.0

    else:

        #
        # Going to use the stored solution.  
        #
        
        logging.info(' Using stored wavelength solution.')

        solution = {'coeffs': wavecalinfo['coeffs'],
                    'covar': wavecalinfo['covar'],
                    'rms': wavecalinfo['rms'],
                    'nlines': wavecalinfo['nlines'],
                    'ngood': wavecalinfo['ngood'],
                    'nbad': wavecalinfo['nbad']}

        stored_solution_offset = -offset
            
    #
    # Creating rectification indices
    #

    indices = []
    for i in range(flatinfo['norders']):
        idxs = wavecal.make_interp_indices_1d(flatinfo['edgecoeffs'][i, :, :],
                                              flatinfo['xranges'][i, :],
                                              flatinfo['slith_arc'],
                                              array_output=True)
              
        indices.append(idxs)
        
    #
    # Write the wavecal file to disk.
    #

    logging.info(' Writing wavecal to disk.')

    if wavecalinfo['wcaltype'] == '1d':

        xdinfo = None

    else:

        xdinfo = {'orderdeg': wavecalinfo['ordrdeg'],
                  'homeorder': wavecalinfo['homeorder']}

    oname = os.path.join(setup.state['cal_path'], output_filename + '.fits')
        
    wavecal.write_wavecal1d_fits(flatinfo['ordermask'],                         
                                 flatinfo['xranges'],
                                 solution['coeffs'],
                                 solution['covar'],
                                 wavecalinfo['dispdeg'],
                                 solution['rms'] * 1e4,
                                 solution['nlines'],
                                 solution['ngood'],
                                 solution['nbad'],
                                 wavecal_pixels,
                                 stored_solution_offset,
                                 spatcal,
                                 indices,
                                 flatinfo['rotation'],
                                 flat_file,
                                 skys_filenames,
                                 oname,
                                 setup.state['version'],
                                 xdinfo=xdinfo,
                                 stored_solution=use_stored_solution)

    logging.info(' Wavecal file ' + output_filename + \
                 '.fits written to disk.\n')


def _make_wavecal_image(arc_files:str | list,
                        sky_files:str | list | None,
                        linearity_correction:bool,
                        detector_info:dict | None,
                        flat_file:str,
                        verbose:bool,
                        qa_plotnumber:int,
                        qa_figuresize:tuple,
                        qa_fontsize:int,
                        qa_show:bool,
                        qa_showscale:bool,
                        qa_showblock:bool,
                        qa_write:bool):
    
    """
    Creates an "arc" image using arc and sky images for wavelength calibration.

    Parameters
    ----------
    arc_files : str or list
        If type is str, then the fullpath to an arc file.
        If type is list, then an (nfiles,) list of fullpaths to arc files.

    sky_files : str or list
        If type is str, then the fullpath to a sky file.
        If type is list, then an (nfiles,) list of fullpaths to sky files.
       
    linearity_correction : {True, False}
        Set to True to correct for non-linearity.
        Set to False to not correct for non-linearity.

    detector_info : dict, deafult None
        A dictionary with any information that needs to be passed to the
        instrument-specific readfits program.  

    flat_file : str
        A string to a pySpextool flat file.
    
    qa_plotnumber : int
        The plot number for the QA plot.
    
    qa_figuresize : tuple
        A (2,) array giving the figure size.

    qa_fontsize : int
        The font size.

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
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].

    Returns
    -------
    ndarray

        An (nrows, ncols) array suitable for use in wavelength calibration.
        If `sky_files` is not None, then 

    """

    #
    # Check parameters and keywords
    #

    check_parameter('_make_wavecal_image', 'arc_files', arc_files, 
                    ['str', 'list'])

    check_parameter('_make_wavecal_image', 'sky_files', sky_files, 
                    ['str', 'list', 'NoneType'])

    check_parameter('_make_wavecal_image', 'linearity_correction',
                    linearity_correction, 'bool')

    check_parameter('_make_wavecal_image', 'detector_info', detector_info,
                    ['dict', 'NoneType'])

    check_parameter('_make_wavecal_image', 'flat_file', flat_file, 'str')

    check_parameter('_make_wavecal_image', 'qa_plotnumber', qa_plotnumber, 
                    'int')

    check_parameter('_make_wavecal_image', 'qa_figuresize', qa_figuresize, 
                    'tuple')

    check_parameter('_make_wavecal_image', 'qa_fontsize', qa_fontsize, 
                    'int')

    check_parameter('_make_wavecal_image', 'verbose', verbose, 
                    ['NoneType','bool'])

    check_parameter('_make_wavecal_image', 'qa_show', qa_show, 
                    ['NoneType','bool'])

    check_parameter('_make_wavecal_image', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])
    
    check_parameter('_make_wavecal_image', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
        
    check_parameter('_make_wavecal_image', 'qa_write', qa_write, 
                    ['NoneType','bool'])

    check_qakeywords(verbose=verbose)
    
    #
    # Load the instrument module for the read_fits program
    #

    module = 'pyspextool.instruments.' + setup.state['instrument'] + \
             '.' + setup.state['instrument']

    instr = importlib.import_module(module)

    #
    # Read the flat and wavecal files
    #

    flatinfo = read_flat_fits(flat_file)

    wavecalfile = os.path.join(setup.state['instrument_path'],
                               flatinfo['mode'] + '_wavecalinfo.fits')

    wavecalinfo = wavecal.read_wavecal_file(wavecalfile)

    #
    # Deal with the arc images first
    #

    # Load the data

    logging.info(' Creating the arc image.')    

    result = instr.read_fits(arc_files,
                             setup.state['linearity_info'],
                             pair_subtract=False,
                             keywords=setup.state['extract_keywords'],
                             rotate=flatinfo['rotation'],
                             linearity_correction=linearity_correction,
                             extra=detector_info,
                             verbose=False)
    
    imgs = result[0]
    masks = result[3]

    #
    # Combine images as necessary and don't worry about the variance
    #

    if np.ndim(imgs) == 3:

        # Scale their intensities to a common flux level

        logging.info(' Scaling the arc images.')

        simgs, svars, scales = scale_data_stack(imgs, None)

        # Now median the scaled images
        
        logging.info(' Medianing the scaled arc images.')

        arc, munc = median_data_stack(simgs)

        # Combine the masks
        
        arc_mask = combine_flag_stack(masks, nbits=1)

    else:
        
        arc = imgs
        arc_mask = masks

    #
    # Now deal with the potential sky frames
    #
    
    if len(wavecalinfo['skyorders']) !=0:

        logging.info(' Creating the sky image.')
        
        # Load the image(s)

        result = instr.read_fits(sky_files,
                                 setup.state['linearity_info'],
                                 pair_subtract=False,
                                 keywords=setup.state['extract_keywords'],
                                 rotate=flatinfo['rotation'],
                                 verbose=verbose)
        imgs = result[0]
        masks = result[3]
        
        # Combine them if necessary

        if np.ndim(imgs) == 3:
            
            logging.info(' Scaling the sky images.')
            
            simgs, svars, scales = scale_data_stack(imgs, None)
            
            # Now median the scaled images
            
            logging.info(' Medianing the scaled sky images.')
            
            sky, munc = median_data_stack(simgs)
 
            # Combine the masks
            
            sky_mask = combine_flag_stack(masks, nbits=1)
            
        else:
            
            sky = imgs
            sky_mask = masks
            
        # Mix the orders

        result = wavecal.mix_orders(arc,
                                    arc_mask,
                                    sky,
                                    sky_mask,
                                    flatinfo['ordermask'],
                                    flatinfo['orders'],
                                    wavecalinfo['arcorders'],
                                    wavecalinfo['skyorders'])
        
        wavecal_image = result[0]
        wavecal_mask = result[1]

    else:
        
        wavecal_image = arc
        wavecal_mask = arc_mask

    #
    # Do the QA plot
    #

    orders_plotinfo = {'edgecoeffs': flatinfo['edgecoeffs'],
                       'xranges': flatinfo['xranges'],
                       'orders': flatinfo['orders']}
        
    if qa_write is True:

        filename = extract.load['wavecal']['output_filename'] + \
            '_image'+setup.state['qa_extension']
        fullpath = os.path.join(setup.state['qa_path'],filename)
        
        plot_image(wavecal_image,
                   mask=wavecal_mask,
                   orders_plotinfo=orders_plotinfo,
                   figure_size=qa_figuresize,
                   font_size=qa_fontsize,
                   output_fullpath=fullpath)
        
    if qa_show is True:
        
        plot_image(wavecal_image,
                   mask=wavecal_mask,
                   plot_number=qa_plotnumber,
                   orders_plotinfo=orders_plotinfo,
                   figure_size=qa_figuresize,
                   font_size=qa_fontsize,
                   showscale=qa_showscale,
                   showblock=qa_showblock)

    #
    # Return the arrays
    #
        
    return flatinfo, wavecalinfo, wavecal_image, wavecal_mask

