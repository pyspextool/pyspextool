from os.path import join, basename
import importlib
import numpy as np
from astropy.io import fits
import logging

from pyspextool import config as setup
from pyspextool.extract import config as extract

from pyspextool.extract import background_subtraction as background
from pyspextool.extract.images import scale_order_background
from pyspextool.io.check import check_parameter, check_qakeywords, check_file
from pyspextool.io.files import make_full_path
from pyspextool.extract.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.io.fitsheader import average_headerinfo
from pyspextool.utils import math
from pyspextool.utils.arrays import idl_unrotate
from pyspextool.utils.split_text import split_text
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.plot.plot_image import plot_image
from pyspextool.io.files import files_to_fullpath
from pyspextool.pyspextoolerror import pySpextoolError


def combine_images(files:str | list,
                   output_filename:str,
                   input_extension:str='.fits*',
                   output_directory:str='proc',
                   correct_nonlinearity:bool=True,
                   beam_mode:str='A',
                   scale_background:bool=False,
                   subtract_background:bool=False,
                   flatfield_filename:str=None,
                   flatfield:bool=False,
                   statistic:str='robust weighted mean',
                   robust_sigma:int | float=8,
                   verbose=None,
                   qa_show=None,
                   qa_showscale:float | int=1.0,
                   qa_showblock:bool=None,
                   qa_write=None):

    """
    To combine (pair-subtracted) images for later extraction.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
        
    output_name : str
        The filename of the combined image to written to disk.      

    input_extension : str, default '.fits*'
        The file extension

    correct_nonlinearity : {True, False}
        Set to True to correct for non-linearity.
        Set to False to not correct for non-linearity.
        Deafult is True

    beam_mode : {'A', 'A-B'}
        Seto to 'A-B' to pair-subtract images before combining.  

    scale_orders : {False, True}
        Set to scale the flux levels in each order to the median value of all
        orders.

    subtract_background : {False, True}
        Set to subtract the background off of each image.  Currently only 
        does the median flux level over the slit.

    flatfield_name : str, optional
        A string given the flat field FITS file associated with the images.  
        Required if `scale_orders`, `background_subtraction` or `flat_field` is 
        set to True.

    flatfield : {False, True}
        Set to flat field the image.

    statistic : {'robust weighted mean', 'robust mean', 'weighted mean', 
                 'mean', 'median'}
        For any of the means, the uncertainty is the standard error on the 
        mean, e.g. s/np.sqrt(n) where s is the sample standard deviation and 
        n is the number of data points.

        For the median, the uncertainty is given by 1.482*MAD/np.sqrt(n) where
        MAD is the median absolute deviation and is given by,

        1.482*median(|x_i-x_med|).

    robust_sigma : float or int, default=8
        The sigma threshold of a robust statistic is requested.   Values are 
        identified as outliers if 

        |x_i - x_med|/MAD > robust_sigma,

        where x_i is the ith data point, x_med is the median of the data, and 
        MAD is the median absolute deviation given by,

        1.482*median(|x_i-x_med|).

    output_directory : {'proc', 'cal'}
        The pyspextool directory into which the file is written.

    qa_show : {None, True, False}
        Set to True/False to override config.state['qa_show'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qa_showsize : tuple, default=(10, 6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}
        Set to True/False to override config.state['qa_write'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.
    
    Returns
    -------
    None
    Writes a FITS file to disk.  


    """

    #
    #  Check parameters and QA keywords
    #

    check_parameter('combine_images', 'files', files, ['str', 'list'],
                    list_types=['str','str'])

    check_parameter('combine_images', 'output_filename', output_filename, 'str')

    check_parameter('combine_images', 'input_extension', input_extension, 'str')

    check_parameter('combine_images', 'output_directory', output_directory,
                    'str', possible_values=['proc', 'cal'])
    
    check_parameter('combine_images', 'correct_nonlinearity',
                    correct_nonlinearity, 'bool')

    check_parameter('combine_images', 'beam_mode', beam_mode, 'str',
                    possible_values=['A', 'A-B'])

    check_parameter('combine_images', 'scale_background', scale_background,
                    'bool')

    check_parameter('combine_images', 'subtract_background',
                    subtract_background, 'bool')

    check_parameter('combine_images', 'flatfield_fiename', flatfield_filename,
                    ['NoneType', 'str'])

    check_parameter('combine_images', 'flatfield', flatfield, 'bool')

    check_parameter('combine_spectra', 'statistic', statistic, 'str')

    check_parameter('combine_spectra', 'robust_sigma', robust_sigma,
                    ['int', 'float'])

    check_parameter('combine_images', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('combine_images', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('combine_images', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('combine_images', 'qa_showscale', qa_showscale,
                    ['int', 'float'])

    check_parameter('combine_images', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Let the user know what you are doing.
    #

    logging.info(' Combining Images \n ---------------------')    
    
#
    # Create the file names

    input_files = files_to_fullpath(setup.state['raw_path'],
                                    files,
                                    setup.state['nint'],
                                    setup.state['suffix'],
                                    input_extension)
    check_file(input_files)

    #
    # Deal with the beam parameter
    #

    if beam_mode == 'A-B':

        # Are there an even number of images?

        if len(input_files) % 2 == 1:

            message = "An even number of images are required for "+ \
            "`beam_mode` to be equal to 'A-B'."
            raise pySpextoolError(message)

        if setup.state['irtf'] is True:

            input_files, indices = reorder_irtf_files(input_files)

        nimages = int(len(input_files) / 2)            
        pair = True

    else:

        nimages = int(len(input_files))
        pair = False

    #
    #  Load the images
    #

    # First we have to figure out whether we are doing background subtraction


    if (scale_background or subtract_background or flatfield) is True:

        # Check if the flat field file is passed.

        if flatfield_filename is None:

            message = 'The parameter `flatfield_filename` is required for '+\
                'flat fielding and/or scaling orders and/or background '+\
                'subtraction.'
            raise pySpextoolError(message)

        # Read the flat field

        file_name = join(setup.state['cal_path'], flatfield_filename)
        check_file(file_name)

        flatinfo = read_flat_fits(file_name)
        rotation = flatinfo['rotation']
        edgecoeffs = flatinfo['edgecoeffs']
        xranges = flatinfo['xranges']
        orders = flatinfo['orders']
        ybuffer = flatinfo['ybuffer']

        unrotate = True

    else:

        rotation = 0
        unrotate = False

    # Check beam mode


            
    
    # Load the instrument module for the read_fits program

    module = 'pyspextool.instrument_data.' + setup.state['instrument'] + \
             '_dir.' + setup.state['instrument']

    instr = importlib.import_module(module)

    logging.info(' Loading images...')    
    
    result = instr.read_fits(input_files,
                             setup.state['linearity_info'],
                             pair_subtract=pair,
                             linearity_correction=correct_nonlinearity,
                             keywords=setup.state['extract_keywords'],
                             rotate=rotation,
                             verbose=qa['verbose'])

    data = result[0]
    var = result[1]
    hdr = result[2]
    mask = result[3]
    
    #
    # Scale the orders to a common intensity level?
    #
    
    if scale_background is True:

            logging.info(' Scaling the orders to a common intensity level.')

            result = scale_order_background(data,
                                            orders,
                                            edgecoeffs,
                                            xranges,
                                            var_stack=var,
                                            ybuffer=ybuffer,
                                            verbose=qa['verbose'])

            data = result[0]
            var = result[1]

    #
    # Subtract the background? 
    #

    if subtract_background is True:

        logging.info(' Subtracting the background.')    
        
        for i in range(nimages):
            
            result = background.median_1dxd(data[i, :, :],
                                            flatinfo['edgecoeffs'],
                                            flatinfo['xranges'],
                                            var=var[i, :, :],
                                            ybuffer=ybuffer)
            data[i, :, :] = result[0]
            var[i, :, :] = result[1]

    #
    # Now combine the images
    #

    logging.info(' Combining images with a ' + statistic + '.')

    if statistic == 'robust weighted mean':
        
        result = math.mean_data_stack(data, robust=robust_sigma,
                                      weights=1 / var, stderr=True)
    elif statistic == 'robust mean':
        
        result = math.mean_data_stack(data, robust=robust_sigma, stderr=True)
        
    elif statistic == 'weighted mean':
        
        result = math.mean_data_stack(data, weights=1 / var, stderr=True)

    elif statistic == 'mean':

        result = math.mean_data_stack(data, stderr=True)

    elif statistic == 'median':

        result = math.median_data_stack(data, stderr=True)

    else:

        message = '`statistic` is unknown.'
        raise pySpextool(message)

    mean = result[0]
    var = result[1] ** 2

    # Now do the mask

    mask = math.combine_flag_stack(mask)

    #
    # Flat field if requested
    #

    if flatfield is True:

        logging.info(' Flat fielding the image.')

        np.divide(mean, extract.state['flat'], out=mean)
        np.divide(var, extract.state['flat'] ** 2, out=var)

    #
    # Unrotate the images if necessary
    #

    if unrotate is True:
        
        mean = idl_unrotate(mean, rotation)
        var = idl_unrotate(var, rotation)
        mask = idl_unrotate(mask, rotation)

    #
    # Write the results to disk
    #

    # Average the headers together

    pair = True if beam_mode == 'A-B' else False

    avehdr = average_headerinfo(hdr, pair=pair)

    # Store the history

    is_history_present = 'HISTORY' in avehdr

    if is_history_present is True:

        old_history = avehdr['HISTORY']

        # remove it from the header_info dictionary

        avehdr.pop('HISTORY')
    
    # Add useful things to header

    avehdr['MODULE'] = ['extract', ' Creation module']

    avehdr['VERSION'] = [setup.state['version'], ' pySpextool version']

    avehdr['IC_STAT'] = [statistic, ' Image combine: combination statistic']

    avehdr['IC_THRSH'] = [robust_sigma,
                          ' Image combine:  robust threshold (if used)']

    avehdr['IC_LNCOR'] = [correct_nonlinearity,
                         ' Image combine: linearity corrected?']

    avehdr['IC_LNMAX'] = [setup.state['linearity_info']['max'],
                          ' Image combine: linearity correction max (DN)']

    avehdr['IC_SCLED'] = [scale_background, ' Image combine: orders scaled?']

    avehdr['IC_BGSUB'] = [subtract_background,
                          ' Imagine combine: background subtraction?']

    avehdr['IC_FLATD'] = [flatfield, ' Imagine combine: flat fielded?']

    avehdr['IC_NIMGS'] = [len(files),
                          ' Image combine: number of images combined']

    avehdr['IC_BEAM'] = [beam_mode, ' Imagine combine: beam mode']

    avehdr['FILENAME'] = [output_filename + '.fits', ' Filename']

    # Create the history

    # Strip the directories.

    file_names = [basename(i) for i in files]

    history = 'This image was created by combining the images ' + \
              ', '.join(file_names) + ' with a ' + statistic

    if statistic[0:6] == 'robust':

        history += ' with a sigma threshold of ' + str(robust_sigma) + '.  '

    else:

        history = history + '.'

    if correct_nonlinearity is True:
        history = history + '  The raw data were corrected for non-linearity.'

    if scale_background is True:
        history = history + '  The orders were scaled to a common flux level.'

    if subtract_background is True:
        history = history + '  Background subtraction was accomplished by ' + \
                  'subtracting the median flux level on a ' + \
                  'column by column basis within each order.'

    if flatfield is True:
        history = history + '  The image was flat fielded.'

    if is_history_present is True:
        
        history = old_history + history


    #        
    # Make the QA plot
    #

    if (scale_background or subtract_background or flatfield) is True:
    
        orders_plotinfo = {'xranges': flatinfo['xranges'],
                           'edgecoeffs': flatinfo['edgecoeffs'],
                           'orders': flatinfo['orders']}

    else:

        orders_plotinfo = None
        
    if qa['show'] is True:
        
        plot_image(mean,
                   mask=mask,
                   orders_plotinfo=orders_plotinfo,
                   figure_size=(setup.plots['square_size'][0]*qa['showscale'],
                                setup.plots['square_size'][1]*qa['showscale']),
                   font_size=setup.plots['font_size']*qa['showscale'],
                   showblock=qa['showblock'],
                   plot_number=setup.plots['flat'])
            
    if qa['write'] is True:

        filename = output_filename + '_combined' + setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

        plot_image(mean,
                   mask=mask,
                   orders_plotinfo=orders_plotinfo,
                   output_fullpath=fullpath,
                   figure_size=setup.plots['square_size'],
                   font_size=setup.plots['font_size'])


        

    #
    # Write the file to disk
    #

    # Create the primary HDU

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the original file header keywords

    keys = list(avehdr.keys())

    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (avehdr[keys[i]][0], avehdr[keys[i]][1])

    # Add the history

    history = split_text(history, length=65)

    for hist in history:
        hdr['HISTORY'] = hist

        # Write the results

    img_hdu = fits.ImageHDU(mean)
    var_hdu = fits.ImageHDU(var)
    flg_hdu = fits.ImageHDU(mask)

    hdu = fits.HDUList([phdu, img_hdu, var_hdu, flg_hdu])
    hdu.writeto(join(setup.state['proc_path'], output_filename + '.fits'),
                overwrite=True)

    logging.info(' Wrote ' + output_filename + '.fits to disk.')
