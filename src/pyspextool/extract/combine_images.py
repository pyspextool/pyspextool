import os
import importlib
import numpy as np
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.extract import config as extract

from pyspextool.extract import background_subtraction as background
from pyspextool.extract.scale_orders import scale_orders as scale_background
from pyspextool.io.check import check_parameter
from pyspextool.io.files import make_full_path
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.io.fitsheader import average_header_info
from pyspextool.utils import math
from pyspextool.utils.arrays import idl_unrotate
from pyspextool.utils.split_text import split_text
from pyspextool.utils.loop_progress import loop_progress


def combine_images(files, output_name, extension='.fits*',
                   linearity_correction=True, beam_mode='A',
                   scale_orders=False, background_subtraction=False,
                   flat_field_name=None, flat_field=False,
                   statistic='robust weighted mean',
                   robust_sigma=8, output_directory='proc', overwrite=True,
                   qa_plot=None, qa_plotsize=(6, 6), qa_file=None,
                   verbose=None):
    """
    To combine (pair-subtracted) images for later extraction

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

    extension : str, default='.fits*'
        The file extension

    linearity_correction : {True, False}
        Set to True to correct for linearity

    beam_mode : {'A', 'A-B'}
        Seto to 'A-B' to pair-subtract images before combining.  

    scale_orders : {False, True}
        Set to scale the flux levels in each order to the median value of all
        orders.

    background_subtraction : {False, True}
        Set to subtract the background off of each image.  Currently only 
        does the median flux level over the slit.

    flat_field_name : str, optional
        A string given the flat field FITS file associated with the images.  
        Required if `scale_orders`, `background_subtraction` or `flat_field` is 
        set to True.

    flat_field : {False, True}
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

    overwrite : {True, False}
        Set to True to overwrite a file on disk with the same name.

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

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.
    
    Returns
    -------
    None
    Writes a FITS file to disk.  


    """

    #
    #  Check parameters
    #

    check_parameter('combine_images', 'files', files, ['str', 'list'])

    if isinstance(files, list):
        check_parameter('combine_images', 'files[0]', files[0], 'str')
        check_parameter('combine_images', 'files[1]', files[1],
                        ['str', 'list', 'int'])

    check_parameter('combine_images', 'output_name', output_name, 'str')

    check_parameter('combine_images', 'output_directory', output_directory,
                    'str', possible_values=['proc', 'cal'])

    check_parameter('combine_images', 'extension', extension, 'str')

    check_parameter('combine_images', 'beam_mode', beam_mode, 'str',
                    possible_values=['A', 'A-B'])

    check_parameter('combine_images', 'linearity_correction',
                    linearity_correction, 'bool')

    check_parameter('combine_images', 'scale_orders', scale_orders, 'bool')

    check_parameter('combine_images', 'background_subtraction',
                    background_subtraction, 'bool')

    check_parameter('combine_images', 'flat_field_name', flat_field_name,
                    ['NoneType', 'str'])

    check_parameter('combine_images', 'flat_field', flat_field, 'bool')

    check_parameter('combine_spectra', 'statistic', statistic, 'str')

    check_parameter('combine_spectra', 'robust_sigma', robust_sigma,
                    ['int', 'float'])

    check_parameter('combine_images', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('combine_images', 'overwrite', overwrite, 'bool')

    check_parameter('combine_images', 'qa_file', qa_file, ['NoneType', 'bool'])

    check_parameter('combine_images', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('combine_images', 'qa_plotsize', qa_plotsize,
                    ['NoneType', 'tuple'])

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
    # Store user inputs
    #

    extract.combine['linearity_correction'] = linearity_correction
    extract.combine['scale_orders'] = scale_orders
    extract.combine['background_subtraction'] = background_subtraction
    extract.combine['flat_field_name'] = flat_field_name
    extract.combine['flat_field'] = flat_field
    extract.combine['overwrite'] = overwrite
    extract.combine['qaplotsize'] = qa_plotsize
    extract.combine['qaplot'] = qa_plot
    extract.combine['qafile'] = qa_file
    extract.combine['verbose'] = verbose

    #
    # Let the user know what you are doing.
    #

    if verbose is True:
        print('Combining Images')
        print('---------------------')

        #
    # Load the instrument module for the read_fits program
    #

    module = 'pyspextool.instrument_data.' + setup.state['instrument'] + \
             '_dir.' + setup.state['instrument']

    instr = importlib.import_module(module)

    #
    # Create the file names
    #

    # Create the file names

    if isinstance(files, str):

        # You are in FILENAME mode

        files = files.replace(" ", "").split(',')
        files = make_full_path(setup.state['raw_path'], files, exist=True)

    else:

        # You are in INDEX mode

        prefix = files[0]
        nums = files[1]

        files = make_full_path(setup.state['raw_path'], nums,
                               indexinfo={'nint': setup.state['nint'],
                                          'prefix': prefix,
                                          'suffix': setup.state['suffix'],
                                          'extension': extension},
                               exist=True)

    #
    # Deal with the beam parameter
    #

    if beam_mode == 'A-B' and setup.state['irtf'] is True:
        files, indices = reorder_irtf_files(files)

    #
    #  Load the images
    #

    # First we have to figure out whether we are doing background subtraction

    if scale_orders is True or background_subtraction is True:

        # Check if the flat field file is passed.

        if flat_field_name is None:
            message = '`flat_field_name` required for scaling/subtraction.'
            raise ValueError(message)

        # Read the flat field

        file_name = os.path.join(setup.state['cal_path'], flat_field_name)
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

    pair = True if beam_mode == 'A-B' else False

    data, var, hdr, mask = instr.read_fits(files,
                                           setup.state['linearity_info'],
                                           pair_subtract=pair,
                                           linearity_correction=linearity_correction,
                                           keywords=setup.state['extract_keywords'],
                                           rotate=rotation,
                                           verbose=verbose)

    if pair is True:

        nimages = int(len(files) / 2)

    else:

        nimages = int(len(files))

    #
    # Scale orders
    #

    if scale_orders is True:

        if verbose is True:
            print('Scaling the orders to a common level...')

        result = scale_background(data, orders, edgecoeffs, xranges,
                                  var_stack=var, ybuffer=ybuffer)

        data = result[0]
        var = result[1]

    #
    # Background subtraction    
    #

    if background_subtraction is True:

        for i in range(nimages):

            if verbose is True:
                loop_progress(i, 0, nimages,
                              message='Subtracting background...')

            result = background.median_1dxd(data[i, :, :], flatinfo['edgecoeffs'],
                                            flatinfo['xranges'], var=var[i, :, :],
                                            ybuffer=flatinfo['ybuffer'])
            data[i, :, :] = result[0]
            var[i, :, :] = result[1]

    #
    # Now combine the images together
    #

    if verbose is True:
        print('Combining images with a ' + statistic + '...')

    if statistic == 'robust weighted mean':

        print(np.shape(data))
        result = math.mean_data_stack(data, robust=robust_sigma, weights=1 / var,
                                      stderr=True)
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
        raise ValueError(message)

    mean = result[0]
    var = result[1] ** 2

    # Now do the mask

    mask = math.combine_flag_stack(mask)

    #
    # Flat field if requested
    #

    if flat_field is True:

        if verbose is True:
            print('Flat fielding the image...')

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

    avehdr = average_header_info(hdr, pair=pair)

    # Store the history

    try:

        old_history = avehdr['HISTORY']

        # remove it from the avehdr

        avehdr.pop('HISTORY')

    except:

        old_history = ''

    # Add useful things to header

    avehdr['MODULE'] = ['extract', ' Creation module']

    avehdr['VERSION'] = [setup.state['version'], ' pySpextool version']

    avehdr['COMBSTAT'] = [statistic, ' Combination statistic']

    avehdr['RBTHRESH'] = [robust_sigma, ' Robust threshold (if used)']

    avehdr['LINCOR'] = [linearity_correction, ' Linearity corrected?']

    avehdr['LINCRMAX'] = [setup.state['linearity_info']['max'],
                          ' Linearity correction maximum (DN)']

    avehdr['ORDSCLED'] = [scale_orders, ' Orders scaled?']

    avehdr['BGSUBED'] = [background_subtraction, ' Background subtraction?']

    avehdr['FLATED'] = [flat_field, ' Flat fielded?']

    avehdr['NIMAGES'] = [len(files), ' Number of images combined']

    avehdr['BEAMMODE'] = [beam_mode, ' Beam mode']

    avehdr['FILENAME'] = [output_name + '.fits', ' Filename']

    # Create the history

    # Strip the directories.

    file_names = [os.path.basename(i) for i in files]

    history = 'This image was created by combining the images ' + \
              ', '.join(file_names) + ' with a ' + statistic

    if statistic[0:6] == 'robust':

        history += ' with a sigma threshold of ' + str(robust_sigma) + '.  '

    else:

        history = history + '.'

    if linearity_correction is True:
        history = history + '  The raw data were corrected for non-linearity.'

    if scale_orders is True:
        history = history + '  The orders were scaled to a common flux level.'

    if background_subtraction is True:
        history = history + '  Background subtraction was accomplished by ' + \
                  'subtracting the median flux level on a ' + \
                  'column by column basis within each order.'

    if flat_field is True:
        history = history + '  The image was flat fielded.'

    history = old_history + history

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

    path = setup.state[output_directory + '_path']

    hdu = fits.HDUList([phdu, img_hdu, var_hdu, flg_hdu])
    hdu.writeto(os.path.join(path, output_name + '.fits'), overwrite=overwrite)

    if verbose is True:
        print('Wrote ' + output_name + '.fits to disk.')
