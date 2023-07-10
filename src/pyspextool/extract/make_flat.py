import importlib
import os
import numpy as np
from pyspextool import config as setup
from pyspextool.extract import config as extract

from pyspextool.extract.locate_orders import locate_orders
from pyspextool.extract.normalize_flat import normalize_flat
from pyspextool.io.check import check_parameter
from pyspextool.io.files import make_full_path
from pyspextool.io.fitsheader import average_header_info
from pyspextool.io.flat import read_flatcal_file
from pyspextool.io.flat import write_flat

from pyspextool.plot.plot_image import plot_image

from pyspextool.utils import math
from pyspextool.utils.split_text import split_text


def make_flat(files, output_name, extension='.fits*', normalize=True,
              overwrite=True, qa_plot=None, qa_plotsize=(8,8),
              qa_file=None, verbose=None):
    """
    To create a (normalized) pyspextool flat field file.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
        
    output_name : str
        The filename of the flat field image to written to disk.      

    extension : str, default='.fits*'
        The file extension

    normalize : {True, False}, optional
        Set to True to normalize the orders.

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.

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

    Returns
    -------
     None
        Writes a FITS file and QA plots to disk.

    """

    #
    # Check the input parmameters
    #

    check_parameter('make_flat', 'files', files, ['str', 'list'])

    if isinstance(files, list):
        check_parameter('make_flat', 'files[0]', files[0], 'str')
        check_parameter('make_flat', 'files[1]', files[1],
                        ['str', 'list', 'int'])

    check_parameter('make_flat', 'files', files, ['str', 'list'])

    check_parameter('make_flat', 'output_name', output_name, 'str')

    check_parameter('make_flat', 'extension', extension, 'str')

    check_parameter('make_flat', 'normalize', normalize, 'bool')

    check_parameter('make_flat', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('make_flat', 'overwrite', overwrite, 'bool')

    check_parameter('make_flat', 'qa_file', qa_file, ['NoneType', 'bool'])

    check_parameter('make_flat', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('make_flat', 'qa_plotsize', qa_plotsize,
                    ['NoneType', 'tuple'])

    #
    # Check the qa variables and set to system default if need be.
    #

    if qa_file is None:
        qa_file = setup.state['qa_file']

    if qa_plot is None:
        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Let the user know what you are doing.
    #
    
    if verbose is True:
        print('Generating Flat Field')
        print('---------------------')        

        
    #
    # Load the instrument module for the read_fits program
    #

    module = 'pyspextool.instrument_data.' + setup.state['instrument'] + \
             '_dir.' + setup.state['instrument']

    instr = importlib.import_module(module)

    #
    # Load the files into memory
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

    # Get the mode name from the first file in order to get the array
    # rotation value.

    img, var, hdr, mask = instr.read_fits(files, setup.state['linearity_info'],
                                          verbose=False)

    mode = hdr[0]['MODE'][0]

    modefile = os.path.join(setup.state['instrument_path'],
                            mode + '_flatinfo.fits')

    modeinfo = read_flatcal_file(modefile)

    # Load the FITS files into memory

    img, var, hdr, mask = instr.read_fits(files,
                                          setup.state['linearity_info'],
                                          rotate=modeinfo['rotation'],
                                    keywords=setup.state['extract_keywords'],
                                          verbose=verbose)

    average_header = average_header_info(hdr)

    # and combine the masks

    flag = math.combine_flag_stack(mask)

    #
    # Combine the images together
    #

    # Now scale their intensities to a common flux level

    if verbose is True:
        print('Scaling images...')

    simgs, svars, scales = math.scale_data_stack(img, None)

    # Now median the scaled images

    if verbose is True:
        print('Medianing the images...')

    med, munc = math.median_data_stack(simgs)

    #
    # Locate the orders
    #

    # Get set up for QA plotting

    if qa_file is True:

        qa_fileinfo = {'figsize': (6,6),
                       'filepath': setup.state['qa_path'],
                       'filename': output_name + '_locateorders',
                       'extension': setup.state['qa_extension']}

    else:

        qa_fileinfo = None
        
    
    if verbose is True:
        print('Locating the orders...')

    edgecoeffs = locate_orders(med, modeinfo['guesspos'],
                                   modeinfo['xranges'], modeinfo['step'],
                                   modeinfo['slith_range'],
                                   modeinfo['edgedeg'], modeinfo['ybuffer'],
                                   modeinfo['flatfrac'], modeinfo['comwidth'],
                                   qa_plot=qa_plot, qa_fileinfo=qa_fileinfo,
                                   qa_plotsize=qa_plotsize)

    #
    # Normalize the flat if requested
    #

    if normalize is True:

        if verbose is True:
            print('Normalizing the median image...')

        nimg, nvar, rms = normalize_flat(med, edgecoeffs, modeinfo['xranges'],
                                         modeinfo['slith_arc'],
                                         modeinfo['nxgrid'],
                                         modeinfo['nygrid'],
                                         var=munc ** 2, oversamp=1,
                                         ybuffer=modeinfo['ybuffer'],
                                         verbose=verbose)

    else:

        nimg = med
        nvar = munc ** 2
        rms = np.full((len(modeinfo['orders'])), np.nan)

    # Protect against zeros

    z = np.where(nimg == 0)
    nimg[z] = 1

    # Make the qaplot

    if qa_plot is True:

        orders_plotinfo = {'edgecoeffs': edgecoeffs,
                           'xranges': modeinfo['xranges'],
                           'orders': modeinfo['orders']}
        
        plot_image(nimg, orders_plotinfo=orders_plotinfo, plot_size=qa_plotsize)
            
    if qa_file is True:

        orders_plotinfo = {'edgecoeffs': edgecoeffs,
                           'xranges': modeinfo['xranges'],
                           'orders': modeinfo['orders']}

        qa_fileinfo = {'figsize': (6,6),
                       'filepath': setup.state['qa_path'],
                       'filename': output_name + '_normalized',
                       'extension': setup.state['qa_extension']}

        plot_image(nimg, orders_plotinfo=orders_plotinfo, file_info=qa_fileinfo)

    #
    # Write the file to disk
    #

    # Create the HISTORY

    if verbose is True:
        print('Writing flat to disk...')

    basenames = []
    for file in files:
        basenames.append(os.path.basename(file))

    history = 'This flat was created by scaling the files ' + \
              ', '.join(str(b) for b in basenames) + ' to a common ' + \
              'median flux level and then medianing the scaled imges.  ' + \
              'The variance is given by (1.4826*MAD)**2/nimages where MAD ' + \
              'is the median absolute deviation.  The zeroth bit of pixels ' + \
              'in the third extension are set if their corresponding ' + \
              'intensity values are greater than LINCORMAX.  User selected ' + \
              'FITS keywords are from the first frame in the series.'

    history = split_text(history)

    # Get the slit widths and resolving power and write to disk

    slitw_arc = float(average_header['SLIT'][0][0:3])
    slitw_pix = slitw_arc / modeinfo['ps']

    resolvingpower = modeinfo['rpppix'] / slitw_pix

    write_flat(nimg, nvar, flag, average_header, modeinfo['rotation'],
               modeinfo['orders'], edgecoeffs, modeinfo['xranges'],
               modeinfo['ybuffer'], modeinfo['ps'], modeinfo['slith_pix'],
               modeinfo['slith_arc'], slitw_pix, slitw_arc, mode, rms,
               resolvingpower, setup.state['version'], history,
               os.path.join(setup.state['cal_path'],
                            output_name + '.fits'), overwrite=overwrite)

    if verbose is True:
        print('Flat field ' + output_name + '.fits written to disk.\n')
