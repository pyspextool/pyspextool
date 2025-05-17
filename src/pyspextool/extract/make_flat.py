import importlib
from os.path import join, basename
import numpy as np
import logging
from decimal import *

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import files_to_fullpath
from pyspextool.io.fitsheader import average_headerinfo
from pyspextool.extract.flat import locate_orders, read_flatcal_file, \
write_flat, normalize_flat
from pyspextool.extract.images import make_ordermask
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils import math
from pyspextool.utils.split_text import split_text
from pyspextool.utils.arrays import idl_rotate


def make_flat(files:list | str,
              output_name:str,
              linearity_correction=False,
              detector_info:dict=None,
              normalize:bool=True,
              verbose:bool=None,
              qa_show:bool=None,
              qa_showscale:int | float=None,
              qa_showblock:bool=None,
              qa_write:bool=None):

    """
    To create a (normalized) pySpextool flat field file.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
        
    output_name : str
        The root of the filename of the flat field image to written to disk.

    linearity_correction : {False, True}
        Set to True to correct for non-linearity.
        Set to False to not correct for non-linearity.

    detector_info : dict, deafult None
        A dictionary with any information that needs to be passed to the
        instrument specific readfits program.  
        
    normalize : {True, False}, optional
        Set to True to normalize the orders.
        Set to False to NOT normalize the orders.

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
     None
        Writes a FITS file and QA plots to disk.

    """

    #
    # Check the input parmameters
    #

    check_parameter('make_flat', 'files', files, ['str', 'list'])

    check_parameter('make_flat', 'output_name', output_name, 'str')

    check_parameter('make_flat', 'linearity_correction', linearity_correction,
                    'bool')

    check_parameter('make_flat', 'detector_info', detector_info,
                    ['dict', 'NoneType'])
    
    check_parameter('make_flat', 'normalize', normalize, 'bool')

    check_parameter('make_flat', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('make_flat', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('make_flat', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('make_flat', 'qa_showscale', qa_showscale,
                    ['NoneType', 'int', 'float'])    

    check_parameter('make_flat', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)        
    #
    # Let the user know what you are doing.
    #

    message = ' Generating Flat Field'
    logging.info(message+'\n'+'-'*(len(message)+5)+'\n')

    #
    # Load the instrument module for the read_fits program
    #

    module = 'pyspextool.instruments.' + setup.state['instrument'] + \
             '.' + setup.state['instrument']

    instr = importlib.import_module(module)

    #
    # Load the files into memory
    #

    # Create the file names

    result = files_to_fullpath(setup.state['raw_path'],
                               files,
                               setup.state['nint'],
                               setup.state['suffix'],
                               setup.state['search_extension'])
        
    files = result[0]
    readmode = result[1]
    
    # Get the mode name from the first file in order to get the array
    # rotation value.

    result = instr.read_fits(files,
                             setup.state['linearity_info'],
                             verbose=False)
    hdr = result[2]
    
    mode = hdr[0]['MODE'][0]
    modefile = join(setup.state['instrument_path'], mode + '_flatinfo.fits')

    modeinfo = read_flatcal_file(modefile)

    # Load the FITS files into memory

    logging.info(' Loading the images.')
    result = instr.read_fits(files,
                             setup.state['linearity_info'],
                             rotate=modeinfo['rotation'],
                             keywords=setup.state['extract_keywords'],
                             linearity_correction=linearity_correction,
                             extra=detector_info,
                             verbose=qa['verbose'])
    
    img = result[0]
    var = result[1]
    hdr = result[2]
    mask = result[3]

    # Average the headers and combine the masks

    average_header = average_headerinfo(hdr)    
    flag = math.combine_flag_stack(mask)

#    # Load the bad pixel mask

#    bad_pixel_mask = idl_rotate(setup.state['raw_bad_pixel_mask'],
#                                modeinfo['rotation'])

    #
    # Combine the images
    #

    # Now scale their intensities to a common flux level

    logging.info(' Scaling the images.')

    simgs, svars, scales = math.scale_data_stack(img, None)

    # Now median the scaled images

    logging.info(' Medianing the images.')

    med, munc = math.median_data_stack(simgs)

    #
    # Locate the orders
    #

    # Get set up for QA plotting

    if qa['write'] is True:

        filename = output_name + '_locateorders' + setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

    else:

        fullpath = None
            
    logging.info(' Locating the orders.')

    result = locate_orders(med,
                           modeinfo['guesspos'],
                           modeinfo['xranges'],
                           modeinfo['step'],
                           modeinfo['slith_range'],
                           modeinfo['edgedeg'],
                           modeinfo['ybuffer'],
                           modeinfo['flatfrac'],
                           modeinfo['comwidth'],
                           qa_plotnumber=setup.plots['locate_orders'],
                           qa_fontsize=setup.plots['font_size'],
                           qa_figuresize=setup.plots['square_size'],
                           qa_show=qa['show'],
                           qa_showblock=qa['showblock'],
                           qa_showscale=qa['showscale'],
                           qa_fullpath=fullpath,
                           debug=False)

    edgecoeffs = result[0]
    xranges = result[1]

    #
    # Normalize the flat if requested
    #

    if normalize is True:

        logging.info(' Normalizing the median flat.')

        nimg, nvar, rms = normalize_flat(med,
                                         edgecoeffs,
                                         xranges,
                                         modeinfo['slith_arc'],
                                         modeinfo['nxgrid'],
                                         modeinfo['nygrid'],
                                         var=munc ** 2,
                                         ybuffer=modeinfo['ybuffer'],
                                         verbose=qa['verbose'])

    else:

        nimg = med
        nvar = munc ** 2
        rms = np.full((len(modeinfo['orders'])), np.nan)
    
    # Protect against zeros

    z = np.where(nimg <= 0)
    nimg[z] = 1

    #
    # Create the order mask
    #

    nrows, ncols = np.shape(med)
    ordermask = make_ordermask(ncols,
                               nrows,
                               edgecoeffs,
                               xranges,
                               modeinfo['orders'])

    #        
    # Make the QA plot
    #
    
    orders_plotinfo = {'edgecoeffs': edgecoeffs,
                       'xranges': xranges,
                       'orders': modeinfo['orders']}
              
    if qa['show'] is True:
        
        plot_image(nimg,
                   mask=flag,
                   orders_plotinfo=orders_plotinfo,
                   figure_size=(setup.plots['square_size'][0]*qa['showscale'],
                                setup.plots['square_size'][1]*qa['showscale']),
                   font_size=setup.plots['font_size']*qa['showscale'],
                   showblock=qa['showblock'],
                   plot_number=setup.plots['flat'])
            
    if qa['write'] is True:

        filename = output_name + '_normalized' + setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

        plot_image(nimg,
                   mask=flag,
                   orders_plotinfo=orders_plotinfo,
                   output_fullpath=fullpath,
                   figure_size=setup.plots['square_size'],
                   font_size=setup.plots['font_size'])

    #
    # Write the file to disk
    #

    # Create the HISTORY

    basenames = []
    for file in files:
        basenames.append(basename(file))

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
    getcontext().prec = 2
    slitw_pix = float(Decimal(slitw_arc) / Decimal(modeinfo['ps']))
    resolvingpower =int(Decimal(modeinfo['rpppix']) / Decimal(slitw_pix))

    # Write it to disk.
    
    output_fullpath =  join(setup.state['cal_path'], output_name + '.fits')

    write_flat(nimg,
               nvar,
               flag,
               ordermask,
               average_header,
               modeinfo['rotation'],
               modeinfo['orders'],
               edgecoeffs,
               xranges,
               modeinfo['ybuffer'],
               modeinfo['ps'],
               modeinfo['slith_pix'],
               modeinfo['slith_arc'],
               slitw_pix,
               slitw_arc,
               mode,
               rms,
               resolvingpower,
               setup.state['version'],
               history,
               output_fullpath,
               overwrite=True)

    logging.info(' Flat field file ' + output_name + '.fits written to disk.\n')
