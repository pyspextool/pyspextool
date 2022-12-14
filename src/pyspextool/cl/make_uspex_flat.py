import numpy as np
import os

from pyspextool.cl import config
from pyspextool.io.files import make_full_path
from pyspextool.io.fitsheader import average_header_info
from pyspextool.io.flat import write_flat
from pyspextool.io.flat import read_flatcal_file
from pyspextool.io.read_uspex_fits import main as readfits
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.math import scale_data_stack
from pyspextool.utils.math import median_data_stack
from pyspextool.utils.split_text import split_text
from pyspextool.calibration.locate_orders import locate_orders
from pyspextool.calibration.normalize_flat import normalize_flat

from astropy.io import fits


def make_uspex_flat(files, output_name, prefix='flat-', suffix='.[ab]',
                    extension='.fits*', input_method='index',
                    normalize=True, clupdate=True,
                    overwrite=True, qafile=False):
    """
    To create a normalized uSpeX flat field file.

    Parameters
    ----------
    files : str or list of str.
        See input_method parameter for more information.
        
    output_name : str
        The filename of the flat field image to written to disk.      

    prefix : str, default='flat-', optional
        The prefix of the FITS file name (i.e. the stuff before the 
        file numbers.)

    suffix : str, default='.[ab]', optional
        The suffix of the FITS file name (i.e. the stuff after the
        file numbers but before the extension.)

    extension : str, default='.fits*', optional
        The file extension

    input_method : {'index', 'filename'}, optional
        `index`: `files` is a str of file numbers, e.g., '1-3,5-10'.

        `filename`: a str or list of strings with full file names, e.g.
            ['flat-00001.a.fits',''flat-00002.a.fits']

    normalize : {True, False}, optional
        Set to True to normalize the orders.

    clupdate : {True, False}, optional
        Set to True for command line updates during execution. 

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.

    qafile : {True, False}, optional
        Set to create the quality assurance plot.

    Returns
    -------
    None
        Writes a FITS file to disk.

    Examples
    --------
    later

    """

    #
    # Load the files into memory
    #

    # Get the user requested keywords

    keywords = config.state['xspextool_keywords']

    # Add GRAT and DIT so that you can determine the mode.  Users will just have
    # to live with them in their files

    if 'GRAT' not in keywords:
        keywords.append('GRAT')

    if 'DIT' not in keywords:
        keywords.append('DIT')

    # Now create the file names

    if input_method == 'index':

        files = make_full_path(config.state['rawpath'], files,
                               indexinfo={'nint': config.state['nint'],
                                          'prefix': prefix,
                                          'suffix': suffix,
                                          'extension': extension},
                               exist=True)

    elif input_method == 'filename':

        files = make_full_path(config.state['rawpath'], files, exist=True)

    else:

        raise ValueError('Unknown input_method.')

    # Load the FITS files into memory

    if clupdate is True:
        print(' ')
        print('Creating the flat field...')
        print('Loading FITS images...')

    img, var, hdr, mask = readfits(files, config.state['linearity_info'],
                                   keywords=keywords, clupdate=clupdate)

    # And then average the headers

    average_header = average_header_info(hdr)

    # and combine the masks

#    print(mask.dtype)
    print('Combine----------------------')
    flag = combine_flag_stack(mask)
#    print(flag.dtype)
    return
    
    #
    # Combine the images together
    #

    # Now scale their intensities to a common flux level

    if clupdate is True:
        print('Scaling images...')

    simgs, svars, scales = scale_data_stack(img, None)

    # Now median the scaled images

    if clupdate is True:
        print('Medianing the images...')

    med, munc = median_data_stack(simgs)

    #
    # Locate the orders
    #

    # Get the mode name and read modefile

    mode = hdr[0]['GRAT'][0]

    modefile = os.path.join(config.state['packagepath'], 'instruments',
                            config.state['instrument'], 'data',
                            mode + '_flatinfo.fits')

    modeinfo = read_flatcal_file(modefile)

    if clupdate is True:
        print('Locating the orders...')

    if qafile is True:

        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
                      'filename': output_name, 'extension': '.pdf'}
    else:
        qafileinfo = None

    edgecoeffs = locate_orders(med, modeinfo['guesspos'], modeinfo['xranges'],
                               modeinfo['step'], modeinfo['slith_range'],
                               modeinfo['edgedeg'], modeinfo['ybuffer'],
                               modeinfo['flatfrac'], modeinfo['comwidth'],
                               qafileinfo=qafileinfo)

    #
    # Normalize the spectrum if requested
    #

    if normalize is True:

        if clupdate is True:
            print('Normalizing the median image...')

        nimg, nvar, rms = normalize_flat(med, edgecoeffs, modeinfo['xranges'],
                                         modeinfo['slith_arc'],
                                         modeinfo['nxgrid'],
                                         modeinfo['nygrid'],
                                         var=munc ** 2, oversamp=1,
                                         ybuffer=modeinfo['ybuffer'],
                                         clupdate=False)

    else:

        nimg = med
        nvar = munc ** 2
        rms = np.full((len(modeinfo['orders'])), np.nan)

    #
    # Write the file to disk
    #

    # Create the HISTORY

    if clupdate is True:
        print('Writing flat to disk...')

    basenames = []
    for file in files:
        basenames.append(os.path.basename(file))

    history = 'This flat was created by scaling the files ' + \
              ', '.join(str(b) for b in basenames) + ' to a common '+ \
              'median flux level and then medianing the scaled imges.  ' + \
              'The variance is given by (1.4826*MAD)**2/nimages where MAD ' + \
              'is the median absolute deviation.  The zeroth bit of pixels ' + \
              'in the third extension are set if their corresponding '+ \
              'intensity values are greater than LINCORMAX.  User selected ' + \
              'FITS keywords are from the first frame in the series.'

    history = split_text(history)

    # Get the slit widths and resolving power and write to disk

    slitw_arc = float(average_header['SLIT'][0][0:3])
    slitw_pix = slitw_arc / modeinfo['ps']

    resolvingpower = modeinfo['rpppix'] / slitw_pix

    write_flat(nimg, nvar, flag, average_header, modeinfo['rotation'],
               modeinfo['orders'], edgecoeffs, modeinfo['xranges'],
               modeinfo['ps'], modeinfo['slith_pix'], modeinfo['slith_arc'],
               slitw_pix, slitw_arc, mode, rms, resolvingpower,
               config.state['version'], history,
               os.path.join(config.state['calpath'], output_name + '.fits'),
               overwrite=overwrite)

    if clupdate is True:
        print('Flat field '+output_name+'.fits written to disk.')

