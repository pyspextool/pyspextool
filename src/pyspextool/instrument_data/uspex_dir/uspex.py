import numpy as np
from astropy.io import fits
import re
import os

from pyspextool.extract.locate_orders import locate_orders
from pyspextool.extract.normalize_flat import normalize_flat
from pyspextool.extract.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.extract.get_line_guess_position import get_line_guess_position
from pyspextool.extract.find_lines_1dxd import find_lines_1dxd
from pyspextool.extract.wavecal_solution_1d import wavecal_solution_1d
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.extract.get_spectral_pixelshift import get_spectral_pixelshift
from pyspextool.extract.make_interp_indices_1d import make_interp_indices_1d
from pyspextool.extract import config
from pyspextool.fit.polyfit import image_poly
from pyspextool.io.check import check_parameter
from pyspextool.io.files import make_full_path
from pyspextool.io.fitsheader import get_header_info
from pyspextool.io.fitsheader import average_header_info
from pyspextool.io.flat import write_flat
from pyspextool.io.flat import read_flatcal_file
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.wavecal import read_line_list
from pyspextool.io.wavecal import read_wavecal_file
from pyspextool.io.wavecal import write_wavecal_1d
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils.for_print import for_print
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.math import scale_data_stack
from pyspextool.utils.math import median_data_stack
from pyspextool.utils.split_text import split_text



try:
    from importlib.resources import files # Python 3.10+
except ImportError:
    from importlib_resources import files # Python <=3.9

def correct_amps(img):

    """
    To correct for bias voltage drift in an uSpeX FITS image

    Parameters
    ----------
    img : numpy array
        An uSpeX image

    Returns
    --------
    numpy.ndarray
        The uSpeX image with the bias variations "corrected".

    Notes
    -----
    There are 32 amplifiers that talk to 64 columns each.  The median
    intensity of the 64 reference pixels at the bottom of image are
    subtracted from all rows in the 64 columns.

    Example
    --------
    later

    """

    for i in range(0, 32):

        xl = 0+64*i
        xr = xl+63

        med = np.median(img[2044:2047+1, xl:xr+1])
        img[:, xl:(xr+1)] -= med

    return img    


def make_flat(files, output_name, prefix='flat-', suffix='.[ab]',
              extension='.fits*', input_method='index', normalize=True,
              verbose=True, overwrite=True, qafile=False):
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

    verbose : {True, False}, optional
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

        files = make_full_path(config.user['setup']['rawpath'], files,
                               indexinfo={'nint': config.state['nint'],
                                          'prefix': prefix,
                                          'suffix': suffix,
                                          'extension': extension},
                               exist=True)

    elif input_method == 'filename':

        files = make_full_path(config.user['setup']['rawpath'], files,
                               exist=True)

    else:

        raise ValueError('Unknown input_method.')

    # Load the FITS files into memory

    if verbose is True:
        print(' ')
        print('Make Flat Field')
        print('---------------')
        print('Loading FITS images...')

    img, var, hdr, mask = read_fits(files, config.state['linearity_info'],
                                    keywords=keywords, verbose=verbose)

    average_header = average_header_info(hdr)

    # and combine the masks

    #print(mask.dtype)
    #print('Combine----------------------')
    flag = combine_flag_stack(mask)
    #print(flag.dtype)
    #return

    #
    # Combine the images together
    #

    # Now scale their intensities to a common flux level

    if verbose is True:
        print('Scaling images...')

    simgs, svars, scales = scale_data_stack(img, None)

    # Now median the scaled images

    if verbose is True:
        print('Medianing the images...')

    med, munc = median_data_stack(simgs)

    #
    # Locate the orders
    #

    # Get the mode name and read modefile

    mode = hdr[0]['GRAT'][0]

    modefile = os.path.join(config.state['instrument_path'],
                            mode + '_flatinfo.fits')

    modeinfo = read_flatcal_file(modefile)

    if verbose is True:
        print('Locating the orders...')

    if qafile is True:

        qafileinfo = {'figsize':(8.5, 11),
                      'filepath':config.user['setup']['qapath'],
                      'filename':output_name, 'extension':'.pdf'}
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

        if verbose is True:
            print('Normalizing the median image...')

        nimg, nvar, rms = normalize_flat(med, edgecoeffs, modeinfo['xranges'],
                                         modeinfo['slith_arc'],
                                         modeinfo['nxgrid'],
                                         modeinfo['nygrid'],
                                         var=munc ** 2, oversamp=1,
                                         ybuffer=modeinfo['ybuffer'],
                                         verbose=False)

        if qafileinfo is not None:

            qafileinfo['filename'] = qafileinfo['filename']+'_normalized'

            orders_plotinfo = {'edgecoeffs':edgecoeffs,
                               'xranges':modeinfo['xranges'],
                               'orders':modeinfo['orders']}
            plot_image(nimg, mask=flag, orders_plotinfo=orders_plotinfo,
                       qafileinfo=qafileinfo)

    else:

        nimg = med
        nvar = munc ** 2
        rms = np.full((len(modeinfo['orders'])), np.nan)

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
               os.path.join(config.user['setup']['calpath'],
               output_name + '.fits'),
               overwrite=overwrite)

    if verbose is True:
        print('Flat field '+output_name+'.fits written to disk.\n')


def make_wavecal(files, flat_file, output_name, prefix='arc-', suffix='.[ab]',
                 extension='.fits*', input_method='index',
                 use_stored_solution=False, verbose=True, qafile_shift=True,
                 qafile_findlines=True, qafile_fitlines=True, overwrite=True):

    #
    # Check the parameters
    #

    check_parameter('make_wavecal', 'files',files,'str')

    check_parameter('make_wavecal', 'flat_file',flat_file,'str')

    check_parameter('make_wavecal', 'output_name',output_name,'str')

    check_parameter('make_wavecal', 'prefix',prefix,'str')

    check_parameter('make_wavecal', 'suffix',suffix,'str')

    check_parameter('make_wavecal', 'extension',extension,'str')

    check_parameter('make_wavecal', 'input_method',input_method,'str')

    check_parameter('make_wavecal', 'use_stored_solution',
                    use_stored_solution,'bool')

    check_parameter('make_wavecal', 'verbose', verbose, 'bool')

    check_parameter('make_wavecal', 'qafile_shift', qafile_shift, 'bool')

    check_parameter('make_wavecal', 'qafile_findlines',
                    qafile_findlines, 'bool')

    check_parameter('make_wavecal', 'qafile_fitlines',
                    qafile_fitlines, 'bool')

    check_parameter('make_wavecal', 'overwrite', overwrite, 'bool')                                
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

    # Create the file names

    if input_method == 'index':

        files = make_full_path(config.user['setup']['rawpath'], files,
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

    if verbose is True:

        print('Make Wavecal File')
        print('-----------------')
        print('Loading FITS images...')

    img, var, hdr, mask = read_fits(files, config.state['linearity_info'],
                                    keywords=keywords, verbose=verbose)

    #
    # Combine images as necessary
    #

    # Now scale their intensities to a common flux level

    if len(files) > 1:

        if verbose is True:
            print('Scaling images...')

        simgs, svars, scales = scale_data_stack(img, None)

    else:
        simgs = img

    # Now median the scaled images

    if len(files) > 1:

        if verbose is True:
            print('Medianing the images...')

        med, munc = median_data_stack(simgs)

    else:
        med = simgs

    #
    # Let's do the extraction of the "arc" spectra
    #

    # Read the flat and wavecal files

    flat_file = os.path.join(config.user['setup']['calpath'], flat_file)
    flatinfo = read_flat_fits(flat_file)

    wavecalfile = os.path.join(config.state['instrument_path'],
                               flatinfo['mode'] + '_wavecalinfo.fits')

    wavecalinfo = read_wavecal_file(wavecalfile)

    # Create wavecal and spatcal images

    wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                             flatinfo['nrows'],
                                             flatinfo['edgecoeffs'],
                                             flatinfo['xranges'],
                                             flatinfo['slith_arc'])

    # Extract the "arc"                                

    appos = np.full((np.size(flatinfo['orders']),1), flatinfo['slith_arc'] / 2)
    spectra = extract_extendedsource_1dxd(med, var, flatinfo['ordermask'],
                                          flatinfo['orders'], wavecal,
                                          spatcal, appos,
                                          wavecalinfo['apradius'],
                                          linmax_bitmask=None,
                                          badpixel_mask=None, bginfo=None,
                                          verbose=verbose)

    #
    # Find the pixel offset between these spectra and the disk spectra
    #

    # Get the anchor order and spectra.

    z = flatinfo['orders'] == wavecalinfo['xcororder']
    z = np.sum(np.where(flatinfo['orders'] == wavecalinfo['xcororder']))

    xanchor = np.arange(int(wavecalinfo['xranges'][z, 0]),
                        int(wavecalinfo['xranges'][z, 1] + 1), dtype=int)
    fanchor = np.squeeze(wavecalinfo['spectra'][z, 1, :])

    # Get the source order

    xsource = np.squeeze(spectra['spectra'][z][0, :])
    fsource = np.squeeze(spectra['spectra'][z][1, :])

    if qafile_shift is True:

        qafileinfo = {'figsize': (8.5, 11),
                      'filepath':config.user['setup']['qapath'],
                      'filename': output_name, 'extension': '.pdf'}
    else:
        qafileinfo = None

    offset = get_spectral_pixelshift(xanchor, fanchor, xsource, fsource,
                                     qafileinfo=qafileinfo)

    #
    # Are we using the stored solution?
    #

    if use_stored_solution is False:

       #
       # Locate the line positions
       #

        # Get the line list to search for lines

        filename = os.path.join(config.state['instrument_path'],
                                wavecalinfo['linelist'])

        lineinfo = read_line_list(filename, delta_to_microns=True)

        #  Determine the guess position and search range for each

        lineinfo = get_line_guess_position(wavecalinfo['spectra'],
                                           wavecalinfo['orders'],
                                           flatinfo['xranges'], lineinfo)

        # Add the shift offset to the results

        lineinfo['xguess'] = lineinfo['xguess'] + offset
        lineinfo['range_min_xguess'] = lineinfo['range_min_xguess'] + offset
        lineinfo['range_max_xguess'] = lineinfo['range_max_xguess'] + offset

        # Now find the lines

        if verbose:
            print('Finding the lines...')

        if qafile_findlines is True:

            qafileinfo = {'figsize': (8.5, 11),
                          'filepath': config.user['setup']['qapath'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the lines

        lineinfo = find_lines_1dxd(spectra['spectra'], wavecalinfo['orders'],
                                   lineinfo, flatinfo['slitw_pix'],
                                   qafileinfo=qafileinfo, verbose=verbose)

        #
        # Let's do the actual calibration
        #

        if verbose:
            print('Determining the wavelength solution...')

        # Get set up for either 1d of 1dxd

        if wavecalinfo['wcaltype'] == '1d':

            figsize = (6,4)
            xd = None

        else:

            figsize = (8.5,11)
            xd = {'homeorder':wavecalinfo['homeorder'],
                  'orderdeg':wavecalinfo['ordrdeg']}

        if qafile_fitlines is True:

            qafileinfo = {'figsize': figsize,
                          'filepath': config.user['setup']['qapath'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the solution

        solution = wavecal_solution_1d(wavecalinfo['orders'], lineinfo,
                                       wavecalinfo['dispdeg'], xd=xd,
                                       qafileinfo=qafileinfo,
                                       verbose=verbose)
        
    else:

#        with open('data.sav', 'rb') as f:
#           spectra, wavecalinfo, lineinfo, flatinfo, offset, wavecal, spatcal,   = pickle.load(f)

        if verbose:
            print('Using stored solution...')
        
        solution = {'coeffs':wavecalinfo['coeffs'],
                    'covar':wavecalinfo['covar'],
                    'rms':wavecalinfo['rms'],
                    'nlines':wavecalinfo['nlines'],
                    'ngood':wavecalinfo['ngood'],
                    'nbad':wavecalinfo['nbad']}

    #
    # Creating rectification indices
    #

    indices = []
    for i in range(flatinfo['norders']):

        idxs = make_interp_indices_1d(flatinfo['edgecoeffs'][i,:,:],
                                      flatinfo['xranges'][i,:],
                                      flatinfo['slith_arc'],
                                      array_output=True)

        indices.append(idxs)
        
        
        
            
    #
    # Write the wavecal file to disk.
    #

    if verbose:
        print('Writing wavecal to disk...')

    if wavecalinfo['wcaltype'] == '1d':

        write_wavecal_1d(flatinfo['ncols'], flatinfo['nrows'],
                         flatinfo['orders'], flatinfo['edgecoeffs'],
                         flatinfo['xranges'], solution['coeffs'],
                         solution['covar'], wavecalinfo['dispdeg'],
                         solution['rms']*1e4, solution['nlines'],
                         solution['ngood'], solution['nbad'],
                         wavecal, spatcal, indices, flatinfo['rotation'],
                         flat_file,
                         os.path.join(config.user['setup']['calpath'],
                                        output_name + '.fits'),
                         config.state['version'],
                         stored_solution=use_stored_solution,
                         overwrite=overwrite)

    elif wavecalinfo['wcaltype'] == '1dxd':

        write_wavecal_1d(flatinfo['ncols'], flatinfo['nrows'],
                         flatinfo['orders'], flatinfo['edgecoeffs'],
                         flatinfo['xranges'], solution['coeffs'],
                         solution['covar'], wavecalinfo['dispdeg'],
                         solution['rms']*1e4, solution['nlines'],
                         solution['ngood'], solution['nbad'],
                         wavecal, spatcal, indices, flatinfo['rotation'],
                         flat_file,
                         os.path.join(config.user['setup']['calpath'],
                                        output_name + '.fits'),
                         config.state['version'],
                         xd={'orderdeg':wavecalinfo['ordrdeg'],
                             'homeorder':wavecalinfo['homeorder']},
                         stored_solution=use_stored_solution,
                         overwrite=overwrite)

    else:
        print('unknown wcaltype.')

    if verbose:
        print('Wavecal '+output_name+'.fits written to disk.\n')                


def read_fits(files, lininfo, keywords=None, pair_subtract=False, rotate=0,
              linearity_correction=True, ampcor=False, verbose=False):

    """
    To read an upgraded SpeX FITS image file.

    Parameters
    ----------
    files : list of str
        A list of fullpaths to FITS files.

    lininfo : dict {'max':int,'bit':int}
        information to identify pixels beyond range of linearity correction

        'max' maximum value in DN
        'bit' the bit to set for pixels beyond `max`

    keywords : list of str, optional
        A list of FITS keyword to retain 

    pair_subtract : {False, True}, optional
        Set to pair subtract the images.  

    rotate : {0,1,2,3,4,5,6,7}, optional 
        Direction  Transpose?  Rotation Counterclockwise
        -------------------------------------------------

        0          No          None
        1          No          90 deg
        2          No          180 deg
        3          No          270 deg
        4          Yes         None
        5          Yes         90 deg
        6          Yes         180 deg
        7          Yes         270 deg

        The directions follow the IDL rotate function convention.
        
    linearity_correction : {True, False}, optional
        Set to correct for non-linearity.

    ampcor : {False, True}, optional
        Set to correct for amplifying drift (see uspexampcor.py)

    Returns
    --------
    tuple 
        The results are returned as (data,var,hdrinfo,bitmask) where
        data = the image(s) in DN/s
        var  = the variance image(s) in (DN/s)**2
        hdrinfo  = a list where element is a dict.  The key is the FITS 
        keyword and the value is a list consiting of the FITS value and FITS 
        comment.

    Example
    --------
    ?

    Modification History
    --------------------
    2022-05-25 - Written by M. Cushing, University of Toledo.
                 Based on the Spextool mc_readuspexfits.pro IDL program.
    """
    # Get setup information

    naxis1 = 2048
    naxis2 = 2048

    nfiles = len(files)

#    dolincor = [0, 1][lincor is not None]

    # Correct for non-linearity?

    if linearity_correction is True:

        linearity_file = os.path.join(config.state['instrument_path'],
                                      'uspex_lincorr.fits')
        lc_coeffs = fits.getdata(linearity_file)

    else:

        lc_coeffs = None

    # Get set up for linearity check

#    lininfo = config.state['linearity_info']

    bias_file = os.path.join(config.state['instrument_path'],'uspex_bias.fits')
    
    hdul = fits.open(bias_file)
    divisor = hdul[0].header['DIVISOR']
    bias = hdul[0].data / divisor
    hdul.close()

    if pair_subtract is True:

        #  Check to make sure the right number of files

        if (nfiles % 2) != 0:

            print('mc_readuspexfits:  Not an even number of images.')
            sys.exit(1)

        else:

            nimages = int(nfiles / 2)

    else:

        nimages = nfiles

    # Make empty arrays

    data = np.empty((nimages, naxis2, naxis1))
    var = np.empty((nimages, naxis2, naxis1))
    hdrinfo = []
    bitmask = np.empty((nimages, naxis2, naxis1), dtype=np.int8)

    # Load the data

    if pair_subtract is True:

        # pair subtraction

        for i in range(0, nimages):
            a = load_data(files[i * 2], lininfo, bias, keywords=keywords,
                          ampcor=ampcor, lccoeffs=lc_coeffs)

            b = load_data(files[i * 2 + 1], lininfo, bias, keywords=keywords,
                          ampcor=ampcor, lccoeffs=lc_coeffs)

            combmask = combine_flag_stack(np.stack((a[3], b[3])),
                                          nbits=lininfo['bit'] + 1)

            data[i, :, :] = idl_rotate(a[0] - b[0], rotate)
            var[i, :, :] = idl_rotate(a[1] + b[1], rotate)
            bitmask[i, :, :] = idl_rotate(combmask, rotate)

            hdrinfo.append(a[2])
            hdrinfo.append(b[2])

    else:

        for i in range(0, nimages):
            im, va, hd, bm = load_data(files[i], lininfo, bias,
                                       keywords=keywords, ampcor=ampcor,
                                       lccoeffs=lc_coeffs)

            data[i, :, :] = idl_rotate(im, rotate)
            var[i, :, :] = idl_rotate(va, rotate)
            bitmask[i, :, :] = idl_rotate(bm, rotate)

            hdrinfo.append(hd)

    return np.squeeze(data), np.squeeze(var), hdrinfo, np.squeeze(bitmask)


def load_data(file, lininfo, bias, keywords=None, ampcor=None, lccoeffs=None):

    readnoise = 12.0  # per single read
    gain = 1.5  # electrons per DN

    hdul = fits.open(file)
    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems

    itime = hdul[0].header['ITIME']
    coadds = hdul[0].header['CO_ADDS']
    ndrs = hdul[0].header['NDR']
    readtime = hdul[0].header['TABLE_SE']
    divisor = hdul[0].header['DIVISOR']

    #  Get set up for error propagation and store total exposure time

    rdvar = (2. * readnoise ** 2) / ndrs / coadds / itime ** 2 / gain ** 2
    crtn = (1.0 - readtime * (ndrs ** 2 - 1.0) / 3. / itime / ndrs)

    #  Read images, get into units of DN.

    img_p = hdul[1].data / divisor
    img_s = hdul[2].data / divisor

    #  Check for linearity maximum

    mskp = ((img_p < (bias - lininfo['max'])) * 2 ** lininfo['bit']).astype(np.uint8)
    msks = ((img_s < (bias - lininfo['max'])) * 2 ** lininfo['bit']).astype(np.uint8)

    #  Combine the masks 

    bitmask = combine_flag_stack(np.stack((mskp, msks)),
                                 nbits=lininfo['bit'] + 1)

    #  Create the image

    img = img_p - img_s

    #  Correct for amplifier offsets

    if ampcor:
        img = correct_uspex_amp(img)

    #  Determine the linearity correction for the image

    if lccoeffs is not None:
        cor = image_poly(img, lccoeffs)
        cor = np.where(cor == 0, 1, cor)

        #  Now set the corrections to unity for pixels > lincormax

        cor = np.where(bitmask == 2 ** lininfo['bit'], 1, cor)

        #  Set black pixel corrections to unity as well.

        cor[:, 0:3 + 1] = 1.0
        cor[:, 2044:2047 + 1] = 1.0
        cor[0:3 + 1, :] = 1.0
        cor[2044:2047 + 1, :] = 1.0

        # Apply the corrections

        img /= cor

        # Delete unecessary files

        del cor, img_p, img_s

    # Create the actual image.
    # Convert image back to total DN for error propagation

    img = img * divisor

    # Compute the variance and the final image

    var = np.absolute(img) * crtn / ndrs / (coadds ** 2) / \
          (itime ** 2) / gain + rdvar 
    img = img / divisor / itime

    # Collect header information

    hdr = get_header(hdul[0].header)

    hdul.close()

    return [img, var, hdr, bitmask]


def get_header(hdr, keywords=None):
    # Grab keywords if requested

    if keywords:

        hdrinfo = get_header_info(hdr, keywords=keywords)

    else:

        hdrinfo = get_header_info(hdr)

    #  Grab require keywords and convert to standard Spextool keywords

    # Airmass 

    hdrinfo['AM'] = [hdr['TCS_AM'], ' Airmass']

    # Hour angle

    val = hdr['TCS_HA']
    m = re.search('[-]', '[' + val + ']')
    if not m:
        val = '+' + val.strip()
    hdrinfo['HA'] = [val, ' Hour angle (hours)']

    # Position Angle

    hdrinfo['PA'] = [hdr['POSANGLE'], ' Position Angle E of N (deg)']

    # Dec 

    val = hdr['TCS_DEC']
    m = re.search('[-]', '[' + val + ']')
    if not m:
        val = '+' + val.strip()
    hdrinfo['DEC'] = [val, ' Declination, FK5 J2000']

    # RA

    hdrinfo['RA'] = [hdr['TCS_RA'].strip(), ' Right Ascension, FK5 J2000']

    # COADDS, ITIME

    coadds = hdr['CO_ADDS']

    itime = hdr['ITIME']
    hdrinfo['ITIME'] = [itime, ' Integration time (sec)']
    hdrinfo['NCOADDS'] = [coadds, ' Number of COADDS']
    hdrinfo['IMGITIME'] = [coadds * itime,
                           ' Image integration time, NCOADDSxITIME (sec)']

    # Time

    hdrinfo['TIME'] = [hdr['TIME_OBS'].strip(), ' Observation time in UTC']

    # Date

    hdrinfo['DATE'] = [hdr['DATE_OBS'].strip(), ' Observation date in UTC']

    # MJD

    hdrinfo['MJD'] = [hdr['MJD_OBS'], ' Modified Julian date OBSDATE+TIME_OBS']

    # FILENAME

    hdrinfo['FILENAME'] = [hdr['IRAFNAME'].strip(), ' Filename']

    # MODE

    hdrinfo['MODE'] = [hdr['GRAT'].strip(), ' Instrument Mode']

    # INSTRUMENT

    hdrinfo['INSTR'] = ['SpeX', ' Instrument']

    # now move the comment key to the end

    comment = hdrinfo.pop('COMMENT')
    hdrinfo['COMMENT'] = comment

    return hdrinfo
