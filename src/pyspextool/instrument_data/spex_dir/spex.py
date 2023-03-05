import numpy as np
from astropy.io import fits
from astropy.time import Time
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
from pyspextool.utils import math
from pyspextool.utils.split_text import split_text



try:
    from importlib.resources import files # Python 3.10+
except ImportError:
    from importlib_resources import files # Python <=3.9



def correct_linearity(image, itime, slowcnts, ndr, coefficients):

    """
    Correct for non-linearity in SpeX images.

    Parameters
    ----------
    image : ndarray
        (nrows, ncols) image in DN (raw/DIVISOR).

    itime : int or float
        The integration time

    slowcnts : int or float
        The slowcnt value

    ndr : int
        The number of none destructive reads

    coefficients : ndarray
        A (nrows, ncols, ncoefficients) array of linearity coefficients.

    Returns
    -------
    A (nrows, ncols) linearity corrected image.

    Notes
    -----
    This implements the non-linearity correction method described in 
    Vacca et al. (2004, Vacca, PASP, 116, 352).  Since the average pedestal 
    and signal images are not stored, they must be estimated as described 
    in S2.3.

    """

    # Store possible values of of slowcnts and corresponding array read times.
    
    slowcnts_possible = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                  14, 15, 16, 17, 18, 19, 20])
        
    tread = np.array([0.10, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32,
                      0.34, 0.36,0.38, 0.40, 0.41, 0.43, 0.45, 0.47, 0.49,
                      0.51])

    # Pull the zeroth coefficient out for use later

    c0 = coefficients[0,:,:]

    # Figure out the slow
    
    z = slowcnts_possible == slowcnts

    #
    # Do the first itteration
    #
    
    # Estimate the mean pedistal and signal images

    mean_pedestal_1 = image*np.squeeze(tread[z])*(ndr+1)/2/itime
    mean_signal_1 = image + mean_pedestal_1
    
    # Do the corrections

    correction = c0/spex_linearity_imagepoly(mean_pedestal_1,coefficients)
    correction = np.clip(correction, 1.0, None)
    corrected_mean_pedestal_1 = mean_pedestal_1*correction

    correction = c0/spex_linearity_imagepoly(mean_signal_1,coefficients)
    correction = np.clip(correction, 1.0, None)
    corrected_mean_signal_1 = mean_signal_1*correction    

    # Create a new estimate of the image

    image_2 = corrected_mean_signal_1 - corrected_mean_pedestal_1 

    #
    # Do the second itteration
    #
    
    # Estimate the mean pedistal and signal images

    mean_pedestal_2 = image_2*np.squeeze(tread[z])*(ndr+1)/2/itime
    mean_signal_2 = image_2 + mean_pedestal_2

    # Do the corrections
    
    correction = c0/spex_linearity_imagepoly(mean_pedestal_2,coefficients)
    correction = np.clip(correction, 1.0, None)
    corrected_mean_pedestal_2 = mean_pedestal_2*correction

    correction = c0/spex_linearity_imagepoly(mean_signal_2,coefficients)
    correction = np.clip(correction, 1.0, None)
    corrected_mean_signal_2 = mean_signal_2*correction    

    # Create a new estimate of the image

    image_3 = corrected_mean_signal_2 - corrected_mean_pedestal_2

    # replace NaNs with zeros

    z = np.isnan(image_3)
    
    if np.sum(z) != 0:

        image_3[z] = 0
        
    #    return image

    return image_3


        
def get_header(header, keywords=None):

    """
    To get header values from the FITS file.

    This routine both grab existing keywords from the FITS header but also 
    creates and populates the Spextool-required ones.

    Parameters
    ----------
    hdr : astropy HDU.header
         see https://docs.astropy.org/en/stable/io/fits/index.html

    keywords : list of str, optional
         A list of keywords to store.

    Returns
    -------
    dict

        A dict where each key is the FITS keyword and each value is a list 
        with [keyword value, keyword comment].

        To access the value and comment of keyword XXX:
        hdrinfo['XXX'][0] = value
        hdrinfo['XXX'][1] = comment

    """

    #
    # Grab keywords if requested, if not grab them all
    #
    
    if keywords:

        hdrinfo = get_header_info(header, keywords=keywords)

    else:

        hdrinfo = get_header_info(header)

    #
    # Now create the required Spextool keywords
    #
        
    # Airmass 

    hdrinfo['AM'] = [header['AIRMASS'], ' Airmass']

    # Hour angle

    val = header['HA']
    m = re.search('[-]', '[' + val + ']')
    if not m:
        val = '+' + val.strip()
    hdrinfo['HA'] = [val, ' Hour angle (hours)']

    # Position Angle

    hdrinfo['PA'] = [header['POSANGLE'], ' Position Angle E of N (deg)']

    # Dec 

    val = header['DEC']
    m = re.search('[-]', '[' + val + ']')
    if not m:
        val = '+' + val.strip()
    hdrinfo['DEC'] = [val, ' Declination']

    # RA

    hdrinfo['RA'] = [header['RA'].strip(), ' Right Ascension']

    # COADDS, ITIME

    coadds = header['CO_ADDS']

    itime = header['ITIME']
    hdrinfo['ITIME'] = [itime, ' Integration time (sec)']
    hdrinfo['NCOADDS'] = [coadds, ' Number of COADDS']
    hdrinfo['IMGITIME'] = [coadds * itime,
                           ' Image integration time, NCOADDSxITIME (sec)']

    # Time

    hdrinfo['TIME'] = [header['TIME_OBS'].strip(), ' Observation time in UTC']

    # Date

    hdrinfo['DATE'] = [header['DATE_OBS'].strip(), ' Observation date in UTC']

    # MJD

    time_string = header['DATE_OBS']+'T'+header['TIME_OBS']
    tobj = Time(time_string)
    
    hdrinfo['MJD'] = [tobj.mjd, ' Modified Julian date OBSDATE+TIME_OBS']

    # FILENAME

    hdrinfo['FILENAME'] = [header['IRAFNAME'].strip(), ' Filename']

    # MODE

    hdrinfo['MODE'] = [header['GRAT'].strip(), ' Instrument Mode']

    # INSTRUMENT

    hdrinfo['INSTR'] = ['SpeX', ' Instrument']

    # now move the comment key to the end

    comment = hdrinfo.pop('COMMENT')
    hdrinfo['COMMENT'] = comment

    return hdrinfo



def load_data(file, linearity_info, keywords=None, linearity_coefficients=None):

    
    """
    To read a single SpeX FITS file, correct for linearity, and grab keywords.

    Parameters
    ----------
    files : str
        The fullpath to a SpeX FITS file.

    linearity_info : dict 

        `'max'` : int
            Value in DN beyond which the linearity correction fails

        `'bit'` : int
            The bit to set in the mask for pixels > dict['max'].

    keywords : list of str, optional
        A list of FITS keyword to retain.  If None, all are stored.


    linear_coefficients : ndarray, optional
        A (ncoefficients, nrows, ncols) array of linearity correction 
        coefficients.

    Returns
    --------
    tuple : (ndarray, ndarray, ndarray, ndarray)

        tuple(0) : The (nrows, ncols) image in DN s-1

        tuple(1) : The (nrows, ncols) variance in (DN s-1)**2

        tuple(2) : A dict where each key is the FITS keyword and each value 
            is a list with [keyword value, keyword comment].

            To access the value and comment of keyword XXX:
            hdrinfo['XXX'][0] = value
            hdrinfo['XXX'][1] = comment

        tuple(3) : A (nrows, ncols) bitset mask.  The zeroth bit of pixels with 
            values > linearity_info['max'] is set.

    """

    #
    # Set the values for the read noise and gain
    #

    readnoise = 50.0  # per single read
    gain = 13  # electrons per DN

    #
    # Open the file and grab important values 
    #
    
    hdul = fits.open(file)
    hdul[0].verify('silentfix')  # this was needed to correct hdr problems

    itime = hdul[0].header['ITIME']
    coadds = hdul[0].header['CO_ADDS']
    ndrs = hdul[0].header['NDR']
    readtime = hdul[0].header['TABLE_MS']/1000.
    divisor = hdul[0].header['DIVISOR']

    #  Get set up for error propagation 

    rdvar = (2. * readnoise**2) / ndrs / coadds / itime**2 / gain**2
    crtn = (1.0 - readtime * (ndrs**2 - 1.0) / 3. / itime / ndrs)

    # Grab the image from the Header Data Unit List
    
    raw_img = hdul[0].data

    #
    #  Identify pixels with DN > linearity_info['max']
    #
    
    bitmask = (raw_img/divisor > linearity_info['max'] *
               2**linearity_info['bit']).astype(np.uint8)

    #  Determine the linearity correction for the image

    if linearity_coefficients is not None:

        slowcnt = hdul[0].header['SLOWCNT']

        # pass in the image divided by divisor
        
        img = correct_linearity(raw_img/divisor, itime, slowcnt, ndrs,
                                linearity_coefficients)

        # convert back to total DN

        img = img*divisor

        # Replace any pixels > linearity_info['max'] with original values

        z = np.where(bitmask == 2**linearity_info['bit'])
        img[z] = raw_img[z]
        
    else:

        img = raw_img
                    
    # Compute the variance and the final image

    var = np.absolute(img) * crtn / ndrs / (coadds**2) / \
          (itime**2) / gain + rdvar 
    img = img / divisor / itime

    # Collect header information

    hdr = get_header(hdul[0].header)

    # Close the file

    hdul.close()

    return [img, var, hdr, bitmask]



def make_flat(files, output_name, prefix='flat', suffix='.[ab]',
              extension='.fits*', input_method='index', normalize=True,
              verbose=True, overwrite=True, qafile=False):
    
    """
    To create a normalized SpeX flat field file.

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

    """

    #
    # Load the files into memory
    #

    # Get the user requested keywords

    keywords = config.state['xspextool_keywords']

    # Add GRAT and DIT so that you can determine the mode.  Users will just
    # have to live with them in their files

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

    if verbose is True:
        print(' ')
        print('Creating the flat field...')
        print('Loading FITS images...')

    img, var, hdr, mask = read_fits(files, config.state['linearity_info'],
                                    rotate=7, keywords=keywords,
                                    clupdate=verbose)

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

    # Get the mode name and read modefile

    mode = hdr[0]['GRAT'][0]

    modefile = os.path.join(config.state['instrument_path'],
                            mode + '_flatinfo.fits')

    modeinfo = read_flatcal_file(modefile)

    if verbose is True:
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

        if verbose is True:
            print('Normalizing the median image...')

        nimg, nvar, rms = normalize_flat(med, edgecoeffs, modeinfo['xranges'],
                                         modeinfo['slith_arc'],
                                         modeinfo['nxgrid'],
                                         modeinfo['nygrid'],
                                         var=munc ** 2, oversamp=1,
                                         ybuffer=modeinfo['ybuffer'],
                                         clupdate=False)

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
               os.path.join(config.state['calpath'], output_name + '.fits'),
               overwrite=overwrite)

    if verbose is True:
        print('Flat field '+output_name+'.fits written to disk.')


def make_wavecal(files, flat_file, output_name, prefix='arc', suffix='.[ab]',
                 extension='.fits*', input_method='index',
                 use_stored_solution=False, verbose=True, qafile_shift=True,
                 qafile_findlines=True, qafile_fitlines=True, overwrite=True):

    """
    Creates a Spextool wavecal file.

    Parameters
    ----------
    files : 

    flat_file : str
        The full path to a Spextool flat file.

    output_name : str
        The output name of the wavecal file sans the extension.

    prefix : 

    suffix :

    extension : 

    input_method : 


    """

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

    check_parameter('make_wavecal', 'clupdate', verbose, 'bool')

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

    if verbose is True:
        print(' ')
        print('Creating the wavecal file...')
        print('Loading FITS images...')

    img, var, hdr, mask = read_fits(files, config.state['linearity_info'],
                                    rotate=7, keywords=keywords,
                                    clupdate=verbose)

    #
    # Combine images as necessary
    #

    # Now scale their intensities to a common flux level

    if len(files) > 1:

        if verbose is True:
            print('Scaling images...')

        simgs, svars, scales = math.scale_data_stack(img, None)

    else:
        simgs = img

    # Now median the scaled images

    if len(files) > 1:

        if verbose is True:
            print('Medianing the images...')

        med, munc = math.median_data_stack(simgs)

    else:
        med = simgs

    #
    # Let's do the extraction of the "arc" spectra
    #

    # Read the flat and wavecal files

    flat_file = os.path.join(config.state['calpath'], flat_file)
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
    spectra = extract_extendedsource_1dxd(med, med, flatinfo['ordermask'],
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

        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
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

        #     Determine the guess position and search range for each

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
                          'filepath': config.state['qapath'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the lines

        lineinfo = find_lines_1dxd(spectra['spectra'], wavecalinfo['orders'],
                                   lineinfo, flatinfo['slitw_pix'],
                                   qafileinfo=qafileinfo, clupdate=verbose)

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
                          'filepath': config.state['qapath'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the solution

        solution = wavecal_solution_1d(wavecalinfo['orders'], lineinfo,
                                       wavecalinfo['dispdeg'], xd=xd,
                                       qafileinfo=qafileinfo,
                                       clupdate=verbose)
        
    else:

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
                         os.path.join(config.state['calpath'],
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
                         os.path.join(config.state['calpath'],
                                        output_name + '.fits'),
                         config.state['version'],
                         xd={'orderdeg':wavecalinfo['ordrdeg'],
                             'homeorder':wavecalinfo['homeorder']},
                         stored_solution=use_stored_solution,
                         overwrite=overwrite)

    else:
        print('unknown wcaltype.')

    if verbose:
        print('Wavecal '+output_name+'.fits written to disk.')


        
def read_fits(files, linearity_info, keywords=None, pair_subtract=False,
              rotate=0, linearity_correction=True, clupdate=False):

    """
    To read a SpeX FITS image file.

    Parameters
    ----------
    files : list of str
        A list of fullpaths to FITS files.

    linearity_info : dict {'max':int,'bit':int}
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

    """
    #
    # Get setup information
    #
    
    naxis1 = 1024
    naxis2 = 1024

    nfiles = len(files)


    #
    # Correct for non-linearity?
    #
    
    if linearity_correction is True:

        linearity_file = os.path.join(config.state['instrument_path'],
                                      'spex_lincorr.fits')
        linearity_coeffs = fits.getdata(linearity_file)

    else:

        linearity_coeffs = None

    #
    # Figure out the number of images that will be returned.
    #
    
    if pair_subtract:

        #  Check to make sure the right number of files

        if (nfiles % 2) != 0:

            print('mc_readuspexfits:  Not an even number of images.')
            sys.exit(1)

        else:

            nimages = int(nfiles / 2)

    else:

        nimages = nfiles

    #
    # Setup output arrays and lists
    #

    data = np.empty((nimages, naxis2, naxis1))
    var = np.empty((nimages, naxis2, naxis1))
    bitmask = np.empty((nimages, naxis2, naxis1), dtype=np.uint8)
    hdrinfo = []

    #
    # Now we load the data
    #
    
    if pair_subtract is True:

        # pair subtraction

        for i in range(0, nimages):

            a = load_data(files[i*2], linearity_info, keywords=keywords,
                          linearity_coefficients=linearity_coeffs)

            b = load_data(files[i*2+1], linearity_info, keywords=keywords,
                          linearity_coefficients=linearity_coeffs)

            combmask = math.combine_flag_stack(np.stack((a[3], b[3])),
                                          nbits=linearity_info['bit'] + 1)

            data[i, :, :] = idl_rotate(a[0] - b[0], rotate)
            var[i, :, :] = idl_rotate(a[1] + b[1], rotate)
            bitmask[i, :, :] = idl_rotate(combmask, rotate)

            hdrinfo.append(a[2])
            hdrinfo.append(b[2])

    else:

        for i in range(0, nimages):

            a = load_data(files[i], linearity_info, keywords=keywords,
                          linearity_coefficients=linearity_coeffs)

            data[i, :, :] = idl_rotate(a[0], rotate)
            var[i, :, :] = idl_rotate(a[1], rotate)
            bitmask[i, :, :] = idl_rotate(a[3], rotate)

            hdrinfo.append(a[2])

    return np.squeeze(data), np.squeeze(var), hdrinfo, np.squeeze(bitmask)



def spex_linearity_imagepoly(image, coefficients):

    """
    To apply the linearity coefficients

    There is a special correction required so the standard image_poly 
    routine will not work.  TBD whether this should really be there...

    Parameters
    ----------
    image : ndarray
        A (nrows, ncols) FITS image 

    coefficients : ndarray
        A (nrows, ncols, ncoefficients) array of coefficients.

    """

    #
    # Perform the correction to the image
    #
   
    c0 = coefficients[0,:,:]
    tread    = 0.36
    corrected_image = image - c0*tread

    #
    # Now apply the coefficients
    #

    result = image_poly(corrected_image, coefficients)

    return result


