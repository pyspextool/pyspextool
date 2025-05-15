import numpy as np
import numpy.typing as npt
from astropy.io import fits
import re
import os
import logging

from pyspextool.fit.polyfit import image_poly
from pyspextool.io.check import check_parameter
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils.math import combine_flag_stack
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.setup_utils import mishu # pooch class for caching files

def correct_uspexbias(img:npt.ArrayLike):

    """
    To correct for bias voltage drift in an uSpeX FITS image

    Parameters
    ----------
    img : ndarray
        An (nrows, ncols) uSpeX image

    Returns
    -------
    ndarray
        The uSpeX image with the bias variations "corrected".

    Notes
    -----
    There are 32 amplifiers that talk to 64 columns each.  The median
    intensity of the 64 reference pixels at the bottom of the image are
    subtracted from all rows in the 64 columns.

    """

    #
    # Check parameter
    #

    check_parameter('correct_uspexbias', 'img', img, 'ndarray')

    #
    # Loop over each amp
    #
    
    corrected_image = img.copy()

    for i in range(0, 32):

        xl = 0+64*i
        xr = xl+63

        med = np.median(img[2044:(2047+1), xl:(xr+1)])    
        corrected_image[:, xl:(xr+1)] = img[:, xl:(xr+1)] - med

    return  corrected_image    


def read_fits(files:list,
              linearity_info:dict,
              keywords:list=None,
              pair_subtract:bool=False,
              rotate:int=0,
              linearity_correction:bool=True,
              extra:dict=None,
              verbose:bool=False):

    """
    To read upgraded SpeX FITS image files.

    Parameters
    ----------
    files : list of str
        A list of fullpaths to uSpeX FITS files.

    linearity_info : dict
        `"max"` : int
            The maximum value in DN beyond which a pixel is deemed non-linear.

        `"bit"` : int8
            The bit to set to identify a pixel that is non-linear in the mask
            that is returned.

    keywords : list of str, default None
        A list of FITS keyword to retain from the raw image. 

    pair_subtract : {False, True}
        Set to True to pair subtract the images.
        Set to False to load the images sequentially.    

    rotate : {0,1,2,3,4,5,6,7}
    
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
        
    linearity_correction : {True, False}
        Set to True to correct for non-linearity.
        Set to False to not correct for non-linearity.    

    extra : dict, default None

        `"correct_bias"` : {True, False}
            Set to True to correct the bias drift.
            Set to False to not correct the bias drift.    
        
    verbose : {False, True}
        Set to True to report progress.
        Set to False to not report progress.    
    
    Returns
    -------
    tuple 

        tuple[0] : ndarray
            An (nimages, nrows, ncols) image in units of DN s-1.

        tuple[1] : ndarray
            An (nimages, nrows, ncols) variance in units of (DN s-1)^2.

        tuple[2] : list
            A (nimages, ) list where each element is a dictionary of
            key:value pairs where `key` is the FITS header keyword and
            `value` is a 2-element list where the first element is the
            FITS value and the second element is the FITS comment.

        tuple[3] : ndarray of uint8
            An (nimages, nrows, ncols)  bitset mask where pixels greater than
            linearity_info['max'] have the bit linearity_info['bit'] set.  


    """
    # Get setup information

    naxis1 = 2048
    naxis2 = 2048

    nfiles = len(files)

    # Correct for non-linearity?
    if linearity_correction is True:
        linearity_file = mishu.fetch('uspex_lincorr.fits')
        lc_coeffs = fits.getdata(linearity_file)
    else:
        lc_coeffs = None

    # Get set up for linearity check
    bias_file = mishu.fetch('uspex_bias.fits')   
    hdul = fits.open(bias_file)
    
    divisor = hdul[0].header['DIVISOR']
    bias = hdul[0].data / divisor
    hdul.close()

    # Check for amplifier correction

    if extra is None:

        extra = {'correct_bias':True}

    #
    # Deal with the logging
    #

    if verbose is True:

        list_filenames = []    
        for file in files:

            list_filenames.append(os.path.basename(file))
            
            file_names = ", ".join(list_filenames)

        
        message = ' Loading images(s) '+file_names

        if extra['correct_bias'] is True:
            
            message = message + ' correcting bias drift,'
            
        else :
            
            message = message + ' not correcting bias drift,'
            
        if linearity_correction is True:
        
            message = message + ' and correcting for non-linearity.'
        
        else:
        
            message = message + ' and not correcting for non-linearity.'

        logging.info(message)                
                   
    if pair_subtract is True:

        #  Check to make sure the right number of files

        if (nfiles % 2) != 0:

            message = ' An even number of images is required if `pair`=True.'
            raise pySpextoolError

        else:

            nimages = int(nfiles / 2)

    else:

        nimages = nfiles

    # Make empty arrays

    data = np.empty((nimages, naxis2, naxis1))
    var = np.empty((nimages, naxis2, naxis1))
    bitmask = np.empty((nimages, naxis2, naxis1), dtype=np.int8)
    hdrinfo = []

    # Load the data

    if pair_subtract is True:

        # pair subtraction

        for i in range(0, nimages):

            if verbose is True:

                loop_progress(i, 0, nimages)

            a = load_uspeximage(files[i * 2],
                                linearity_info,
                                bias,
                                keywords=keywords,
                                correct_bias=extra['correct_bias'],
                                linearity_coeffs=lc_coeffs)
            
            b = load_uspeximage(files[i * 2 + 1],
                                linearity_info,
                                bias,
                                keywords=keywords,
                                correct_bias=extra['correct_bias'],
                                linearity_coeffs=lc_coeffs)

            combmask = combine_flag_stack(np.stack((a[3], b[3])),
                                          nbits=linearity_info['bit'] + 1)

            data[i, :, :] = idl_rotate(a[0] - b[0], rotate)
            var[i, :, :] = idl_rotate(a[1] + b[1], rotate)
            bitmask[i, :, :] = idl_rotate(combmask, rotate)

            hdrinfo.append(a[2])
            hdrinfo.append(b[2])

    else:

        for i in range(0, nimages):

            if verbose is True:

                loop_progress(i, 0, nimages)            

            result = load_uspeximage(files[i],
                                     linearity_info,
                                     bias,
                                     keywords=keywords,
                                     correct_bias=extra['correct_bias'],
                                     linearity_coeffs=lc_coeffs)

            data[i, :, :] = idl_rotate(result[0], rotate)
            var[i, :, :] = idl_rotate(result[1], rotate)
            bitmask[i, :, :] = idl_rotate(result[3], rotate)

            hdrinfo.append(result[2])

    return np.squeeze(data), np.squeeze(var), hdrinfo, np.squeeze(bitmask)


def load_uspeximage(file:list,
                    linearity_info:dict,
                    bias:npt.ArrayLike,
                    keywords=None,
                    correct_bias=True,
                    linearity_coeffs=None):

    """
    To load a single uSpeX image into memory.

    Parameters
    ----------
    files : str
        A fullpath to a uSpeX FITS file.

    linearity_info : dict
        `"max"` : int
            The maximum value in DN beyond which a pixel is deemed non-linear.

        `"bit"` : int8
            The bit to set to identify a pixel that is non-linear in the mask
            that is returned.
        
    bias : ndarray
        An (nrows, ncols) uSpeX bias frame.

    keywords : list of str, default None
        A list of FITS keyword to retain from the raw image. 

    correct_bias : {True, False}
        Set to True to correct the bias drift.
        Set to False to not correct the bias drift.    

    linearity_coeffs : ndarray, default None
        An (ncoeffs, nrows, ncols) array of linearity correction coefficients.
    
    Returns
    -------
    list

        list[0] : ndarray
            An (nrows, ncols) image in units of DN s-1.

        list[1] : ndarray
            An (nrows, ncols) variance in units of (DN s-1)^2.

        list[2] : dict
            A dictionary of key:value pairs where `key` is the FITS header
            keyword and `value` is a 2-element list where the first element
            is the FITS value and the second element is the FITS comment.

        list[3] : ndarray of uint8
            An (nrows, ncols)  bitset mask where pixels greater than
            linearity_info['max'] have the bit linearity_info['bit'] set.  
    
    """

    #
    # Check parameters
    #

    check_parameter('load_uspeximage', 'file', file, 'str')

    check_parameter('load_uspeximage', 'linearity_info', linearity_info,
                    'dict')    

    check_parameter('load_uspeximage', 'bias', bias, 'ndarray')

    check_parameter('load_uspeximage', 'keywords', keywords,
                    ['list','NoneType'])        

    check_parameter('load_uspeximage', 'correct_bias', correct_bias, 'bool')

    check_parameter('load_uspeximage', 'correct_bias', correct_bias, 'bool')

    check_parameter('load_uspeximage', 'linearity_coeffs', linearity_coeffs,
                    ['ndarray', 'NoneType'])        
    
            
    readnoise = 12.0  # per single read
    gain = 1.5  # electrons per DN

    hdul = fits.open(file, ignore_missing_data=True)
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
            
    mskp = (img_p < (bias - linearity_info['max'])) * 2 ** linearity_info['bit']
    msks = (img_s < (bias - linearity_info['max'])) * 2 ** linearity_info['bit']

    #  Combine the masks 

    stack = np.stack((mskp.astype(np.uint8), msks.astype(np.uint8)))

    bitmask = combine_flag_stack(stack,
                                 nbits=linearity_info['bit'] + 1)

    #  Create the image

    img = img_p - img_s
    
    #  Correct for amplifier offsets

    if correct_bias:

        img = correct_uspexbias(img)
        
    #  Determine the linearity correction for the image

    if linearity_coeffs is not None:

        linearity_correction = True
        
        cor = image_poly(img, linearity_coeffs)
        cor = np.where(cor == 0, 1, cor)

        #  Now set the corrections to unity for pixels > lincormax

        cor = np.where(bitmask == 2 ** linearity_info['bit'], 1, cor)

        #  Set black pixel corrections to unity as well.

        cor[:, 0:3 + 1] = 1.0
        cor[:, 2044:2047 + 1] = 1.0
        cor[0:3 + 1, :] = 1.0
        cor[2044:2047 + 1, :] = 1.0

        # Apply the corrections

        img /= cor

        # Delete unecessary files

        del cor, img_p, img_s

    else:

        linearity_correction = False
        
    # Create the actual image.
    # Convert image back to total DN for error propagation
    #

    img = img * divisor

    
    # Compute the variance and the final image

    var = np.absolute(img) * crtn / ndrs / (coadds ** 2) / \
          (itime ** 2) / gain + rdvar 
    img = img / divisor / itime
    
    # Collect header information

    hdr = get_uspexheader(hdul[0].header,
                          keywords=keywords)

    # Add on the linearity correction and bias drift

    hdr['LINMAX'] = (linearity_info['max'], ' Maximum linear value (DN)')
    hdr['LINCOR'] = (linearity_correction, ' Linearity correction?')
    hdr['BIASCOR'] = (correct_bias, ' Bias drift correction?')        

       
    hdul.close()

    return [img, var, hdr, bitmask]



def get_uspexheader(hdr,
                    keywords=None):

    """
    To obtain header values from the FITS header

    Parameters:
    ----------
    hdr : Header
        An astropy HDU object with an attribute of .header 

    keywords : list of str, deafult None
        A list of string keywords to store.   

    Returns
    -------
    dict

        A dictionary with the following keys along with the keys for user
        passed `keywords`.

        `"AM"` : float
            The airmass of the observation.

        `"HA"` : str
            The hour angle of the observations (+-hh:mm:ss.ss).

        `"PA"` : float
            The position angle of the observation (degrees)

        `"DEC"` : str
            The declination of the observation (+-dd:mm:ss.s)
    
        `"RA"` : str
            The right ascension of the observation (hh:mm:ss.s)    

        `"ITIME"` : float
            The exposure time of the observation (sec)

        `"NCOADDS"` : int
            The number of coadds of the observations.

        `"IMGITIME"` : float
            The total exposure time of the observations, ITIME*NCOADDS (sec)

        `"TIME"` : str
            The UTC time of the observations (hh:mm:ss.ssssss)
    
        `"DATE"` : str
            The UTC date of the observation (yyy-mm-dd)
    
        `"MJD"` : float                               
            The modified Julian date of the observation 
    
        `"FILENAME"` : str
            The filename.  

        `"MODE"` : str
            The instrument mode.

        `"INSTR"` : str
            The instrument.
    
    """

    #
    # Check parameters
    #

    check_parameter('get_uspexheader', 'hdr', hdr, 'Header')

    check_parameter('get_uspexheader', 'keywords', keywords,
                    ['list','NoneType'])    

    #
    # Grab specific keywords if requested.  Otherwise, grab all of them
    #
    
    hdrinfo = get_headerinfo(hdr,keywords=keywords)

    #
    #  Grab require keywords and convert to standard Spextool keywords
    #

    
    # Airmass 

    try:

        airmass = hdr['TCS_AM']

    except KeyError as e:

        airmass = np.nan

    hdrinfo['AM'] = [airmass, ' Airmass']

    # Hour angle

    try:

        hourangle = hdr['TCS_HA']

        m = re.search('[-]', '[' + hourangle + ']')
        if not m:
            hourangle = '+' + hourangle.strip()
            
    except KeyError as e:

        hourangle = 'nan'
        
    hdrinfo['HA'] = [hourangle, ' Hour angle (hours)']

    # Position Angle

    hdrinfo['PA'] = [hdr['POSANGLE'], ' Position Angle E of N (deg)']


    # Declination

    try:

        declination = hdr['TCS_DEC']

        m = re.search('[-]', '[' + declination + ']')
        if not m:
            declination = '+' + declination.strip()
            

    except KeyError as e:

        declination = 'nan'

    hdrinfo['DEC'] = [declination, ' Declination, FK5 J2000']

    # Right Ascension

    try:

        ra = hdr['TCS_RA'].strip()

    except KeyError as e:

        declination = 'nan'

    hdrinfo['RA'] = [declination, ' Right Ascension, FK5 J2000']

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

    # Now grab any user COMMENT

#    print(hdr['COMMENT'])
#    comment = str(hdr['COMMENT'])
#    comment = comment.split('=')[1]    
#    hdrinfo['USERCOM'] = [comment[2:-1], ' User comment']

    return hdrinfo
