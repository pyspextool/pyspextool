import numpy as np
import numpy.typing as npt
from astropy.io import fits
from astropy.time import Time
import re
import sys

from pyspextool.fit.polyfit import image_poly
from pyspextool.io.check import check_parameter
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils import math
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.setup_utils import mishu


def correct_linearity(image:npt.ArrayLike,
                      itime:int | float,
                      slowcnts:int | float,
                      ndr:int,
                      coefficients:npt.ArrayLike):

    """
    Correct for non-linearity in SpeX images.

    Parameters
    ----------
    image : ndarray
        (nrows, ncols) image in DN (raw/DIVISOR).

    itime : int or float
        The integration time (sec)

    slowcnts : int or float
        The slowcnt value

    ndr : int
        The number of none destructive reads

    coefficients : ndarray
        An (nrows, ncols, ncoefficients) array of linearity coefficients.

    Returns
    -------
    An (nrows, ncols) linearity corrected image.

    Notes
    -----
    This implements the non-linearity correction method described in 
    Vacca et al. (2004, Vacca, PASP, 116, 352).  Since the average pedestal 
    and signal images are not stored, they must be estimated as described 
    in S2.3.

    """

    #
    # Check parameters
    #
    
    check_parameter('correct_linearity', 'image', image, 'ndarray')

    check_parameter('correct_linearity', 'itime', itime, ['int','float'])

    check_parameter('correct_linearity', 'slowcnts', slowcnts, ['int','float'])

    check_parameter('correct_linearity', 'ndr', ndr, 'int')

    check_parameter('correct_linearity', 'coeffiients', coefficients, 'ndarray')
            
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
    z_nan = np.isnan(correction)
    correction[z_nan] = 1.0

    corrected_mean_pedestal_1 = mean_pedestal_1*correction

    correction = c0/spex_linearity_imagepoly(mean_signal_1,coefficients)
    correction = np.clip(correction, 1.0, None)
    z_nan = np.isnan(correction)
    correction[z_nan] = 1.0

    corrected_mean_signal_1 = mean_signal_1*correction    

    # Create a new estimate of the image

    image_2 = corrected_mean_signal_1 - corrected_mean_pedestal_1 

    #
    # Do the second itteration
    #
    
    # Estimate the mean pedistal and signal images

    mean_pedestal_2 = image_2*np.squeeze(tread[z])*(ndr+1)/2/itime
    mean_signal_2 = image + mean_pedestal_2

    # Do the corrections
    
    correction = c0/spex_linearity_imagepoly(mean_pedestal_2,coefficients)
    correction = np.clip(correction, 1.0, None)
    z_nan = np.isnan(correction)
    correction[z_nan] = 1.0

    corrected_mean_pedestal_2 = mean_pedestal_2*correction

    correction = c0/spex_linearity_imagepoly(mean_signal_2,coefficients)
    correction = np.clip(correction, 1.0, None)
    z_nan = np.isnan(correction)
    correction[z_nan] = 1.0

    corrected_mean_signal_2 = mean_signal_2*correction    

    # Create a new estimate of the image

    image_3 = corrected_mean_signal_2 - corrected_mean_pedestal_2

    # replace NaNs with zeros

    z = np.isnan(image_3)
    
    if np.sum(z) != 0:

        image_3[z] = 0
        
    #    return image

    return image_3


        
def get_header(header,
               keywords=None):

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

        hdrinfo = get_headerinfo(header, keywords=keywords)

    else:

        hdrinfo = get_headerinfo(header)
    
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

    # RA

    hdrinfo['RA'] = [header['RA'].strip(), ' Right Ascension']
    
    # Dec 

    val = header['DEC']
    m = re.search('[-]', '[' + val + ']')
    if not m:
        val = '+' + val.strip()
    hdrinfo['DEC'] = [val, ' Declination']

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

    # Now grab any user COMMENT

    comment = str(header['COMMENT'])
    comment = comment.split('=')[1]    
    hdrinfo['USERCOM'] = [comment[2:-1], ' User comment']

    return hdrinfo



def load_data(file,
              linearity_info,
              keywords=None,
              linearity_coefficients=None):

    
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
    
    hdul = fits.open(file, ignore_missing_end=True)
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

    hdr = get_header(hdul[0].header, keywords=keywords)

    # Close the file

    hdul.close()

    return [img, var, hdr, bitmask]



def read_fits(files,
              linearity_info,
              keywords=None,
              pair_subtract=False,
              rotate=0,
              linearity_correction=True,
              verbose=False,
              extra=None):

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
        linearity_file = mishu.fetch("spex_lincorr.fits")
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

            if verbose is True:
                loop_progress(i, 0, nimages)
                
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

            if verbose is True:
                loop_progress(i, 0, nimages)
            
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
