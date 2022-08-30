import numpy as np
from astropy.io import fits
import re
from pyspextool.utils.arrays import idl_rotate
from pyspextool.fit.polyfit import image_poly
from pyspextool.utils.math import combine_flag_stack
from pyspextool.io.fitsheader import get_header_info
#import re
#import sys


def read_uspex_fits(files, lininfo, keywords=None, pair=False, rotate=0,
                    lincor=None, ampcor=False, clupdate=False):
    """
    To read an upgraded SpeX FITS image file.

    Parameters
    ----------
    files : list of str
        A list of fullpaths to FITS files.

    lininfo : dict {'bias':str,'max':int,'bit':int}
        information to identify pixels beyond range of linearity correction

        'bias' is the fullpath to the bias frame
        'max' maximum value in DN
        'bit' the bit to set for pixels beyond `max`

    keywords : list of str, optional
        A list of FITS keyword to retain 

    pair : {False, True}, optional
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
        
    lincor : str, optional
        the fullpath to the FITS file of linearity correction coefficients

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

    Procedure
    ---------
    ?


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

    dolincor = [0, 1][lincor is not None]

    # Correct for non-linearity?

    if dolincor:

        lc_coeffs = fits.getdata(lincor)

    else:

        lc_coeffs = None

    # Get set up for lineary check

    hdul = fits.open(lininfo['bias'])
    divisor = hdul[0].header['DIVISOR']
    bias = hdul[0].data / divisor
    hdul.close()

    if pair:

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

    if pair is True:

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

    if not pair:

        for i in range(0, nimages):
            im, va, hd, bm = load_data(files[i], lininfo, bias, keywords=keywords,
                                       ampcor=ampcor, lccoeffs=lc_coeffs)

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

    mskp = (img_p < (bias - lininfo['max'])) * 2 ** lininfo['bit']
    msks = (img_s < (bias - lininfo['max'])) * 2 ** lininfo['bit']

    #  Combine the masks 

    bitmask = combine_flag_stack(np.stack((mskp, msks)), nbits=lininfo['bit'] + 1)

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

    var = np.absolute(img) * crtn / ndrs / (coadds ** 2) / (itime ** 2) / gain + rdvar
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
