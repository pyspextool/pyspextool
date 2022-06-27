from numpy import empty as npempty
from numpy import int8 as npint8
from numpy import stack as npstack
from numpy import squeeze as npsqueeze
from numpy import where as npwhere
from numpy import absolute as npabsolute
from numpy import shape as npshape
from astropy.io import fits
import re
import sys
from gethdrinfo import gethdrinfo
from combflagstack import combflagstack
from imgpoly import imgpoly
from uspexampcor import uspexampcor
from idlrotate import idlrotate


#
#=============================================================================
#
def readuspexfits(files,lininfo,keywords=None,pair=False,rotate=0,\
                  lincor=None,ampcor=False,clupdate=False):

    """
    To read an (upgraded) SpeX FITS image file.

    Parameters
    ----------------
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

    ampor : {False, True}, optional 
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

#
# Get setup information
#

    NAXIS1=2048
    NAXIS2=2048
    
    nfiles = len(files)
    
    dolincor = [0,1][lincor is not None]

# Correct for non-linearity?

    if dolincor:

        lc_coeffs = fits.getdata(lincor)
                
    else:

        lc_coeffs = None
        
# Get set up for lineary check

    hdul = fits.open(lininfo['bias'])            
    DIVISOR = hdul[0].header['DIVISOR']
    bias = (hdul[0].data)/DIVISOR
    hdul.close()

    if pair:

#  Check to make sure the right number of files

       if (nfiles % 2) != 0:

           print('mc_readuspexfits:  Not an even number of images.')
           sys.exit(1)
           
       else:

           nimages = int(nfiles/2)

    else:

        nimages = nfiles

# Make empty arrays

    data    = npempty((nimages,NAXIS2,NAXIS1))
    var     = npempty((nimages,NAXIS2,NAXIS1))
    hdrinfo = []
    bitmask = npempty((nimages,NAXIS2,NAXIS1),dtype=npint8)

#
# Load the data
#        
    if pair is True:

# pair subtraction

        for i in range(0,nimages):

            A = loaddata(files[i*2],lininfo,bias,\
                         keywords=keywords,ampcor=ampcor,lccoeffs=lc_coeffs)

            B = loaddata(files[i*2+1],lininfo,bias,\
                         keywords=keywords,ampcor=ampcor,lccoeffs=lc_coeffs)

            combmask=combflagstack(npstack((A[3],B[3])),nbits=lininfo['bit']+1)
            
            data[i,:,:]    = idlrotate(A[0]-B[0],rotate)
            var[i,:,:]     = idlrotate(A[1]+B[1],rotate)
            bitmask[i,:,:] = idlrotate(combmask,rotate)

            hdrinfo.append(A[2])
            hdrinfo.append(B[2])
            
    if not pair:

        for i in range(0,nimages):

            im,va,hd,bm = loaddata(files[i],lininfo,bias,keywords=keywords,\
                                   ampcor=ampcor,lccoeffs=lc_coeffs)

            data[i,:,:]    = idlrotate(im,rotate)
            var[i,:,:]     = idlrotate(va,rotate)
            bitmask[i,:,:] = idlrotate(bm,rotate)

            hdrinfo.append(hd)

    return(npsqueeze(data),npsqueeze(var),hdrinfo,npsqueeze(bitmask))

#
#=============================================================================
#
def loaddata(file,lininfo,bias,keywords=None,ampcor=None,lccoeffs=None):
    

    readnoise = 12.0  #  per single read
    gain = 1.5        #  electrons per DN

    hdul = fits.open(file)
    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems
    
    ITIME    = hdul[0].header['ITIME']
    COADDS   = hdul[0].header['CO_ADDS']
    NDRS     = hdul[0].header['NDR']
    READTIME = hdul[0].header['TABLE_SE']
    DIVISOR  = hdul[0].header['DIVISOR']
    
#  Get set up for error propagation and store total exposure time

    rdvar   = (2.*readnoise**2)/NDRS/COADDS/ITIME**2/gain**2
    crtn    = (1.0 - READTIME*(NDRS**2 -1.0)/3./ITIME/NDRS)

#  Read images, get into units of DN.
    
    img_P = (hdul[1].data)/DIVISOR
    img_S = (hdul[2].data)/DIVISOR
            
#  Check for linearity maximum

    mskP = (img_P < (bias-lininfo['max']))*2**lininfo['bit']
    mskS = (img_S < (bias-lininfo['max']))*2**lininfo['bit']              

#  Combine the masks 
                
    bitmask=combflagstack(npstack((mskP,mskS)),nbits=lininfo['bit']+1)
        
#  Create the image

    img = img_P-img_S
        
#  Correct for amplifier offsets

    if ampcor:

        img = uspexampcor(img)

#  Determine the linearity correction for the image

    if lccoeffs is not None:

        cor = imgpoly(img,lccoeffs)
        cor = npwhere(cor == 0,1,cor)        
                
#  Now set the corrections to unity for pixels > lincormax

        cor = npwhere(bitmask == 2**lininfo['bit'],1,cor)
                
#  Set black pixel corrections to unity as well.

        cor[:,0:3+1] = 1.0
        cor[:,2044:2047+1] = 1.0
        cor[0:3+1,:] = 1.0
        cor[2044:2047+1,:] = 1.0                

# Apply the corrections

        img/=cor

# Delete unecessary files

        del cor,img_P,img_S

# Create the actual image.
# Convert image back to total DN for error propagation

    img = img*DIVISOR

# Compute the variance and the final image

    var=npabsolute(img)*crtn/NDRS/(COADDS**2)/(ITIME**2)/gain + rdvar
    img = img/DIVISOR/ITIME
    
# Collect header information

    hdr = gethdr(hdul[0].header)
                
    hdul.close()

    return[img,var,hdr,bitmask]

#
#=============================================================================
#
def gethdr(hdr,keywords=None):

    
# Grab keywords if requested    

    if keywords:

        hdrinfo = gethdrinfo(hdr,keywords=keywords)
    
    else:

        hdrinfo = gethdrinfo(hdr)        

#  Grab require keywords and convert to standard Spextool keywords

# Airmass 
        
    hdrinfo['AM'] = [hdr['TCS_AM'],' Airmass']

# Hour angle
    
    val = hdr['TCS_HA']
    m = re.search('[-]','['+val+']')
    if not m: val = '+'+val.strip()
    hdrinfo['HA'] = [val,' Hour angle (hours)']

# Position Angle
    
    hdrinfo['PA'] = [hdr['POSANGLE'],' Position Angle E of N (deg)']

# Dec 
    
    val = hdr['TCS_DEC']
    m = re.search('[-]','['+val+']')
    if not m: val = '+'+val.strip()
    hdrinfo['DEC'] = [val,' Declination, FK5 J2000']

# RA
    
    hdrinfo['RA'] = [hdr['TCS_RA'].strip(),' Right Ascension, FK5 J2000']

# COADDS, ITIME
    
    coadds = hdr['CO_ADDS']


    itime = hdr['ITIME']
    hdrinfo['ITIME'] = [itime,' Integration time (sec)']
    hdrinfo['NCOADDS'] = [coadds,' Number of COADDS']    
    hdrinfo['IMGITIME'] = [coadds*itime,\
                           ' Image integration time, NCOADDSxITIME (sec)']

# Time

    
    hdrinfo['TIME'] = [hdr['TIME_OBS'].strip(),' Observation time in UTC']

# Date
    
    hdrinfo['DATE'] = [hdr['DATE_OBS'].strip(),' Observation date in UTC']        

# MJD
    
    hdrinfo['MJD'] = [hdr['MJD_OBS'],' Modified Julian date OBSDATE+TIME_OBS']

# FILENAME
    
    hdrinfo['FILENAME'] = [hdr['IRAFNAME'].strip(),' Filename']

# MODE
    
    hdrinfo['MODE'] = [hdr['GRAT'].strip(),' Instrument Mode']

# INSTRUMENT
    
    hdrinfo['INSTR'] = ['SpeX',' Instrument']

# now move the comment key to the end
    
    comment = hdrinfo.pop('COMMENT')
    hdrinfo['COMMENT'] = comment
    
    return(hdrinfo)
