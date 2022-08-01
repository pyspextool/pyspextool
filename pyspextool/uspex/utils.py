"""Functions for SpeX."""


import os.path
import re, sys
import numpy as np
from numpy import median as median
from astropy.io.fits import getheader
from numpy import empty as npempty
from numpy import int8 as npint8
from numpy import stack as npstack
from numpy import squeeze as npsqueeze
from numpy import where as npwhere
from numpy import absolute as npabsolute
from numpy import shape as npshape
from astropy.io import fits

from pyspextool.utils.ioutils import check_file, check_dir, mkfullpath
from pyspextool.utils.fitsutils import medcomb, avehdrs, gethdrinfo, idlrotate
from pyspextool.utils.image import scalestack
from pyspextool.utils.spec import findorders, combflagstack
from pyspextool.utils.text import splittext, wherelist
from pyspextool.io.fits import readinstrfile, readflatcalfile, writeflat
from pyspextool.calibration.cal import normspecflat
from pyspextool.fit.image import imgpoly


def make_uspex_flat(files,instrfilepath,modefilepath,biasfilepath,oname,\
                    rawpath=None,calpath=None,prefix='flat-',\
                    suffix='.[ab].fits*',input_method='index',\
                    normalize=True,qafilename=None,clupdate=True,\
                    overwrite=True):

    '''
    To create a normalized uSpeX flat field file.


    Input Parameters
    ----------------
    files : str or list of str
        

    instrfilepath : str
        The directory to the instrument file.  Will be removed later.


    modefilepath : str
        The directory to the mode file.  Will be removed later.


    biasfilepath : str
        The directory to the bias file.  Will be removed later.


    oname : str
        The filename of the flat field image to written to disk.      


    rawpath : str, optional
        The path to the raw data.


    calpath : str, optional
        The path to the directory to write the flat.


    prefix : str, default='flat-', optional
        The prefix of the FITS file name (i.e. the stuff before the 
        file numbers.)


    suffix : str, default='.[ab].fits*', optional
        The prefix of the FITS file name (i.e. the stuff after the 
        file numbers.)


    input_method : {'index', 'filename'}, optional
        `index`: `files` is a str of file numbers, e.g., '1-3,5-10'.

        `filename`: a str or list of strings with full file names, e.g.
            ['flat-00001.a.fits',''flat-00002.a.fits']


    normalize : {True, False}, optional
        Set to True to normalize the orders.


    qafilename : str, optional
        The full path and filename for the finderorders quality 
        assurance plot.


    clupdate : {True, False}, optional
        Set to True for command line updates during execution. 
        

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.


    Returns
    -------
    None
        Writes a FITS file to disk.


    Notes
    -----
    None


    Examples
    --------
    later


    Modification History
    --------------------
    2022-06-28 - Written by M. Cushing, University of Toledo.

    '''
    
    # Construct file names 
    
    instrfile = os.path.join(instrfilepath,'uspex.dat')

    biasfile = os.path.join(biasfilepath,'uSpeX_bias.fits')
        
    # Check that things things exist

    test = check_dir([rawpath,calpath])

    test = check_file([instrfile,biasfile])

    # Read instrument file.  Grab user requested keywords.

    instrinfo = readinstrfile(instrfile)
    keywords = instrinfo['XSPEXTOOL_KEYWORDS']

    # Add GRAT and DIT so that you can determine the mode.  Users will just have
    # to live with them in their files

    z = wherelist(keywords,'GRAT')
    if bool(z) is False: keywords.append('GRAT')

    z = wherelist(keywords,'DIT')
    if bool(z) is False: keywords.append('DIT')        

    # Now create the file names
    
    if input_method == 'index':

        files = mkfullpath(rawpath,files,\
                           indexinfo={'nint':instrinfo['NINT'],\
                                      'prefix':prefix,'suffix':suffix},\
                           exist=True)

    elif input_method == 'filename':

        files = mkfullpath(rawpath,files,exist=True)
        
    else:

        raise ValueError('Unknown input_method')

    # Load the FITS files into memory

    if clupdate is True: print('Loading FITS images...')
        
    lininfo = {'bias':biasfile,'max':instrinfo['LINCORMAX'],'bit':0}
    img,var,hdr,mask = readuspexfits(files,lininfo,keywords=keywords,\
                                    clupdate=clupdate)

    # Average the headers

    avehdr = avehdrs(hdr)
    
    # Combine the masks

    flag = combflagstack(mask)
    
    # Now scale their intensities to a common flux level

    if clupdate is True: print('Scaling images...')
    
    simgs,svars,scales = scalestack(img)

    # Now median the scaled images

    if clupdate is True: print('Medianing the images...')

    med,munc = medcomb(simgs)

    # Get the mode name and read modefile

    mode = hdr[0]['GRAT'][0]

    modefile = modefilepath+mode+'_flatinfo.fits'

    test = check_file(modefile)    

    modeinfo = readflatcalfile(modefile)

    # Locate the orders
    
    if clupdate is True: print('Locating the orders...')
    
    edgecoeffs = findorders(med,modeinfo['guesspos'],modeinfo['xranges'],\
                            modeinfo['step'],modeinfo['slith_range'],\
                            modeinfo['edgedeg'],modeinfo['ybuffer'],\
                            modeinfo['flatfrac'],modeinfo['comwidth'],\
                            qafig=qafilename)

    # Normalize the spectrum if requested

    if normalize is True:
    
        if clupdate is True: print('Normalizing the median image...')

        nimg, nvar, rms = normspecflat(med,edgecoeffs,\
                                        modeinfo['xranges'],\
                                        modeinfo['slith_arc'],\
                                        modeinfo['nxgrid'],\
                                        modeinfo['nygrid'],\
                                        var=munc**2,oversamp=1,\
                                        ybuffer=modeinfo['ybuffer'],\
                                        clupdate=False)

    else:

        nimg = med
        nvar = munc**2
        rms = np.full((len(modeinfo['orders'])),np.nan)

    # Create the HISTORY

    basenames = []
    for file in files:

        basenames.append(os.path.basename(file))

    history = 'This flat was created by scaling the files '+\
      ', '.join(str(b) for b in basenames)+' to a common median flux '+\
      'level and then medianing the scaled imges.  The variance is '+\
      'given by (1.4826*MAD)**2/nimages where MAD is the median absolute '+\
      'deviation.  The zeroth bit of pixels in the third extension are '+\
      'set if their corresponding intensity values are greater than '+\
      'LINCORMAX.  User selected FITS keywords are from the first frame '+\
      'in the series.'

    history = splittext(history)
    
        
    # Get the slit widths and resolving power and write to disk

    slitw_arc = float(avehdr['SLIT'][0][0:3])
    slitw_pix = slitw_arc/modeinfo['ps']

    resolvingpower = modeinfo['rpppix']/slitw_pix

    writeflat(nimg,nvar,flag,avehdr,modeinfo['rotation'],modeinfo['orders'],\
                edgecoeffs,modeinfo['xranges'],modeinfo['ps'],\
                modeinfo['slith_pix'],modeinfo['slith_arc'],\
                slitw_pix,slitw_arc,mode,rms,resolvingpower,'1.0beta',\
                history,os.path.join(calpath,oname),overwrite=overwrite) 
    

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

    # Get setup information

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

    # Load the data
     
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


def uspexampcor(img):

    '''
    To correct for bias voltage drift in an uSpeX FITS image

    Input Parameters
    ----------------
    img : numpy array
        An uSpeX image

    Returns
    --------
    numpy.ndarray
        The uSpeX image with the bias variations "corrected".

    Procedure
    ---------
    There are 32 amplifiers that talk to 64 columns each.  The median 
    intensity of the 64 reference pixels at the bottom of image are 
    subtracted from all rows in the 64 columns.

    Example
    --------
    result = uspexampcor(img) (how do we do this properly?)


    Modification History
    --------------------
    2022-05-25 - Written by M. Cushing, University of Toledo.
                Based on the Spextool mc_uspexampcor.pro IDL program.
    '''

    for i in range(0,32):

        xl = 0+64*i
        xr = xl+63

        med = median(img[2044:2047+1,xl:xr+1])
        img[:,xl:(xr+1)]-=med

    return(img)