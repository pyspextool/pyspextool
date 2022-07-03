import os.path
import numpy as np
from astropy.io.fits import getheader
from check_file import check_file
from check_dir import check_dir
from readinstrfile import readinstrfile
from readflatcalfile import readflatcalfile
from mkfullpath import mkfullpath
from readuspexfits import readuspexfits
from scalestack import scalestack
from medstack import medstack
from findorders import findorders
from normspecflat import normspecflat
from avehdrs import avehdrs
from combflagstack import combflagstack
from writeflat import writeflat
from splittext import splittext
from wherelist import wherelist

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

    med,munc = medstack(simgs)

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
    
