import os
import numpy as np
from pyspextool.cl import config
from pyspextool.io.files import check_file
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.files import make_full_path
from pyspextool.io.read_uspex_fits import read_uspex_fits
from pyspextool.utils.math import scale_data_stack
from pyspextool.utils.math import median_data_stack
from pyspextool.io.flat import read_flatcal_file
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.wavecal import read_wavecal_file
from pyspextool.calibration.simulate_wavecal_1dxd import simulate_wavecal_1dxd 
from pyspextool.spectroscopy.get_spectral_pixelshift import get_spectral_pixelshift
from pyspextool.spectroscopy.extract_extendedsource_1dxd import extract_extendedsource_1dxd 

def make_uspex_wavecal(files,flatfile,oname,prefix='arc-',\
                    suffix='.[ab].fits*',input_method='index',\
                    normalize=True,qafilename=None,clupdate=True,\
                    overwrite=True):

# Construct file names 

    instrfile = os.path.join(config.state['packagepath'],'instruments',
                             config.state['instrument'],'data',
                             config.state['instrument']+'.dat')
    
    biasfile = os.path.join(config.state['packagepath'],'instruments',
                             config.state['instrument'],'data',
                             config.state['instrument']+'_bias.fits')

# Check that things things exist

    test = check_file([instrfile,biasfile])

# Read instrument file.  Grab user requested keywords.

    instrinfo = read_instrument_file(instrfile)
    keywords = instrinfo['XSPEXTOOL_KEYWORDS']

# Add GRAT and DIT so that you can determine the mode.  Users will just have
# to live with them in their files

    if 'GRAT' not in keywords:
        keywords.append('GRAT')

    if 'DIT' not in keywords:
        keywords.append('DIT')            

# Now create the file names
    
    if input_method == 'index':

        files = make_full_path(config.state['rawpath'],files,\
                           indexinfo={'nint':instrinfo['NINT'],\
                                      'prefix':prefix,'suffix':suffix},\
                           exist=True)

    elif input_method == 'filename':

        files = make_full_path(rawpath,files,exist=True)
        
    else:

        print('make_flat: Unknown input_method.')
        sys.exit(1)

# Load the FITS files into memory

    if clupdate is True: print('Loading FITS images...')
        
    lininfo = {'bias':biasfile,'max':instrinfo['LINCORMAX'],'bit':0}
    img,var,hdr,mask = read_uspex_fits(files,lininfo,keywords=keywords,\
                                        clupdate=clupdate)

# Now scale their intensities to a common flux level

    if clupdate is True and len(files) > 1:

        print('Scaling images...')
        simgs,svars,scales = scale_data_stack(img,None)

    else: simgs = img
        
# Now median the scaled images

    if clupdate is True and len(files) >1:

        print('Medianing the images...')

        med,munc = median_data_stack(simgs)

    else: med = simgs

# Read the flat and wavecal files
        
    flatinfo = read_flat_fits(os.path.join(config.state['calpath'],flatfile))

    wavecalfile = os.path.join(config.state['packagepath'],'instruments',
                               config.state['instrument'],'data',
                            flatinfo['mode']+'_wavecalinfo.fits')

    wavecalinfo = read_wavecal_file(wavecalfile)    

    
# Create wavecal and spatcal images
    
    wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                             flatinfo['nrows'],
                                             flatinfo['edgecoeffs'],
                                             flatinfo['xranges'],
                                             flatinfo['slith_arc'])


    
# Extract the "arc"                                

    spec = extract_extendedsource_1dxd(med, var, flatinfo['ordermask'],
                                       flatinfo['orders'], wavecal, spatcal,
                                       flatinfo['slith_arc']/2,
                                       wavecalinfo['apradius'],
                                       linmax_bitmask=None, badpixel_mask=None,
                                       bginfo=None,clupdate=True)

##### Do the cross correlation    
        
# Get the anchor order

    z = flatinfo['orders'] == wavecalinfo['xcororder']


    xanchor = np.arange(int(wavecalinfo['xranges'][z,0]),\
                        int(wavecalinfo['xranges'][z,1]+1),dtype=int)
    fanchor = np.squeeze(wavecalinfo['spectra'][z,1,:])

# Get the source order

    key = 'OR'+str(int(flatinfo['orders'][z])).zfill(3)+'_AP01'
    xsource = np.squeeze(spec[key][0,:])
    fsource = np.squeeze(spec[key][1,:])

    qafileinfo = {'figsize':(8.5,11), 'filepath':config.state['qapath'],
                  'filename':oname+'_findpixelshift','extension':'.pdf'}
        
    offset = get_spectral_pixelshift(xanchor, fanchor, xsource, fsource,
                                     qafileinfo=qafileinfo)

    
    





# Get the mode name and read modefile

#    mode = hdr[0]['GRAT'][0]
#
#    modefile = modefilepath+mode+'_flatinfo.fits'
#
#    test = check_file(modefile)    
#
#    modeinfo = readmodefile(modefile)
#
## Locate the orders
#    
#    if clupdate is True: print('Locating the orders...')
#    
#    edgecoeffs = findorders(med,modeinfo['guesspos'],modeinfo['xranges'],\
#                            modeinfo['step'],modeinfo['slith_range'],\
#                            modeinfo['edgedeg'],modeinfo['ybuffer'],\
#                            modeinfo['flatfrac'],modeinfo['comwidth'],\
#                            qafig=qafilename)
#
## Normalize the spectrum if requested
#
#    if normalize is True:
#    
#        if clupdate is True: print('Normalizing the median image...')
#
#        nimg, nvar, rms = normspecflat(med,edgecoeffs,\
#                                        modeinfo['xranges'],\
#                                        modeinfo['slith_arc'],\
#                                        modeinfo['nxgrid'],\
#                                        modeinfo['nygrid'],\
#                                        var=munc**2,oversamp=1,\
#                                        ybuffer=modeinfo['ybuffer'],\
#                                        clupdate=False)
#
#    else:
#
#        nimg = med
#        nvar = munc**2
#        rms = np.full((len(modeinfo['orders'])),np.nan)
#
## Create the HISTORY
#
#    basenames = []
#    for file in files:
#
#        basenames.append(os.path.basename(file))
#
#    history = 'This flat was created by scaling the files '+\
#      ', '.join(str(b) for b in basenames)+' to a common median flux '+\
#      'level and then medianing the scaled imges.  The variance is '+\
#      'given by (1.4826*MAD)**2/nimages where MAD is the median absolute '+\
#      'deviation.  The zeroth bit of pixels in the third extension are '+\
#      'set if their corresponding intensity values are greater than '+\
#      'LINCORMAX.  User selected FITS keywords are from the first frame '+\
#      'in the series.'
#
#    history = splittext(history)
#    
#        
## Get the slit widths and resolving power and write to disk
#
#    slitw_arc = float(avehdr['SLIT'][0][0:3])
#    slitw_pix = slitw_arc/modeinfo['ps']
#
#    resolvingpower = modeinfo['rpppix']/slitw_pix
#
#    writeflat(nimg,nvar,flag,avehdr,modeinfo['rotation'],modeinfo['orders'],\
#                edgecoeffs,modeinfo['xranges'],modeinfo['ps'],\
#                modeinfo['slith_pix'],modeinfo['slith_arc'],\
#                slitw_pix,slitw_arc,mode,rms,resolvingpower,'1.0beta',\
#                history,oname,overwrite=overwrite) 
