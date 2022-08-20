import os
import numpy as np
import pickle
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
from pyspextool.io.wavecal import read_line_list
from pyspextool.calibration.get_line_guess_position import get_line_guess_position
from pyspextool.utils.for_print import for_print
from pyspextool.calibration.find_lines_1dxd import find_lines_1dxd

#from pyspextool.calibration.find_lines_1dxd import find_lines_1dxd

def make_uspex_wavecal(files,flatfile,oname,prefix='arc-',\
                    suffix='.[ab].fits*',input_method='index',\
                    normalize=True,qafilename=None,clupdate=True,\
                    overwrite=True):

## Construct file names 
#
#    instrfile = os.path.join(config.state['packagepath'],'instruments',
#                             config.state['instrument'],'data',
#                             config.state['instrument']+'.dat')
#    
#    biasfile = os.path.join(config.state['packagepath'],'instruments',
#                             config.state['instrument'],'data',
#                             config.state['instrument']+'_bias.fits')
#
## Check that things things exist
#
#    test = check_file([instrfile,biasfile])
#
## Read instrument file.  Grab user requested keywords.
#
#    instrinfo = read_instrument_file(instrfile)
#    keywords = instrinfo['XSPEXTOOL_KEYWORDS']
#
## Add GRAT and DIT so that you can determine the mode.  Users will just have
## to live with them in their files
#
#    if 'GRAT' not in keywords:
#        keywords.append('GRAT')
#
#    if 'DIT' not in keywords:
#        keywords.append('DIT')            
#
## Now create the file names
#    
#    if input_method == 'index':
#
#        files = make_full_path(config.state['rawpath'],files,\
#                           indexinfo={'nint':instrinfo['NINT'],\
#                                      'prefix':prefix,'suffix':suffix},\
#                           exist=True)
#
#    elif input_method == 'filename':
#
#        files = make_full_path(rawpath,files,exist=True)
#        
#    else:
#
#        print('make_flat: Unknown input_method.')
#        sys.exit(1)
#
## Load the FITS files into memory
#
#    if clupdate is True: print('Loading FITS images...')
#        
#    lininfo = {'bias':biasfile,'max':instrinfo['LINCORMAX'],'bit':0}
#    img,var,hdr,mask = read_uspex_fits(files,lininfo,keywords=keywords,\
#                                        clupdate=clupdate)
#
## Now scale their intensities to a common flux level
#
#    if clupdate is True and len(files) > 1:
#
#        print('Scaling images...')
#        simgs,svars,scales = scale_data_stack(img,None)
#
#    else: simgs = img
#        
## Now median the scaled images
#
#    if clupdate is True and len(files) >1:
#
#        print('Medianing the images...')
#
#        med,munc = median_data_stack(simgs)
#
#    else: med = simgs
#
## Read the flat and wavecal files
#        
#    flatinfo = read_flat_fits(os.path.join(config.state['calpath'],flatfile))
#
#    wavecalfile = os.path.join(config.state['packagepath'],'instruments',
#                               config.state['instrument'],'data',
#                               flatinfo['mode']+'_wavecalinfo.fits')
#
#    wavecalinfo = read_wavecal_file(wavecalfile)    
#
#    
## Create wavecal and spatcal images
#    
#    wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
#                                             flatinfo['nrows'],
#                                             flatinfo['edgecoeffs'],
#                                             flatinfo['xranges'],
#                                             flatinfo['slith_arc'])
#
#
#    
## Extract the "arc"                                
#
#    spec = extract_extendedsource_1dxd(med, var, flatinfo['ordermask'],
#                                       flatinfo['orders'], wavecal, spatcal,
#                                       flatinfo['slith_arc']/2,
#                                       wavecalinfo['apradius'],
#                                       linmax_bitmask=None, badpixel_mask=None,
#                                       bginfo=None,clupdate=True)
#
###### Do the cross correlation    
#        
## Get the anchor order
#
#    z = flatinfo['orders'] == wavecalinfo['xcororder']
#
#
#    xanchor = np.arange(int(wavecalinfo['xranges'][z,0]),\
#                        int(wavecalinfo['xranges'][z,1]+1),dtype=int)
#    fanchor = np.squeeze(wavecalinfo['spectra'][z,1,:])
#
## Get the source order
#
#    key = 'OR'+str(int(flatinfo['orders'][z])).zfill(3)+'_AP01'
#    xsource = np.squeeze(spec[key][0,:])
#    fsource = np.squeeze(spec[key][1,:])
#
#    qafileinfo = {'figsize':(8.5,11), 'filepath':config.state['qapath'],
#                  'filename':oname+'_findpixelshift','extension':'.pdf'}
#        
#    offset = get_spectral_pixelshift(xanchor, fanchor, xsource, fsource,
#                                     qafileinfo=qafileinfo)
#
#    
## Get the line list
#
#
#    filename = os.path.join(config.state['packagepath'],'instruments',
#                            config.state['instrument'],'data',
#                            wavecalinfo['linelist'])
#
#    lineinfo = read_line_list(filename)
#
#
#    with open('data.sav', 'wb') as f:
#        pickle.dump([spec, wavecalinfo, lineinfo, flatinfo, offset], f)
#
#    return
        
    with open('data.sav', 'rb') as f:
        spec, wavecalinfo, lineinfo, flatinfo, offset = pickle.load(f)

        
    lineinfo['leftwin'] = lineinfo['leftwin']/1e4
    lineinfo['rghtwin'] = lineinfo['rghtwin']/1e4     
    lineinfo = get_line_guess_position(wavecalinfo['spectra'],
                                       wavecalinfo['orders'],
                                       flatinfo['xranges'], lineinfo)

    
    lineinfo['xguess'] = lineinfo['xguess'] + offset
    lineinfo['xleft'] = lineinfo['xleft'] + offset
    lineinfo['xright'] = lineinfo['xright'] + offset        

    lineinfo = find_lines_1dxd(spec, wavecalinfo['orders'], lineinfo,
                               flatinfo['slitw_pix'],
                               qafileinfo={'figsize':(8.5,11)})

    


