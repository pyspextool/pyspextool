import numpy as np
import importlib
from astropy.io import fits

from pyspextool.calibration.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.cl import config
from pyspextool.io.check_parameter import check_parameter
from pyspextool.io.files import *
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.io.wavecal import read_wavecal_fits
from pyspextool.utils.arrays import idl_rotate



def load_image(files, flat_name, *wavecal_name, reduction_mode='A-B',
               directory='raw',suffix=None, flat_field=True,
               linearity_correction=True, clupdate=True):

    
    #
    # Check the parameters
    #

    check_parameter('load_image', 'files', files, ['str', 'list'])
    check_parameter('load_image', 'flat_name', flat_name, 'str')
    if len(wavecal_name) !=0:
        check_parameter('load_image', 'wavecal_name', wavecal_name[0], 'str')
    check_parameter('load_image', 'reduction_mode', reduction_mode, 'str',
                    possible_values=['A', 'A-B', 'A-Sky'])
    check_parameter('load_image', 'directory', directory, 'str',
                    possible_values=['raw', 'cal', 'proc'])
    check_parameter('load_image', 'suffix', suffix, ['str', 'NoneType'])
    check_parameter('load_image', 'flat_field', flat_field, 'bool')
    check_parameter('load_image', 'linearity_correction',
                    linearity_correction, 'bool')        
    check_parameter('load_image', 'clupdate', clupdate, 'bool')    

    # Get the readfits module
    
    module =  'pyspextool.io.'+config.state['readfits']
    readfits = importlib.import_module(module)
    
    # Get the path

    if directory == 'raw':

        path = config.state['rawpath']

    elif directory == 'cal':

        path = config.state['calpath']

    elif directory == 'proc':

        path = config.state['procpath']                
        
    #
    # Check for file existence
    #

    full_flat_name = os.path.join(config.state['calpath'], flat_name)
    check_file(full_flat_name)

    if len(wavecal_name) !=0:

        dowavecal = True
        full_wavecal_name = os.path.join(config.state['calpath'],
                                         wavecal_name[0])
        
    else:

        dowavecal = False

    #
    # Create the file names
    #

    if isinstance(files, str):

        # You are in FILENAME mode

        files = make_full_path(path, files, exist=True)
        
    else:

        # You are in INDEX mode

        nums = files[0]
        prefix = files[1]

        files = make_full_path(path, nums,
                               indexinfo={'nint': config.state['nint'],
                                          'prefix': prefix,
                                          'suffix': config.state['suffix'],
                                          'extension': '.fits*'},
                               exist=True)

    # Got the right number?
    
    nfiles = len(files)

    if reduction_mode == 'A' and nfiles != 1:

        message = 'The A reduction mode requires 1 image.'
        raise ValueError(message)        

    if reduction_mode == 'A-Sky' and nfiles != 1:

        message = 'The A-Sky reduction mode requires 1 image.'
        raise ValueError(message)        

    
    if reduction_mode == 'A-B' and nfiles != 2:

        message = 'The A-B reduction mode requires 2 images.'
        raise ValueError(message)

        
    if config.state['irtf'] is True:
        # Reorder files to abab
            
        files = reorder_irtf_files(files)

    # load the files into the state variable

    config.state['workfiles'] = files

    #
    # Load the flat field image
    #

    if clupdate is True:
        print('Loading the flat...')
    
    flatinfo = read_flat_fits(full_flat_name)

    config.state['flat'] = flatinfo['flat']
    config.state['omask'] = flatinfo['ordermask']
    config.state['edgecoeffs'] = flatinfo['edgecoeffs']
    config.state['orders'] = flatinfo['orders']
    config.state['norders'] = flatinfo['norders']
    config.state['plate_scale'] = flatinfo['ps']
    config.state['slith_arc'] = flatinfo['slith_arc']
    config.state['slith_pix'] = flatinfo['slith_pix']
    config.state['slitw_arc'] = flatinfo['slitw_arc']
    config.state['slitw_pix'] = flatinfo['slitw_pix']
    config.state['resolvingpower'] = flatinfo['rp']
    config.state['xranges'] = flatinfo['xranges']
    config.state['rotation'] = flatinfo['rotation']
    config.state['ncols'] = flatinfo['ncols']
    config.state['nrows'] = flatinfo['nrows']
    config.state['badpixelmask'] = idl_rotate(config.state['rawbadpixelmask'],
                                              flatinfo['rotation'])    

    # Let's deal with the mode

    if config.state['modename'] != flatinfo['mode']:
 
        config.state['modename'] = flatinfo['mode']
        config.state['psdoorders'] = np.ones(flatinfo['norders'])
        config.state['xsdoorders'] = np.ones(flatinfo['norders'])

    #
    # Load the wavcal image
    #
    
    if clupdate is True:
        print('Loading the wavecal...')

    if dowavecal is True:
        
        wavecalinfo = read_wavecal_fits(full_wavecal_name,rotate=True)
        wavecal = wavecalinfo['wavecal']
        spatcal = wavecalinfo['spatcal']        
            
    else:

        wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                                 flatinfo['nrows'],
                                                 flatinfo['edgecoeffs'],
                                                 flatinfo['xranges'],
                                                 flatinfo['slith_arc'])

    #
    # Load the data
    #

    if clupdate is True:

        if linearity_correction is True:
           
            print('Loading the image and correcting for non-linearity...') 

        else:

            print('Loading the image and not correcting for non-linearity...') 

        
    if reduction_mode == 'A':

        if directory == 'raw':
            img, var, hdr, mask = readfits.main(files,
                                            config.state['linearity_info'],
                                keywords=config.state['xspextool_keywords'],
                                linearity_correction=linearity_correction,
                                            clupdate=clupdate)


        else:
            print('Do later')


    elif reduction_mode =='A-B':

        img, var, hdr, mask = readfits.main(files,
                                            config.state['linearity_info'],
                                            pair=True, 
                                    keywords=config.state['xspextool_keywords'],
                                    linearity_correction=linearity_correction,
                                                clupdate=clupdate)
    else:

        print('do later.')

    #
    # Flat field the data
    #

    if flat_field is True:

        if clupdate is True:
            print('Flat fielding the image...')

        np.divide(img, flatinfo['flat'], out=img)
        np.divide(var, flatinfo['flat']**2, out=var)
            
    else:

        if clupdate is True:
            print('Image not flat fielded...')

    #
    # Rectify the orders
    #
    
    
