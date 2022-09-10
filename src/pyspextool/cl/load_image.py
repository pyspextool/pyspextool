import numpy as np

from pyspextool.cl import config
from pyspextool.io.check_parameter import check_parameter
from pyspextool.io.files import *
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files

def load_image(files, flat_name, *wavecal_name, reduction_mode='A-B',
               directory='raw',suffix=None):


    #
    # Check the parameters
    #

    check_parameter('load_images', 'files', files, ['str', 'list'])
    check_parameter('load_images', 'flat_name', flat_name, 'str')
    check_parameter('load_images', 'wavecal_name', wavecal_name[0], 'str')
    check_parameter('load_images', 'reduction_mode', reduction_mode, 'str',
                    possible_values=['A', 'A-B', 'A-Sky'])
    check_parameter('load_images', 'directory', directory, 'str',
                    possible_values=['raw', 'cal', 'proc'])
    check_parameter('load_images', 'suffix', suffix, ['str', 'NoneType'])

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
        
    # Load the flat field image

    info = read_flat_fits(full_flat_name)

    config.state['flat'] = info['flat']
    config.state['omask'] = info['omask']
    
    
