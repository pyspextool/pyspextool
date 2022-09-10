import pyspextool.io.check_parameter import check_parameter

def load_images(fullpaths, flat_file, reduction_mode, *wavecal_name):

    check_parameter('load_images', 'fullpaths', fullpaths, ['list'])
    check_parameter('load_images', 'flat_file', flat_file, 'str')
    check_parameter('load_images', 'reduction_mode', reduction_mode, 'str')    
    
