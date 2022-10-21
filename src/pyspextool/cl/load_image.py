import numpy as np
import importlib
from astropy.io import fits
import glob
import os

from pyspextool.calibration.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import *
from pyspextool.io.files import *
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.io.wavecal import read_wavecal_fits
from pyspextool.plot.plot_image import plot_image
from pyspextool.spectroscopy.rectify_order import rectify_order
from pyspextool.utils.arrays import idl_rotate
from pyspextool.fit.polyfit import poly_1d


def load_image(files, flat_name, *wavecal_name, reduction_mode='A-B',
               directory='raw',suffix=None, flat_field=True,
               linearity_correction=True, clupdate=True, iplot=False,
               qafile=False):

    '''
    To load an (pair-subtracted, sky/dark subrtracted) image into memory.

    Parameters
    ----------

    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the index numbers of the files and 
        files[1] is a string giving the perfix, e.g.
        ['1-2', 'spc'].
    
    flat_name : str
        The full name of a Spextool flat file.

    wavecal_name : str, optional
        The full name of a Spextool wavecal file.

    reduction_mode : {'A-B', 'A', 'A-Sky/Dark'}, optional
        The image reduction mode.

    directory : {'raw', 'cal', 'proc'}, optional
        The directory containing the file(s) to load.

    suffix : str, NoneType, optional
        An optional suffix if in index mode.

    flat_field : {True, False}, optional
        Set to False to not flat field the image.

    linearity_correction : {True, False}, optional
        Set to False to not correct for linearity.

    clupdate : {True, False}, optional
        Set to True for command line updates during execution. 

    iplot : {False, True}, optional
        Set to True for an interactive plot to appear.

    qafile : {False, True}, optional
        Set to True for a QA plot to be generated.  

    '''
    
    #
    # Check the parameters
    #

    check_parameter('load_image', 'files', files, ['str', 'list'])

    check_parameter('load_image', 'flat_name', flat_name, 'str')

    if len(wavecal_name) !=0:
        check_parameter('load_image', 'wavecal_name', wavecal_name[0], 'str')

    check_parameter('load_image', 'reduction_mode', reduction_mode, 'str',
                    possible_values=['A', 'A-B', 'A-Sky/Dark'])

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
    config.state['flatfile'] = flat_name

    if len(wavecal_name) !=0:

        dowavecal = True
        full_wavecal_name = os.path.join(config.state['calpath'],
                                         wavecal_name[0])
        check_file(full_wavecal_name)
        config.state['wavecalfile'] = wavecal_name[0]
        
    else:

        dowavecal = False
        config.state['wavecalfile'] = None

    #
    # Create the file names
    #

    if isinstance(files, str):

        # You are in FILENAME mode

        config.state['filereadmode'] = 'filename'
        files = files.replace(" ", "").split(',')
        input_files = make_full_path(path, files, exist=True)
        
    else:

        # You are in INDEX mode

        config.state['filereadmode'] = 'index'
        
        nums = files[0]
        prefix = files[1]

        # Create the files to read into memory
        
        indexinfo = {'nint': config.state['nint'], 'prefix': prefix,
                     'suffix':config.state['suffix'], 'extension': '.fits*'}
        
        input_files = make_full_path(path, nums, indexinfo=indexinfo,
                                     exist=True)

        # Create the corresponding output file names.
        
        indexinfo = {'nint': config.state['nint'],
                     'prefix': config.state['output_prefix'],
                     'suffix':'', 'extension': '.fits'}
        
        output_files = make_full_path(config.state['procpath'], nums,
                                      indexinfo=indexinfo)        

    # Got the right number?
    
    nfiles = len(input_files)

    if reduction_mode == 'A' and nfiles != 1:

        message = 'The A reduction mode requires 1 image.'
        raise ValueError(message)        

    if reduction_mode == 'A-Sky/Dark' and nfiles != 1:

        message = 'The A-Sky/Dark reduction mode requires 1 image.'
        raise ValueError(message)        

    if reduction_mode == 'A-B' and nfiles != 2:

        message = 'The A-B reduction mode requires 2 images.'
        raise ValueError(message)

    # Do we need to reorder because we are IRTF?
            
    if config.state['irtf'] is True:
            
        input_files = reorder_irtf_files(input_files)

        if config.state['filereadmode'] == 'index':
        
            output_files = reorder_irtf_files(output_files)        
            config.state['output_files'] = output_files
        
    # Create the qafilename root
        
    basenames = []
    for file in input_files:

        basename = os.path.basename(file)
        root = os.path.splitext(basename)
        if root[1] == '.gz':

            root = os.path.splittext(root[0])

        basenames.append(root[0])

    qafilename = '_'.join(basenames)
    config.state['qafilename'] = qafilename


    #
    # Load the flat field image
    #

    if clupdate is True:
        print('Loading the flat...')
    
    flatinfo = read_flat_fits(full_flat_name)

    config.state['reductionmode'] = reduction_mode
    config.state['flat'] = flatinfo['flat']
    config.state['ordermask'] = flatinfo['ordermask']
    config.state['edgecoeffs'] = flatinfo['edgecoeffs']
    config.state['edgedeg'] = flatinfo['edgedeg']    
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
        config.state['psdoorders'] = np.ones(flatinfo['norders'], dtype=int)
        config.state['xsdoorders'] = np.ones(flatinfo['norders'], dtype=int)

    #
    # Load the wavcal image
    #
    
    if clupdate is True:
        print('Loading the wavecal...')

    if dowavecal is True:
        
        wavecalinfo = read_wavecal_fits(full_wavecal_name,rotate=True)
        wavecal = wavecalinfo['wavecal']
        spatcal = wavecalinfo['spatcal']
        wavecaltype = wavecalinfo['wctype']

        #
        # Get the atmospheric transmission 
        #

        # First we have to get the possible file names
        
        fullpath = glob.glob(os.path.join(config.state['packagepath'],
                                          'data','atran*.fits') ) 

        # Then strip the paths off

        basenames = [os.path.basename(x) for x in fullpath]

        # Now get the resolving powers
        
        rps = np.array([int(x[5:x.find('.')]) for x in basenames])

        # Find the closest one

        deltas = rps - flatinfo['rp']
        z = deltas == np.min(deltas)

        # Load that file

        array = fits.getdata(np.array(fullpath)[z][0])
        config.state['atrans_wave'] = array[0,:]
        config.state['atrans_trans'] = array[1,:]
        
        
    else:

        wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                                 flatinfo['nrows'],
                                                 flatinfo['edgecoeffs'],
                                                 flatinfo['xranges'],
                                                 flatinfo['slith_arc'])

        config.state['atrans_wave'] = np.nan
        config.state['atrans_trans'] = np.nan

    config.state['wavecal'] = wavecal
    config.state['spatcal'] = spatcal
        
        
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
            img, var, hdr, mask = readfits.main(input_files,
                                            config.state['linearity_info'],
                                keywords=config.state['xspextool_keywords'],
                                linearity_correction=linearity_correction,
                                            clupdate=clupdate)


        else:
            print('Do later')


    elif reduction_mode =='A-B':

        img, var, hdr, mask = readfits.main(input_files,
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

    config.state['workimage'] = img
    config.state['varimage'] = var
    config.state['hdrinfo'] = hdr
    
    #
    # Rotate the bad pixel mask
    #

    config.state['badpixelmask'] = idl_rotate(config.state['rawbadpixelmask'],
                                              flatinfo['rotation'])
    #
    # Rectify the orders
    #
    rectorders = []
    indices = wavecalinfo['rectindices']
    for i in range(config.state['norders']):

        order = rectify_order(img, indices[i]['xidx'], indices[i]['yidx'])
#                              var=var, bpmask=config.state['badpixelmask'],
#                              bsmask = mask)


        # Now get the wavelength solution to tack on


        bot = np.ceil(poly_1d(indices[i]['x'],
                      flatinfo['edgecoeffs'][i,0,:])).astype(int)

        w = wavecal[bot,np.fix(indices[i]['x']).astype(int)]

        order.update({'w':w})
        order.update({'y':indices[i]['y']})        
        
        rectorders.append(order)
        
    # Store the results

    config.state['rectorders'] = rectorders

    #
    # Do the plotting
    #

    order_plotinfo = {'xranges':config.state['xranges'],
                      'edgecoeffs':config.state['edgecoeffs'],
                      'orders':config.state['orders']}
    
    if iplot is True:

        plot_image(config.state['workimage'],orders_plotinfo=order_plotinfo)

    if qafile is True:

        qafileinfo = {'figsize': (7,7), 'filepath':config.state['qapath'],
                      'filename':qafilename+'_image',
                      'extension':config.state['qaextension']}
       
        plot_image(config.state['workimage'],
                   orders_plotinfo=order_plotinfo,
                   qafileinfo=qafileinfo)        
    
    # Set the continue flags

    config.state['continue'] = 1
