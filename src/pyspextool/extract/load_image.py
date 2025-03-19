import importlib
import numpy as np
from astropy.io import fits
from os.path import join, basename as osbasename, splitext
import logging
import glob

from pyspextool import config as setup
from pyspextool.extract import config as extract

from pyspextool.extract.wavecal import simulate_wavecal_1dxd
from pyspextool.extract.images import rectify_order
from pyspextool.io.check import check_parameter, check_qakeywords, check_file
from pyspextool.io.files import files_to_fullpath
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.extract.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.extract.wavecal import read_wavecal_fits
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.utils.math import combine_flag_stack
from pyspextool.pyspextoolerror import pySpextoolError


def load_image(files:str | list,
               flat_name:str,
               wavecal_name:str,
               output_filenames:str=None,
               output_prefix:str='spectra',
               input_extension:str='.fits*',
               load_directory='raw',
               flat_field=True,
               linearity_correction=True,
               detector_info:dict=None,
               write_rectified_orders:bool=False,
               do_all_steps:bool=False,
               verbose:bool=None,
               qa_show:bool=None,
               qa_showscale:float | int=None,
               qa_showblock:bool=None,
               qa_write:bool=None):

    """
    To load a (pair-subtracted, sky/dark subrtracted) image into memory.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.

        e.g. ['spc', '1-2']
    
    flat_name : str
        The full name of a pySpextool flat file.

    wavecal_name : str or None
        The full name of a pySpextool wavecal file.
    
    output_files : str or list of str, optional
        A str or list of str of output files names.  Only required if `files` 
        gives file names intead of a prefix and index numbers.

    directory : {'raw', 'cal', 'proc'}, optional
        The directory containing the file(s) to load.

    suffix : str, NoneType, optional
        An optional suffix if in index mode.

    flat_field : {True, False}, optional
        Set to False to not flat field the image.

    linearity_correction : {True, False}, optional
        Set to False to not correct for linearity.

    detector_info : dict, deafult None
        A dictionary with any information that needs to be passed to the
        instrument specific readfits program.  
    
    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does NOT affect plots written
        to disk.  Set to None to default to setup.state['qa_scale'].    
   
    do_all_steps : {False, True}, optional
        Set to True to skip loading the flat and wavecals.

    Returns
    -------
    None.  Loads data into the config.extract variable.

    """

    #
    # Check the parameters and QA keywords
    #
    
    check_parameter('load_image', 'files', files, ['str', 'list'],
                    list_types=['str','str'])

    check_parameter('load_image', 'flat_name', flat_name, 'str')

    check_parameter('load_image', 'wavecal_name', wavecal_name,
                    ['NoneType','str'])
    
    check_parameter('load_image', 'output_filenames', output_filenames,
                    ['NoneType','str','list'])
        
    check_parameter('load_image', 'load_directory', load_directory, 'str',
                    possible_values=['raw', 'proc'])

    check_parameter('load_image', 'flat_field', flat_field, 'bool')

    check_parameter('load_image', 'linearity_correction',
                    linearity_correction, 'bool')

    check_parameter('load_image', 'detector_info', detector_info,
                    ['NoneType','dict'])
    
    check_parameter('load_image', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('load_image', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Create and ensure that the flat field and wavecal file existence
    #

    full_flat_name = join(setup.state['cal_path'], flat_name)
    check_file(full_flat_name)

    if isinstance(wavecal_name, str):

        dowavecal = True
        full_wavecal_name = join(setup.state['cal_path'], wavecal_name)
        check_file(full_wavecal_name)
        extract.state['wavecalfile'] = wavecal_name        

    else:

        dowavecal = False
        extract.state['wavecalfile'] = None
    
    #
    # Determine the path in which `files` is located.  
    #
    
    if load_directory == 'raw':

        load_path = setup.state['raw_path']
        image_type = 'raw'
        
    elif load_directory == 'proc':

        load_path = setup.state['proc_path']
        image_type = 'combined'
        
        
    #
    # Create the full input file names
    #

    results = files_to_fullpath(load_path,
                                files,
                                setup.state['nint'],
                                setup.state['suffix'],
                                input_extension)

    input_files = results[0]
    file_readmode = results[1]

    check_file(input_files)

    n_inputfiles = len(input_files)    
    
    #
    # Create the full output file names
    #

    if output_filenames is not None:

        # The user wants to use their own file names

        # Convert to a list if a str
        
        if isinstance(output_filenames,str) is True:
            
            output_filenames = [output_filenames]
            
        # Check to make sure the number of files giveb equals the number of
        # files we are loading
            
        if n_inputfiles != len(output_filenames):

            message = 'The number of files in the keyword '+ \
                '`output_filenames` does not match the number of files in '+ \
                'the parameter `files`.'
            raise pySpextoolError(message)

    else:

        # The user wants them created by us.  Do it differently depending
        # on what readmode was selected.
        
        if file_readmode == 'index':
        
            files[0] = output_prefix
            
            result = files_to_fullpath('',
                                       files,
                                       setup.state['nint'],
                                       '',
                                       '',
                                       exist=False)

            output_filenames = result[0]        

        else:

            output_filenames = []
            for file in input_files:

                root = splitext(osbasename(file))
                if root[1] == '.gz':
                    root = splitext(root[0])
                    
                output_filenames.append(output_prefix+'_'+root[0]+'.fits')

    #
    # Determine the reduction mode
    #

    if n_inputfiles == 1:

       reduction_mode = 'A'

    elif n_inputfiles == 2:

        reduction_mode = 'A-B'
    
    else:

        message = 'More than two files cannot be passed.'
        raise pySpextoolError(message)

    if do_all_steps is False:
        logging.info(' Setting reduction mode to '+reduction_mode+'.')
        
    # Do we need to reorder because we are IRTF?

    if reduction_mode == 'A-B' and setup.state['irtf'] is True:

        input_files, indices = reorder_irtf_files(input_files)

        tmp = np.array(output_filenames)
        tmp = tmp[indices]
        output_filenames = tmp.tolist()
                           
    extract.state['output_files'] = output_filenames
    
    # Create the qafilename root

    basenames = []
    for file in input_files:

        basename = osbasename(file)
        root = splitext(basename)
        if root[1] == '.gz':
            root = splitext(root[0])

        basenames.append(root[0])

    qafilename = '_'.join(basenames)
    extract.state['qafilename'] = qafilename
    
    #
    # Is this call a part of a do_all_steps call?  If so, load just the data
    #
    
    if do_all_steps is False:

        #
        # Load the flat field image and store important things
        #

        logging.info(' Loading the flat file '+flat_name+'.')

        flatinfo = read_flat_fits(full_flat_name)

        extract.state['flat_filename'] = flat_name
        extract.state['reductionmode'] = reduction_mode
        extract.state['flat'] = flatinfo['flat']
        extract.state['flat_bitmask'] = flatinfo['bitmask']        
        extract.state['ordermask'] = flatinfo['ordermask']
        extract.state['edgecoeffs'] = flatinfo['edgecoeffs']
        extract.state['edgedeg'] = flatinfo['edgedeg']
        extract.state['orders'] = flatinfo['orders']
        extract.state['norders'] = flatinfo['norders']
        extract.state['plate_scale'] = flatinfo['ps']
        extract.state['slith_arc'] = flatinfo['slith_arc']
        extract.state['slith_pix'] = flatinfo['slith_pix']
        extract.state['slitw_arc'] = flatinfo['slitw_arc']
        extract.state['slitw_pix'] = flatinfo['slitw_pix']
        extract.state['resolvingpower'] = flatinfo['rp']
        extract.state['xranges'] = flatinfo['xranges']
        extract.state['rotation'] = flatinfo['rotation']
        extract.state['ncols'] = flatinfo['ncols']
        extract.state['nrows'] = flatinfo['nrows']

        # Let's deal with doorders arrays.
        # If the mode is the same, do nothing.
        # If the mode is different, or this is the first load, create doorders
        # arrays
        
#        if extract.state['modename'] != flatinfo['mode']:

        extract.state['modename'] = flatinfo['mode']
        extract.state['doorders'] = np.ones(flatinfo['norders'], dtype=int)

        #
        # Load the wavcal image
        #
        
        if dowavecal is True:

            logging.info(' Loading the wavecal file '+wavecal_name+'.')
            
            wavecalinfo = read_wavecal_fits(full_wavecal_name,
                                            rotate=True)
            
            wavecal = wavecalinfo['wavecal']
            spatcal = wavecalinfo['spatcal']
            indices = wavecalinfo['rectindices']
            wavecaltype = wavecalinfo['wctype']
            
            #
            # Get the atmospheric transmission 
            #

            # First we have to get the possible file names

            fullpath = glob.glob(join(setup.state['package_path'],
                                      'data', 'atran*.fits'))

            # Then strip the paths off

            basenames = [osbasename(x) for x in fullpath]

            # Now get the resolving powers

            rps = np.array([int(x[5:x.find('.')]) for x in basenames])

            # Find the closest one

            deltas = rps - flatinfo['rp']
#            z = deltas == np.min(np.abs(deltas))
            zind = np.argmin(np.abs(deltas))
                             
            logging.info(' Loading the atmospheric tranmission at R='+\
                         str(rps[zind])+'.')
            
            # Load that file

#            array = fits.getdata(np.array(fullpath)[z][0])
            array = fits.getdata(np.array(fullpath)[zind])
            atmosphere = {'wavelength':array[0, :], 'transmission':array[1, :]}
            
            # Get units

            xunits = 'um'
            latex_xunits = r'$\mu$m'
            latex_xlabel = r'Wavelength ($\mu$m)'
                        
        else:

            result = simulate_wavecal_1dxd(flatinfo['ncols'],
                                           flatinfo['nrows'],
                                           flatinfo['edgecoeffs'],
                                           flatinfo['xranges'],
                                           flatinfo['slith_arc'])

            wavecal = result[0]
            spatcal = result[1]            
            indices = result[2]

            for i in range(extract.state['norders']):

                indices[i]['w'] = indices[i]['x']
                indices[i]['a'] = indices[i]['y']                
            
            atmosphere  = None
            xunits = 'pixel'
            latex_xunits = 'pixel'
            latex_xlabel = 'Wavelength (pixel)'
            wavecaltype = None

        # Store useful things

        extract.state['wavecal'] = wavecal
        extract.state['spatcal'] = spatcal
        extract.state['rectindices'] = indices
        extract.state['atmosphere'] = atmosphere
        extract.state['xunits'] = xunits
        extract.state['latex_xunits'] = latex_xunits
        extract.state['latex_xlabel'] = latex_xlabel
        extract.state['latex_yunits'] = 'DN s$^{-1}$'
        extract.state['latex_ylabel'] = 'Count Rate (DN s$^{-1}$)'
        extract.state['latex_ulabel'] = r'$\sigma$ (DN s$^{-1}$)'        
        extract.state['wavecaltype'] = wavecaltype
        
    #
    # Load the data
    #

    # Get the readfits module

    module = 'pyspextool.instruments.' + setup.state['instrument'] + \
             '.' + setup.state['instrument']

    instr = importlib.import_module(module)
        
    # Now go reduction mode by reduction mode

    if reduction_mode == 'A':
        
        # This could either be a standard A or a combined image.

        if image_type == 'raw':

            result = instr.read_fits(input_files,
                                     setup.state['linearity_info'],
                                     rotate=extract.state['rotation'],
                                     keywords=setup.state['extract_keywords'],
                                     linearity_correction=linearity_correction,
                                     extra=detector_info,
                                     verbose=qa['verbose'])

            img = result[0]
            var = result[1]
            hdrinfo = result[2]
            linearity_mask = result[3]
            
        else:

            # This is a combined image so just read in.  

            hdul = fits.open(input_files[0])
            hdul[0].verify('silentfix')

            hdr = hdul[0].header

            hdrinfo = [get_headerinfo(hdr)]
           
            img = idl_rotate(hdul[1].data, extract.state['rotation'])
            var = idl_rotate(hdul[2].data, extract.state['rotation'])
            linearity_mask = idl_rotate(hdul[3].data, extract.state['rotation'])
            
            hdul.close()

    elif reduction_mode == 'A-B':

        result = instr.read_fits(input_files,
                                 setup.state['linearity_info'],
                                 pair_subtract=True,
                                 rotate=extract.state['rotation'],
                                 keywords=setup.state['extract_keywords'],
                                 linearity_correction=linearity_correction,
                                 extra=detector_info,
                                 verbose=qa['verbose'])

        img = result[0]
        var = result[1]
        hdrinfo = result[2]
        linearity_mask = result[3]
        
    else:

        print('do later.')
    
    #
    # Flat field the data
    #

    if flat_field is True:

        logging.info(' Flat fielding the image.')                

        extract.state['flat_fielded'] = True
        
        np.divide(img, extract.state['flat'], out=img)
        np.divide(var, extract.state['flat'] ** 2, out=var)
        
        # Combine the masks

        flag_mask = combine_flag_stack(np.stack((extract.state['flat_bitmask'],
                                                 linearity_mask)))
        
    else:

        logging.info(' Image not flat field.')
        extract.state['flat_fielded'] = False        
        flag_mask = linearity_mask
              
    # Store the results

    extract.state['workimage'] = img
    extract.state['varimage'] = var
    extract.state['maskimage'] = flag_mask    
    extract.state['hdrinfo'] = hdrinfo
 
    #
    # Rotate the bad pixel mask
    #

    bad_pixel_mask = idl_rotate(setup.state['raw_bad_pixel_mask'],
                                extract.state['rotation'])

    extract.state['bad_pixel_mask'] = bad_pixel_mask
        
    #
    # Rectify the orders
    #

    logging.info(' Rectifying the orders.')                
    
    rectorders = []
    indices = extract.state['rectindices']
    for i in range(extract.state['norders']):

        loop_progress(i,0,extract.state['norders'])
        result = rectify_order(img,
                               indices[i]['xidx'],
                               indices[i]['yidx'],
                               variance=var,
                               bad_pixel_mask=bad_pixel_mask,
                               flag_mask=flag_mask)

        rectorder = {'wavelengths':indices[i]['w'],
                     'angles':indices[i]['a'],
                     'image':result['image'],
                     'variance':result['variance'],
                     'badpixel_mask':result['bpmask'],
                     'flag_mask':result['flagmask']}

        rectorders.append(rectorder)        


    # Store the results

    extract.state['rectorders'] = rectorders

    #        
    # Make the QA plot
    #

    orders_plotinfo = {'xranges': extract.state['xranges'],
                       'edgecoeffs': extract.state['edgecoeffs'],
                       'orders': extract.state['orders']}
              
    if qa['show'] is True:
        
        plot_image(img,
                   mask=linearity_mask,
                   orders_plotinfo=orders_plotinfo,
                   figure_size=(setup.plots['square_size'][0]*qa['showscale'],
                                setup.plots['square_size'][1]*qa['showscale']),
                   font_size=setup.plots['font_size']*qa['showscale'],
                   showblock=qa['showblock'],
                   plot_number=setup.plots['combine_image'])
            
    if qa['write'] is True:

        filename = qafilename + '_image' + setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

        plot_image(img,
                   mask=linearity_mask,
                   orders_plotinfo=orders_plotinfo,
                   output_fullpath=fullpath,
                   figure_size=setup.plots['square_size'],
                   font_size=setup.plots['font_size'])

    #
    # Write rectified orders
    #

    if write_rectified_orders is True:

        # Get output file name 
        
        basename = qafilename+'_rectifiedorders'+'.fits'
        fullpath = join(setup.state['proc_path'], basename)
        
        logging.info(' Writing rectified orders to file '+basename+'.')

        # Create the primary HDU

        phdu = fits.PrimaryHDU()
        hdr = phdu.header

        hdr['MODE'] = (flatinfo['mode'], ' Instrument Mode')
        hdr['NORDERS'] = (flatinfo['norders'], ' Number of orders identified')
        orders = ','.join(str(o) for o in flatinfo['orders'])
        hdr['ORDERS'] = (orders, ' Orders')
        
        list_hdu = [phdu]
        
        # Loop over each order and create the array of images
        
        
        for i in range(extract.state['norders']):

            ny, nx = np.shape(rectorders[i]['image'])
        
            wavemap = np.tile(rectorders[i]['wavelengths'],(ny,1))
            spatmap = np.rot90(np.tile(rectorders[i]['angles'], (nx, 1)), k=3)

            array = np.full((6, ny,nx), np.nan)
            array[0,:,:] = wavemap
            array[1,:,:] = spatmap
            array[2,:,:] = rectorders[i]['image']
            array[3,:,:] = rectorders[i]['variance']
            array[4,:,:] = rectorders[i]['badpixel_mask']
            array[5,:,:] = rectorders[i]['flag_mask']
            
            idx_hdu = fits.ImageHDU(np.float32(array))        
            list_hdu.append(idx_hdu)

        # Now make the FITS file
            
        hdu = fits.HDUList(list_hdu)
        hdu.writeto(fullpath, overwrite=True)

    #
    # Set the done variables
    #

    extract.state['load_done'] = True
    extract.state['profiles_done'] = False
    extract.state['apertures_done'] = False
    extract.state['orders_done'] = False
    extract.state['trace_done'] = False
    extract.state['parameters_done'] = False
    extract.state['extract_done'] = False
