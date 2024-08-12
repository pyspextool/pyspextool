import importlib
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.extract import config as extract

from pyspextool.extract.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.extract.rectify_order import rectify_order
from pyspextool.io.check import *
from pyspextool.io.files import *
from pyspextool.io.fitsheader import get_header_info
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.reorder_irtf_files import reorder_irtf_files
from pyspextool.io.wavecal import read_wavecal_fits
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils.arrays import idl_rotate
from pyspextool.utils.math import combine_flag_stack
from pyspextool.fit.polyfit import poly_1d


def load_image(files, flat_name, wavecal_name, *output_files,
               reduction_mode='A-B', directory='raw', suffix=None,
               flat_field=True, linearity_correction=True, verbose=None,
               do_all_steps=False, qa_show=None, qa_write=None,
               qa_showsize=(6, 6)):
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
        The full name of a Spextool flat file.

    wavecal_name : str or NoneType, optional
        The full name of a Spextool wavecal file.

    output_files : str or list of str, optional
        A str or list of str of output files names.  Only required if `files` 
        gives file names intead of a prefix and index numbers.

    reduction_mode : {'A-B', 'A', 'A-Sky'}, optional
        The image reduction mode.

    directory : {'raw', 'cal', 'proc'}, optional
        The directory containing the file(s) to load.

    suffix : str, NoneType, optional
        An optional suffix if in index mode.

    flat_field : {True, False}, optional
        Set to False to not flat field the image.

    linearity_correction : {True, False}, optional
        Set to False to not correct for linearity.

    verbose : {None, True, False}, optional
        Set to True/False to override config.setup['verbose'].

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the 
        Pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qa_showsize : tuple, default=(6,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    do_all_steps : {False, True}, optional
        Set to True to skip loading the flat and wavecals.

    Returns
    -------
    None.  Loads data into the config.extract variable.

    """

    #
    # Check the parameters
    #

    check_parameter('load_image', 'files', files, ['str', 'list'])

    if isinstance(files, list):
        check_parameter('load_image', 'files[0]', files[0], 'str')
        check_parameter('load_image', 'files[1]', files[1],
                        ['str', 'list', 'int'])

    check_parameter('load_image', 'flat_name', flat_name, 'str')

    check_parameter('load_image', 'wavecal_name', wavecal_name,
                    ['NoneType','str'])

    if len(output_files) != 0:
        check_parameter('load_image', 'output_files', output_files[0], 'str')                
    check_parameter('load_image', 'reduction_mode', reduction_mode, 'str',
                    possible_values=['A', 'A-B', 'A-Sky'])

    check_parameter('load_image', 'directory', directory, 'str',
                    possible_values=['raw', 'cal', 'proc'])

    check_parameter('load_image', 'suffix', suffix, ['str', 'NoneType'])

    check_parameter('load_image', 'flat_field', flat_field, 'bool')

    check_parameter('load_image', 'linearity_correction',
                    linearity_correction, 'bool')

    check_parameter('load_image', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('load_image', 'qa_showsize', qa_showsize,
                    ['NoneType', 'tuple'])

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_write is None:
        qa_write = setup.state['qa_write']

    if qa_show is None:
        qa_show = setup.state['qa_show']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Store user inputs
    #

    extract.load['doflat'] = flat_field
    extract.load['dolinearity'] = linearity_correction
    extract.load['qaplotsize'] = qa_showsize
    extract.load['flatfile'] = flat_name
    extract.load['wavecalfile'] = wavecal_name
    extract.load['qaplot'] = qa_show
    extract.load['qafile'] = qa_write
    extract.load['verbose'] = verbose

    # Get the readfits module

    module = 'pyspextool.instrument_data.' + setup.state['instrument'] + \
             '_dir.' + setup.state['instrument']

    instr = importlib.import_module(module)

    # Get the path

    if directory == 'raw':

        path = setup.state['raw_path']

    elif directory == 'cal':

        path = setup.state['cal_path']

    elif directory == 'proc':

        path = setup.state['proc_path']

    #
    # Check for file existence
    #

    full_flat_name = os.path.join(setup.state['cal_path'], flat_name)
    check_file(full_flat_name)

    if isinstance(wavecal_name, str):

        dowavecal = True
        full_wavecal_name = os.path.join(setup.state['cal_path'], wavecal_name)
        check_file(full_wavecal_name)

    else:

        dowavecal = False
        extract.load['wavecalfile'] = None

    #
    # Create the file names
    #

    if isinstance(files, str):

        # You are in FILENAME mode

        extract.state['filereadmode'] = 'filename'
        files = files.replace(" ", "").split(',')
        input_files = make_full_path(path, files, exist=True)
        output_files = make_full_path(setup.state['proc_path'],
                                      [output_files[0]])

    else:

        # You are in INDEX mode

        extract.state['filereadmode'] = 'index'

        prefix = files[0]
        nums = files[1]

        extract.state['prefix'] = prefix

        # Create the files to read into memory

        indexinfo = {'nint': setup.state['nint'], 'prefix': prefix,
                     'suffix': setup.state['suffix'], 'extension': '.fits*'}

        input_files = make_full_path(path, nums, indexinfo=indexinfo,
                                     exist=True)

        # Create the corresponding output file names.

        indexinfo = {'nint': setup.state['nint'],
                     'prefix': extract.state['output_prefix'],
                     'suffix': '', 'extension':''}

        output_files = make_full_path(setup.state['proc_path'], nums,
                                      indexinfo=indexinfo)

    # Got the right number?

    nfiles = len(input_files)

    if reduction_mode == 'A' and nfiles != 1:
        message = 'The A reduction mode requires 1 image.'
        raise ValueError(message)

    if reduction_mode == 'A-Sky' and nfiles != 1:
        message = 'The A-Sky reduction mode requires 1 image.'
        raise ValueError(message)

    if reduction_mode == 'A-B' and nfiles != 2:
        message = 'The A-B reduction mode requires 2 images.'
        raise ValueError(message)

        # Do we need to reorder because we are IRTF?

    if reduction_mode == 'A-B' and setup.state['irtf'] is True:

        input_files, indices = reorder_irtf_files(input_files)
        if extract.state['filereadmode'] == 'index':

            tmp = np.array(output_files)
            tmp = tmp[indices]
            output_files = tmp.tolist()
            
    extract.state['output_files'] = output_files

    # Create the qafilename root

    basenames = []
    for file in input_files:

        basename = os.path.basename(file)
        root = os.path.splitext(basename)
        if root[1] == '.gz':
            root = os.path.splitext(root[0])

        basenames.append(root[0])

    qafilename = '_'.join(basenames)
    extract.state['qafilename'] = qafilename

    if do_all_steps is False:

        #
        # Load the flat field image
        #

        if verbose is True:
            print('Loading the flat...')

        flatinfo = read_flat_fits(full_flat_name)

        extract.state['reductionmode'] = reduction_mode
        extract.state['flat'] = flatinfo['flat']
        extract.state['bitmask'] = flatinfo['bitmask']        
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

        # Let's deal with the mode

        if extract.state['modename'] != flatinfo['mode']:
            extract.state['modename'] = flatinfo['mode']
            extract.state['psdoorders'] = \
                np.ones(flatinfo['norders'], dtype=int)
            extract.state['xsdoorders'] = \
                np.ones(flatinfo['norders'], dtype=int)

        #
        # Load the wavcal image
        #

        if verbose is True:
            print('Loading the wavecal...')

        if dowavecal is True:

            wavecalinfo = read_wavecal_fits(full_wavecal_name, rotate=True)
            wavecal = wavecalinfo['wavecal']
            spatcal = wavecalinfo['spatcal']
            indices = wavecalinfo['rectindices']
#            dispersions = wavecalinfo['dispersions']            
            
            #
            # Get the atmospheric transmission 
            #

            # First we have to get the possible file names

            fullpath = glob.glob(os.path.join(setup.state['package_path'],
                                              'data', 'atran*.fits'))
            # Then strip the paths off

            basenames = [os.path.basename(x) for x in fullpath]

            # Now get the resolving powers

            rps = np.array([int(x[5:x.find('.')]) for x in basenames])

            # Find the closest one

            deltas = rps - flatinfo['rp']
            z = deltas == np.min(deltas)

            # Load that file

            array = fits.getdata(np.array(fullpath)[z][0])
            atmosphere = {'wavelength':array[0, :], 'transmission':array[1, :]}

            # Get units

            xunits = 'um'
            latex_xunits = r'$\mu$m'
            latex_xlabel = r'Wavelength ($\mu$m)'
            
            
        else:

            wavecal, spatcal, indices = simulate_wavecal_1dxd(flatinfo['ncols'],
                                                             flatinfo['nrows'],
                                                        flatinfo['edgecoeffs'],
                                                            flatinfo['xranges'],
                                                        flatinfo['slith_arc'])

#            dispersions = None
            atmosphere  = None

            xunits = 'pixel'
            latex_xunits = 'pixel'
            latex_xlabel = 'Wavelength (pixel)'

            

        extract.state['wavecal'] = wavecal
        extract.state['spatcal'] = spatcal
        extract.state['rectindices'] = indices
#        extract.state['dispersions'] = dispersions
        extract.state['atmosphere'] = atmosphere
        extract.state['xunits'] = xunits
        extract.state['latex_xunits'] = latex_xunits
        extract.state['latex_xlabel'] = latex_xlabel
        
    #
    # Load the data
    #
    
    if verbose is True:

        # Get the file names

        list_filenames = []
        for file in input_files:
            list_filenames.append(os.path.basename(file))

        file_names = ', '
        file_names = file_names.join(list_filenames)

        if linearity_correction is True:

            message = 'Loading ' + file_names + \
                      ' and correcting for non-linearity...'
            print(message)

        else:

            message = 'Loading ' + file_names + \
              ' and not correcting for non-linearity...'
            print(message)

    if reduction_mode == 'A':

        if directory == 'raw':
            img, var, hdrinfo, mask = instr.read_fits(input_files,
                                                  setup.state['linearity_info'],
                                            rotate=extract.state['rotation'],
                                    keywords=setup.state['extract_keywords'],
                                    linearity_correction=linearity_correction,
                                                  verbose=verbose)

        else:

            # Read the data 

            hdul = fits.open(input_files[0])
            hdul[0].verify('silentfix')

            hdr = hdul[0].header

            hdrinfo = [get_header_info(hdr)]
           
            img = idl_rotate(hdul[1].data, extract.state['rotation'])
            var = idl_rotate(hdul[2].data, extract.state['rotation'])
        
            hdul.close()

    elif reduction_mode == 'A-B':

        img, var, hdrinfo, mask = instr.read_fits(input_files,
                                              setup.state['linearity_info'],
                                              pair_subtract=True,
                                              rotate=extract.state['rotation'],
                                    keywords=setup.state['extract_keywords'],
                                    linearity_correction=linearity_correction,
                                              verbose=verbose)
    else:

        print('do later.')
        
    #
    # Flat field the data
    #

    if flat_field is True:

        if verbose is True:
            print('Flat fielding the image...')

        np.divide(img, extract.state['flat'], out=img)
        np.divide(var, extract.state['flat'] ** 2, out=var)

        # Combine the masks

        mask = combine_flag_stack(np.stack((extract.state['bitmask'],mask)))
    
    else:

        if verbose is True:
            print('Image not flat fielded...')

    # Store the results

    extract.state['workimage'] = img
    extract.state['varimage'] = var
    extract.state['maskimage'] = mask    
    extract.state['hdrinfo'] = hdrinfo

    #
    # Rotate the bad pixel mask
    #

    extract.state['bad_pixel_mask'] = \
        idl_rotate(setup.state['raw_bad_pixel_mask'],
                   extract.state['rotation'])

    #
    # Rectify the orders
    #

    rectorders = []
    indices = extract.state['rectindices']
    for i in range(extract.state['norders']):
        order = rectify_order(img, indices[i]['xidx'], indices[i]['yidx'])
        
        # Now get the wavelength solution to tack on

        bot = np.ceil(poly_1d(indices[i]['x'],
                              extract.state['edgecoeffs'][i, 0, :])).astype(int)

        w = extract.state['wavecal'][bot, np.fix(indices[i]['x']).astype(int)]

        order.update({'wavelength': w})
        order.update({'angle': indices[i]['y']})

        rectorders.append(order)

    # Store the results

    extract.state['rectorders'] = rectorders
    
    #
    # Do the plotting
    #

    order_plotinfo = {'xranges': extract.state['xranges'],
                      'edgecoeffs': extract.state['edgecoeffs'],
                      'orders': extract.state['orders']}

    if qa_show is True:
        number = plot_image(extract.state['workimage'],
                            orders_plotinfo=order_plotinfo,
                            plot_size=qa_showsize,
                            plot_number=extract.state['image_plotnum'])
        extract.state['image_plotnum'] = number

    if qa_write is True:
        qafileinfo = {'figsize': (7, 7),
                      'filepath': setup.state['qa_path'],
                      'filename': qafilename + '_image',
                      'extension': setup.state['qa_extension']}

        plot_image(extract.state['workimage'], orders_plotinfo=order_plotinfo,
                   file_info=qafileinfo)

    #
    # Set the done variables
    #

    extract.state['load_done'] = True
    extract.state['type_done'] = False
    extract.state['profiles_done'] = False
    extract.state['apertures_done'] = False
    extract.state['orders_done'] = False
    extract.state['trace_done'] = False
    extract.state['parameters_done'] = False
    extract.state['extract_done'] = False
