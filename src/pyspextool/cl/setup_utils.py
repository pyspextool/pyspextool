import os

from pyspextool.cl import config # sets initial state dictionary
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter
from pyspextool.io.check import check_path
try:
    from importlib.resources import files # Python 3.10+
except ImportError:
    from importlib_resources import files # Python <=3.9


def setup_paths(instrument_name=None, raw_path=None, cal_path=None, proc_path=None,
          qa_path=None, verbose=True, qaextension='.pdf'):

    """
    Set up basic information for pyspextool to run.

    Parameters
    ----------
    instrument : str, default = 'uspex', optional
        The name of the instrument.

    raw_path : str, optional
        The path to the raw directory.

    cal_path : str, optional
        The path to the calibration directory.

    proc_path : str, optional
        The path to the processed directory.

    raw_path : str, optional
        The path to the quality assurance directory.

    verbose : bool, default = True
        Set to report the setup results.

    qaextension : {'pdf', 'png'}, optional
        Set the file type for all QA plots.

    Returns
    -------
    None

    Examples
    --------
    later

    """

    #
    # Check parameters
    #

    if instrument_name is None:
        instrument_name = config.state['instruments'][0]

    # if not isinstance(instrument_name, str):
    #    message = f'Instrument name, {instrument_name} should be a string not {type(instrument_name)}.' \
    #             f'Possible instruments are: ' \
    #             f"{str.join(', ', config.state['instruments'])}."
    #    raise TypeError(message)
    # check_parameter('setup', 'instrument_name', instrument_name, ['str', 'NoneType'])

    check_parameter('setup', 'rawpath', raw_path, ['str', 'NoneType'])

    check_parameter('setup', 'calpath', cal_path, ['str', 'NoneType'])

    check_parameter('setup', 'procpath', proc_path, ['str', 'NoneType'])

    check_parameter('setup', 'qapath', qa_path, ['str', 'NoneType'])

    check_parameter('setup', 'qaextension', qaextension, 'str',
                    possible_values=['.pdf','.png'])

    #
    # Set the package path
    #

    config.state['package_path'] = str(files('pyspextool'))

    #
    # Set the instrument and associated things
    #
    
    set_instrument_state(instrument_name)

    #
    # Now deal with the reduction paths.
    #

    set_reduction_paths(raw_path, cal_path, proc_path, qa_path)

    #
    # Store and set things
    #

    # the qa extension filetype

    config.state['qaextension'] = qaextension

    # Set the pscontinue and xscontinue variables

    config.state['pscontinue'] = 0
    config.state['xscontinue'] = 0

    if verbose is True:

        print('Pyspextool Setup')
        print('----------------')
        print('Instrument: ', config.state['instrument_name'])
        print()
        print('Rawpath: ', config.state['rawpath'])
        print('Calpath: ', config.state['calpath'])
        print('Procpath: ', config.state['procpath'])
        print('Qapath: ', config.state['qapath'])
        print()
        print('QA Extension:', config.state['qaextension'])


def set_instrument_state(instrument_name: str):
    
    if not isinstance(instrument_name, str):
        message = f"Instrument name is not a string: {instrument_name}, {type(instrument_name)}"
        return TypeError(message)

    instrument_data_path = os.path.join(config.state['package_path'],
                                        'instrument_data',
                                        instrument_name+'_dir')

    if os.path.isdir(instrument_data_path) is True:
        config.state['instrument_name'] = instrument_name
        config.state['instrument_path'] = instrument_data_path
    else:
        message = f"Needed instrument data for {instrument_name} is not found. Expecting it here: \n" \
                  f"{instrument_data_path}."
        raise FileNotFoundError(message)

    instrument_info_file = os.path.join(instrument_data_path, instrument_name + '.dat')
    #instrument_info_file = instrument_data_path.joinpath(instrument_name + '.dat')

    if os.path.isfile(instrument_info_file):
        instrument_info = read_instrument_file(instrument_info_file)
    else:
        raise FileNotFoundError(instrument_info_file)

    if instrument_name in ['uspex', 'spex']:
        config.state['irtf'] = True

    # Fill out the state variables
    
    config.state['readfits'] = instrument_info['READFITS']

    config.state['suffix'] = instrument_info['SUFFIX']

    config.state['nint'] = instrument_info['NINT']

    config.state['xspextool_keywords'] = instrument_info['XSPEXTOOL_KEYWORDS']

    bias_file = os.path.join(instrument_data_path,
                             config.state['instrument_name'] + '_bias.fits')
    #bias_file = instrument_data_path.joinpath(config.state['instrument_name'] + '_bias.fits')

    if os.path.isfile(bias_file):
        config.state['biasfile'] = bias_file
        config.state['lincormax'] = instrument_info['LINCORMAX']
        config.state['linearity_info'] = {'bias': bias_file,
                                          'max': config.state['lincormax'],
                                          'bit': 0}
    else:
        raise FileNotFoundError(bias_file)

    bad_pixel_mask_file = os.path.join(instrument_data_path,
                                       config.state['instrument_name'] + '_bdpxmk.fits')
    #bad_pixel_mask_file = instrument_data_path.joinpath(config.state['instrument_name'] + '_bdpxmk.fits')
    if os.path.isfile(bad_pixel_mask_file):
        config.state['bad_pixel_mask_file'] = bad_pixel_mask_file
    else:
        raise FileNotFoundError(bad_pixel_mask_file)


def set_reduction_paths(raw_path, cal_path, proc_path, qa_path):

    #
    # Load the .pyspextool file if it exists.
    #

    home_path = os.path.expanduser('~')
    file_name = '.pyspextool_' + config.state['instrument_name'] + '.dat'
    full_path = os.path.join(home_path, file_name)

    # Does the file exist?
    
    if os.path.isfile(full_path) is True:

        # Yes.  Load the previous paths.
        
        f = open(full_path, 'r')
        paths = []

        for line in f:
            paths.append(line.strip())

        config.state['rawpath'] = paths[0]
        config.state['calpath'] = paths[1]
        config.state['procpath'] = paths[2]
        config.state['qapath'] = paths[3]

    else:

        # No.  Use the current working directory
        
        cwd = os.path.abspath(os.getcwd())
        config.state['rawpath'] = cwd
        config.state['calpath'] = cwd
        config.state['procpath'] = cwd
        config.state['qapath'] = cwd

    #
    # Now let's modify the paths based on the user requests.
    #

    if raw_path is not None:

        raw_path = check_path(raw_path, make_absolute=True)
        config.state['rawpath'] = raw_path

    if cal_path is not None:

        cal_path = check_path(cal_path, make_absolute=True)
        config.state['calpath'] = cal_path

    if proc_path is not None:

        proc_path = check_path(proc_path, make_absolute=True)
        config.state['procpath'] = proc_path

    if qa_path is not None:

        qa_path = check_path(qa_path, make_absolute=True)
        config.state['qapath'] = qa_path

    #
    # Now write the paths to the user home directory
    #

    f = open(os.path.join(home_path, '.pyspextool_' + \
                          config.state['instrument_name'] + '.dat'), 'w')
    f.write('%s \n' % config.state['rawpath'])
    f.write('%s \n' % config.state['calpath'])
    f.write('%s \n' % config.state['procpath'])
    f.write('%s \n' % config.state['qapath'])
    f.close()