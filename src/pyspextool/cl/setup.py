import os
from astropy.io import fits
from pyspextool.cl import config # sets initial state dictionary
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter
try:
    from importlib.resources import files # Python 3.10+
except ImportError:
    from importlib_resources import files # Python <=3.9


def setup(instrument_name='uspex', rawpath=None, calpath=None, procpath=None, qapath=None,
          clupdate=False, qaextension='.pdf'):

    """
    Set up basic information for pyspextool to run.

    Parameters
    ----------
    instrument : str, default = 'uspex', optional
        The name of the instrument.
    rawpath : str, optional
        The path to the raw directory.
    calpath : str, optional
        The path to the calibration directory.
    procpath : str, optional
        The path to the processed directory.
    rawpath : str, optional
        The path to the quality assurance directory.
    clupdate : bool, default = False
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

    # if not isinstance(instrument_name, str):
    #    message = f'Instrument name, {instrument_name} should be a string not {type(instrument_name)}.' \
    #              f'Possible instruments are: ' \
    #              f"{str.join(', ', config.state['instruments'])}."
    #    raise TypeError(message)
    # check_parameter('setup', 'instrument_name', instrument_name, ['str', 'NoneType'])

    check_parameter('setup', 'rawpath', rawpath, ['str', 'NoneType'])

    check_parameter('setup', 'calpath', calpath, ['str', 'NoneType'])

    check_parameter('setup', 'procpath', procpath, ['str', 'NoneType'])

    check_parameter('setup', 'qapath', qapath, ['str', 'NoneType'])

    check_parameter('setup', 'qaextension', qaextension, 'str')

    set_instrument_state(instrument_name)

    # Now do the paths.  First we check for the .path file in the user's home
    # directory

    homedir = os.path.expanduser('~')
    filename = os.path.join(homedir, '.pyspextool_' + \
                            config.state['instrument_name'] + '.dat')
    filename = os.path.join(homedir, filename)

    if os.path.isfile(filename) is True:
        f = open(filename, 'r')
        paths = []
        for line in f:
            paths.append(line.strip())
        config.state['rawpath'] = paths[0]
        config.state['calpath'] = paths[1]
        config.state['procpath'] = paths[2]
        config.state['qapath'] = paths[3]
    else:
        cwd = os.path.abspath(os.getcwd())
        config.state['rawpath'] = cwd
        config.state['calpath'] = cwd
        config.state['procpath'] = cwd
        config.state['qapath'] = cwd

    # Now let's modify the paths based on the user

    # Is rawpath passed?

    if rawpath is not None:

        # deal with ~/

        rawpath = os.path.expanduser(rawpath)

        # is the path real?

        if os.path.isdir(rawpath) is True:

            # get the absolute path 

            config.state['rawpath'] = os.path.abspath(rawpath)

        else:

            # Path is not real.

            message = '`rawpath` not a directory.'
            raise ValueError(message)

    # Is calpath passed?

    if calpath is not None:

        # deal with ~/

        calpath = os.path.expanduser(calpath)

        # is the path real?

        if os.path.isdir(calpath) is True:

            # get the absolute path 

            config.state['calpath'] = os.path.abspath(calpath)

        else:

            # Path is not real.

            message = '`calpath` not a directory.'
            raise ValueError(message)

    # Is procpath passed?

    if procpath is not None:

        # deal with ~/

        procpath = os.path.expanduser(procpath)

        # is the path real?

        if os.path.isdir(procpath) is True:

            # get the absolute path 

            config.state['procpath'] = os.path.abspath(procpath)

        else:

            # Path is not real.

            message = '`procpath` not a directory.'
            raise ValueError(message)

    # Is qapath passed?

    if qapath is not None:

        # deal with ~/

        qapath = os.path.expanduser(qapath)

        # is the path real?

        if os.path.isdir(qapath) is True:

            # get the absolute path 

            config.state['qapath'] = os.path.abspath(qapath)

        else:

            # Path is not real.

            message = '`qapath` not a directory.'
            raise ValueError(message)

    # Now write the paths to the user home directory

    f = open(os.path.join(homedir, '.pyspextool_' + \
                          config.state['instrument_name'] + '.dat'), 'w')
    f.write('%s \n' % config.state['rawpath'])
    f.write('%s \n' % config.state['calpath'])
    f.write('%s \n' % config.state['procpath'])
    f.write('%s \n' % config.state['qapath'])
    f.close()

    if clupdate is True:
        print('Pyspextool Setup')
        print('----------------')
        print('Instrument: ', config.state['instrument_name'])
        print()
        print('rawpath: ', config.state['rawpath'])
        print('calpath: ', config.state['calpath'])
        print('procpath: ', config.state['procpath'])
        print('qapath: ', config.state['qapath'])

    # Now store things
    config.state['qaextension'] = qaextension
        
    # Set the continue variables
    config.state['pscontinue'] = 0
    config.state['xscontinue'] = 0    


def set_instrument_state(instrument_name: str):
    if not isinstance(instrument_name, str):
        message = f"{instrument_name}, {type(instrument_name)} boo"
        return TypeError(message)

    instrument_data_path = files('pyspextool.instrument_data').\
        joinpath(instrument_name + '_dir')

    if os.path.isdir(instrument_data_path) is True:
        config.state['instrument_name'] = instrument_name
    else:
        print(instrument_data_path)
        print(config.state['instrument_name'])
        print(instrument_name)
        message = f"Unknown instrument.  Possible instruments are:" \
                  f"{str.join(', ', config.state['instruments'])}."
        raise ValueError(message)

    instrument_info_file = instrument_data_path.joinpath(config.state['instrument_name'] + '.dat')

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

    bias_file = instrument_data_path.joinpath(config.state['instrument_name'] + '_bias.fits')
    if os.path.isfile(bias_file):
        config.state['biasfile'] = bias_file
        config.state['lincormax'] = instrument_info['LINCORMAX']
        config.state['linearity_info'] = {'bias': bias_file,
                                          'max': config.state['lincormax'],
                                          'bit': 0}
    else:
        raise FileNotFoundError(bias_file)

    bad_pixel_mask_file = instrument_data_path.joinpath(config.state['instrument_name'] + '_bdpxmk.fits')
    if os.path.isfile(bad_pixel_mask_file):
        config.state['rawbadpixelmask'] = fits.getdata(bad_pixel_mask_file)
    else:
        raise FileNotFoundError(bad_pixel_mask_file)

    print(config.state['lincormax'])