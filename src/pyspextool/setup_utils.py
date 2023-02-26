import os

from pyspextool.extract import config # sets initial state dictionary
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter
from pyspextool.io.check import check_path
from pyspextool.io.check import check_file
from pyspextool.utils.split_text import split_text
try:
    from importlib.resources import files # Python 3.10+
except ImportError:
    from importlib_resources import files # Python <=3.9


def pyspextool_setup(instrument=config.state['instruments'][0],
                     raw_path=None, cal_path=None, proc_path=None,
                     qa_path=None, verbose=True, qaextension='.pdf'):

    
    
    """
    Set the pyspextool instrument and paths

    Parameters
    ----------
    instrument : str, optional
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
    # We are not going to check parameters because the two routines we
    # call do so.
    #

    if instrument is not None:

        set_instrument(instrument)

    set_paths(raw_path=raw_path, cal_path=cal_path, proc_path=proc_path,
              qa_path=qa_path, verbose=verbose, qaextension=qaextension)


def set_paths(raw_path=None, cal_path=None, proc_path=None, qa_path=None,
              verbose=True, qaextension='.pdf'):

    """
    Set the pyspextool paths

    Parameters
    ----------
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

    check_parameter('set_paths', 'raw_path', raw_path, ['str', 'NoneType'])

    check_parameter('set_paths', 'cal_path', cal_path, ['str', 'NoneType'])

    check_parameter('set_paths', 'proc_path', proc_path, ['str', 'NoneType'])

    check_parameter('set_paths', 'qa_path', qa_path, ['str', 'NoneType'])

    check_parameter('set_paths', 'verbose', verbose, 'bool')    

    check_parameter('set_paths', 'qaextension', qaextension, 'str',
                    possible_values=config.state['qaextensions'])

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

    # Set the qa extension filetype

    config.state['qaextension'] = qaextension

    # Set the pscontinue and xscontinue variables

    config.state['pscontinue'] = 0
    config.state['xscontinue'] = 0

    if verbose is True:

        print()
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
        print('----------------')
        print()

def set_instrument(instrument_name):

    """
    Set the instrument.

    Parameters
    ----------
    instrument_name : str
        The name of the instrument.

    Returns
    -------
    None

    Examples
    --------
    later

    """

    if instrument_name is None:

        insturment_name = config.state['instruments'][0]

    #
    # Check parameter and store results
    #
        
    check_parameter('set_instrument', 'instrument_name', instrument_name,
                    'str')

    config.state['instrument_name'] = instrument_name
    
    #
    # Set the package path
    #

    config.state['package_path'] = str(files('pyspextool'))

    #
    # Check to make sure the instrument path exists
    #
    
    instrument_data_path = os.path.join(config.state['package_path'],
                                        'instrument_data',
                                        instrument_name+'_dir')

    check_path(instrument_data_path)

    config.state['instrument_path'] = instrument_data_path

    #
    # Now get the instrument file and load
    #

    instrument_info_file = os.path.join(instrument_data_path,
                                        instrument_name + '.dat')
    
    check_file(instrument_info_file)

    instrument_info = read_instrument_file(instrument_info_file)

    if instrument_name in ['uspex', 'spex']:
        config.state['irtf'] = True

    # Fill out the state variables
    
    config.state['suffix'] = instrument_info['SUFFIX']

    config.state['nint'] = instrument_info['NINT']

    config.state['xspextool_keywords'] = instrument_info['XSPEXTOOL_KEYWORDS']

    # Now store linearity numbers

    config.state['lincormax'] = instrument_info['LINCORMAX']
    config.state['linearity_info'] = {'max': config.state['lincormax'],
                                      'bit': 0}

    # Get the bad pixel mask
    
    bad_pixel_mask_file = os.path.join(instrument_data_path,
                                       config.state['instrument_name'] + \
                                       '_bdpxmk.fits')

    check_file(bad_pixel_mask_file)
    
    config.state['bad_pixel_mask_file'] = bad_pixel_mask_file
