import os
from astropy.io import fits
import numpy as np

from pyspextool import config as setup
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter, check_path, check_file

try:
    from importlib.resources import files  # Python 3.10+
except ImportError:
    from importlib_resources import files  # Python <=3.9


def pyspextool_setup(instrument=setup.state['instruments'][0],
                     raw_path=None, cal_path=None, proc_path=None,
                     qa_path=None, verbose=True, qa_extension=None,
                     qa_file=None, qa_plot=None):
    """
    Set the pyspextool instrument and paths

    Parameters
    ----------
    instrument : str, optional
        The name of the instrument.  Must be one of 
        config.setup['instruments'].
    
    raw_path : str, optional
        The path to the raw directory.

    cal_path : str, optional
        The path to the calibration directory.

    proc_path : str, optional
        The path to the processed directory.

    qa_path : str, optional
        The path to the quality assurance directory.

    raw_path : str, optional
        The path to the quality assurance directory.

    verbose : bool, default = True
        Set to report the setup results.

    qa_extension : {None, True, False}, optional
        Set True/False to override setup.state['qa_extension']

    qa_file : {None, True, False}, optional
        Set True/False to override setup.state['qa_file']

    qa_plot : {None, True, False}, optional
        Set True/False to override setup.state['qa_plot']

    Returns
    -------
    None

    """

    #
    # We are not going to check parameters because the two routines we
    # call do so.
    #

    if instrument is not None:
        set_instrument(instrument)

    set_parameters(raw_path=raw_path, cal_path=cal_path, proc_path=proc_path,
                   qa_path=qa_path, verbose=verbose, qa_extension=qa_extension,
                   qa_file=qa_file, qa_plot=qa_plot)


def set_parameters(raw_path=None, cal_path=None, proc_path=None, qa_path=None,
                   verbose=True, qa_extension=None, qa_plot=None, qa_file=None):
    """
    Set the pyspextool parameters

    Parameters
    ----------
    raw_path : str, optional
        The path to the raw directory.

    cal_path : str, optional
        The path to the calibration directory.

    proc_path : str, optional
        The path to the processed directory.

    qa_path : str, optional
        The path to the quality assurance directory.

    verbose : {True, False}, optional
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file. 

    qa_extension : {None, 'pdf', 'png'}, optional
        Set to override setup.state['qa_extension']

    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.state['qa_file'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter('set_parameters', 'raw_path', raw_path, ['str', 'NoneType'])

    check_parameter('set_parameters', 'cal_path', cal_path, ['str', 'NoneType'])

    check_parameter('set_parameters', 'proc_path', proc_path,
                    ['str', 'NoneType'])

    check_parameter('set_parameters', 'qa_path', qa_path, ['str', 'NoneType'])

    check_parameter('set_parameters', 'verbose', verbose, 'bool')

    check_parameter('set_parameters', 'qa_extension', qa_extension,
                    ['NoneType', 'str'],
                    possible_values=setup.state['qa_extensions'])

    check_parameter('set_parameters', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('set_parameters', 'qa_file', qa_file, ['NoneType', 'bool'])
    #
    # Load the .pyspextool file if it exists.
    #

    home_path = os.path.expanduser('~')
    file_name = '.pyspextool_' + setup.state['instrument'] + '.dat'
    full_path = os.path.join(home_path, file_name)

    # Does the file exist?

    if os.path.isfile(full_path) is True:

        # Yes.  Load the previous paths.

        f = open(full_path, 'r')
        paths = []

        for line in f:
            paths.append(line.strip())

        setup.state['raw_path'] = paths[0]
        setup.state['cal_path'] = paths[1]
        setup.state['proc_path'] = paths[2]
        setup.state['qa_path'] = paths[3]

    else:

        # No.  Use the current working directory

        cwd = os.path.abspath(os.getcwd())
        setup.state['raw_path'] = cwd
        setup.state['cal_path'] = cwd
        setup.state['proc_path'] = cwd
        setup.state['qa_path'] = cwd

    #
    # Now let's modify the paths based on the user requests.
    #

    if raw_path is not None:
        raw_path = check_path(raw_path, make_absolute=True)
        setup.state['raw_path'] = raw_path

    if cal_path is not None:
        cal_path = check_path(cal_path, make_absolute=True)
        setup.state['cal_path'] = cal_path

    if proc_path is not None:
        proc_path = check_path(proc_path, make_absolute=True)
        setup.state['proc_path'] = proc_path

    if qa_path is not None:
        qa_path = check_path(qa_path, make_absolute=True)
        setup.state['qa_path'] = qa_path

    #
    # Now write the paths to the user home directory
    #

    f = open(os.path.join(home_path, '.pyspextool_' + \
                          setup.state['instrument'] + '.dat'), 'w')
    f.write('%s \n' % setup.state['raw_path'])
    f.write('%s \n' % setup.state['cal_path'])
    f.write('%s \n' % setup.state['proc_path'])
    f.write('%s \n' % setup.state['qa_path'])
    f.close()

    # Set the qa extension filetype

    setup.state['qa_extension'] = qa_extension

    #
    # Check defaults
    #

    if qa_file is not None:
        setup.state['qa_file'] = qa_file

    if qa_plot is not None:
        setup.state['qa_plot'] = qa_plot

    if qa_extension is not None:

        setup.state['qa_extension'] = qa_extension

    else:

        setup.state['qa_extension'] = setup.state['qa_extensions'][0]
            
    # Set the verbose

    setup.state['verbose'] = verbose

    if verbose is True:
        print()
        print('Pyspextool Setup')
        print('----------------')
        print('Instrument: ', setup.state['instrument'])
        print()
        print('Rawpath: ', setup.state['raw_path'])
        print('Calpath: ', setup.state['cal_path'])
        print('Procpath: ', setup.state['proc_path'])
        print('Qapath: ', setup.state['qa_path'])
        print()
        print('QA Extension:', setup.state['qa_extension'])
        print('QA Plot:', setup.state['qa_plot'])
        print('QA File:', setup.state['qa_file'], '\n')


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

    """

    if instrument_name is None:
        instrument_name = setup.state['instruments'][0]

    #
    # Check parameter and store results
    #

    check_parameter('set_instrument', 'instrument_name', instrument_name,
                    'str')

    setup.state['instrument'] = instrument_name

    #
    # Set the package path
    #

    setup.state['package_path'] = str(files('pyspextool'))

    #
    # Check to make sure the instrument path exists
    #

    instrument_data_path = os.path.join(setup.state['package_path'],
                                        'instrument_data',
                                        instrument_name + '_dir')
    data_path = os.path.join(setup.state['package_path'],
                                        'data')

    check_path(instrument_data_path)
    check_path(data_path)

    setup.state['instrument_path'] = instrument_data_path

    #
    # Now get the instrument file and load
    #

    instrument_info_file = os.path.join(instrument_data_path,
                                        instrument_name + '.dat')

    check_file(instrument_info_file)

    instrument_info = read_instrument_file(instrument_info_file)

    if instrument_name in ['uspex', 'spex']:
        setup.state['irtf'] = True

    # Fill out the state variables

    setup.state['suffix'] = instrument_info['SUFFIX']

    setup.state['nint'] = instrument_info['NINT']

    setup.state['extract_keywords'] = instrument_info['XSPEXTOOL_KEYWORDS']

    setup.state['combine_keywords'] = instrument_info['COMBINE_KEYWORDS']

    setup.state['telluric_keywords'] = instrument_info['TELLURIC_KEYWORDS']   

    # Now store linearity numbers

    setup.state['lincormax'] = instrument_info['LINCORMAX']
    setup.state['linearity_info'] = {'max': setup.state['lincormax'],
                                     'bit': 0}

    # Get the bad pixel mask

    bad_pixel_mask_file = os.path.join(instrument_data_path,
                                       setup.state['instrument'] + \
                                       '_bdpxmk.fits')

    check_file(bad_pixel_mask_file)

    # setup.state['raw_bad_pixel_mask'] = fits.getdata(bad_pixel_mask_file)
    # accessing the file breaks remote tests

    #
    # Grab the Spextool keywords
    #

    keywords_path = os.path.join(data_path, 'pyspextool_keywords.dat')

    keywords = np.loadtxt(keywords_path, comments='#', dtype='str').tolist()

    setup.state['pyspextool_keywords'] = keywords
