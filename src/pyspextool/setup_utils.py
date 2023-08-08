import os
from astropy.io import fits
import numpy as np
import logging
from pyspextool import config as setup
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter, check_path, check_file
from pyspextool.plot.plot_image import plot_image
from importlib.resources import files  # Python 3.10+



logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


def pyspextool_setup(instrument=setup.state['instruments'][0], paths=None,
                      verbose=True, qa_extension=None,
                     qa_file=None, qa_plot=None):
    """
    Set the pyspextool instrument and paths

    Parameters
    ----------
    instrument : str, optional
        The name of the instrument.  Must be one of 
        config.setup['instruments'].
    
    paths : dict, optional
        A dictionary of paths.  Should contain the following keys:
            raw_path : str, optional.
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


    # Set up verbose scale and logging

    setup.state['verbose'] = verbose

    if verbose is True:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)


    # TODO: change to set_paths
    state = set_parameters(paths, verbose=verbose, qa_extension=qa_extension,
                   qa_file=qa_file, qa_plot=qa_plot)

    
    # TODO: make separate functions for qa state parameters.

    


    msg = f"""
    Pyspextool Setup
    ----------------
    Instrument: {setup.state['instrument']}

    Rawpath: {setup.state['raw_path']}
    Calpath: {setup.state['cal_path']}
    Procpath: {setup.state['proc_path']}
    Qapath: {setup.state['qa_path']}

    QA Extension: {setup.state['qa_extension']}
    QA Plot: {setup.state['qa_plot']}
    QA File: {setup.state['qa_file']}
    """

    logging.debug(msg)


def set_parameters(paths=None,
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

    check_parameter('set_parameters', 'raw_path', paths['raw_path'], ['str', 'NoneType'])

    check_parameter('set_parameters', 'cal_path', paths['cal_path'], ['str', 'NoneType'])

    check_parameter('set_parameters', 'proc_path', paths['proc_path'],
                    ['str', 'NoneType'])

    check_parameter('set_parameters', 'qa_path', paths['qa_path'], ['str', 'NoneType'])

    check_parameter('set_parameters', 'verbose', verbose, 'bool')

    check_parameter('set_parameters', 'qa_extension', qa_extension,
                    ['NoneType', 'str'],
                    possible_values=setup.state['qa_extensions'])

    check_parameter('set_parameters', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('set_parameters', 'qa_file', qa_file, ['NoneType', 'bool'])
    #
    # Load the .pyspextool file if it exists.
    #

    # Change to use the current working directory rather than HOME
    # home_path = os.path.expanduser('~')
    # file_name = '.pyspextool_' + setup.state['instrument'] + '.dat'
    # full_path = os.path.join(home_path, file_name)

    # Does the file exist?

    # if os.path.isfile(full_path) is True:

    #     # Yes.  Load the previous paths.

    #     f = open(full_path, 'r')
    #     paths = []

    #     for line in f:
    #         paths.append(line.strip())

    #     setup.state['raw_path'] = paths[0]
    #     setup.state['cal_path'] = paths[1]
    #     setup.state['proc_path'] = paths[2]
    #     setup.state['qa_path'] = paths[3]

    # else:

    #     # No.  Use the current working directory

    # Should only use the current working directory if the user has not defined the paths
    cwd = os.path.abspath(os.getcwd())

    # #
    # Now let's modify the paths based on the user requests.
    #

    if paths['raw_path'] is not None:
        paths['raw_path'] = check_path(paths['raw_path'], make_absolute=True)
        setup.state['raw_path'] = paths['raw_path']
        logging.debug(f"Set raw_path to {paths['raw_path']}")   
    else:
        setup.state['raw_path'] = cwd

    if paths['cal_path'] is not None:
        try:
            paths['cal_path'] = check_path(paths['cal_path'], make_absolute=True)
            logging.debug(f"Set cal_path to {paths['cal_path']}")
        except ValueError:
            os.mkdir(paths['cal_path'])
            paths['cal_path'] = check_path(paths['cal_path'], make_absolute=True)
            logging.info(f"Created cal_path directory {paths['cal_path']}")

        setup.state['cal_path'] = paths['cal_path']
    else:
        setup.state['cal_path'] = cwd

    if paths['proc_path'] is not None:
        paths['proc_path'] = check_path(paths['proc_path'], make_absolute=True)
        setup.state['proc_path'] = paths['proc_path']
    else:
        setup.state['proc_path'] = cwd

    if paths['qa_path'] is not None:
        try:
            paths['qa_path'] = check_path(paths['qa_path'], make_absolute=True)
            logging.debug(f"Set qa_path to {paths['qa_path']}")
        except ValueError:
            os.mkdir(paths['qa_path'])
            paths['qa_path'] = check_path(paths['qa_path'], make_absolute=True)
            logging.info(f"Created qa_path directory {paths['qa_path']}")    
            
        setup.state['qa_path'] = paths['qa_path']
    else:
        setup.state['qa_path'] = cwd

    #
    # Now write the paths to the user home directory
    #

    # dat_file_name = f".pyspextool_{setup.state['instrument']}.dat"
    # f = open(os.path.join(home_path, dat_file_name), 'w')
    # f.write('%s \n' % setup.state['raw_path'])
    # f.write('%s \n' % setup.state['cal_path'])
    # f.write('%s \n' % setup.state['proc_path'])
    # f.write('%s \n' % setup.state['qa_path'])
    # f.close()

    # logging.info(f'Created {dat_file_name} in {home_path}')

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

    return setup.state


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

    setup.state['raw_bad_pixel_mask'] = fits.getdata(bad_pixel_mask_file)

    #
    # Grab the Spextool keywords
    #

    keywords_path = os.path.join(data_path, 'pyspextool_keywords.dat')

    keywords = np.loadtxt(keywords_path, comments='#', dtype='str').tolist()

    setup.state['pyspextool_keywords'] = keywords

    logging.info("Instrument state set to somethnig else")

    return setup.state

