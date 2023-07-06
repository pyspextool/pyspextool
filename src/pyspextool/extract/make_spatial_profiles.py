import numpy as np
from scipy import interpolate

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.extract.profiles import make_1d_profile


def make_spatial_profiles(qa_plot=None, qa_file=None, qa_plotsize=(6, 10),
                          verbose=None):
    """
    To create 1D "average" spatial profiles of the orders.

    Parameters 
    ----------
    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(6, 10)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.state['qa_file'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    Returns 
    -------
    None
    Fills the config.extract['profiles'] variable and optionally creates a 
    QA file.

    """

    #
    # Check to make sure we can proceed.
    #

    if extract.state['type_done'] is False:
        message = 'Previous steps not completed.'
        print(message)
        return

    #
    # Check the parameters
    #

    check_parameter('make_spatial_profiles', 'qa_plot', qa_plot,
                    ['NoneType', 'bool'])

    check_parameter('make_spatial_profiles', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('make_spatial_profiles', 'qa_file', qa_file,
                    ['NoneType', 'bool'])

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_file is None:
        qa_file = setup.state['qa_file']

    if qa_plot is None:
        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Save user inputs
    #

    extract.profiles['qaplot'] = qa_plot
    extract.profiles['qafile'] = qa_file
    extract.profiles['verbose'] = verbose

    #
    # Build the profiles
    #

    if verbose is True:
        print('Creating the spatial profiles...')

    profiles = []
    for i in range(extract.state['norders']):

        atmos = extract.state['atmosphere']
        
        # Create the profile
        
        y, mean = make_1d_profile(extract.state['rectorders'][i],
                                  atmospheric_transmission=atmos)

        profiles.append({'order': extract.state['orders'][i], 'angle': y,
                         'profile': mean})

    # Store the results
        
    extract.state['profiles'] = profiles

    # Do the QA plotting
    
    if qa_plot is True:
        number = plot_profiles(extract.state['profiles'],
                               extract.state['slith_arc'],
                               np.ones(extract.state['norders'], dtype=int),
                               plot_size=qa_plotsize,
                               plot_number=extract.state['profiles_plotnum'])
        extract.state['profiles_plotnum'] = number
        
        
    if qa_file is True:
        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': setup.state['qa_path'],
                      'filename': extract.state['qafilename'] + '_profiles',
                      'extension': setup.state['qa_extension']}

        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      file_info=qafileinfo)

    # Set the done variable

    extract.state['profile_done'] = True
