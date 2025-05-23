import numpy as np
import logging
from os.path import join

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.extract.profiles import make_1d_profile

def make_profiles(verbose:bool=None,
                  qa_show:bool=None,
                  qa_showscale:float | int=None,
                  qa_showblock:bool=None,
                  qa_write:bool=None):

    """
    To create 1D "average" spatial profiles of the orders.

    Parameters 
    ----------

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot to the screen.
        Set to False to not show a QA plot to the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_showblock'].
    
    qa_showscale : float or int, default None, 
        The scale factor by which to increase or decrease the default size.
        If None, then it defaults to setup.state['qa_showscale'].

    
    Returns 
    -------
    None
    Fills the config.extract['profiles'] variable and optionally creates a 
    QA file.

    """

    #
    # Check to make sure we can proceed.
    #

    if extract.state['load_done'] is False:

        message = 'Previous steps complete.  Please run extract.load_images.'
        raise ValueError(message)

    #
    # Check the parameters and QA parameters
    #

    
    check_parameter('make_profiles', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('make_profiles', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('make_profiles', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('make_profiles', 'qa_showscale', qa_showscale,
                    ['NoneType', 'int', 'float'])

    check_parameter('make_profiles', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    #
    # Build the profiles
    #

    logging.info(' Creating the 1D spatial profiles.')

    profiles = []
    for i in range(extract.state['norders']):

        atmos = extract.state['atmosphere']
        
        # Create the profile

        result = make_1d_profile(extract.state['rectorders'][i],
                                 atmospheric_transmission=atmos)
        
        angles = result[0]
        profile = result[1]
       
        profiles.append({'order': extract.state['orders'][i], 'angles': angles,
                         'profile': profile})

    # Store the results
        
    extract.state['profiles'] = profiles

    #
    # Do the QA plotting
    #
    
    if qa['show'] is True:
               
        plot_profiles(extract.state['profiles'],
                      extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      plot_number=setup.plots['profiles'],
                      profilestack_max=setup.plots['profilestack_max'],
                      profile_size=setup.plots['profile_size'],
                      font_size=setup.plots['font_size'],
                      showscale=qa['showscale'],
                      showblock=qa['showblock'])

    if qa['write'] is True:

        filename = extract.state['qafilename'] + '_profiles' + \
            setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

        plot_profiles(extract.state['profiles'],
                      extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      profilestack_max=setup.plots['profilestack_max'],
                      profile_size=setup.plots['profile_size'],
                      font_size=setup.plots['font_size'],
                      output_fullpath=fullpath)

    # Set the done variable

    extract.state['profile_done'] = True

