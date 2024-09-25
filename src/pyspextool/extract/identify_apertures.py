import numpy as np
import numpy.typing as npt
import logging
from os.path import join

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.profiles import find_peaks
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.pyspextoolerror import pySpextoolError


def identify_apertures(find_naps:int=None,
                       guess_positions:int | float | list=None,
                       aperture_positions:int | float | list=None,
                       fwhm:int | float=0.8,
                       verbose:bool=None,
                       qa_show:bool=None,
                       qa_showscale:float | int=None,
                       qa_showblock:bool=None,
                       qa_write:bool=None):
    
    """
    To determine the locations of spectral extraction apertures


    Parameters
    ----------
    find_naps : {1,2,3,4} default None
        If given, the function will automatically attempt to identify
        `find_naps` peaks in the profile.

    guess_positions : int, float, list
        If given, the function will use these positions in an attempt to
        identify peaks in the profile.

    aperture_positions : int, float, list
        If given, the function will use store these apertures.

    fwhm: float, default 0.8 (arcseconds).
        The approximate FWHM of the peak to be identified.  Only used 
        if `method` is 'auto' or 'guess'.

    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, default None
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    Returns
    -------
    None
        Stores results in:

        extract.state['aperture_positions']
        extract.state['aperture_signs']
        extract.state['naps']
        extract.state['apertures_done']
    
    """
    
    #
    # Check the load_done variable
    #

    if extract.state['profile_done'] is False:

        message = "Previous steps not complete.  Please run make_profiles.py."
        raise pySpextoolError(message)

    #
    # Check parameters and QA keywords 
    #

    check_parameter('identify_apertures', 'naps', find_naps, ['int','NoneType'],
                    possible_values=[1,2,3,4, None])

    check_parameter('identify_apertures', 'guess_positions', guess_positions,
                    ['int', 'float', 'list', 'NoneType'])    

    check_parameter('identify_apertures', 'aperture_positions',
                    aperture_positions, ['int', 'float', 'list', 'NoneType'])
    
    check_parameter('identify_apertures', 'fwhm', fwhm, ['int', 'float'])

    check_parameter('identify_apertures', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('identify_apertures', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('identify_apertures', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('identify_apertures', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('identify_apertures', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])


    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Make sure only one of `find_naps`, `guess_positions`, and
    # `aperture_positions` was passed.
    #

    test = int(np.sum([find_naps is None, guess_positions is None,
            aperture_positions is None]))

    if test == 3:

        message = 'Must use one of `find_naps`, `guess_positions`, and '+ \
            '`aperture_positions`. when calling identify_apertures.'
        raise pySpextoolError(message)
        
    if test < 2:

        message = 'Can only use one of `find_naps`, `guess_positions`, and '+ \
            '`aperture_positions` when calling identify_apertures.'
        raise pySpextoolError(message)

    #
    # Update the command line if requested
    #

    logging.info(' Locating the apertures.')

    #
    # Now do the search based on the user request
    #
    
    norders = len(extract.state['profiles'])

    if find_naps is not None:

        method = 'auto'        
        peaks = find_naps
        naps = find_naps
        extract.state['aperture_type'] = 'auto'
        
    if guess_positions is not None:

        method = 'guess'
        peaks = np.tile(guess_positions, (norders, 1))
        naps = np.shape(peaks)[1]
        extract.state['aperture_type'] = 'guess'        

    if aperture_positions is not None:

        method = 'fixed'
        peaks = np.tile(aperture_positions, (norders, 1))
        naps = np.shape(peaks)[1]
        extract.state['aperture_type'] = 'fixed'        
               
    apertures, apsigns = find_peaks(extract.state['profiles'],
                                    {'method': method, 'peaks': peaks},
                                    fwhm=fwhm)

    #
    # Determine the average apsign
    #

    average_apsign = np.sum(apsigns, axis=0) / np.sum(np.abs(apsigns), axis=0)
    apsigns = np.empty(naps, dtype=int)

    for i in range(naps):
        apsigns[i] = 1 if average_apsign[i] > 0 else -1

    #
    # Store the results into the config variable
    #

    extract.state['aperture_positions'] = apertures
    extract.state['aperture_signs'] = apsigns
    extract.state['naps'] = naps

    # Report the aperture signs
    
    signs = ', '.join(list(apsigns.astype(str)))
    signs = signs.replace('-1', '-')
    signs = signs.replace('1', '+')
    message = ' Aperture signs are (' + signs + ').'
    logging.info(message)
        
    #
    # Do the QA plotting
    #
    
    if qa['show'] is True:
               
        plot_profiles(extract.state['profiles'],
                      extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      aperture_positions=apertures,
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
                      aperture_positions=apertures,
                      profilestack_max=setup.plots['profilestack_max'],
                      profile_size=setup.plots['profile_size'],
                      font_size=setup.plots['font_size'],
                      output_fullpath=fullpath)

    #
    # Set continue variable
    #

    extract.state['apertures_done'] = True
