import numpy as np

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.find_peaks import find_peaks
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_profiles import plot_profiles


def locate_aperture_positions(apertures, method='auto', fwhm=0.8, qa_show=None,
                              qa_write=None, qa_showsize=(6, 10), verbose=None):
    """
    To determine the locations of spectral extraction apertures


    Parameters
    ----------
    apertures : int, float, list, or numpy.ndarray
        The exact type depends on the extraction type and method:

        Extended Source Extraction:
        An int, float, list, or numpy.ndarray of aperture positions.

        Point Source Extraction:
        `method` = 'auto'
            The int number of apertures to search
        `method` = 'guess'
            An int, float, list, numpy.darray of apertures positions.
        `method` = 'fixed'
            An int, float, list, numpy.darray of apertures positions.

    method : {'auto', 'guess', fixed'}, optional 
        The method by which point source apertures are identified.  It has 
        no impact if the extraction is of an extended source.

    fwhm: float, default 0.8 (arcseconds).
        The approximate FWHM of the peak to be identified.  Only used 
        if `method` is 'auto' or 'guess'.

    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, default=(6,6)
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
    Sets config.extract['apertures'], config.extract['apsigns'], and 
    config.extract['naps'] and if config.extract['exttype'] ='xs' then 
    also calls trace_apertures.
    
    """

    #
    # Check the load_done variable
    #

    if extract.state['profile_done'] is False:

        message = "extract.state['profile_done']=False.  "+\
          "Previous steps not complete."        
        raise ValueError(message)

    #
    # Check parameters that don't depend on extraction type.
    #

    check_parameter('define_aperture_positions', 'apertures', apertures,
                    ['int', 'float', 'list', 'ndarray'])

    check_parameter('define_aperture_positions', 'method', method, 'str',
                    possible_values=['auto', 'fixed', 'guess'])

    check_parameter('define_aperture_positions', 'fwhm', fwhm, ['int', 'float'])

    check_parameter('define_aperture_positions', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('define_aperture_positions', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('define_aperture_positions', 'qa_showsize',
                    qa_showsize, 'tuple')

    check_parameter('load_image', 'verbose', verbose, ['NoneType', 'bool'])

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
    # Store the user passed parameters 
    #

    extract.apertures['apertures'] = apertures
    extract.apertures['method'] = method
    extract.apertures['fwhm'] = fwhm
    extract.apertures['qaplot'] = qa_show
    extract.apertures['qafile'] = qa_write
    extract.apertures['qaplotsize'] = qa_showsize
    extract.apertures['verbose'] = verbose

    #
    # Update the command line if requested
    #

    if verbose is True:
        print('Locating the apertures...')

    #
    # Get useful numbers
    #

    norders = len(extract.state['profiles'])

    #
    # Check the extraction type
    #

    if extract.state['type'] == 'ps':

        #
        # This is a point source extraction
        #

        #
        # Which method was requested?
        #

        if method == 'auto':
            check_parameter('define_aperture_positions', 'apertures',
                            apertures, 'int')

            naps = apertures
            apertures, apsigns = find_peaks(extract.state['profiles'],
                                            {'method': 'auto', 'peaks': naps},
                                            fwhm=fwhm)

        if method == 'guess':
            check_parameter('define_aperture_positions', 'apertures',
                            apertures, ['int', 'float', 'list', 'ndarray'])

            apertures = np.tile(apertures, (norders, 1))

            apertures, apsigns = find_peaks(extract.state['profiles'],
                                            {'method': 'guess',
                                             'peaks': apertures},
                                            fwhm=fwhm)

            naps = int(np.size(apertures) / norders)

        if method == 'fixed':
            check_parameter('define_aperture_positions', 'apertures',
                            apertures, ['int', 'float', 'list', 'ndarray'])

            apertures = np.tile(apertures, (norders, 1))

            apertures, apsigns = find_peaks(extract.state['profiles'],
                                            {'method': 'fixed',
                                             'peaks': apertures},
                                            fwhm=fwhm)

            naps = int(np.size(apertures) / norders)

    else:

        #
        # This is an extended source extraction
        #

        try:

            naps = len(apertures)

        except TypeError:

            naps = 1

        apertures = np.tile(apertures, (norders, 1))
        apsigns = np.full_like(apertures, 1, dtype=int)

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

    extract.state['apertures'] = apertures
    extract.state['apsigns'] = apsigns
    extract.state['naps'] = naps

    if verbose is True:
        signs = ', '.join(list(apsigns.astype(str)))
        signs = signs.replace('-1', '-')
        signs = signs.replace('1', '+')
        message = 'Aperture signs are (' + signs + ')...'
        print(message)

    #
    # Plot the results
    #

    if qa_show is True:
        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      apertures=extract.state['apertures'],
                      plot_number=extract.state['profiles_plotnum'],
                      plot_size=qa_showsize)

    if qa_write is True:
        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': setup.state['qa_path'],
                      'filename': extract.state['qafilename'] +
                                  '_aperturepositions',
                      'extension': setup.state['qa_extension']}

        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      apertures=extract.state['apertures'],
                      file_info=qafileinfo)

    #
    # Set continue variable
    #

    extract.state['apertures_done'] = True

    #
    # Now run the trace if the source is extended
    #

    # if extract.state['type'] == 'xs':
    #    trace_apertures(verbose=verbose, qaplot=qaplot, qafile=qafile)
