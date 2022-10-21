import numpy as np

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.cl.trace_apertures import trace_apertures
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.spectroscopy.find_peaks import find_peaks


def locate_aperture_positions(apertures, method='auto', fwhm=0.8, iplot=True,
                              qafile=False, clupdate=True):
    """
    To determine the locations of spectral extraction apertures


    Parameters
    ----------
    apertures : int, float, list, or numpy.ndarray
        The exact type depends on the extraction type and method.

        Extended Source Extraction:
        An int, float, list, of numpy.ndarray of aperture positions.

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

    iplot : {True, False}, optional

    qafile : {False, True}, optional

    clupdate: {True, False}, optional
        Set to True to report the operation to the command line.

    Returns
    -------
    None

    Notes
    -----
    Just does some organizing and calls find_peaks.  

    Examples
    --------
    later
    
    """

    #
    # Update the command line if requested
    #

    if clupdate is True:
        print('Locating the apertures...')

    
    #
    # Check continue variable
    #

    check_continue(3)

    #
    # Check parameters that don't depend on extraction type.
    #

    check_parameter('define_aperture_positions', 'method', method, 'str',
                    possible_values=['auto', 'fixed', 'guess'])

    check_parameter('define_aperture_positions', 'fwhm', fwhm, ['int', 'float'])

    check_parameter('define_aperture_positions', 'iplot', iplot, 'bool')

    check_parameter('define_aperture_positions', 'qafile', qafile, 'bool')


    #
    # Get useful numbers
    #

    norders = len(config.state['profiles'])

    #
    # Check the extraction type
    #

    if config.state['exttype'] == 'ps':

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
            apertures, apsigns = find_peaks(config.state['profiles'],
                                            {'method': 'auto', 'peaks': naps},
                                            fwhm=fwhm)

        if method == 'guess':
            check_parameter('define_aperture_positions', 'apertures',
                            apertures, ['int', 'float', 'list', 'ndarray'])

            apertures = np.tile(apertures, (norders, 1))

            apertures, apsigns = find_peaks(config.state['profiles'],
                                            {'method': 'guess', 'peaks': apertures},
                                            fwhm=fwhm)

            naps = int(np.size(apertures) / norders)

        if method == 'fixed':
            check_parameter('define_aperture_positions', 'apertures',
                            apertures, ['int', 'float', 'list', 'ndarray'])

            apertures = np.tile(apertures, (norders, 1))

            apertures, apsigns = find_peaks(config.state['profiles'],
                                            {'method': 'fixed', 'peaks': apertures},
                                            fwhm=fwhm)

            naps = int(np.size(apertures) / norders)

            #
        # Determine the average apsign
        #

        if norders > 1:

            average_apsign = np.sum(apsigns, axis=0) / np.sum(np.abs(apsigns),
                                                              axis=0)
            apsigns = np.empty(naps, dtype=int)

            for i in range(naps):
                apsigns[i] = 1 if average_apsign[i] > 0 else -1

        else:

            apsigns = np.squeeze(apsigns)

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
    # Store the results into the config variable
    #

    config.state['apertures'] = apertures
    config.state['apsigns'] = apsigns
    config.state['naps'] = naps

    #
    # Set continue variable
    #

    if config.state['exttype'] == 'ps':

        config.state['continue'] = 4

    else:

        config.state['continue'] = 4
        
    #
    # Plot the results
    #

    if iplot is True:
        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      np.ones(config.state['norders'], dtype=int),
                      apertures=config.state['apertures'])

    if qafile is True:
        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
                      'filename': config.state['qafilename'] + \
                                  '_aperturepositions', 'extension': '.pdf'}

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      np.ones(config.state['norders'], dtype=int),
                      apertures=config.state['apertures'],
                      qafileinfo=qafileinfo)


    #
    # Now run the trace if the source is extended
    #
        
    if config.state['exttype'] == 'xs':

        trace_apertures(clupdate=clupdate, iplot=iplot, qafile=qafile)
        
