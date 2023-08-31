import numpy as np

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.trace_spectrum_1dxd import trace_spectrum_1dxd
from pyspextool.extract.trace_to_xy import trace_to_xy
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_image import plot_image


def trace_apertures(fit_degree=2, step_size=5, summation_width=5,
                    centroid_threshold=2, fwhm=0.8, verbose=None,
                    qa_show=None, qa_write=None, qa_showsize=(6, 6)):
    """
    Command line call for tracing apertures.

    Parameters
    ----------

    fit_degree : int, optional, default 2
        The polynomial degree for the fit.

    step_size : int, optional, default 5
        The step size as the function moves across the array identifying
        the peaks in the spatial dimension.

    summation_width : int, optional, default 5
        The number of columns to combine in order to increase the S/N
        of the data fit to identify the peaks in the spatial dimension.

    centroid_threshold : int, optional, default 2
        If (fit-guess) > centroid_threshold to peak identification is found
        to fail.

    fwhm : float, optional, default 0.8
        The approximate FWHM of the peaks in arcseconds.

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
    Updates config.state['tracecoeffs'] variable.  


    Notes
    -----
    Just collects information from the config.state variable, calls
    trace_spectrum_1dxd, and then optional plots the qa plot.


    Examples
    --------
    later

    """

    #
    # Check if we can proceed
    #

    if extract.state['select_done'] is False:

        message = "extract.state['select_done']=False.  "+\
          "Previous steps not complete."        
        raise ValueError(message)

    #
    # Check parameters
    #

    check_parameter('trace_apertures', 'fit_degree', fit_degree, 'int')

    check_parameter('trace_apertures', 'step_size', step_size, 'int')

    check_parameter('trace_apertures', 'summation_width', summation_width,
                    'int')

    check_parameter('trace_apertures', 'centroid_threshold',
                    centroid_threshold, 'int')

    check_parameter('trace_apertures', 'fwhm', fwhm, 'float')

    check_parameter('trace_apertures', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_showsize', qa_showsize, 'tuple')

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
    # Store user inputs
    #

    extract.trace['fitdeg'] = fit_degree
    extract.trace['step'] = step_size
    extract.trace['sumwidth'] = summation_width
    extract.trace['centresh'] = centroid_threshold
    extract.trace['fwhm'] = fwhm
    extract.trace['qafile'] = qa_write
    extract.trace['qaplot'] = qa_show
    extract.trace['qaplotsize'] = qa_showsize
    extract.trace['verbose'] = verbose

    #
    # Run the trace
    #

    if extract.state['type'] == 'ps':

        # Point source extraction.  Must actually trace the apertures

        doorders = extract.state['psdoorders']

        z = doorders == 1
        trace = trace_spectrum_1dxd(extract.state['workimage'],
                                    extract.state['ordermask'],
                                    extract.state['orders'][z],
                                    extract.state['wavecal'],
                                    extract.state['spatcal'],
                                    extract.state['xranges'][z, :],
                                    extract.state['apertures'][z, :],
                                    fit_degree=fit_degree,
                                    step_size=step_size,
                                    centroid_threshold=centroid_threshold,
                                    fwhm=fwhm, verbose=verbose)
        
        tracecoeffs = trace['coeffs']

        if qa_show is True or qa_write is True:
            norders, naps = np.shape(extract.state['apertures'])
            info = trace_to_xy(extract.state['ordermask'],
                               extract.state['wavecal'],
                               extract.state['spatcal'],
                               extract.state['xranges'],
                               extract.state['orders'], doorders, naps,
                               trace['coeffs'], verbose=verbose)

            plotinfo = {'x': trace['x'], 'y': trace['y'],
                        'goodbad': trace['goodbad'], 'fits': info}

    else:

        # Must be an extended source extraction

        # Build the trace coefficients from the aperture positions

        norders, naps = np.shape(extract.state['apertures'])

        doorders = extract.state['xsdoorders']
        z = doorders == 1
        ndoorders = np.sum(doorders)

        tracecoeffs = np.full((ndoorders * naps, 2), 0)
        tracecoeffs[:, 0] = np.ravel(extract.state['apertures'][z, :])

        # Do the plotting stuff

        if qa_show is True or qa_write is True:
            norders, naps = np.shape(extract.state['apertures'])
            info = trace_to_xy(extract.state['ordermask'],
                               extract.state['wavecal'],
                               extract.state['spatcal'],
                               extract.state['xranges'],
                               extract.state['orders'], doorders, naps,
                               tracecoeffs)

            plotinfo = {'fits': info}

    #
    # Store the results
    #

    extract.state['tracecoeffs'] = tracecoeffs

    #
    # Plot stuff up if requested
    #

    if qa_show is True:

        number = plot_image(extract.state['workimage'], trace_plotinfo=plotinfo,
                   plot_number=extract.state['image_plotnum'],
                   plot_size=qa_showsize)
        extract.state['image_plotnum'] = number

        
    if qa_write is True:
        qafileinfo = {'figsize': (7, 7),
                      'filepath': setup.state['qa_path'],
                      'filename': extract.state['qafilename'] + '_trace',
                      'extension': setup.state['qa_extension']}

        plot_image(extract.state['workimage'], trace_plotinfo=plotinfo,
                   file_info=qafileinfo)

    #
    # Set the done variable
    #

    extract.state['trace_done'] = True
