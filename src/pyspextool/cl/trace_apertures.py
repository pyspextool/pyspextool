import numpy as np

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_image import plot_image
from pyspextool.spectroscopy.trace_spectrum_1dxd import trace_spectrum_1dxd
from pyspextool.spectroscopy.trace_to_xy import trace_to_xy


def trace_apertures(fit_degree=2, step_size=5, summation_width=5,
                    centroid_threshold=2, fwhm=0.8, clupdate=True,
                    iplot=True, qafile=False):
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

    clupdate : {True, False}, optional
        Set to True for command line updates during execution.

    iplot : {True, False}, optional
        Set to True to plot the qa plot interactively.

    qafile : {False, True}, optional
        Set to True to write the qa plot to disk.

    Returns
    -------
    None

    Notes
    -----
    Just collects information from the config.state variable, calls
    trace_spectrum_1dxd, and then optional plots the qa plot.


    Examples
    --------
    later

    """

    #
    # Check continue variable
    #
    check_continue(4)

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

    check_parameter('trace_apertures', 'clupdate', clupdate, 'bool')

    check_parameter('trace_apertures', 'iplot', iplot, 'bool')

    check_parameter('trace_apertures', 'qafile', qafile, 'bool')

    #
    # Run the trace
    #

    if config.state['exttype'] == 'ps':

        # Point source extraction.  Must actually trace the apertures

        doorders = config.state['psdoorders']

        z = doorders == 1
        trace = trace_spectrum_1dxd(config.state['workimage'],
                                    config.state['ordermask'],
                                    config.state['orders'][z],
                                    config.state['wavecal'],
                                    config.state['spatcal'],
                                    config.state['xranges'][z, :],
                                    config.state['apertures'][z, :],
                                    fit_degree=fit_degree,
                                    step_size=step_size,
                                    centroid_threshold=centroid_threshold,
                                    fwhm=fwhm, clupdate=clupdate)

        if iplot is True or qafile is True:
            norders, naps = np.shape(config.state['apertures'])
            info = trace_to_xy(config.state['ordermask'],
                               config.state['wavecal'],
                               config.state['spatcal'],
                               config.state['xranges'],
                               config.state['orders'], doorders, naps,
                               trace['coeffs'])

            plotinfo = {'x': trace['x'], 'y': trace['y'],
                        'goodbad': trace['goodbad'], 'fits': info}

    else:

        # Must be an extended source extraction

        # Build the trace coefficients from the aperture positions
                
        norders, naps = np.shape(config.state['apertures'])        
        
        doorders = config.state['xsdoorders']
        z = doorders == 1
        ndoorders = np.sum(doorders)        

        coeffs = np.full((ndoorders*naps,2),0)
        coeffs[:,0] = np.ravel(config.state['apertures'][z,:])

        # Do do the plotting stuff
            
        if iplot is True or qafile is True:
            norders, naps = np.shape(config.state['apertures'])
            info = trace_to_xy(config.state['ordermask'],
                               config.state['wavecal'],
                               config.state['spatcal'],
                               config.state['xranges'],
                               config.state['orders'], doorders, naps,
                               coeffs)

            plotinfo = {'fits': info}

    #
    # Plot stuff up if requested
    #

            
    if iplot is True:
        plot_image(config.state['workimage'], trace_plotinfo=plotinfo)

    if qafile is True:
        qafileinfo = {'figsize': (7, 7), 'filepath': config.state['qapath'],
                      'filename': config.state['qafilename'] + '_trace',
                      'extension': '.pdf'}

        plot_image(config.state['workimage'], trace_plotinfo=plotinfo,
                   qafileinfo=qafileinfo)

    #
    # Set the continue flag
    #

    config.state['continue'] = 5
