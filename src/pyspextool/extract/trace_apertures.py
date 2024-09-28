import numpy as np
import logging
from os.path import join

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.trace import trace_spectrum_1dxd, trace_to_xy
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.plot.plot_image import plot_image
from pyspextool.pyspextoolerror import pySpextoolError

def trace_apertures(fit_degree:int=2,
                    step_size:int=5,
                    summation_width:int=5,
                    centroid_threshold:int=2,
                    seeing_fwhm:int | float=0.8,
                    verbose:bool=None,
                    qa_show:bool=None,
                    qa_showscale:float | int=None,
                    qa_showblock:bool=None,
                    qa_write:bool=None):

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

    seeing_fwhm : float, optional, default 0.8
        The approximate FWHM of the peaks in arcseconds.

    Returns
    -------
    None
        Updates config.state['tracecoeffs'], extract.state['trace_done'].

    """

    #
    # Check if we can proceed
    #

    if extract.state['select_done'] is False:

        message = "Previous steps not complete.  Please run locate_apertures.py"
        raise pySpextoolError(message)

    #
    # Check parameters and QA parameters
    #

    check_parameter('trace_apertures', 'fit_degree', fit_degree, 'int')

    check_parameter('trace_apertures', 'step_size', step_size, 'int')

    check_parameter('trace_apertures', 'summation_width', summation_width,
                    'int')

    check_parameter('trace_apertures', 'centroid_threshold',
                    centroid_threshold, 'int')

    check_parameter('trace_apertures', 'seeing_fwhm', seeing_fwhm,
                    ['int','float'])

    check_parameter('trace_apertures', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('trace_apertures', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('trace_apertures', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)    
    #
    # Run the trace
    #

    if extract.state['aperture_type'] == 'fixed':

        # Build the trace coefficients from the aperture positions

        norders, naps = np.shape(extract.state['aperture_positions'])

        doorders = extract.state['doorders']
        z = doorders == 1
        ndoorders = np.sum(doorders)

        tracecoeffs = np.full((ndoorders * naps, 2), 0)
        tracecoeffs[:, 0] = np.ravel(extract.state['aperture_positions'][z, :])

    else:
        
        # Must actually trace the apertures

        z = extract.state['doorders'] == 1

        logging.info(' Tracing apertures.')
        trace = trace_spectrum_1dxd(extract.state['workimage'],
                                    extract.state['ordermask'],
                                    extract.state['orders'][z],
                                    extract.state['wavecal'],
                                    extract.state['spatcal'],
                                    extract.state['xranges'][z, :],
                                    extract.state['aperture_positions'][z, :],
                                    fit_degree=fit_degree,
                                    step_size=step_size,
                                    centroid_threshold=centroid_threshold,
                                    fwhm=seeing_fwhm,
                                    verbose=qa['verbose'])
        
        tracecoeffs = trace['coeffs']

    #
    # Store the results
    #

    extract.state['tracecoeffs'] = tracecoeffs

    #
    # Do the QA plots
    #

    if qa['show'] is True or qa['write'] is True:

            norders, naps = np.shape(extract.state['aperture_positions'])

            info = trace_to_xy(extract.state['ordermask'],
                               extract.state['wavecal'],
                               extract.state['spatcal'],
                               extract.state['xranges'],
                               extract.state['doorders'],
                               naps,
                               tracecoeffs,
                               verbose=verbose)

            if extract.state['aperture_type'] != 'fixed':

                trace_info = {'x': trace['x'], 'y': trace['y'],
                              'goodbad': trace['goodbad'], 'fits': info}
                
            else:

                trace_info = {'fits': info}                
                  
    if qa['show'] is True:
        
        plot_image(extract.state['workimage'],
                   mask=extract.state['maskimage'],
                   figure_size=(setup.plots['square_size'][0]*qa['showscale'],
                                setup.plots['square_size'][1]*qa['showscale']),
                   font_size=setup.plots['font_size']*qa['showscale'],
                   showblock=qa['showblock'],
                   plot_number=setup.plots['combine_image'],
                   trace_plotinfo=trace_info)
            
    if qa['write'] is True:

        filename = extract.state['qafilename'] + '_trace' + \
            setup.state['qa_extension']
        fullpath = join(setup.state['qa_path'],filename)

        plot_image(extract.state['workimage'],
                   mask=extract.state['maskimage'],
                   output_fullpath=fullpath,
                   figure_size=setup.plots['square_size'],
                   font_size=setup.plots['font_size'],
                   trace_plotinfo=trace_info)

    #
    # Set the done variable
    #

    extract.state['trace_done'] = True

