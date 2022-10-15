import numpy as np

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_image import plot_image
from pyspextool.spectroscopy.trace_spectrum_1dxd import trace_spectrum_1dxd
from pyspextool.spectroscopy.trace_to_xy import trace_to_xy
from pyspextool.utils.for_print import for_print


def trace_apertures(fit_degree=2, step_size=5, summation_width=5,
                    centroid_threshold=2, fwhm=0.8, clupdate=True,
                    iplot=True, qafile=False):

    

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
                                config.state['xranges'][z,:],
                                config.state['apertures'][z,:],
                                fit_degree=fit_degree, step_size=step_size,
                                centroid_threshold=centroid_threshold,
                                fwhm=fwhm, clupdate=clupdate)

        if iplot is True or qafile is True:
    
            plotinfo = {'x':trace['x'], 'y':trace['y'],
                        'goodbad':trace['goodbad']}
                
    else:

        # Must be an extended source extraction
        
        doorders = config.state['xsdoorders']

        
    norders, naps = np.shape(config.state['apertures'])
    trace_to_xy(config.state['ordermask'], config.state['wavecal'],
                config.state['spatcal'], config.state['xranges'],
                config.state['orders'], doorders, naps, trace['coeffs'])
                

    if iplot is True:
   
        plot_image(config.state['workimage'], trace_plotinfo=plotinfo)

    if qafile is True:

        qafileinfo = {'figsize': (7,7), 'filepath':config.state['qapath'],
                      'filename':config.state['qafilename']+'_trace',
                      'extension':'.pdf'}

        plot_image(config.state['workimage'], trace_plotinfo=plotinfo,
                   qafileinfo=qafileinfo)            
    
