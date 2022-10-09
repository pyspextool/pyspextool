from pyspextool.cl import config
from pyspextool.io.check import check_parameter


def trace_apertures(fit_degree=2, step_size=5, summation_width=5,
                    centroid_threshold=2, fwhm=0.8):

    #
    # Check parameters
    #

    check_parameter('trace_apertures', 'fit_degree', fit_degree, 'int')

    check_parameter('trace_apertures', 'set_size', set_size, 'int')

    check_parameter('trace_apertures', 'summation_width', summation_width,
                    'int')

    check_parameter('trace_apertures', 'centroid_threshold', centroid_threshold,
                    'float')            

    #
    # 
    
    
