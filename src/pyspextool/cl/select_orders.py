import numpy as np

from pyspextool.cl import config
from pyspextool.cl.trace_apertures import trace_apertures
from pyspextool.io.check import check_parameter
from pyspextool.io.files import extract_filestring
from pyspextool.plot.plot_profiles import plot_profiles


def select_orders(include=None, exclude=None, include_all=False, clupdate=True,
                  iplot=False, qafile=True):

    """
    To set which orders are to be traced and extracted

    Parameters
    ----------
    include : int, list, str, optional
        If the type is int, the single order to include.
        If the type is list, a list of integer orders to include.
        If the type is str, a str giving the orders, e.g. '1-3,4,5'.

    exclude : int, list, str, optional
        If the type is int, the single order to include.
        If the type is list, a list of integer orders to include.
        If the type is str, a str giving the orders, e.g. '1-3,4,5'.

    include_all : {False, True}, optional
        Set to include all orders.

    iplot : {False, True}, optional
        Set to True to plot the profiles interactively.

    Examples
    --------
    later

    """

    #
    # Update command line if requested.
    #
    if clupdate is True:
        print('Updating order selection...')
    
    #
    # Ensure only one optional argument is passed
    #

    if include is not None and exclude is not None:
        message = 'Cannot use both parameters `include` and `remove`.'
        raise ValueError(message)

    #
    # Check parameters
    #

    check_parameter('select_orders', 'include', include,
                    ['NoneType', 'int', 'list', 'str'])

    check_parameter('select_orders', 'exclude', exclude,
                    ['NoneType', 'int', 'list', 'str'])

    check_parameter('select_orders', 'include_all', include_all, 'bool')

    check_parameter('select_orders', 'iplot', iplot, 'bool')

    #
    # Do the checks
    #

    if include is not None:

        if include.__class__.__name__ == 'int':
            include = np.array(include, dtype=int)

        if include.__class__.__name__ == 'list':
            include = np.array(include, dtype=int)

        if include.__class__.__name__ == 'str':
            include = np.array(extract_filestring(include, 'index'), dtype=int)

        # Find the overlap 

        test = np.isin(config.state['orders'], include)

        # Test to make sure they are orders you are allowed work with

        if np.sum(test) != np.size(include):
            message = 'A requested order does not exist.'
            raise ValueError(message)

    if exclude is not None:

        if exclude.__class__.__name__ == 'int':
            exclude = np.array(exclude, dtype=int)

        if exclude.__class__.__name__ == 'list':
            exclude = np.array(exclude, dtype=int)

        if exclude.__class__.__name__ == 'str':
            exclude = np.array(extract_filestring(exclude, 'index'), dtype=int)

        # Find the overlap 

        test = ~np.isin(config.state['orders'], exclude)

        # Test to make sure they are orders you are allowed work with

        if np.sum(~test) != np.size(exclude):
            message = 'A requested order does not exist.'
            raise ValueError(message)

    if include_all is True:
        test = np.full(config.state['norders'], 1)

    #
    # Set the correct doorders variable
    #

    if config.state['exttype'] == 'xs':

        config.state['xsdoorders'] = test
        doorders = test

    else:

        config.state['psdoorders'] = test
        doorders = test

    if iplot is True:
        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'])

    if qafile is True:
        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': config.state['qapath'],
                      'filename': config.state['qafilename'] + \
                                  '_aperturepositions', 'extension': '.pdf'}

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'],
                      qafileinfo=qafileinfo)

    #
    # Do the trace if the extraction is extended source
    #
        
    if config.state['exttype'] == 'xs':

        trace_apertures(clupdate=clupdate, iplot=iplot, qafile=qafile)        

    
        
