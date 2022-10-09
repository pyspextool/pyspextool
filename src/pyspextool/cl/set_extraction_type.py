from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter

def set_extraction_type(type,clupdate=True):

    '''
    Sets the extraction type

    Parameters
    ----------
    type : str
    The extraction type, 'point source' or 'extended source' (or 'ps','xs')

    clupdate : {True, False} optional
        Set to report the extraction type to the command line.

    Returns
    -------
    None

    Notes
    -----
    None

    Examples
    --------
    later

    '''

    #
    # Check parameter
    #

    check_parameter('set_extraction_type', 'type', type, 'str',
                    possible_values=['ps', 'xs', 'point source',
                                     'extended source'])

    #
    # Check the continue variable
    #
    
    check_continue(1)


    # Set the extraction type
    
    config.state['exttype'] = type

    # Update if requested
    
    if clupdate is True:

        label = 'point source' if type == 'ps' else 'extended source'
        print('Setting extraction type to '+label+'...')

    # Set the continue flags

    config.state['continue'] = 2 
