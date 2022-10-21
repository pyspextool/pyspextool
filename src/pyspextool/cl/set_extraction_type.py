from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter


def set_extraction_type(extraction_type, clupdate=True):

    """
    Sets the extraction type

    Parameters
    ----------
    extraction_type : str
    The extraction type, 'point source' or 'extended source' (or 'ps', 'xs')

    clupdate : {True, False} optional
        Set to report the extraction type to the command line.

    Returns
    -------
    None
    Sets the config.state['exttype'] variable.

    Notes
    -----
    None

    Examples
    --------
    ---To be updated---
    set_extraction_type('ps')
    set_extraction_type('xs')

    """

    #
    # Check the continue variable
    #
    
    check_continue(1)
    
    #
    # Check parameter
    #

    check_parameter('set_extraction_type', 'extraction_type',
                    extraction_type, 'str',
                    possible_values=['ps', 'xs', 'point source',
                                     'extended source'])

    # Set the extraction type

    if extraction_type == 'point source' or extraction_type == 'ps':

        extraction_type = 'ps'

    if extraction_type == 'extended source' or extraction_type == 'xs':

        extraction_type = 'xs'        
    
    config.state['exttype'] = extraction_type

    # Update if requested
    
    if clupdate is True:

        label = 'point source' if extraction_type == 'ps' else 'extended source'
        print('Setting extraction type to '+label+'...')

    # Set the continue flags

    config.state['continue'] = 2
