from pyspextool.cl import config
from pyspextool.io.check_parameter import check_parameter

def set_extraction_type(type,clupdate=True):

    '''
    Sets the extraction type

    Parameters
    ----------
    type : str
        The extraction type, `ps` or `xs`.

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
                    possible_values=['ps','xs'])

    config.state['exttype'] = type

    if clupdate is True:

        label = 'point source' if type == 'ps' else 'extended source'
        
        print()
        print('Extraction type set to '+label+'.')
