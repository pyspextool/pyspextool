from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter


def set_extraction_type(extraction_type, verbose=None):

    """
    Sets the extraction type

    Parameters
    ----------
    extraction_type : str
    The extraction type, 'point source' or 'extended source' (or 'ps', 'xs')

    verbose : {None, True, False} optional
        Set True/False to override setup.state['verbose'].

    Returns
    -------
    None
    Sets the extract.type['type'] variable.

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
    # Check the load_done variable
    #
    
    if extract.state['load_done'] is False:

        message = "extract.state['load_done']=False.  "+\
          "Previous steps not complete."        
        raise ValueError(message)
    
    #
    # Check parameter
    #

    check_parameter('set_extraction_type', 'extraction_type',
                    extraction_type, 'str',
                    possible_values=['ps', 'xs', 'point source',
                                     'extended source'])

    check_parameter('set_extraction_type', 'verbose', verbose,
                    ['NoneType', 'bool'])    

    # Set the extraction type

    if extraction_type == 'point source' or extraction_type == 'ps':

        extraction_type = 'ps'

    if extraction_type == 'extended source' or extraction_type == 'xs':

        extraction_type = 'xs'        
    
    extract.state['type'] = extraction_type
    extract.type['type'] = extraction_type    

    # Check verbose
    
    if verbose is None:

        verbose = setup.state['verbose']

    extract.type['verbose'] = verbose
    
    # Update if requested
    
    if verbose is True:

        label = 'point source' if extraction_type == 'ps' else 'extended source'
        print('Setting extraction type to '+label+'...')

    # Set the done variable

    extract.state['type_done'] = True
