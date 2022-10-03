import numpy as np

from pyspextool.cl import config
from pyspextool.io.check_parameter import check_parameter
from pyspextool.io.files import extract_filestring

def select_orders(include=None, exclude=None, plot_profiles=True):

    # Ensure only one optional argument is passed
    
    if include is not None and exclude is not None:

        message = 'Cannot use both include and remove'
        raise ValueError(message)
    
    #
    # Check parameters
    #
    check_parameter('select_orders','include',include,
                    ['NoneType', 'int', 'list', 'str'])

    check_parameter('select_orders','exclude',exclude,
                    ['NoneType', 'int', 'list', 'str'])

    check_parameter('select_orders','plot_profiles',plot_profiles,
                    'bool')        
    
    
    if include is not None:

        if include.__class__.__name__ == 'int':

            include = np.array(include,dtype=int)

        if include.__class__.__name__ == 'list':

            include = np.array(include,dtype=int)            

        if include.__class__.__name__ == 'str':

            include = np.array(extract_filestring(include,'index'),dtype=int)

        # Test to make sure they are orders you are allowed work with

        test = np.isin(include, config.state['orders'])

        if np.sum(test) != np.size(include):

            message = 'At least one order requested is not available.'
            raise ValueError(message)

        else:
            test = ~np.isin(config.state['orders'], include)

            
    if exclude is not None:

        if exclude.__class__.__name__ == 'int':

            exclude = np.array(exclude,dtype=int)

        if exclude.__class__.__name__ == 'list':

            exclude = np.array(exclude,dtype=int)            

        if exclude.__class__.__name__ == 'str':

            exclude = np.array(extract_filestring(exclude, 'index'), dtype=int)

        # Test to make sure they are orders you are allowed work with

        test = np.isin(exclude, config.state['orders'])

        if np.sum(test) != len(orders):

            message = 'At least one order requested is not available.'
            raise ValueError(message)

        else:

            test = np.isin(config.state['orders'], np.array(exclude, dtype=int))
            
    doorders = np.full(config.state['norders'],1,dtype=int)
    doorders[test] = 0
    config.state['doorders'] = test
    print(doorders)
    

    
    
    
    
