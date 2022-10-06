from pyspextool.cl import config
from pyspextool.io.check import check_parameter

def check_continue(number):

    #
    # Check parameter
    #
    
    check_parameter('check_continue', 'number', number, 'int')

    # Do the check

    if config.state['exttype'] == 'ps':

        if number < config.state['pscontinue']:

            message = 'Previous steps not completed.'
            raise ValueError(message)
        

    if config.state['exttype'] == 'xs':
    
        if number < config.state['xscontinue']:

            message = 'Previous steps not completed.'
            raise ValueError(message)
