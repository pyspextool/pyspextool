import sys

from pyspextool.cl import config
from pyspextool.io.check import check_parameter

def check_continue(number):

    #
    # Check parameter
    #
    
    check_parameter('check_continue', 'number', number, 'int')

    # Do the check

#    if config.state['exttype'] == 'ps':
        
#        scontinue = config.state['pscontinue']

#    else:

#        scontinue = config.state['xscontinue']

    if config.state['continue'] < number:

        message = 'Previous steps not completed.  Please run '+\
                  config.state['steps'][number-1]+'.'
        print(message)
        sys.exit(1)

        
