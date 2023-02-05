import sys

from pyspextool.extract import config
from pyspextool.io.check import check_parameter

def check_continue(number):

    #
    # Check parameter
    #
    
    check_parameter('check_continue', 'number', number, 'int')

    if config.state['continue'] < number:

        message = 'Previous steps not completed.'
        print()
        print(message)
        print()
        sys.exit(1)

        
