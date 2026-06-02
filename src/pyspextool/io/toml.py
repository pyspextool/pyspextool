from typing import Any

from pyspextool.io.check import check_parameter
from pyspextool.pyspextoolerror import pySpextoolError

def get_tomlparameter(
    parameter_name:list | str,
    default_toml:dict,
    user_toml:dict | None,
    user_value):

    """
    To return the value from the user, user toml file, or default toml file.

    Parameters
    ----------
    parameter_name : list or str
       
    default_toml : dict
    
    user_toml : dict

    user_value : 
        A value passed by the user for this parameter.

    """
    
    #
    # Check the parameters
    # 

    check_parameter('get_tomlparameter', 'parameter_name',
                    parameter_name, ['list', 'str'])

    check_parameter('get_tomlparameter', 'default_toml',
                    default_toml, 'dict')

    check_parameter('get_tomlparameter', 'user_toml',
                    user_toml, ['dict','NoneType'])

    #
    # First check to see if the parameter exists in the default toml
    #

    default_value = default_toml

    try:
        for key in parameter_name:
            default_value = default_value[key]

    except KeyError:

        message = 'The parameter '+".".join(parameter_name)+' does not '+\
            'exist in '+default_toml['toml_fullpath']+'.'
        raise pySpextoolError(message)
        
    #
    # Now start the process 
    #

    if user_value is not None:

        # Return the user value

        return user_value

    if user_toml is not None:

        # Check to see if it is exists

        value = user_toml

        try:
            for key in parameter_name:
                value = value[key]

            # Return the user toml value

            return value

        except KeyError:

            # Return the default toml value

            return default_value

        
