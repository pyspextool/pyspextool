import numpy as np
import os
import glob
import logging

from pyspextool import config as setup
from pyspextool.pyspextoolerror import pySpextoolError



def check_path(path:str,
               make_absolute:bool=False):

    """
    To check whether a path exists.

    Parameters
    ----------
    path : str
        A path (can be a relative path).

    make_absolute : False, True, optional
        Set to True to replace the user path with its absolute path.

    Returns
    -------
    str
        `path`, if the path exists.  If `make_absolute` is set to True,
        `path` is converted to an absolute path before being returned.

    Examples
    --------
    later

    """

    #
    # Check parameter
    #

    check_parameter('check_path', 'path', path, 'str')

    # Expand user path just in case.

    path = os.path.expanduser(path)

    # Now check if the path exists.

    result = os.path.exists(path)

    cwd = os.path.abspath(os.getcwd())
    
    if result is False:

        message = f'The path {path} does not exist. '\
            f'The current working directory is {cwd}'
        raise pySpextoolError(message)

    else:

        # it is good, now convert to absolute if requested.

        if make_absolute is True:

            path = os.path.abspath(path)

    return path


def check_file(files:str,
               raise_error:bool=True):

    """
    To check whether a file exists, and to resolve wildcards in the name.

    Parameters
    ----------------
    files : str, list of str
        The path to a file or a list of paths to files.

    raise_error : {True, False}
        Set to True to raise a pySpextoolError if the file does not exist.
        Set to False to not raise a pySpextoolError if the file does not exist.
        
    Returns
    -------
    str, list of str
        `files`, if the file or files exist, otherwise None.

    Notes
    -----
    The program is capable of using unix wildcards, see glob.
    Returning the file names instead of True might seem odd, but
    it allows glob to find the proper file name using wildcards.

    """

    # Make it a list just in case

    files = [files] if type(files) is str else files

    # Now loop and check

    i = 0
    for file in files:  

        test = glob.glob(file)
        if not test:

            if raise_error is True:

                message = 'File '+file+' not found.'
                raise pySpextoolError(message)

            else:

                return None

            
        else:

            if len(test) > 1:

                if raise_error is True:
                
                    message = 'More than one file matches '+file+'.'
                    raise pySpextoolError(message)

                else:

                    return None
                
            else:

                files[i] = test[0]
                i += 1

    # Now return the results properly

    if len(files) == 1:

        return files[0]

    else:

        return files


def check_parameter(caller_name:str,
                    parameter_name:str,
                    parameter,
                    types,
                    *dimens,
                    possible_values=None,
                    list_types=None,
                    ndarray_type=None):

    """
    check_parameter


    Parameters
    ----------
    caller_name : str
        The name of the function calling check_parameter

    parameter_name : str
        The name of the parameter

    parameter : any
        The parameter to be checked

    types : str or list 
        The possinble types to test against the type of `parameter`.

    dimens: int or list, optional
        The possible dimensions to test against `parameter` if the type of 
        `parameter` is 'ndarray'.

    Returns
    -------
    bool

    Raises
    ------
    TBD


    Notes
    -----

    Examples
    --------
    later

    """
    
    # convert user types to a list just in case

    if isinstance(types, str):
        types = [types]

    # Get the type of the parameter

    parameter_type = type(parameter).__name__

    # Check whether `parameter_type` is `types`.  
    
    if parameter_type in types:

    # Has *dimens been passed and is the type of `parameter` ndarray?

        if len(dimens) != 0 and type(parameter).__name__ == 'ndarray':

    # Convert the dimens to a list if need be
            
            if isinstance(dimens[0], int):
                dimens = [dimens[0]]

            else: 
                dimens = dimens[0]  # because *dimens is a tuple

            # Get the dimension of the array
                
            dimen = len(parameter.shape)
            
            if (dimen not in dimens):

                dimens_str = [str(x) for x in dimens]
                
                message = 'Parameter `'+str(parameter_name)+'` of function ' + \
                caller_name+' has dimension '+str(dimen)+ \
                '.  Acceptable dimension are '+', '.join(dimens_str)+'.'

                raise ValueError(message)
                return False

    else:
        
        message = 'Parameter `'+str(parameter_name)+'` of function '+ \
          caller_name+' has type '+parameter_type+'.  Acceptable types are '+ \
          ', '.join(types)+'.'

        raise TypeError(message)
        return False
        

    if possible_values is not None:
        
        if (parameter not in possible_values):        

            values_str = ['`'+str(x)+'`' for x in possible_values]

            message = 'Parameter `'+str(parameter_name)+'` of function '+ \
              caller_name+' has a value of `'+str(parameter)+ \
              '`.  Acceptable values are, '+', '.join(values_str)+'.'
            
            raise ValueError(message)

    if list_types is not None:

        if parameter_type == 'list':

            for i in range(len(list_types)):

                if type(parameter[i]).__name__ != list_types[i]:

                    message = 'Parameter `'+str(parameter_name) + \
                        '['+str(i)+']` of function ' + \
                    caller_name+' has type ' + \
                    type(parameter[i]).__name__+ \
                    '.  Acceptable type is '+list_types[i]+'.'

                    raise ValueError(message)
                    return False

#    if ndarray_type is not None:
#        
#        if parameter_type == 'ndarray':
#
#            if parameter[i].dtype != ndarray_type:
#            
#                    message = 'Parameter `'+str(parameter_name) + \
#                        '['+str(i)+']` of function ' + \
#                    caller_name+' has dtype ' + \
#                    parameter[i].dtype+ \
#                    '.  Acceptable type is '+ndarray_type+'.'
#
#                    raise ValueError(message)
#                    return False

        
def check_range(values,
                value_range,
                test,
                variable_name=None):

    """
    To check whether a set of numbers or number is in a given range

    Parameters
    ----------
    values : int or float or array_like
        the set of numbers or number under question

    value_range : int or float or array_like
        a single value of a (2,) array_like giving the range.

    test : str {'gt', 'ge', 'lt', 'le', 'gtlt', 'gtle', 'gelt', 'gele'}
        the test requested.

        if 'gt', 'ge', 'lt', or 'le' `values` must be a single value.
        if 'gtlt', 'gtle', 'gelt', or 'gele' `values` must be a (2,) array_like

    variable_name : str, optional
        The name of the variable being tested.  Useful if called from a 
        function.  

    Returns
    -------
    None

    Examples
    --------
    later


    """

    #
    # Check parameter
    #

    check_parameter('check_value', 'values', values,
                    ['int', 'float', 'list', 'ndarray'])

    check_parameter('check_value', 'values', value_range,
                    ['int', 'float', 'list', 'ndarray'])

    check_parameter('check_value', 'values', test, 'str',
                    possible_values=['gt', 'ge', 'lt', 'le', 'gtlt', 'gtle',
                                     'gelt', 'gele'])

    check_parameter('check_value', 'variable_name', variable_name,
                    ['str', 'NoneType'])    

    #
    # Get basic values
    #

    nrange = np.size(value_range)
    values = np.array(values)
    ndat = np.size(values)
    
    #
    # Now make sure values has the right form given the requested test
    #

    name = '`values`' if variable_name is None else '`'+variable_name+'`'
    
    if test in ['gt', 'ge', 'lt', 'le', 'odd', 'even'] and nrange != 1:

        message = name+' should be a single value.'
        raise ValueError(message)

    if test in ['gtlt', 'gtle', 'gelt', 'gele'] and nrange != 2:
        
        message = name+' should be two elements.'
        raise ValueError(message)    

    #
    # Create error messages
    #

    message_1 = ''
    message_2 = ''    

    if test in ['gt', 'gtlt', 'gtle']:

        value = value_range if nrange == 1 else value_range[0]
        message_1 = str(value)+' <'

    if test in ['ge', 'gelt', 'gele']:

        value = value_range if nrange == 1 else value_range[0]
        message_1 = str(value)+' <='

    if test in ['lt', 'gtlt', 'gelt']:

        value = value_range if nrange == 1 else value_range[1]
        message_2 = '< '+str(value)

    if test in ['le', 'gele', 'gtle']:

        value = value_range if nrange == 1 else value_range[1]
        message_2 = '<= '+str(value)

    message = name+' is out of range.  '+message_1+' '+name+' '+message_2+'.'

    #
    # Start the tests
    #

    if test == 'gt':

        test_value = value_range if nrange == 1 else value_range[0]
        z = values > float(test_value)
        if np.sum(z) != ndat:

            raise ValueError(message)

    if test == 'ge':

        test_value = value_range if nrange == 1 else value_range[0]
        z = values >= float(test_value)
        if np.sum(z) != ndat:

            raise ValueError(message)

    if test == 'lt':

        test_value = value_range if nrange == 1 else value_range[0]
        z = values < float(test_value)
        if np.sum(z) != ndat:

            raise ValueError(message)

    if test == 'le':

        test_value = value_range if nrange == 1 else value_range[0]
        z = values <= float(test_value)
        if np.sum(z) != ndat:

            raise ValueError(message)                        

    if test == 'gtlt':

        z1 = values > float(value_range[0])
        z2 = values < float(value_range[1])        
        if np.sum(z1*z2) != ndat:

            raise ValueError(message)

    if test == 'gtle':

        z1 = values > float(value_range[0])
        z2 = values <= float(value_range[1])        
        if np.sum(z1*z2) != ndat:

            raise ValueError(message)

    if test == 'gelt':

        z1 = values >= float(value_range[0])
        z2 = values < float(value_range[1])        
        if np.sum(z1*z2) != ndat:

            raise ValueError(message)

    if test == 'gele':

        z1 = values >= float(value_range[0])
        z2 = values <= float(value_range[1])        
        if np.sum(z1*z2) != ndat:

            raise ValueError(message)                


def check_qakeywords(**kwargs):

    """
    Checks user input pySpextool QA keywords against those set in the setup.

    

    """

    keywords = list(kwargs.keys())

    output = {}
    
    for keyword in keywords:
        
        value = kwargs.get(keyword)

        if keyword == 'verbose':
            
            if value is None:

                if setup.state["verbose"] is True:
                    logging.getLogger().setLevel(logging.INFO)
                    output['verbose'] = True
                    
                if setup.state["verbose"] is False:
                    logging.getLogger().setLevel(logging.ERROR)
                    output['verbose'] = False
                    
            if value is True:
                logging.getLogger().setLevel(logging.INFO)
                output['verbose'] = True
                
            if value is False:
                logging.getLogger().setLevel(logging.ERROR)
                output['verbose'] = False
                
        elif keyword == 'show':

            if value is None:
        
                qa_show = setup.state['qa_show']

            else:

                qa_show = value

            output['show'] = qa_show

            if qa_show is True:

                if setup.state['qa_path'] is None:

                    message = 'The `qa_path` cannot be None if `qa_show` '+\
                        'is True.'
                    raise pySpextoolError(message)
                                       
        elif keyword == 'showscale':

            if value is None:
        
                qa_showscale = setup.state['qa_showscale']

            else:

                qa_showscale = value

            output['showscale'] = qa_showscale

        elif keyword == 'showblock':

            if value is None:
        
                qa_showblock = setup.state['qa_showblock']

            else:

                qa_showblock = value

            output['showblock'] = qa_showblock                

        elif keyword == 'write':

            if value is None:
        
                qa_write = setup.state['qa_write']

            else:

                qa_write = value

            output['write'] = qa_write

            if qa_write is True:

                if setup.state['qa_path'] is None:

                    message = 'The `qa_path` cannot be None if `qa_write` '+\
                        'is True.'
                    raise pySpextoolError(message)
            
            
            
        else:

            message = "Keyword '"+keyword+"'"+" is unknown."
            raise pySpextoolError(message)
            
    return output



def check_sansfits(file:str,
                   variable_name:str):

    """
    Determines whether the string ends in '.fits'

    Parameters
    ----------
    file : str
        A file with a potential suffix of '.fits'.

    Returns
    -------
    None
    
    """

    #
    # Check parameters
    #

    check_parameter('check_sansfits', 'file', file, 'str')

    check_parameter('check_sansfits', 'variable_name', variable_name, 'str')

    if file[-5:] == '.fits':

        message = "The variable `"+variable_name+"` cannot have a '.fits' "+\
            "suffix."
        raise pySpextoolError(message)
        
