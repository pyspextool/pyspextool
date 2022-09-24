def check_parameter(caller_name, parameter_name, parameter, types, *dimens,
                    possible_values=None):

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
    
    if (parameter_type in types):

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
                
                message = 'Parameter `'+str(parameter_name)+'` of ' + \
                caller_name+' has dimension '+str(dimen)+ \
                '.  Acceptable dimension are '+', '.join(dimens_str)+'.'

                raise ValueError(message)
                return False

    else:
        
        message = 'Parameter `'+str(parameter_name)+'` of '+ \
          caller_name+' has type '+parameter_type+'.  Acceptable types are '+ \
          ', '.join(types)+'.'

        raise TypeError(message)
        return False
        

    if possible_values is not None:

        if (parameter not in possible_values):        

            values_str = ['`'+str(x)+'`' for x in possible_values]

            message = 'Parameter `'+str(parameter_name)+'` of '+ \
              caller_name+' has a value of `'+str(parameter)+ \
              '`.  Acceptable values are, '+', '.join(values_str)+'.'            
            raise ValueError(message)

    return True
