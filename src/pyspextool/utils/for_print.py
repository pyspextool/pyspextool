import numpy as np
from pyspextool.io.check import check_parameter


def for_print(*args):

    """
    Prints lists/arrays of data as columns to the command line.

    Input Parameters
    ----------------
    *args : lists or 1D numpy.ndarrays
        "Arrays" of values to be printed as columns together.

    Returns
    -------
    None

    Notes
    -----
    Inputs must have the same length.

    Examples
    --------
    > x = [1,2,3]
    > y = [4,5,6]
    > forprint(x,y)
      1 4
      2 5
      3 6

    """

    #
    # How many arguents?
    #
    
    nargs = len(args)

    # Convert the tuple to a list so that you can modify each argument

    args = list(args)

    #
    # Do some massaging of each argument
    #
    
    for i in range(nargs):

        # check the parameter to make sure it is a list or ndarray
        
        label = 'a['+str(i)+']'
        check_parameter('for_print', label, args[i], ['list','ndarray'])
        
        if isinstance(args[i],np.ndarray):

            # if it is an ndarray, squeeze it to ensure it has only 1 dimension
            
            args[i] = np.squeeze(args[i])            

        else:

            # otherwise, convert it to a ndarray
            
            args[i] = np.array(args[i])

    #
    # Check to make sure all inputs have the same length
    #
    
    lengths = []
    for i in range(nargs):

        lengths.append(np.size(args[i]))

    if lengths.count(lengths[0]) != nargs:

        raise Exception('Inputs are not the same length.')

    #
    # Now print things
    #
    
    args = tuple(args)

    for i in range(np.size(args[0])):

        for j in range(nargs):

            end = ' ' if j != nargs-1 else '\n'
            print(args[j][i], end=end)
