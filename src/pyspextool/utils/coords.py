import numpy as np
import numpy.typing as npt

from pyspextool.io.check import check_parameter


def ten(val:str | list | npt.ArrayLike,
        toradians=False):

    """
    Converts a sexigesimal number to a decimal

    Parameters
    ----------
    val : str, list or ndarray
          A sexigesimal number that is a colon-delimited string or 
          3-element list or ndarray

    toradians : {False, True}
          Set to True to multiply by np.pi/180.
          Set to False to not multiply by np.pi/180.    
    
    Returns
    -------
    float
        The decimal number

    Notes
    -----
    Based on the IDL Astronomy User's Library sixty program.  
    Basic manipulation and formatting
    
    Examples
    --------
    > x = '-00:00:40.04424'
    > ten(x)

    -0.0111234

    > x = [-0.0,0.0,40.04424]
    > ten(x)

    -0.0111234

    > x = [0.0,0.0,-40.04424]
    > ten(x)

    -0.0111234

    > import numpy as np
    > x = np.array([-0.0,0.0,40.04424])
    > ten(x)

    -0.0111234

    > import numpy as np
    > x = np.array([0.0,0.0,-40.04424])
    > ten(x)

    -0.0111234

    """

    #
    # Check parameter
    #

    check_parameter('ten', 'val', val, ['str', 'list', 'ndarray'])

    check_parameter('ten', 'toradians', toradians, 'bool')    

    #
    # Figure out what the name is 
    #
    
    typ = type(val).__name__

    # String input

    if typ == 'str':

        val0 = val.replace('++','+')
        hms = (val0.split(':'))

        decimal = abs(float(hms[0])) + float(hms[1]) / 60. + \
            float(hms[2]) / 3600.

        # Grab the first element of the string to test for positivity 

        posneg = hms[0][0]

        if posneg == '+':

            pass 
            
        elif posneg == '-':

            decimal *= -1
                
        # A list of numpy array

    elif typ == 'list' or typ == 'ndarray':

        # Convert to positive (to deal with things like [-0.0,0.0.40]

        decimal = abs(float(val[0]) + float(val[1]) / 60. + \
                      float(val[2]) / 3600.)

        # Check for negative

        prod = val[0] * val[1] * val[2]
        if prod == -0.0:
            
            decimal *= -1

    if toradians is True:

        decimal *= np.pi/180.

            
    return decimal


def sixty(val:int | float,
          colons:dict=None,
          trailsign:bool=False):

    """
    Converts a decimal number to sexigesmal

    Parameters
    ----------
    val : int or float 
        A decimal value

    colons: dict of {'dec':int, 'plus':int}, optional
        If given, the output is a colon-separated string.  'dec' gives 
        the number of decimal places on the third value set 'plus' 
        for plus sign if the first value is positive.  

    trailsign: {True, False}
        If `colons` is given, this has no affect.  By default, (False),
        the first non-zero value has the negative sign.  Setting 
        `trailsign` forces the first element to have the negative sign.  

    Returns
    -------
    list, str
        list: a 3-element list giving the sexigesimal values
        str: a string of the form hh:mm:ss 

    Note
    ----
    Based on the IDL Astronomy User's Library sixty program.  
    Basic manipulation and formating

    Examples
    --------
    > x = -0.0111234
    > sixty(x)
    
    [0.0, 0.0, -40.04424]

    > x = -0.0111234
    > sixty(x,trainsign=True)
    
    [-0.0, 0.0, 40.04424]

    > x = -0.0111234
    > sixty(x,colons={'dec':2,'plus':1},trailsign=True)

    -00:00:40.04

    > x = 0.0111234
    > sixty(x)

    [0.0, 0.0, 40.04424]

    > x = 0.0111234
    > sixty(x,trailsign=True)

    [0.0, 0.0, 40.04424]

    > x = 0.0111234
    > sixty(x,colons={'dec':2,'plus':1},trailsign=True)

    +00:00:40.04

    > x = 0.0111234
    > sixty(x,colons={'dec':3,'plus':0},trailsign=True)

    00:00:40.044

    """

    #
    # Check parameters
    #

    check_parameter('sixty', 'val', val, ['int', 'float'])

    check_parameter('sixty', 'colons', colons, ['NoneType', 'dict'])

    check_parameter('sixty', 'trailsign', trailsign, 'bool')        
    
    # First check to see if the value is negative

    neg = [0, 1][val < 0]

    # Convert the value to positive degrees, minutes, and seconds

    ss = abs(3600. * val)
    mm = abs(60. * val)
    dd = abs(val)

    # Now determine the positive sexigesimal values

    sexg = [int(dd)]
    sexg.append(int(mm - 60. * sexg[0]))
    sexg.append(ss - 3600. * sexg[0] - 60 * sexg[1])

    # Check to see whether a colon-separated string is requested.

    if colons is not None:

        # Create the positive colon-separated string.

        fmt = '{:.' + str(colons['dec']) + 'f}'

        result = str(sexg[0]).zfill(2) + ':' + \
                 str(sexg[1]).zfill(2) + ':' + \
                 fmt.format(sexg[2]).zfill(3 + colons['dec'])

        # Now deal with the positive/negative issue

        if neg == 1:

            result = '-' + result

        else:

            result = [result, '+' + result][colons['plus'] == 1]

    if colons is None:

        sexg = [float(x) for x in sexg]
        if neg == 1:

            if trailsign is True:

                sexg[0] *= -1

            else:

                if sexg[0] != 0:

                    sexg[0] *= -1

                elif sexg[1] != 0:

                    sexg[1] *= -1

                elif sexg[2] != 0:

                    sexg[2] *= -1

        result = sexg

    return result
