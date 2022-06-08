from astropy.visualization import (PercentileInterval,ZScaleInterval,MinMaxInterval)

def img_range(arr,info):

    """
    a wrapper for astropy to get image ranges

    Input Parameters
    ----------------
    arr : numpy.npdarray
        an array (typically a 2D image)

    info : float or str
        float - the fraction of pixels to keep in the range (0,1)
        str - 'zscale' or 'minmax'

    Returns
    --------
    tuple
         interval based on user request

    Procedure
    ---------
    calls astropy.visualization routines and uses all default settings


    Examples
    --------
    later
    

    Modification History
    --------------------
    2022-06-07 - Written by M. Cushing, University of Toledo.

    """
    
    if isinstance(info,(float)) is True:

# This is if the users requests a fraction

        interval = PercentileInterval(info)
        range = interval.get_limits(arr)
        return range

    elif info == 'zscale':

## This is if the users requests the zscale
        
        interval = ZScaleInterval()
        range = interval.get_limits(arr)
        return range

    elif info == 'minmax':

## This is if the users requests the zscale

        interval = MinMaxInterval()
        range = interval.get_limits(arr)
        return range

    else:

        return None

