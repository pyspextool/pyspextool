from pyspextool.io.check_parameter import check_parameter

def reorder_irtf_files(fullpaths):

    """
    reorder_irtf_files


    Parameters
    ----------
    fullpaths : str or list
        The fullpaths of IRTF FITS files


    Returns
    -------
    list
         A list fullpaths names ordered such that the beams are ababab...

    Notes
    -----

    Examples
    --------
    later

    """

    check_parameter('abba_to_abab', 'fullpaths', fullpaths, ['str','list'])

    if isinstance(fullpaths, str):
        fullpaths = [fullpaths]
    
    nfiles = len(fullpaths)

    if nfiles % 2 !=0:

        message = '`fullpaths` must contain an even number of images.'
        raise ValueError(message)

    for i in range((nfiles-1)//2):

        # Pull the beam label from the file name
        
        position = fullpaths[i*2].rfind('.fits')
        beam = fullpaths[i*2][position-1]

        # If it is 'b', switch the filesnames
        
        if beam == 'b':

            tmp = fullpaths[i*2] 
            fullpaths[i*2] = fullpaths[i*2+1]
            fullpaths[i*2+1] = tmp


    return fullpaths
