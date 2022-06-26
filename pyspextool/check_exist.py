import os.path
import sys

def check_exist(path):

    '''

    Input Parameters
    ----------------
    path : str
        The (full) path of a directory or a file.


    Returns
    --------
    str
        `path`, if the path or file exists, otherwise None.


    Notes
    -----
    Just checks to see if the file exists.


    Examples
    --------
    later

    '''
    
    result = os.path.exists(path)

    if result == True:

        return(path)

    else:

        print('check_exist: The file or path "',path,\
              '" does not exist.',sep='')
        sys.exit(1)
