import os.path
import sys

def check_exist(paths):

    '''

    Input Parameters
    ----------------
    paths : str, list of str
        The (full) path(s) or a directory or a file.


    Returns
    --------
    str, list of str
        `paths`, if the paths or files exists, otherwise None.


    Notes
    -----
    Just checks to see if the file exists.


    Examples
    --------
    later

    '''

# Make it a list just in case

    paths = [paths] if type(paths) is str else paths

# Now loop and check
    
    for path in paths:  

        result = os.path.exists(path)

        if result == False:
            
            print('check_exist: The file or path "',path,\
                '" does not exist.',sep='')
            sys.exit(1)

# Now return the results properly
            
    if len(paths) == 1:

        return(paths[0])

    else:

        return(paths)        
