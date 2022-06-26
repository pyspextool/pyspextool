import os.path
import sys

def check_dir(paths):

    '''
    To check whether a directory exists.


    Input Parameters
    ----------------
    paths : str, list of str
        The path to a directory or a list of paths to directories.


    Returns
    -------
    str, list of str
        `paths`, if the directory or directories exist, otherwise None.


    Notes
    -----
    Just checks to see if the directory exists.


    Examples
    --------
    later when the package stuff is worked out.

    '''

# Make it a list just in case

    paths = [paths] if type(paths) is str else paths

# Now loop and check
    
    for path in paths:  

        result = os.path.exists(path)

        if result == False:
            
            print('check_dir: The directory "',path,\
                '" does not exist.',sep='')
            sys.exit(1)

# Now return the results properly
            
    if len(paths) == 1:

        return(paths[0])

    else:

        return(paths)        
