import os.path
import glob
import sys

def check_file(files):

    '''
    To check whether a file exists, and to resolve wildcards in the name.


    Input Parameters
    ----------------
    files : str, list of str
        The path to a file or a list of paths to files.


    Returns
    -------
    str, list of str
        `files`, if the file or files exist, otherwise None.


    Notes
    -----
    The program is capable of using unix wildcards, see glob.


    Examples
    --------
    later when the package stuff is worked out

    '''

# Make it a list just in case

    files = [files] if type(files) is str else files

# Now loop and check

    i = 0
    for file in files:  

        test = glob.glob(file)
        if not test:

            print('check_file: The file "',file,'" does not exist.',sep='')
            sys.exit(1)
            
        else:

            if len(test) > 1:

                print('check_file: More than one files matches "',\
                          file,'"',sep='')
                sys.exit(1)                
                
            else:
                
                files[i] = test[0]
                i +=1
                
# Now return the results properly
            
    if len(files) == 1:

        return(files[0])

    else:

        return(files)        
