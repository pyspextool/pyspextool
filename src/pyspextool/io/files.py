import os
import sys

def check_directory(paths):

    """
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

    """

    # Make it a list just in case

    paths = [paths] if type(paths) is str else paths

    # Now loop and check
    
    for path in paths:  

        result = os.path.exists(path)

        if result is False:
            
            print('check_dir: The directory "', path,
                  '" does not exist.', sep='')
            sys.exit(1)

    # Now return the results properly
            
    if len(paths) == 1:

        return paths[0]

    else:

        return paths


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


def extract_filestring(string, method):

    """
    Extracts the indices or filenames from a comma-separated string

    Parameters
    ----------
    string : str
        a comma separated string of either file names or file index numbers.

    method : {'index','filename'}
        'index' if the values passed are index values and 'filename' if the 
        values passed our file names

    Returns
    --------
    list
        `method`='index'
         A list of integers giving the individual file numbers

        `method`='filename'
         A list of strings giving the individual file names

    Procedure
    ---------
    `method` = 'index'
    1.  separate into groups based on the comma.
    2.  loop over each group, and separate based on dash.
    3.  if no dash detected, append group (which is a number) to output list.
    3.  if dash detected, generate sequential numbers between limits and add to
        output list.
    
    'method` = 'filename'
    1.  separate into groups based on the comma.

    Examples
    --------
    > extract_filestring('1-3,5,7,10-12','index')
    [1, 2, 3, 5, 7, 10, 11, 12]

    > fsextract('spc00001.a.fits,spc00002.a.fits','filename')
    ['spc00001.a.fits', 'spc00002.a.fits']


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_fsextract.pro

    """

    #  Check whether you are in index or filename mode
    
    if method == 'index':

        #  Split on the comma
        
        groups = string.split(',')

        #  Now loop and split on the dash, full in the missing numbers, and convert
        #  to a string
        
        oarr = []
        for group in groups:
            lowupp = group.split('-')
            if len(lowupp) == 1:

                # no dash, just add to output list               

                oarr.append(int(lowupp[0]))

            else:

                # dash dectected, generate sequential numbers and add to output list
                
                arr = list(range(int(lowupp[0]), int(lowupp[1])+1))
                oarr+=arr

        return oarr

    elif method == 'filename':

        # Just split on the comma and return
        
        return string.split(',')
        
    else:

        print('method unknown.')
        return


def make_fullpath(dir, files, indexinfo: None, exist=False):
    """
    constructs fullpath strings for files

    Parameters
    ----------
    dir : str
        the directory where the files are located

    files : list 
        a list of strings that either contain the index numbers of the 
        files or the file names

    indexinfo : dict of {'nint':int,'prefix':str,'suffix':str}, optional
        a dictionary giving the information necessary to create the file 
        names from the index numbers.
        
        nint : int
            the length of the number.  zeros fill unused spaces.

        prefix : str
            the prefix of the file

        suffix : str
            the suffix of the file.

    exist : {False, True}, optional
        set to test whether the file exists

    Returns
    --------
    list
        A list of strings giving the fullpath of the files

    Procedure
    ---------
    Lots of paying attend to details and testing

    Examples
    --------

    > files = '1-5'
    > dir = '../../uSpeXdata/raw/'
    > mkfullpath(dir,files,indexinfo={'nint':5,'prefix':'spc-',
                 'suffix':'.[ab].fits'},exist=True)

    ['../../uSpeXdata/raw/spc-00001.a.fits', 
     '../../uSpeXdata/raw/spc-00002.b.fits', 
     '../../uSpeXdata/raw/spc-00003.b.fits', 
     '../../uSpeXdata/raw/spc-00004.a.fits', 
     '../../uSpeXdata/raw/spc-00005.a.fits']                 

     since the files exists locally.

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_mkfullpath.pro

    """

    #  Check whether you are in filename or index mode

    if indexinfo:

        #  Parse the index numbers

        files = extract_filestring(files, 'index')

        #  Check to see if any of the numbers are too large.

        for test in files:

            if test > 10 ** indexinfo['nint']:
                print('File numbers >=',
                      10 ** indexinfo['nint'], 'are not allowed.')
                return

                # Now create the file names

        output = [dir + indexinfo['prefix'] +
                  str(root).zfill(indexinfo['nint']) +
                  indexinfo['suffix'] for root in files]

    else:

        output = [dir + root for root in files]

        #  Now let's check to see if the file actually exists

    if exist is True:
        test = check_file(output)

    return output    
    
