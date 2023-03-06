import os

from pyspextool.io.check import check_file
from pyspextool.io.check import check_parameter

def extract_filestring(string, method):

    """
    Extracts the indices or filenames from a comma-separated string.

    Parameters
    ----------
    string : str
        a comma separated string of either file names or file index numbers.

    method : {'index', 'filename'}
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

    > extract_filestring('spc00001.a.fits,spc00002.a.fits','filename')
    ['spc00001.a.fits', 'spc00002.a.fits']

    """

    #
    # Check parameters
    #
    check_parameter('extract_filestring', 'string',  string, 'str')

    check_parameter('extract_filestring', 'method',  method, 'str')    

    #
    #  Check whether you are in index or filename mode
    #
    
    if method == 'index':

        #  Split on the comma
        
        groups = string.split(',')

        #  Now loop and split on the dash, full in the missing numbers,
        #  and convert to a string
        
        oarr = []
        for group in groups:
            lowupp = group.split('-')
            if len(lowupp) == 1:

                # no dash, just add to output list               

                oarr.append(int(lowupp[0]))

            else:

                # dash detected, generate sequential numbers and add to
                # output list
                
                arr = list(range(int(lowupp[0]), int(lowupp[1])+1))
                oarr+=arr

        return oarr

    elif method == 'filename':

        # Just split on the comma and return
        
        return string.split(',')
        
    else:

        message = 'Unknown method: `index` or `filename`.'
        raise ValueError(message)



def make_full_path(dir, files, indexinfo=None, exist=False):
    
    """
    Constructs fullpath strings for files.

    Parameters
    ----------
    dir : str
        the directory where the files are located

    files : list, str, int 
        a list of strings that either contain the index numbers of the 
        files or the file names

    indexinfo : dict, optional
        A dictionary giving the information necessary to create the file 
        names from the index numbers.

        `'nint'` : int
            the length of the number.  zeros fill unused spaces.

        `'prefix'` : str
            the prefix of the file

        `'suffix'` : str
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
                 'suffix':'.[ab].fits'})

    ['../../uSpeXdata/raw/spc-00001.a.fits', 
     '../../uSpeXdata/raw/spc-00002.b.fits', 
     '../../uSpeXdata/raw/spc-00003.b.fits', 
     '../../uSpeXdata/raw/spc-00004.a.fits', 
     '../../uSpeXdata/raw/spc-00005.a.fits']                 

    """

    #
    # Check parameters
    #

    check_parameter('make_full_path', 'dir', dir, 'str')

    check_parameter('make_full_path', 'files', files, ['list', 'str', 'int'])

    check_parameter('make_full_path', 'indexinfo', indexinfo, 'dict')

    check_parameter('make_full_path', 'exist', exist, 'bool')            

    #
    #  Check whether you are in filename or index mode
    #
    
    if indexinfo:

        #  Parse the index numbers if required

        if isinstance(files, str):

            files = extract_filestring(files, 'index')

        # Take integers and make them lists.

        if isinstance(files, int):

            files = [files]
            
        #  Check to see if any of the numbers are too large.

        for test in files:

            if test > 10 ** indexinfo['nint']:

                message = 'File numbers >='+str(10**indexinfo['nint'])+\
                          'are not allowed.'
                raise ValueError(message)

                # Now create the file names

        output = [os.path.join(dir,indexinfo['prefix'] +
                  str(root).zfill(indexinfo['nint']) +
                  indexinfo['suffix']+indexinfo['extension'])
                  for root in files]

    else:

        output = [os.path.join(dir,root) for root in files]

        #  Now let's check to see if the file actually exists

    if exist is True:
        test = check_file(output)

    return output    
    
