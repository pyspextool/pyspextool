from fsextract import fsextract
import glob

def mkfullpath(dir,files,indexinfo:None,exist=False):

    """
    constructs fullpath strings for files

    Input Parameters
    ----------------
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

        files = fsextract(files,'index')

#  Check to see if any of the numbers are too large.

        for test in files:

            if test > 10**indexinfo['nint']:

                print('File numbers >=',\
                      10**indexinfo['nint'] ,'are not allowed.')
                return 
        
# Now create the file names

        output = [dir+indexinfo['prefix']+\
                  str(root).zfill(indexinfo['nint'])+\
                  indexinfo['suffix'] for root in files]

    else:

        output = [dir+root for root in files] 


#  Now let's check to see if the file actually exists
        
    if exist == True:

        i = 0        
        for name in output:

            test = glob.glob(name)
            if not test:

                print('File ',name, 'does not exist.')
                return

            else:

                if len(test) > 1:

                    print('More than one file matches ',name)
                    return
                    
                else:
                
                    output[i] = test[0]

            i+=1

    return(output)
