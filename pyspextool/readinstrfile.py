from numpy import loadtxt as nploadtxt
from numpy import array as nparray
from numpy import where as npwhere
from numpy import size as npsize
from sys import exit as exit

def readinstrfile(filename):

    """
    To read a Spextool instrument configuration file

   Input Parameters
    ----------------
    filename : str
        The name of a Spextool instrument file.

    Returns
    --------
    dict
        Later, when you finalize the pySpextool cal file

    Procedure
    ---------
    Just a bunch of searching for values.

    Examples
    --------
    > readinstrfile('uspex.dat')

    Modification History
    --------------------
    2022-06-01 - Written by M. Cushing, University of Toledo.
        Based on the Spextool IDL program mc_readinstrfile.pro.

    """

# Read the file into to string arrays
    
    labels,vals = nploadtxt(filename,comments='#',delimiter='=',\
                            unpack=True,dtype='str')
    
# Strip any edge white spaces from poor formatting of the file
    
    labels = [value.strip() for value in labels]
    labels = nparray(labels)
    
    vals = [value.strip() for value in vals]
    vals = nparray(vals)

# Start the search for each keyword.

    output = {}
    
# INSTRUMENT   
    
    keyword = 'INSTRUMENT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# NROWS
    
    keyword = 'NROWS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# NCOLS
    
    keyword = 'NCOLS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())


# STDIMAGE
    
    keyword = 'STDIMAGE'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# PLOTWINSIZE
    
    keyword = 'PLOTWINSIZE'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [float(value) for value in tmp]    
    output[keyword] = tmp
    
# NINT
    
    keyword = 'NINT'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# NSUFFIX
    
    keyword = 'NSUFFIX'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())    

# BADPIXMASK
    
    keyword = 'BADPIXMASK'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# CALMODULE
    
    keyword = 'CALMODULE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# FILEREADMODE
    
    keyword = 'FILEREADMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# OPREFIX
    
    keyword = 'OPREFIX'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()        

# SUFFIX
    
    keyword = 'SUFFIX'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# FITSREADPROGRAM
    
    keyword = 'FITSREADPROGRAM'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()            


# REDUCTIONMODE
    
    keyword = 'REDUCTIONMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# COMBMODE
    
    keyword = 'COMBMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# COMBSTAT
    
    keyword = 'COMBSTAT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# COMBTHRESH
    
    keyword = 'COMBTHRESH'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()                        

# COMBODIR
    
    keyword = 'COMBODIR'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()                        

# PSNAPS
    
    keyword = 'PSNAPS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# PSNAPS
    
    keyword = 'PSNAPS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())            
    
# OPTEXTRACT
    
    keyword = 'OPTEXTRACT'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# AVEPROF
    
    keyword = 'AVEPROF'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# PSPSFRAD
    
    keyword = 'PSPSFRAD'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

# PSBGSUB
    
    keyword = 'PSBGSUB'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# PSBGSTART
    
    keyword = 'PSBGSTART'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

# PSBGWIDTH
    
    keyword = 'PSBGWIDTH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

# PSBGDEG
    
    keyword = 'PSBGDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())            

# XSBGSUB
    
    keyword = 'XSBGSUB'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# XSBGREG
    
    keyword = 'COMBSTAT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

# XSBGDEG
    
    keyword = 'XSBGDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# TRACEDEG
    
    keyword = 'TRACEDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# TRACESTEP
    
    keyword = 'TRACESTEP'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# TRACESUMAP
    
    keyword = 'TRACESUMAP'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# TRACESIGTHRESH
    
    keyword = 'TRACESIGTHRESH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

# TRACEWINTHRESH
    
    keyword = 'TRACEWINTHRESH'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

# BADPIXELTHRESH
    
    keyword = 'BADPIXELTHRESH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

# LINCORMAX
    
    keyword = 'LINCORMAX'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())                                

# AMPCOR
    
    keyword = 'AMPCOR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# LINCOR
    
    keyword = 'LINCOR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# FLATFIELD
    
    keyword = 'FLATFIELD'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# PLOTXCORR
    
    keyword = 'PLOTXCORR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# RECTMETHOD
    
    keyword = 'RECTMETHOD'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# FIXBADPIXELS
    
    keyword = 'FIXBADPIXELS'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

# XSPEXTOOL KEYWORDS
    
    keyword = 'XSPEXTOOL_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XSPEXTOOL_KEYWORDS'] = tmp

# XCOMBSPEC KEYWORDS
    
    keyword = 'XCOMBSPEC_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XCOMBSPEC_KEYWORDS'] = tmp

# XTELLCOR KEYWORDS
    
    keyword = 'XTELLCOR_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XTELLCOR_KEYWORDS'] = tmp                                
    
    return output
    
#
#==========================================================================
#
def find(labels,keyword):

    z = npwhere(labels == keyword)
    if npsize(z):

        return z
        
    else:

        print('Cannot find keyword ',keyword,'.',sep='')
        exit(1)


    return(output)
