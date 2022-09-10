import numpy as np


def read_instrument_file(filename):
    """
    To read a Spextool instrument configuration file

    Parameters
    ----------
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
    > read_instrument_file('uspex.dat')

    Modification History
    --------------------
    2022-06-01 - Written by M. Cushing, University of Toledo.
        Based on the Spextool IDL program mc_readinstrfile.pro.

    """

    # Read the file into to string arrays

    labels, vals = np.loadtxt(filename, comments='#', delimiter='=',
                              unpack=True, dtype='str')

    # Strip any edge white spaces from poor formatting of the file

    labels = [value.strip() for value in labels]
    labels = np.array(labels)

    vals = [value.strip() for value in vals]
    vals = np.array(vals)

    # Start the search for each keyword.

    output = {}

    # INSTRUMENT   

    keyword = 'INSTRUMENT'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # NROWS

    keyword = 'NROWS'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # NCOLS

    keyword = 'NCOLS'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # STDIMAGE

    keyword = 'STDIMAGE'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # PLOTWINSIZE

    keyword = 'PLOTWINSIZE'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [float(value) for value in tmp]
    output[keyword] = tmp

    # NINT

    keyword = 'NINT'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # NSUFFIX

    keyword = 'NSUFFIX'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # BADPIXMASK

    keyword = 'BADPIXMASK'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # CALMODULE

    keyword = 'CALMODULE'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # FILEREADMODE

    keyword = 'FILEREADMODE'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # OPREFIX

    keyword = 'OPREFIX'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # SUFFIX

    keyword = 'SUFFIX'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # EXTENSION

    keyword = 'EXTENSION'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()    

    # FITSREADPROGRAM

    keyword = 'FITSREADPROGRAM'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # REDUCTIONMODE

    keyword = 'REDUCTIONMODE'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # COMBMODE

    keyword = 'COMBMODE'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # COMBSTAT

    keyword = 'COMBSTAT'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # COMBTHRESH

    keyword = 'COMBTHRESH'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # COMBODIR

    keyword = 'COMBODIR'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # PSNAPS

    keyword = 'PSNAPS'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # PSNAPS

    keyword = 'PSNAPS'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # OPTEXTRACT

    keyword = 'OPTEXTRACT'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # AVEPROF

    keyword = 'AVEPROF'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # PSPSFRAD

    keyword = 'PSPSFRAD'
    z = find_keyword(labels, keyword)
    output[keyword] = float(vals[z].item())

    # PSBGSUB

    keyword = 'PSBGSUB'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # PSBGSTART

    keyword = 'PSBGSTART'
    z = find_keyword(labels, keyword)
    output[keyword] = float(vals[z].item())

    # PSBGWIDTH

    keyword = 'PSBGWIDTH'
    z = find_keyword(labels, keyword)
    output[keyword] = float(vals[z].item())

    # PSBGDEG

    keyword = 'PSBGDEG'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # XSBGSUB

    keyword = 'XSBGSUB'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # XSBGREG

    keyword = 'COMBSTAT'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # XSBGDEG

    keyword = 'XSBGDEG'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # TRACEDEG

    keyword = 'TRACEDEG'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # TRACESTEP

    keyword = 'TRACESTEP'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # TRACESUMAP

    keyword = 'TRACESUMAP'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # TRACESIGTHRESH

    keyword = 'TRACESIGTHRESH'
    z = find_keyword(labels, keyword)
    output[keyword] = float(vals[z].item())

    # TRACEWINTHRESH

    keyword = 'TRACEWINTHRESH'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # BADPIXELTHRESH

    keyword = 'BADPIXELTHRESH'
    z = find_keyword(labels, keyword)
    output[keyword] = float(vals[z].item())

    # LINCORMAX

    keyword = 'LINCORMAX'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # AMPCOR

    keyword = 'AMPCOR'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # LINCOR

    keyword = 'LINCOR'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # FLATFIELD

    keyword = 'FLATFIELD'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # PLOTXCORR

    keyword = 'PLOTXCORR'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # RECTMETHOD

    keyword = 'RECTMETHOD'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # FIXBADPIXELS

    keyword = 'FIXBADPIXELS'
    z = find_keyword(labels, keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]
    output[keyword] = tmp

    # XSPEXTOOL KEYWORDS

    keyword = 'XSPEXTOOL_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XSPEXTOOL_KEYWORDS'] = tmp

    # XCOMBSPEC KEYWORDS

    keyword = 'XCOMBSPEC_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XCOMBSPEC_KEYWORDS'] = tmp

    # XTELLCOR KEYWORDS

    keyword = 'XTELLCOR_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XTELLCOR_KEYWORDS'] = tmp

    return output


def find_keyword(labels, keyword):
    z = np.where(labels == keyword)
    if np.size(z):

        return z

    else:

        print('Cannot find keyword ', keyword, '.', sep='')
        exit(1)

    return output
