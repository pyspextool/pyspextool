import numpy as np


def read_instrument_file(filename:str):

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

    """

    # Read the file into to string arrays

    labels, vals = np.loadtxt(filename,
                              comments='#',
                              delimiter='=',
                              unpack=True,
                              dtype='str')

    # Strip any edge white spaces from poor formatting of the file

    labels = [value.strip() for value in labels]
    labels = np.array(labels)

    vals = [value.strip() for value in vals]
    vals = np.array(vals)

    # Start the search for each keyword.

    output = {}

    # NINT

    keyword = 'NINT'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # BADPIXMASK

    keyword = 'BADPIXMASK'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # SUFFIX

    keyword = 'SUFFIX'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # FITSREADPROGRAM

    keyword = 'READFITS'
    z = find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # LINCORMAX

    keyword = 'LINCORMAX'
    z = find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # XSPEXTOOL KEYWORDS

    keyword = 'EXTRACT_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['EXTRACT_KEYWORDS'] = tmp

    # COMBINE KEYWORDS

    keyword = 'COMBINE_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['COMBINE_KEYWORDS'] = tmp

    # XTELLCOR KEYWORDS

    keyword = 'TELLURIC_KEYWORD'
    z = find_keyword(labels, keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['TELLURIC_KEYWORDS'] = tmp

    return output


def find_keyword(labels,
                 keyword):

    z = np.where(labels == keyword)
    if np.size(z):

        return z

    else:

        print('Cannot find keyword ', keyword, '.', sep='')
        exit(1)

    return output
