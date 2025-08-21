import numpy as np


def read_instrument_file(
    filename:str):

    """
    To read a Spextool instrument configuration file.  

    Parameters
    ----------
    filename : str
        The fullname path to the Spextool instrument file.

    Returns
    --------
    dict
        Later, when you finalize the pySpextool cal file

    """

    # Read the file into to string arrays

    labels, vals = np.loadtxt(
        filename,
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
    z = _find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # BADPIXMASK

    keyword = 'BADPIXMASK'
    z = _find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # SUFFIX

    keyword = 'SUFFIX'
    z = _find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # FITSREADPROGRAM

    keyword = 'READFITS'
    z = _find_keyword(labels, keyword)
    output[keyword] = vals[z].item()

    # LINCORMAX

    keyword = 'LINCORMAX'
    z = _find_keyword(labels, keyword)
    output[keyword] = int(vals[z].item())

    # KEYWORDS

    keyword = 'KEYWORD'
    z = _find_keyword(labels, keyword)

    output['KEYWORDS'] = []
    output['COMBINE_IGNORE_KEYWORDS'] = []
    output['TELLURIC_IGNORE_KEYWORDS'] = []
    output['MERGE_IGNORE_KEYWORDS'] = []
    for value in vals[z]:

        list = value.split(' ')

        output['KEYWORDS'].append(list[0].strip())
        if 'combine' not in list[1::]:

            output['COMBINE_IGNORE_KEYWORDS'].append(list[0].strip())

        if 'telluric' not in list[1::]:

            output['TELLURIC_IGNORE_KEYWORDS'].append(list[0].strip())

        if 'merge' not in list[1::]:

            output['MERGE_IGNORE_KEYWORDS'].append(list[0].strip())


    return output


def _find_keyword(
    labels,
    keyword):

    z = np.where(labels == keyword)
    if np.size(z):

        return z

    else:

        print('Cannot find keyword ', keyword, '.', sep='')
        exit(1)

    return output
