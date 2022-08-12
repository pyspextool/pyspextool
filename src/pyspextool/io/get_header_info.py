import re


def get_header_info(hdr, keywords=None):

    """
    Pulls (user requested) hdr values and comments from a FITS file


    Parameters
    ----------
    hdr : astropy HDU.header
         see https://docs.astropy.org/en/stable/io/fits/index.html


    keywords : list of str, optional
         A list of keywords to pull.


    Returns
    -------
    dict
        A dict where each key is the FITS keyword and each value is a list 
        with [keyword value,keyword comment].

        To access the value and comment of keyword XXX:
        hdrinfo['XXX'][0] = value
        hdrinfo['XXX'][1] = comment


    Notes
    -----
    The program can use the unix wildcard * for keywords with a common
    basename, e.g. COEFF_1, COEFF_2, COEFF_3 can be obtained as COEFF*.

    
    Examples
    --------
    later


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_gethdrinfo.pro

    """

    # Open an empty dict

    hdrinfo = {}

    # Must do comments and history separately, so they can be at the end
    # in order

    docomment = 0
    dohistory = 0

    # Check to see whether the users passed a set of keywords to grab.

    if keywords is not None:

        # Fill in element with the keyword name and then a list of the
        # value and comment 

        for keyword in keywords:

            # We need to check to see if the request has a wild card.

            m = re.search('[*]', '[' + keyword + ']')

            if m:

                for eight in list(hdr.keys()):

                    n = re.search(keyword[0:-1], eight)
                    if n: hdrinfo[eight] = [hdr[eight], hdr.comments[eight]]

            else:

                if keyword == 'COMMENT':
                    docomment = 1
                    continue

                if keyword == 'HISTORY':
                    dohistory = 1
                    continue

                hdrinfo[keyword] = [hdr[keyword], hdr.comments[keyword]]

    else:

        for name in list(hdr.keys()):

            if name == 'COMMENT':
                docomment = 1
                continue

            if name == 'HISTORY':
                dohistory = 1
                continue

            hdrinfo[name] = [hdr[name], hdr.comments[name]]

            # Now do the comments if need be

    if docomment:

        comments = []
        for line in hdr['COMMENT']:
            comments.append(str(line.replace('=', '')))

        hdrinfo['COMMENT'] = comments

    # Now do the history if need be

    if dohistory:

        history = []
        for line in hdr['HISTORY']:
            history.append(str(line))

        hdrinfo['HISTORY'] = history

    return (hdrinfo)
