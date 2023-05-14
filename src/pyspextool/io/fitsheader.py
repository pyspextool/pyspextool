import re
from pyspextool.utils import coords
from pyspextool.utils.add_entry import add_entry

def average_header_info(hdrs, pair=False):

    """
    Averages a pySpextool header lists

    Parameters
    ----------
    hdrs : list
        list of dictionaries where each element is of the form:
        {key:[val,com]} where key is the FITS keyword, "val" is the 
        value, and "com" is the comment.

    pair : {False, True}, optional 
        Set to True if the Spextool header list was generated for 
        pair subtraced images

    Returns
    --------
    dict
        a pySpextool header list with the values for TIME,
        DATE, HA, MJD, AIRMASS replaced with start, average, and end 
        values for each keyword, e.g. SRT_HA, AVE_HA, END_HA.  A new 
        keyword, ITOTIMAGE gives the total exposure time.

    Procedure
    ---------
    Mostly a bunch of string parsing and watching for oddities the signs 
    of sexigesimal numbers.

    Examples
    --------
    Later

    """
    # Store the first hdrinfo

    hdrinfo = hdrs[0]
    nfiles = len(hdrs)

    # Now star the averaging

    # ============================ TIME and DATE ==============================

    # Grab the values

    times = [item["TIME"][0] for item in hdrs]

    dates = [item["DATE"][0] for item in hdrs]

    # Determine if you are observing over a day or not

    days = list(set(dates))
    ndays = len(days)

    # Convert to decimal hours

    dvals = []
    for i in range(0, nfiles):
        dvals.append(coords.ten(times[i]))

    if ndays == 1:

        ave = sum(dvals) / nfiles

        # Convert back to sexigesimal

        avetime = coords.sixty(ave, colons={'dec': 2, 'plus': 1})
        avedate = dates[0]

    elif ndays == 2:

        zday2 = wherelist(dates, days[1])

        # Add 24 to the second day times

        for i in zday2:
            dvals[i] += 24

        # Compute the average 

        ave = sum(dvals) / nfiles

        # Check to see what date the average falls on and adjust 

        if ave >= 24:

            avetime = coords.sixty(ave - 24, colons={'dec': 2, 'plus': 0})
            avedate = dates[1]

        else:

            avetime = coords.sixty(ave, colons={'dec': 2, 'plus': 0})
            avedate = dates[0]

    else:
        print('Observations over three days not allowed.')

    # Add the new keywords

    hdrinfo = add_entry(hdrinfo, 'TIME', 'after', 'SRT_TIME',
                        [times[0], 'Start observations time in UTC'])

    hdrinfo = add_entry(hdrinfo, 'SRT_TIME', 'after', 'AVE_TIME',
                        [avetime, 'Average observation time in UTC'])

    hdrinfo = add_entry(hdrinfo, 'AVE_TIME', 'after', 'END_TIME',
                        [times[-1], 'End observation time in UTC'])

    # Remove the current keyword

    hdrinfo.pop('TIME')

    # Add the new keywords

    hdrinfo = add_entry(hdrinfo, 'DATE', 'after', 'SRT_DATE',
                        [dates[0], 'Start observations date in UTC'])

    hdrinfo = add_entry(hdrinfo, 'SRT_DATE', 'after', 'AVE_DATE',
                        [avedate, 'Average observation date in UTC'])

    hdrinfo = add_entry(hdrinfo, 'AVE_DATE', 'after', 'END_DATE',
                        [dates[-1], 'End observation date in UTC'])

    # Remove the current keyword

    hdrinfo.pop('DATE')

    # ============================== HOUR ANGLE ===============================

    # Grab the values

    vals = [item["HA"][0] for item in hdrs]

    # Convert to decimal hours

    dvals = []
    for i in range(0, nfiles):
        dvals.append(coords.ten(vals[i]))

    ave = sum(dvals) / nfiles

    # Convert back to sexigesimal

    avetime = coords.sixty(ave, colons={'dec': 2, 'plus': 1})

    # Add the new keywords

    hdrinfo = add_entry(hdrinfo, 'HA', 'after', 'SRT_HA',
                        [vals[0], 'Start hour angle (hours)'])

    hdrinfo = add_entry(hdrinfo, 'SRT_HA', 'after', 'AVE_HA',
                        [avetime, 'Average hour angle (hours)'])

    hdrinfo = add_entry(hdrinfo, 'AVE_HA', 'after', 'END_HA',
                        [vals[-1], 'End hour angle (hours)'])

    # Remove the current keyword

    hdrinfo.pop('HA')

    # ================================== MJD ==================================

    # Grab the values

    vals = [item["MJD"][0] for item in hdrs]

    # Get the number of significant digits

    junk = str(vals[0]).split('.')
    ndigits = len(junk[1])
    fmt = '{:.' + str(ndigits) + 'f}'

    # Add the new keywords

    hdrinfo = add_entry(hdrinfo, 'MJD', 'after', 'SRT_MJD',
                        [vals[0], 'Start modified Julian date'])

    hdrinfo = add_entry(hdrinfo, 'SRT_MJD', 'after', 'AVE_MJD',
                        [fmt.format(sum(vals) / nfiles),
                        'Average modified Julian date'])

    hdrinfo = add_entry(hdrinfo, 'AVE_MJD', 'after', 'END_MJD',
                        [vals[-1], 'End airmass'])

    # Remove the current keyword

    hdrinfo.pop('MJD')

    # ================================= AIRMASS ===============================

    # Grab the values

    vals = [item["AM"][0] for item in hdrs]

    # Get the number of significant digits

    junk = str(vals[0]).split('.')
    ndigits = len(junk[1])
    fmt = '{:.' + str(ndigits) + 'f}'

    hdrinfo = add_entry(hdrinfo, 'AM', 'after', 'SRT_AM',
                        [vals[0], 'Start airmass'])

    hdrinfo = add_entry(hdrinfo, 'SRT_AM', 'after', 'AVE_AM',
                        [fmt.format(sum(vals) / nfiles), 'Average airmass'])

    hdrinfo = add_entry(hdrinfo, 'AVE_AM', 'after', 'END_AM',
                        [vals[-1], 'End airmass'])

    # Remove the current keyword

    hdrinfo.pop('AM')

    # ================================= EXPTOT ===============================

    # Grab the values

    vals = [item["IMGITIME"][0] for item in hdrs]

    if pair is True:

        tot = vals[0] * (nfiles / 2)
        com = 'Total integration time PER BEAM (sec)'

    else:

        tot = sum(vals)
        com = 'Total integration time (sec)'

    hdrinfo = add_entry(hdrinfo, 'IMGITIME', 'after', 'TOTITIME', [tot, com])

    return hdrinfo


def get_header_info(hdr, keywords=None, ignore_missing_keywords=False):

    """
    Pulls (user requested) keyword values and comments from a FITS file


    Parameters
    ----------
    hdr : astropy HDU.header
        see https://docs.astropy.org/en/stable/io/fits/index.html


    keywords : list of str, optional
        A list of keywords to pull.

    ignore_missing_keywords : {False, True}, optional
        Set to True to ignore keywords not present in the hdr. 


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
                    if n:

                        # Test if it exits

                        test = keyword in hdr
                        if test is True:
                
                            hdrinfo[eight] = [hdr[eight], hdr.comments[eight]]

            else:

                if keyword == 'COMMENT':
                    docomment = 1
                    continue

                if keyword == 'HISTORY':
                    dohistory = 1
                    continue

                # Test if it exits

                test = keyword in hdr
                if test is True:
                
                    hdrinfo[keyword] = [hdr[keyword], hdr.comments[keyword]]

    else:

        for name in list(hdr.keys()):

            if name == '':
                continue
            
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
