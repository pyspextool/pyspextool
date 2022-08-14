from pyspextool.utils import coords
from pyspextool.utils import add_entry_to_dictionary


def average_header_info(hdrs, pair=False):

    """
    Averages a pySpextool headers lists

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

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """
    # Store the first hdrinfo

    hdrinfo = hdrs[0]
    nfiles = len(hdrs)

    # Now star the averaging

    # ============================ TIME and DATE ==================================

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

    # ============================== HOUR ANGLE ===================================

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

    # ================================== MJD ======================================

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

    # ================================= AIRMASS ===================================

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

    # ================================= EXPTOT ===================================

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
