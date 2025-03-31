import numpy as np
import fnmatch


from pyspextool.io.check import check_parameter
from pyspextool.utils import coords
from pyspextool.utils.add_entry import add_entry

def average_headerinfo(hdrs:list,
                       pair:bool=False):

    """
    Averages a pySpextool header lists

    Parameters
    ----------
    hdrs : list
        (nheaders,) list of dictionaries where each dictionary element is of 
        the form: {key:[val,com]} where key is the FITS keyword, "val" is the 
        value, and "com" is the comment.

    pair : {False, True}
        Set to True if the Spextool header list was generated for 
        pair subtraced images.

    Returns
    -------
    dict
        a pySpextool header list with the values for TIME,
        DATE, HA, MJD, AIRMASS replaced with start, average, and end 
        values for each keyword, e.g. SRT_HA, AVE_HA, END_HA.  A new 
        keyword, ITOTIMAGE gives the total exposure time.


    """

    #
    # Check parameters
    #

    check_parameter('average_headerinfo', 'hdrs', hdrs, 'list')

    check_parameter('average_headerinfo', 'pair', pair, 'bool')


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

        avetime = coords.sixty(ave, colons={'dec': 2, 'plus': 0})
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
    nhourangles = 0
    for i in range(0, nfiles):
    
        if vals[i].strip() == 'nan':

            continue

        else:
            
            dvals.append(coords.ten(vals[i]))
            nhourangles += 1

    if nhourangles == 0:

        # All nans

        start = 'nan'
        ave = 'nan'
        end = 'nan'

    else:

        ave = sum(dvals) / nhourangles

        # Convert back to sexigesimal

        ave = coords.sixty(ave, colons={'dec': 2, 'plus': 1})

        start = vals[0]
        end = vals[-1]



    # Add the new keywords

    hdrinfo = add_entry(hdrinfo, 'HA', 'after', 'SRT_HA',
                        [start, 'Start hour angle (hours)'])

    hdrinfo = add_entry(hdrinfo, 'SRT_HA', 'after', 'AVE_HA',
                        [ave, 'Average hour angle (hours)'])

    hdrinfo = add_entry(hdrinfo, 'AVE_HA', 'after', 'END_HA',
                        [end, 'End hour angle (hours)'])

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
                        [vals[0], 'Start Modified Julian Date'])

    hdrinfo = add_entry(hdrinfo, 'SRT_MJD', 'after', 'AVE_MJD',
                        [float(fmt.format(sum(vals) / nfiles)),
                        'Average Modified Julian Date'])

    hdrinfo = add_entry(hdrinfo, 'AVE_MJD', 'after', 'END_MJD',
                        [vals[-1], 'End Modified Julian Date'])

    # Remove the current keyword

    hdrinfo.pop('MJD')

    # ================================= AIRMASS ===============================

    # Grab the values

    vals = [item["AM"][0] for item in hdrs]
    nvals = len(vals)
    nanmask = np.isnan(vals)

    if nvals == np.sum(nanmask):

        # They are all NaNs 

        start = np.nan
        ave = np.nan
        end = np.nan
        ndigits = 1

    else:

        # At least 1 isn't a NaN
        
        start = vals[0]
        ave = np.nanmean(vals)
        end = vals[-1]

        junk = str(vals[int(~nanmask[0])]).split('.')
        ndigits = len(junk[1])
        
    # Add the results 
    
    fmt = '{:.' + str(ndigits) + 'f}'

    hdrinfo = add_entry(hdrinfo, 'AM', 'after', 'SRT_AM',
                        [start, 'Start airmass'])

    hdrinfo = add_entry(hdrinfo, 'SRT_AM', 'after', 'AVE_AM',
                        [float(fmt.format(ave)),
                         'Average airmass'])

    hdrinfo = add_entry(hdrinfo, 'AVE_AM', 'after', 'END_AM',
                        [end, 'End airmass'])

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


def get_headerinfo(hdr,
                   keywords=None,
                   ignore_missing_keywords=False):

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

    #
    # Check parameters
    #

    check_parameter('get_header_info', 'hdr', hdr, 'Header')

    check_parameter('get_header_info', 'keywords', keywords,
                    ['list','NoneType'])

    check_parameter('get_header_info', 'ignore_missing_keywords',
                    ignore_missing_keywords, 'bool')    

    #
    # Get setup
    #
    
    # Must do comments and history separately, so they can be at the end
    # in order

    docomment = False
    dohistory = False
    
    # Open an empty dict

    hdrinfo = {}

    # Grab the headers' keywords

    header_keywords = list(hdr.keys())

    # Does the user pass a list of keywords?  
    
    if keywords is not None:

        # Loop over each keyword passed by the user.

        for keyword in keywords:

            # Do the basic checks
            
            if keyword == '':
                continue
            
            if keyword == 'COMMENT':
                docomment = True
                continue

            if keyword == 'HISTORY':
                dohistory = True
                continue

            # Now start the search of that keyword against the header

            matches = []
            for header_keyword in header_keywords:

                test = fnmatch.fnmatch(header_keyword, keyword)
                if test is True:

                    matches.append(header_keyword)

            # Does it find any?
                    
            if len(matches) == 0:

                # It does not.  Now proceed as the user requests.
                
                if ignore_missing_keywords is False:

                    hdrinfo[keyword] = [None,'']

            else: 

                # It does.  Store the results
                
                for match in matches:

                    hdrinfo[match] = [hdr[match], hdr.comments[match]]
                
    else:

        # Just store each keyword in the header
        
        for name in header_keywords:

            # Do basic tests
            
            if name == '':
                continue
            
            if name == 'COMMENT':
                docomment = 1
                continue

            if name == 'HISTORY':
                dohistory = 1
                continue

            # Store the value
            
            hdrinfo[name] = [hdr[name], hdr.comments[name]]

    #        
    # Now do the history and comments if need be
    #
    
    if docomment is True:

        comments = []
        for line in hdr['COMMENT']:
            comments.append(str(line.replace('=', '')))

        hdrinfo['COMMENT'] = comments

    # Now do the history if need be

    if dohistory is True:

        history = []
        for line in hdr['HISTORY']:
            history.append(str(line))

        hdrinfo['HISTORY'] = history

    return (hdrinfo)
