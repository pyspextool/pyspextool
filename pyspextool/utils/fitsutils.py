"""Functions for FITS utlitiy management."""


from astropy.io import fits
import re
from sys import exit
from astropy.visualization import PercentileInterval, ZScaleInterval, MinMaxInterval
import numpy as np
from math import nan as mnan
from numpy import rot90 as nprot90
from numpy import flipud as npflipud
from numpy import fliplr as npfliplr
from numpy import transpose as nptranspose

from pyspextool.utils.text import dictaddentry, wherelist
from pyspextool.utils.image import ten, sixty


def avehdrs(hdrs,pair=False):

    """
    Averages a pySpextool headers lists

    Input Parameters
    ----------------
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
        a pySpextool header list with the values for values for TIME, 
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
    
    #============================ TIME and DATE ==================================

    # Grab the values
    
    times = [item["TIME"][0] for item in hdrs]

    dates = [item["DATE"][0] for item in hdrs]

    # Determine if you are observing over a day or not

    days = list(set(dates))
    ndays = len(days)

    # Convert to decimal hours
    
    dvals = []
    for i in range(0,nfiles):
            
        dvals.append(ten(times[i]))
            
    if ndays == 1:
        
        ave = sum(dvals)/nfiles

        # Convert back to sexigesimal

        avetime = sixty(ave,colons={'dec':2,'plus':1})
        avedate = dates[0]
        
    elif ndays == 2:

        zday2 = wherelist(dates,days[1])

        # Add 24 to the second day times

        for i in zday2:

            dvals[i] += 24

        # Compute the average 
            
        ave = sum(dvals)/nfiles

        # Check to see what date the average falls on and adjust 

        if ave >= 24:

            avetime = sixty(ave-24,colons={'dec':2,'plus':0})
            avedate = dates[1]
                        
        else:

            avetime = sixty(ave,colons={'dec':2,'plus':0})
            avedate = dates[0]        
        
    else: print('Observations over three days not allowed.')
    
    # Add the new keywords

    hdrinfo = dictaddentry(hdrinfo,'TIME','after','SRT_TIME',\
                           [times[0],'Start observations time in UTC'])

    hdrinfo = dictaddentry(hdrinfo,'SRT_TIME','after','AVE_TIME',\
                           [avetime,'Average observation time in UTC'])

    hdrinfo = dictaddentry(hdrinfo,'AVE_TIME','after','END_TIME',\
                           [times[-1],'End observation time in UTC'])

    # Remove the current keyword

    hdrinfo.pop('TIME')

    # Add the new keywords

    hdrinfo = dictaddentry(hdrinfo,'DATE','after','SRT_DATE',\
                           [dates[0],'Start observations date in UTC'])

    hdrinfo = dictaddentry(hdrinfo,'SRT_DATE','after','AVE_DATE',\
                           [avedate,'Average observation date in UTC'])

    hdrinfo = dictaddentry(hdrinfo,'AVE_DATE','after','END_DATE',\
                           [dates[-1],'End observation date in UTC'])

    # Remove the current keyword

    hdrinfo.pop('DATE')            
    
    #============================== HOUR ANGLE ===================================

    # Grab the values

    vals = [item["HA"][0] for item in hdrs]

    # Convert to decimal hours
    
    dvals = []
    for i in range(0,nfiles):

        dvals.append(ten(vals[i]))

    ave = sum(dvals)/nfiles

    # Convert back to sexigesimal

    avetime = sixty(ave,colons={'dec':2,'plus':1})
    
    # Add the new keywords

    hdrinfo = dictaddentry(hdrinfo,'HA','after','SRT_HA',\
                           [vals[0],'Start hour angle (hours)'])

    hdrinfo = dictaddentry(hdrinfo,'SRT_HA','after','AVE_HA',\
                           [avetime,'Average hour angle (hours)'])

    hdrinfo = dictaddentry(hdrinfo,'AVE_HA','after','END_HA',\
                           [vals[-1],'End hour angle (hours)'])

    # Remove the current keyword

    hdrinfo.pop('HA')    
 
    #================================== MJD ======================================     

    # Grab the values

    vals = [item["MJD"][0] for item in hdrs]    
    
    # Get the number of significant digits

    junk = str(vals[0]).split('.')
    ndigits = len(junk[1])
    fmt = '{:.'+str(ndigits)+'f}'

    # Add the new keywords

    hdrinfo = dictaddentry(hdrinfo,'MJD','after','SRT_MJD',\
                           [vals[0],'Start modified Julian date'])

    hdrinfo = dictaddentry(hdrinfo,'SRT_MJD','after','AVE_MJD',\
                           [fmt.format(sum(vals)/nfiles),\
                            'Average modified Julian date'])

    hdrinfo = dictaddentry(hdrinfo,'AVE_MJD','after','END_MJD',\
                           [vals[-1],'End airmass'])

    # Remove the current keyword

    hdrinfo.pop('MJD')    

    #================================= AIRMASS ===================================      

    # Grab the values

    vals = [item["AM"][0] for item in hdrs]

    # Get the number of significant digits

    junk = str(vals[0]).split('.')
    ndigits = len(junk[1])
    fmt = '{:.'+str(ndigits)+'f}'

    hdrinfo = dictaddentry(hdrinfo,'AM','after','SRT_AM',\
                           [vals[0],'Start airmass'])

    hdrinfo = dictaddentry(hdrinfo,'SRT_AM','after','AVE_AM',\
                           [fmt.format(sum(vals)/nfiles), 'Average airmass'])

    hdrinfo = dictaddentry(hdrinfo,'AVE_AM','after','END_AM',\
                           [vals[-1],'End airmass'])
                           
    # Remove the current keyword

    hdrinfo.pop('AM')

    #================================= EXPTOT ===================================

    # Grab the values

    vals = [item["IMGITIME"][0] for item in hdrs]
    
    if pair is True:

        tot = vals[0]*(nfiles/2)
        com = 'Total integration time PER BEAM (sec)'
        
    else:

        tot = sum(vals)
        com = 'Total integration time (sec)'    


    hdrinfo = dictaddentry(hdrinfo,'IMGITIME','after','TOTITIME',[tot, com])     
    
    return(hdrinfo)


def gethdrinfo(hdr,keywords=None):

    '''
    Pulls (user requested) hdr values and comments from a FITS file


    Input Parameters
    ----------------
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

    '''

    # Open an empty dict
    
    hdrinfo = {}

    # Must do comments and history separately so they can be at the end
    # in order

    docomment = 0
    dohistory = 0

    # Check to see whether the users passed a set of keywords to grab.

    if keywords is not None:
    
        # Fill in element with the keyword name and then a list of the
        # value and comment 

       for keyword in keywords:

            # We need to check to see if the request has a wild card.

           m = re.search('[*]','['+keyword+']')
           
           if m:

               for eight in list(hdr.keys()):

                   n = re.search(keyword[0:-1],eight)
                   if n: hdrinfo[eight] = [hdr[eight],hdr.comments[eight]]
               
           else:

               if keyword == 'COMMENT':

                   docomment = 1
                   continue

               if keyword == 'HISTORY':

                   dohistory = 1
                   continue               
               
               hdrinfo[keyword] = [hdr[keyword],hdr.comments[keyword]]

    else:

        for name in list(hdr.keys()):

            if name == 'COMMENT':

                docomment = 1
                continue

            if name == 'HISTORY':

                dohistory = 1
                continue                           
        
            hdrinfo[name] = [hdr[name],hdr.comments[name]]            

    # Now do the comments if need be
            
    if docomment:

        comments = []
        for line in hdr['COMMENT']:

            comments.append(str(line.replace('=','')))

        hdrinfo['COMMENT'] = comments

    # Now do the history if need be
        
    if dohistory:

        history = []        
        for line in hdr['HISTORY']:

            history.append(str(line))

        hdrinfo['HISTORY'] = history
        
    return(hdrinfo)


def getimgrange(arr,info):

    """
    a wrapper for astropy to get image ranges

    Input Parameters
    ----------------
    arr : numpy.npdarray
        an array (typically a 2D image)

    info : float or str
        float - the fraction of pixels to keep in the range (0,1)
        str - 'zscale' or 'minmax'

    Returns
    --------
    tuple
         interval based on user request

    Procedure
    ---------
    calls astropy.visualization routines and uses all default settings


    Examples
    --------
    later
    

    Modification History
    --------------------
    2022-06-07 - Written by M. Cushing, University of Toledo.

    """
    
    if isinstance(info,(float)) is True:

        # This is if the users requests a fraction

        interval = PercentileInterval(info)
        range = interval.get_limits(arr)
        return range

    elif info == 'zscale':

        ## This is if the users requests the zscale
        
        interval = ZScaleInterval()
        range = interval.get_limits(arr)
        return range

    elif info == 'minmax':

        ## This is if the users requests the zscale

        interval = MinMaxInterval()
        range = interval.get_limits(arr)
        return range

    else:

        return None


def idlrotate(img,direction):

    """
    Rotates and/or tranposes an image (rotate is in multiples of 90deg)

    Input Parameters
    ----------------
    img : numpy.ndarray
        an image

    direction : int

        Direction  Transpose?  Rotation Counterclockwise
        -------------------------------------------------

        0          No          None
        1          No          90 deg
        2          No          180 deg
        3          No          270 deg
        4          Yes         None
        5          Yes         90 deg
        6          Yes         180 deg
        7          Yes         270 deg

    The directions follow the IDL rotate function convention.


    Returns
    --------
    numpy.ndarray
        the rotated/transposed image

    Procedure
    ---------

    uses numpy.rot90, flipud, fliplr, ant transpose to rotate and/or 
    transpose the image as requested.

    Examples
    --------

    > import numpy as np
    > img = np.array([[1,2,3],[4,5,6],[7,8,9]])
    > idlrotate(img,6)

    [[9 6 3]
    [8 5 2]
    [7 4 1]]

    Modification History
    --------------------
    2022-06-02 - Written by M. Cushing, University of Toledo.
    Based on the IDL program rotate.pro

    """
    
    if direction == 0:

        return(img)        

    elif direction == 1:

        return(nprot90(img,3))

    elif direction == 2:

        return(nprot90(img,2))

    elif direction == 3:

        return(nprot90(img,1))
        
    elif direction == 4:

        return(nptranspose(img))
        
    elif direction == 5:

        return(npfliplr(img))        
        
    elif direction == 6:

        return(npfliplr(nprot90(img,1)))
    
    elif direction == 7:               

        return(npflipud(img))        
        
    else:

        print('idlrotate:  Unknown direction.')
        exit(1)


def medcomb(data,mask=None,stderr=True):

    """
    Median a spectral or image stack with optional mask


    Input Parameters
    ----------------
    data : numpy.ndarray
        either a stack of spectra [nspec,npoints] or a stack of images 
        [nimgs,nrows,ncols].  

    mask : numpy.ndarray, optional
        a mask array with the same shape as `data`.  
        0 = bad, 1=good

    stderr : {True, False}, optional
        Set to return 1.4826*MAD/sqrt(n) instead of just 1.4826*MAD 
        (see Procedure)

    Returns
    --------
    list
        list[0] : numpy.ndarray 
            the median of the spectral or image stack

        list[1] : numpy.ndarray
            the uncertainty of the spectral or image stack (see Procedure)

    Procedure
    ---------
    Spectral stack:

        in this case, the data have the shape [nspec,npoints].  The 
        median of the stack is computed producing an array of size 
        [npoints]. At each spectral point, the median absolute deviation 
        (MAD=median(|data-med}) is computed.  

    Image stack:

        in this case, the data have the shape [nspec,nrows,ncols].  The 
        median of the stack is computed producing an array of size 
        [nrows,ncols]. At each image point, the median absolute deviation 
        (MAD=median(|data-med}) is computed.  


    The estimate of the standard deviation assuming the data arises from 
    a gaussian is given by 1.4826*MAD.  Finally, if stderr is set, then the 
    standard error is computed as 1.4826*MAD/root(n), where n is the number 
    of data points at a given spectral point.

    Note:  Points excluded by the mask are ignore in all calculations.

    Examples
    --------

    > import numpy as np
    > ss = np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
    > msk = np.ones((3,4),dtype=int)
    > msk[0,0] = 0
    > print(ss)
    > print(msk)
    > med,unc=medcomb(ss,mask=msk)
    > print(med)
    > print(unc)

      [[1 2 3 7]
       [0 5 8 2]
       [2 9 4 6]]
      [[0 1 1 1]
       [1 1 1 1]
       [1 1 1 1]]
      [1. 5. 4. 6.]
      [1.04835651 2.56793853 0.85597951 0.85597951]

    > istack = np.array([[[1,2,3],[4,5,6],[7,8,9]],\
                         [[6,3,1],[9,2,4],[1,5,0]],\
                         [[3,4,9],[5,7,7],[3,9,1]],\
                         [[1,6,5],[2,1,9],[5,2,7]]])              
    > msk = np.ones((4,3,3),dtype=int)
    > msk[0,0,0] = 0
    > print('Image Stack Test')
    > print(istack)
    > print(msk)
    > med,unc=medcomb(istack,mask=msk)
    > print(med)
    > print(unc)

      [[[1 2 3]
        [4 5 6]
        [7 8 9]]

       [[6 3 1]
        [9 2 4]
        [1 5 0]]

       [[3 4 9]
        [5 7 7]
        [3 9 1]]

       [[1 6 5]
        [2 1 9]
        [5 2 7]]]
      [[[0 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]

       [[1 1 1]
        [1 1 1]
        [1 1 1]]]
      [[3.  3.5 4. ]
       [4.5 3.5 6.5]
       [4.  6.5 4. ]]
      [[1.71195902 0.7413     1.4826    ]
       [1.11195    1.4826     1.11195   ]
       [1.4826     1.4826     2.59455   ]]

    Modification History
    --------------------
    2022-06-01 - Written by M. Cushing, University of Toledo.
        Based on the Spextool IDL program mc_medcomb.pro.

    """

    # Get array dimensions

    ndimen = np.ndim(data)
    shape = np.shape(data)

    # If no mask passed, create one.

    if mask is None:

        mask = np.ones(shape,dtype=int)

    # Now search and replace any masked pixels with NaNs
        
    data = np.where(mask != 0, data,mnan)

    # Spectral or image stack?

    if ndimen == 2:

        tileshape = (shape[0],1) # spectral stack

    elif ndimen == 3:

        tileshape = (shape[0],1,1) # image stack

    else:

        print('Unknown data shape.')
        return
    
    # Compute the median and median absolute deviation

    med = np.nanmedian(data,axis=0)

    mad = np.nanmedian(np.abs(data-np.tile(med,tileshape)),axis=0)

    if stderr is not None:
            
        mad *= 1.4826   # assume gaussian distribution
        unc = mad/np.sqrt(np.sum(mask, axis=0))
            
    else:

        unc = mad / 1.4826   # assume gaussian distribution            

    return[med,unc]