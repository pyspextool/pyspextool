from astropy.io import fits
import re

def mc_gethdrinfo(hdr,keywords=None):

    '''
    Pulls (user requested) hdr values and comments from a FITS file

    Input Parameters:
        hdr - an astropy header of an HDU
              see https://docs.astropy.org/en/stable/io/fits/index.html

    Optional Parameters:
        keywords - a list of strings of keywords to obtain.

    Output Parameters:
        A dict where each key is the FITS keyword and each value is a list 
        with [keyword value,keyword comment].

        To access the value and comment of keyword XXX:
        hdrinfo['XXX'][0] = value
        hdrinfo['XXX'][1] = comment

    Procedure:

    Example:
        NA

    Modification History:
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the gethdrinfo.pro IDL program.
    '''

# Open an empty dict
    
    hdrinfo = {}

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

               hdrinfo[keyword] = [hdr[keyword],hdr.comments[keyword]]


    else:

        for name in list(hdr.keys()):

           hdrinfo[name] = [hdr[name],hdr.comments[name]]            
           
    return(hdrinfo)
