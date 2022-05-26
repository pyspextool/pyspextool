from astropy.io import fits
from numpy import empty as npempty
from numpy.polynomial.polynomial import polyval as nppolyval

def readflatinfo(file):

    """
    To read a Spextool flatinfo calibration file.


    Input Parameters
    ----------------
    file : str
        A Spextool flatinfo file, e.g. ShortXD_flatinfo.fits

    Returns
    --------
    dict
       A dictionary with the following keywords:

       rotation : int
                  IDL rotation command for the order numbers to increase 
                  upwards and wavelength increase to the right

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

       slith_arc  : int
                   slit height in arcseconds
       slith_pix  : int
                   nominal slit height in pixels
       slith_rng  : list 
                   [int,int] range of slit heights in pixels
       orders     : list
                   order numbers
       rppix      : float
                   resolving power per pixel
       ps         : float
                   plate scale in arcseconds per pixel
       fixed      : {0,1}
                   set if the top and bottom of the slit align with rows
       step       : int
                   step size in pixels for tracing
       flatfrac   : float
                   see findorders.py
       comwin     : int 
                   see findorders.py
       edgedeg    : int
                   polynomial degree for the edges of the orders
       norm_nxg   : int
                   see normspecflat.py and fiterpolate.py
       norm_nyg   : int
                   see normspecflat.py and fiterpolate.py
       ybuffer    : int
                   see normspecflat.py
       ycororder  : int
                   see adjustguesspos.py
       xranges    : numpy.ndarray, int, [norders,2] 
                   the columns over which extraction can occur
       edgecoeffs : numpy.ndarray, float, [norders,2,edgedeg+1] 
       guesspos   : numpy.ndarray, int, [norders,2] 
                   location of where to start the order search, 
                   see findorders.py

    Procedure
    ---------
    Just reads FITS header information, and calculates the guess positions.

    Example
    --------
    result = readflatinfo('ShortXD_flatinfo.fits')

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """
    
# Open the file, grab the mask
    
    hdul = fits.open(file)
    omask = hdul[0].data

# Clean the header and grab important keywords
    
    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems

    val = hdul[0].header['ROTATION']
    result = {'rotation':val}

    val = hdul[0].header['SLTH_ARC']
    result.update({'slith_arc':val})

    val = hdul[0].header['SLTH_PIX']
    result.update({'slith_pix':val})

    val = hdul[0].header['SLTH_RNG'].split(',')
    val = [int(x) for x in val]
    result.update({'slith_rng':val})

    val = hdul[0].header['ORDERS'].split(',')
    orders = [int(x) for x in val]
    norders = len(orders)
    result.update({'orders':orders})            

    val = hdul[0].header['RPPIX']
    result.update({'rpppix':val})

    val = hdul[0].header['PLTSCALE']
    result.update({'ps':val})

    val = hdul[0].header['FIXED']
    result.update({'fixed':val})        

    if not val:

        val = hdul[0].header['STEP']
        result.update({'step':val})

        val = hdul[0].header['FLATFRAC']
        result.update({'flatfrac':val})

        val = hdul[0].header['COMWIN']
        result.update({'comwin':val})                                

    deg = int(hdul[0].header['EDGEDEG'])
    result.update({'edgedeg':deg})

    val = hdul[0].header['NORM_NXG']
    result.update({'norm_nxg':val})

    val = hdul[0].header['NORM_NYG']
    result.update({'norm_nyg':val})

#    val = hdul[0].header['OVERSAMP']
#    result.update({'oversamp':val})

    val = hdul[0].header['YBUFFER']
    result.update({'ybuffer':val})

    val = hdul[0].header['YBUFFER']
    result.update({'ybuffer':val})

    val = hdul[0].header['YCORORDR']
    result.update({'ycororder':val})                                    

# Now get the edge coefficients and the xranges

    xranges = npempty((norders,2),dtype=int)
    edgecoeffs = npempty((norders,2,deg+1))
    guesspos = npempty((norders,2),dtype=int)
    
    for i in range(0,norders):

# Get the xrange and guess position x position
        
        val = hdul[0].header['OR'+str(orders[i]).zfill(3)+'_XR'].split(',')
        val = [int(x) for x in val]
        xranges[i,:] = val

        guesspos[i,0] = sum(xranges[i,:])/2

# Now grab the edgecoefficients
        
        for j in range(0,deg+1):

            keyt = 'OR'+str(orders[i]).zfill(3)+'_T'+str(j+1).zfill(1)
            keyb = 'OR'+str(orders[i]).zfill(3)+'_B'+str(j+1).zfill(1)

            edgecoeffs[i,0,j] = hdul[0].header[keyb]
            edgecoeffs[i,1,j] = hdul[0].header[keyt]

# Now determine the guess position y position 
            
        bot = nppolyval(guesspos[i,0],edgecoeffs[i,0,:])
        top = nppolyval(guesspos[i,0],edgecoeffs[i,1,:])

        guesspos[i,1] = int((bot+top)/2)


    result.update({'xranges':xranges})
    result.update({'edgecoeffs':edgecoeffs})
    result.update({'guesspos':guesspos})
    hdul.close()
    

    return(result)
