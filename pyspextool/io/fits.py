"""Functions for reading/writing FITS files."""


from pyspextool.utils.text import dictaddentry
from astropy.io import fits
import numpy as np
from numpy import empty as npempty
from numpy.polynomial.polynomial import polyval as nppolyval
from numpy import loadtxt as nploadtxt
from numpy import array as nparray
from numpy import where as npwhere
from numpy import size as npsize
from sys import exit as exit


def readflat(flatfile):

    '''
    To read a pySpextool FITS flat image into memory.


    Input Parameters
    ----------------
    flatfile : str
        The full path to a pySpextool FITS flat file.


    Returns
    -------
    dict 
        ``"flat"``
            (numpy.ndarray) -> (nrows,cols) flat image

        ``"var"``
            (numpy.ndarray) -> (nrows,cols) variance image

        ``"bitmask"``
            (numpy.ndarray of int) -> (nrows,cols) bitmask image.  First 
            bit is set if a pixel in the original stack was over the 
            linmax limit.

        ``"ordermask"``
            (numpy.ndarray of int) -> (nrows,cols) order mask image.  
            A pixel is set to its order number.

        ``"ncols"``
            float -> The number of columns.

        ``"nrows"``
            float -> The number of rows.

        ``"mode"``
            str -> The instrument mode name.

        ``"norders"``
            int -> The number of orders on the image.

        ``"orders"``
            numpy.ndarray of int -> The orders numbers.  orders[0] is 
            the order number of the order nearest the bottom of the 
            image after rotation. 

        ``"edgedeg"``
            int -> The polynomial degree for the fits to the edge of 
            the order.

        ``"ps"``
            float -> The plate scale (arcseconds per pixel).

        ``"slith_pix"``
            float -> The slit height (pixels).

        ``"slith_arc"``
            float -> The nominal slit height (arcseconds).

        ``"slitw_arc"``
            float -> The slit width (arcseconds).

        ``"slitw_pix"`` 
            float -> The slit width (pixels).

        ``"rp"``
            float -> The nominal resolving power.

        ``"rotation"``
            {0, 1, 2, 3, 4, 5, 6, 7} -> Direction to rotate a raw image 
            so that the dispersion direction roughly aligns with the 
            columns of the array and wavelength increases to the right, 
            and the spatial axis roughly aligns with the rows of the array.  
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

        ``"edgecoeffs"``
            array_like of float -> (norders,`edgedeg`+1,2) array giving 
            the polynomial coefficients delineating the top and bottom 
            of each order.  edgecoeffs[0,0,:] gives the coefficients for 
            the bottom of the order closest to the bottom of the image 
            and edgecoeffs[0,1,:] gives the coefficients for the top of 
            said order.  

        ``"xranges"``
            array_like of float -> An (norders,2) array giving the 
            column numbers over which to operate.  xranges[0,0] gives 
            the starting column number for the order nearest the bottom 
            of the image and xranges[0,1] gives the end column number 
            for said order.

        ``"rms"``
            list of float -> (norders,) list of RMS values for each order.



    Notes
    -----
    None


    Examples
    --------
    later

    Modification History
    --------------------
    2022-06-28 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_readflat.pro.

    '''

    # Read the data 
    
    hdul = fits.open(flatfile)
    hdul[0].verify('silentfix')

    hdr = hdul[0].header

    flat = hdul[1].data
    var = hdul[2].data
    mask = hdul[3].data

    hdul.close()

    # create flatinfo dictionary

    flatinfo = {'flat':flat}
    flatinfo.update({'var':var})
    flatinfo.update({'bitmask':mask})    

    shape = np.shape(flat)    
    flatinfo.update({'ncols':shape[1]})
    flatinfo.update({'nrows':shape[0]})
    
    flatinfo.update({'mode':hdr['MODE']})

    norders = hdr['NORDERS']
    flatinfo.update({'norders':norders})

    orders = hdr['ORDERS'].split(',')
    orders = [int(o) for o in orders]
    flatinfo.update({'orders':np.array(orders,dtype=int)})

    edgedeg = hdr['EDGEDEG']
    flatinfo.update({'edgedeg':edgedeg})

    flatinfo.update({'ps':hdr['PLTSCALE']})
    flatinfo.update({'slith_arc':hdr['SLTH_ARC']})
    flatinfo.update({'slith_pix':hdr['SLTH_PIX']})
    flatinfo.update({'slitw_arc':hdr['SLTW_ARC']})
    flatinfo.update({'slitw_pix':hdr['SLTW_PIX']})
    flatinfo.update({'rp':hdr['RP']})
    flatinfo.update({'rotation':hdr['ROTATION']})

    # Grab the edge coeffiecients, xranges, and rms values

    ordermask = np.zeros([flatinfo['nrows'],flatinfo['ncols']],dtype=int)
    edgecoeffs = np.empty([norders,2,edgedeg])
    xranges = np.empty([norders,2],dtype=int)
    rms = np.empty([norders])

    for i in range(norders):

        root = 'OR'+str(orders[i]).zfill(3)

        for j in range(edgedeg):

            edgecoeffs[i,0,j] = hdr[root+'_B*'][j]
            edgecoeffs[i,1,j] = hdr[root+'_T*'][j]

        xranges[i,:]= [int(x) for x in hdr[root+'_XR'].split(',')]

        #  May not have an RMS if it wasn't normalized
        
        try:

            rms[i] = hdr[root+'RMS']

        except KeyError as err:

            rms[i] = np.nan

        # Create order mask            

        x = np.arange(xranges[i,0],xranges[i,1]+1,1,dtype=int)
        bot = np.polynomial.polynomial.polyval(x,edgecoeffs[i,0,:])
        top = np.polynomial.polynomial.polyval(x,edgecoeffs[i,1,:])

        for j in range(len(x)):

            ordermask[np.floor(bot[j]).astype('int'):\
                      np.ceil(top[j]).astype('int'),x[j]]=orders[i]
                    
    flatinfo.update({'edgecoeffs':edgecoeffs})
    flatinfo.update({'xranges':xranges})
    flatinfo.update({'rms':rms})
    flatinfo = dictaddentry(flatinfo,'bitmask','after','ordermask',\
                            ordermask)

    return(flatinfo)


def readflatcalfile(file):

    """
    Reads a Spextool flatinfo calibration file.


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

       slith_arc : int
           The slit height in arcseconds.

       slith_pix : int
           The nominal slit height in pixels.

       slith_range  : list of int
           (2,) list range of slit heights in pixels

       orders : list of int
           The order numbers.

       rppix : float
           resolving power per pixel

       ps : float
           plate scale in arcseconds per pixel

       step : int
           step size in pixels for tracing

       flatfrac : float
           see findorders.py

       comwidth : int 
           The window in units of pixels used to compute the 
           center-of-mass (COM) (see findorders.py)

       edgedeg : int
           Polynomial fit degree for the edges of the orders
           (see findorders.py)

       norm_nxg : int
           See normspecflat.py and fiterpolate.py

       norm_nyg : int
           See normspecflat.py and fiterpolate.py

       oversamp : float
           See normspecflat.py 

       ybuffer  : int
           See normspecflat.py

       ycororder : int
           See adjustguesspos.py

       xranges : array_like of int, [norders,2] 
           An (norders,2) array giving the column numbers over which to 
           search.  sranges[0,0] gives the starting column number for 
           the first order and sranges[0,1] gives the end column number 
           for the first order.

       edgecoeffs : array_like of float
           (norders,2,ncoeffs) array giving the polynomial 
           coefficients delineating the top and bottom of each order.  
           edgecoeffs[0,0,:] gives the coefficients for the bottom of 
           the order closests to the bottom of `img` and 
           edgecoeffs[0,1,:] gives the coefficients for the top of said 
           order.  

       guesspos : array_like of int
           An (norders,2) array giving the positions to start the 
           search.  guesspos[0,0] gives the column number for the 
           first order and guesspos[0,1] gives the row number for the 
           first order.  Typically the positions are near the center 
           of the image and center of the slit.

    Procedure
    ---------
    Just reads FITS header information, and calculates the guess 
    positions based on the edgecoeffs and the xranges

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
    result.update({'slith_range':val})

    val = hdul[0].header['ORDERS'].split(',')
    orders = [int(x) for x in val]
    norders = len(orders)
    result.update({'orders':orders})            

    val = hdul[0].header['RPPIX']
    result.update({'rpppix':val})

    val = hdul[0].header['PLTSCALE']
    result.update({'ps':val})

    #    val = hdul[0].header['FIXED']
    #    result.update({'fixed':val})        

    #    if not val:

    val = hdul[0].header['STEP']
    result.update({'step':val})
    
    val = hdul[0].header['FLATFRAC']
    result.update({'flatfrac':val})
    
    val = hdul[0].header['COMWIN']
    result.update({'comwidth':val})                                

    deg = int(hdul[0].header['EDGEDEG'])
    result.update({'edgedeg':deg})

    val = hdul[0].header['NORM_NXG']
    result.update({'nxgrid':int(val)})

    val = hdul[0].header['NORM_NYG']
    result.update({'nygrid':int(val)})

    val = hdul[0].header['OVERSAMP']
    result.update({'oversamp':val})

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


def readinstrfile(filename):

    """
    To read a Spextool instrument configuration file

   Input Parameters
    ----------------
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
    > readinstrfile('uspex.dat')

    Modification History
    --------------------
    2022-06-01 - Written by M. Cushing, University of Toledo.
        Based on the Spextool IDL program mc_readinstrfile.pro.

    """

    # Read the file into to string arrays
    
    labels,vals = nploadtxt(filename,comments='#',delimiter='=',\
                            unpack=True,dtype='str')
    
    # Strip any edge white spaces from poor formatting of the file
    
    labels = [value.strip() for value in labels]
    labels = nparray(labels)
    
    vals = [value.strip() for value in vals]
    vals = nparray(vals)

    # Start the search for each keyword.

    output = {}
    
    # INSTRUMENT   
    
    keyword = 'INSTRUMENT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # NROWS
    
    keyword = 'NROWS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # NCOLS
    
    keyword = 'NCOLS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # STDIMAGE
    
    keyword = 'STDIMAGE'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # PLOTWINSIZE
    
    keyword = 'PLOTWINSIZE'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [float(value) for value in tmp]    
    output[keyword] = tmp
    
    # NINT
    
    keyword = 'NINT'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # NSUFFIX
    
    keyword = 'NSUFFIX'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())    

    # BADPIXMASK
    
    keyword = 'BADPIXMASK'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # CALMODULE
    
    keyword = 'CALMODULE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # FILEREADMODE
    
    keyword = 'FILEREADMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # OPREFIX
    
    keyword = 'OPREFIX'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()        

    # SUFFIX
    
    keyword = 'SUFFIX'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # FITSREADPROGRAM
    
    keyword = 'FITSREADPROGRAM'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()            


    # REDUCTIONMODE
    
    keyword = 'REDUCTIONMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # COMBMODE
    
    keyword = 'COMBMODE'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # COMBSTAT
    
    keyword = 'COMBSTAT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # COMBTHRESH
    
    keyword = 'COMBTHRESH'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()                        

    # COMBODIR
    
    keyword = 'COMBODIR'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()                        

    # PSNAPS
    
    keyword = 'PSNAPS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # PSNAPS
    
    keyword = 'PSNAPS'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())            
    
    # OPTEXTRACT
    
    keyword = 'OPTEXTRACT'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # AVEPROF
    
    keyword = 'AVEPROF'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # PSPSFRAD
    
    keyword = 'PSPSFRAD'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

    # PSBGSUB
    
    keyword = 'PSBGSUB'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # PSBGSTART
    
    keyword = 'PSBGSTART'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

    # PSBGWIDTH
    
    keyword = 'PSBGWIDTH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

    # PSBGDEG
    
    keyword = 'PSBGDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())            

    # XSBGSUB
    
    keyword = 'XSBGSUB'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # XSBGREG
    
    keyword = 'COMBSTAT'
    z = find(labels,keyword)
    output[keyword] = vals[z].item()

    # XSBGDEG
    
    keyword = 'XSBGDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # TRACEDEG
    
    keyword = 'TRACEDEG'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # TRACESTEP
    
    keyword = 'TRACESTEP'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # TRACESUMAP
    
    keyword = 'TRACESUMAP'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # TRACESIGTHRESH
    
    keyword = 'TRACESIGTHRESH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

    # TRACEWINTHRESH
    
    keyword = 'TRACEWINTHRESH'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())

    # BADPIXELTHRESH
    
    keyword = 'BADPIXELTHRESH'
    z = find(labels,keyword)
    output[keyword] = float(vals[z].item())

    # LINCORMAX
    
    keyword = 'LINCORMAX'
    z = find(labels,keyword)
    output[keyword] = int(vals[z].item())                                

    # AMPCOR
    
    keyword = 'AMPCOR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # LINCOR
    
    keyword = 'LINCOR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # FLATFIELD
    
    keyword = 'FLATFIELD'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # PLOTXCORR
    
    keyword = 'PLOTXCORR'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # RECTMETHOD
    
    keyword = 'RECTMETHOD'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # FIXBADPIXELS
    
    keyword = 'FIXBADPIXELS'
    z = find(labels,keyword)
    tmp = vals[z].item().split(' ')
    tmp = [int(value) for value in tmp]    
    output[keyword] = tmp

    # XSPEXTOOL KEYWORDS
    
    keyword = 'XSPEXTOOL_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XSPEXTOOL_KEYWORDS'] = tmp

    # XCOMBSPEC KEYWORDS
    
    keyword = 'XCOMBSPEC_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XCOMBSPEC_KEYWORDS'] = tmp

    # XTELLCOR KEYWORDS
    
    keyword = 'XTELLCOR_KEYWORD'
    z = find(labels,keyword)
    tmp = vals[z]
    tmp = [value.strip() for value in tmp]
    output['XTELLCOR_KEYWORDS'] = tmp                                
    
    return output
    

def find(labels,keyword):

    z = npwhere(labels == keyword)
    if npsize(z):

        return z
        
    else:

        print('Cannot find keyword ',keyword,'.',sep='')
        exit(1)

    return(output)


def readwavecalfile(file):

    '''
    To read a pySpextool wavecal calibration file.

    Input Parameters
    ----------------
    file : str
        The fullpath to a calibration file.

        Typicallt the file will end in "_wavecalinfo.fits".

    Returns
    -------
    dict
        A dictionary with the following keywords:

        spectra : numpy.ndarray of float
            (4,n_wave) array 
            spectra[0,:] = wave
            spectra[1,:] = ''flux''
            spectra[2,:] = uncertainty
            spectra[3,:] = bit-set mask

        norders : int
            The number of orders

        orders : numpy.ndarray of int
            The order numbers.  orders[0] is the order closest to the bottom 
            of the array.

        wcaltype : str
            '1D' - A single spectrum where wavelength aligns with the columns 
                   of the array.
            '1DXD' - Cross-dispered spectra where wavelength aligns with the 
                   columns of the array.
            '2D' - A single spectrum where wavelength does not aligns with 
                   the columns of the array.
            '2DXD' - Cross-dispered spectra where wavelength does not aligns 
                   with the columns of the array.

        linelist : str
            The name of the file with the line list.

        xranges : array_like of float
            An (norders,2) array giving the column numbers over which to 
            operate.  xranges[0,0] gives the starting column number for the 
            order nearest the bottom of the image and xranges[0,1] gives 
            the end column number for said order.  

        apradius : float
            The extraction aperture radius (in arcseconds).

        xcororder : int
            The order number to be used for the cross-correlation.

        dispdeg : int
            The polynomial order for the pixel-to-wavelength coefficients.

        if `wcaltype`== '1DXD'

        ordrdeg : int
            The polynomial order for the 2D pixel-to-wavelength coefficients.         
        p2wcoeffs : numpy.ndarray of float
            The polynomial coefficients to convert a pixel number to a 
            wavelength.
      
        homeorder : int
            The order the other orders are scaled to for the fit.

        wavefmt : str
            The format string for wavelengths.

        spatfmt : str
            The format string for spatial angles.            

    Notes
    -----
    

    Examples
    --------
    later

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_readwavecalinfo.pro

    '''

    hdul = fits.open(file)
    spectra = hdul[0].data

    # Clean the header and grab important keywords
    
    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems

    result = {'spectra':spectra}
    
    #    naps = hdul[0].header['NAPS']
    #    result.update({'naps':naps})

    norders = hdul[0].header['NORDERS']
    result.update({'norders':norders})

    val = hdul[0].header['ORDERS'].split(',')
    orders = np.array([int(x) for x in val])
    result.update({'orders':orders})          

    wcaltype = hdul[0].header['WCALTYPE']
    result.update({'wcaltype':wcaltype})

    linelist = hdul[0].header['LINELIST']
    result.update({'wcaltype':linelist})        

    xranges = np.empty((norders,2),dtype=int)
    for i in range(norders):

        val = hdul[0].header['OR'+str(orders[i]).zfill(3)+'_XR'].split(',')
        val = [int(x) for x in val]
        xranges[i,:] = val

    result.update({'xranges':xranges})

    val = hdul[0].header['EXTAP']
    result.update({'apradius':val})

    # Get the cross correlation spectrum
    
    xcororder = hdul[0].header['XCORORDR']
    result.update({'xcororder':xcororder})                

    z = orders == xcororder
    xcorspec = np.squeeze(spectra[z,:,:])
    xcorspec[0,:] = np.arange(xranges[z,0],xranges[z,1]+1)

    result.update({'xcorspec':xcorspec})                    

    if wcaltype == '1DXD': 

        dispdeg = hdul[0].header['DISPDEG']
        result.update({'dispdeg':dispdeg})
                
        ordrdeg = hdul[0].header['ORDRDEG']
        result.update({'ordrdeg':ordrdeg})        

        # Now get the pixel to wavelength coefficients

        ncoeffs = (dispdeg+1)*(ordrdeg+1)
        p2wcoeffs = np.empty((ncoeffs))

        for i in range(ncoeffs):

            key = 'P2W_C'+str(i).zfill(2)
            p2wcoeffs[i] = hdul[0].header[key]

        result.update({'p2wcoeffs':pswcoeffs})
            
        val = hdul[0].header['HOMEORDR']
        result.update({'homeorder':val})

    val = hdul[0].header['WAVEFMT']
    result.update({'wavefmt':val})

    val = hdul[0].header['SPATFMT']
    result.update({'spatfmt':val})            

    hdul.close()


    return result


def writeflat(flat,var,flag,hdrlist,rotate,orders,edgecoeffs,xranges,\
              ps,slith_pix,slith_arc,slitw_pix,slitw_arc,modename,rms,rp,\
              version,history,oname,linmax=None,overwrite=True):

    '''
    To write a Spextool flat FITS file to disk.


    Input Parameters
    ----------------
    flat : numpy.ndarray of float
        The flat field image.


    var : numpy.ndarray of float
        The variance image.


    flag : numpy.ndarray of int
        A biset image.  


    rotate : {0, 1, 2, 3, 4, 5, 6, 7}, optional 
        Direction to rotate a raw image so that the dispersion direction 
        roughly aligns with the columns of the array and wavelength 
        increases to the right, and the spatial axis roughly aligns with 
        the rows of the array.  


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


    orders : list of int
        The order numbers.  orders[0] is the order number of the 
        order nearest the bottom of the image after rotation. 


    edgecoeffs : array_like of float
        (norders,ncoeffs+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  


    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives the 
        end column number for said order.


    ps : float
        The plate scale (arcseconds per pixel).


    slith_pix : float
        The slit height (pixels).


    slith_arc : float
        The nominal slit height (arcseconds).


    slitw_arc : float
        The slit width (arcseconds).


    slitw_pix : float
        The slit width (pixels).


    modename : str
        The name of the instrument mode.


    rms : list of float
        (norders,) list of RMS values for each order.


    rp: float
        The nominal resolving power.


    version : str
        The version of the pySpextool software used to create the flat.


    oname : str
        The filename of the flat field image to written to disk.  


    linmax : int, optional
        The linearity maximum used to identified pixels beyond the 
        linearity limit.


    history : list of str, optional
        The history string spliced to fit in a FITS file.


    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.
        

    Returns
    -------
    None
        Writes a FITS file to PROCPATH.


    Notes
    -----
    None


    Examples
    --------
    Later


    Modification History
    --------------------
    2022-06-29 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program writeflat.pro.

    '''

    # Get basic things

    norders = len(orders)
    edgedeg = edgecoeffs.shape[-1]
    
    # Create the primary HDU
    
    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the hdrlist keywords and values
    
    keys = hdrlist.keys()
    for key in keys:

        if key == 'COMMENT':

            comments = hdrlist['COMMENT']
            
        elif key != 'HISTORY':  # probably not necessary, just in case
            
            hdr[key] = tuple(hdrlist[key])
    
        # Now add new ones

        hdr['MODE'] = (modename,' Instrument Mode')
        hdr['NORDERS'] = (norders, ' Number of orders identified')
        hdr['ORDERS']=(','.join(str(o) for o in orders),'Orders identified')
        hdr['PLTSCALE'] = (ps, 'Plate scale (arcseconds per pixel)')
        hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit length (pixels)')
        hdr['SLTH_ARC'] = (slith_arc, ' Slit length (arcseconds)')
        hdr['SLTW_PIX'] = (slitw_pix, ' Slit width (pixels)')
        hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
        hdr['RP'] = (rp, ' Nominal resolving power')
        hdr['ROTATION'] = (rotate, ' IDL rotate value')
        hdr['VERSION'] = (version, ' Spextool version')

    # Record linearity maximum if given 
        
    if linmax is not None:

        hdr['LINMAX'] = (linmax,' Linearity maximum')

    # Add the RMS values.  Check to make sure not NaN

    if np.sum(np.isnan(rms)) == 0:
        
        for i in range(norders):

            name = 'OR'+str(orders[i]).zfill(3)+'RMS'
            comment = ' RMS of normalized order '+str(orders[i]).zfill(3)
            hdr[name] = (rms[i],comment)

    # Add the xranges

    for i in range(norders):

            name = 'OR'+str(orders[i]).zfill(3)+'_XR'
            comment = ' Extraction range for order '+str(orders[i]).zfill(3)
            hdr[name] = (','.join(str(x) for x in xranges[i,:]),comment)
        
    # Add the edgecoeffs

    hdr['EDGEDEG']=(edgedeg,' Degree of the polynomial fit to order edges')
    for i in range(norders):

        for j in range(2):

            for k in range(edgedeg):

                if j == 0:

                    name = 'OR'+str(orders[i]).zfill(3)+'_B'+str(k+1)
                    comment = ' a'+str(k)+\
                      ' edge coefficient for bottom of order '+\
                      str(orders[i]).zfill(3)

                    hdr[name] = (edgecoeffs[i,j,k], comment)

                if j == 1:

                    name = 'OR'+str(orders[i]).zfill(3)+'_T'+str(k+1)
                    comment = ' a'+str(k)+\
                      ' edge coefficient for top of order '+\
                      str(orders[i]).zfill(3)

                    hdr[name] = (edgecoeffs[i,j,k], comment)

    # Now add the comments 

    if 'comments' in locals():
                    
        for com in comments:

            hdr['COMMENT'] = com

    # and then the history
            
    label = '         ============ Spextool History ============'
    hdr['HISTORY'] = label
    for hist in history:
            
        hdr['HISTORY'] = hist

    # Write the results

    img_hdu = fits.ImageHDU(flat)
    var_hdu = fits.ImageHDU(var)
    flg_hdu = fits.ImageHDU(flag)
            
    hdu = fits.HDUList([phdu,img_hdu,var_hdu,flg_hdu])
    hdu.writeto(oname,overwrite=overwrite)
            
