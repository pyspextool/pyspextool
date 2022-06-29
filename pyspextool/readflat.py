from dictaddentry import dictaddentry
from astropy.io import fits
import numpy as np

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
            list of int -> The orders numbers.  orders[0] is the order 
            number of the order nearest the bottom of the image after 
            rotation. 

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
    flatinfo.update({'orders':orders})

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


    
