import numpy as np
from scipy import interpolate

def rectifyorder1d(img,edgecoeffs,xranges,slith_arc,oversamp=1,ybuffer=0):

    """
    Rectifies a spectral order 

    The function "straightens" a spectral order onto a uniform rectangular 
    grid


    Input Parameters
    ----------------

    img : array_like of float
        An (nrows,ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    edgecoeffs : array_like of float
        (norders,2,ncoeffs+1) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for 
        the top of said order.  

    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of `img` and xranges[0,1] gives the end 
        column number for said order.

    slith_arc : float
        The slit height in arcseconds

    oversamp : float, optional
        The factor by which to oversample the slit during rectification. 
        The new slit length in pixels is given by 
        round(oversamp*min-slit-length-in-native-pixels).  

    ybuffer : int, optional
        The number of native pixels from the top and bottom of the slit to 
        avoid during the operation.  Useful to account for the fact that 
        the drop-off in intensity at the edge of the slit is not a 
        heaviside function but rather occurs over a few pixels.  

    
    Returns
    -------
    spatmap, rimg : numpy.ndarray of float, numpy.ndarray of float
    
        spatmap - each pixel is set to its angular position on the sky

        rimg - the resampled order

    Notes
    -----
    The function assumnes that the slit image is aligned with the columns of 
    the detector.  The length (height) of the slit in native pixels can vary
    across an order due to anamorphic magnification and so the plate scale 
    must be computed at every column.  


    Examples
    --------
    later?


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_fsextract.pro

    """

# Get basic info and create basic things

    nrows,ncols = img.shape

    startcol = xranges[0]
    stopcol  = xranges[1]

    order_ncols = stopcol-startcol+1
    cols = np.arange(order_ncols,dtype=int)+startcol
    rows = np.arange(nrows,dtype=int)

# Get the bottom and top of the slit 

    botedge = np.polynomial.polynomial.polyval(cols,edgecoeffs[0,:])
    topedge = np.polynomial.polynomial.polyval(cols,edgecoeffs[1,:])

# Get the conversion between pixels and arcseconds

    dif = topedge-botedge
    arctopix = dif/slith_arc
    
# Create a uniform grid of arcsecond values for the new order

    mindif = np.floor(min(dif))

    nrslit = round(oversamp*mindif) # number of pixels in resampled slit
    
    rslit_arc = np.arange(nrslit)*slith_arc/(nrslit-1)

# Set the pixels ybuffer from the ends to the values at ybuffer from the edge

    ybuffer = round(ybuffer*oversamp)
    
    rslit_arc[0:ybuffer-1] = rslit_arc[ybuffer-1]
    rslit_arc[(nrslit-ybuffer):nrslit] = rslit_arc[nrslit-ybuffer]

# Now do the interpolation one column at a time

# NOTE: this is not the fast way, but I don't see an IDL equivalent of
# interpolate in python.  Looking through SOFIA's repo they seem to have
# written a bunch of their own so it may not be possible with standard
# scipy packages.

    order = np.empty((nrslit,order_ncols))
    
    for i in range(0,order_ncols):

# Get the bottom and top of the slit positions at cols[i].   Use floor and
# ceil to ensure the interpolation has no extrapolation.
        
        slitbot_pix = np.floor(botedge[i]).astype(int)
        slittop_pix = np.ceil(topedge[i]).astype(int)

# Yank out the row values and convert to arcseconds
        
        slit_pix = rows[slitbot_pix:(slittop_pix+1)]-botedge[i]
        slit_arc = slit_pix/arctopix[i]

# Do the interpolation and store
        
        f = interpolate.interp1d(slit_arc,img[slitbot_pix:slittop_pix+1,i])
        order[:,i] = f(rslit_arc)

# Create the spatial map using rslit_arc

    spatmap = np.rot90(np.tile(rslit_arc,(order_ncols,1)),k=3)
        
    return (spatmap,order)
    
