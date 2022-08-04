"""Functions for calibration."""


import numpy as np
from pyspextool.fit.interpolate import fiterpolate
from scipy import interpolate
from scipy.signal import medfilt2d
from pyspextool.utils.math import loopprogress, moments


def normspecflat(img,edgecoeffs,xranges,slith_arc,nxgrid,nygrid,var=None,\
                 oversamp=1,ybuffer=0,clupdate=False):

    '''
    Normalize spectral flat field image


    Input Parameters
    ----------------
    img : array_like of float
        An (nrows,ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 


    edgecoeffs : array_like of float
        (norders,ncoeffs+1,2) array giving the polynomial coefficients 
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


    nxgrid : int


    nygrid : int


    var : numpy.ndarray 
        

    oversamp : float, default 0, optional
        The factor by which to oversample the slit during rectification. 
        The new slit length in pixels is given by 
        round(oversamp*min-slit-length-in-native-pixels). 


    ybuffer : int, default 1, optional
        The number of native pixels from the top and bottom of the slit to 
        avoid during the operation.  Useful to account for the fact that 
        the drop-off in intensity at the edge of the slit is not a 
        heaviside function but rather occurs over a few pixels.  


    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray

        nimg :  An image of the same size as `img` with the spectral ordered 
        normalized to unity, and all other pixels set to unity.


        nvar :  A variance image of the same size as `img` with the spectral 
            ordered normalized by the same model as `nimg`.  All other pixels
            are set to NaN.


        rms : (norders,) array giving the RMS value of the pixels


    Notes
    -----
    Resampling each order onto a rectangular grid.  Uses fiterpolate 
    to create a surface model and then resamples the model back onto the 
    original pixels.  


    Examples
    --------
    later


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_normspecflat.pro

    '''    

    # Get basic info and do basic things

    nrows, ncols = img.shape

    ndimen = edgecoeffs.ndim
    if ndimen == 2: norders = 1
    if ndimen == 3: norders = edgecoeffs.shape[0]

    nimg = np.full((nrows,ncols),1.0)
    nvar = np.full((nrows,ncols),np.nan) if var is not None else None
    rms  = np.empty(norders)    

    # Start the loop

    for i in range(norders):        

        # Rectify the order

        spatmap, order= rectifyorder1d(img,edgecoeffs[i,:,:],xranges[i,:],\
                                       slith_arc,oversamp=oversamp,\
                                       ybuffer=ybuffer)

        # Fiterpolate the results after median smoothing to minimize bad pixels

        model = fiterpolate(medfilt2d(order,kernel_size=(5,5)),nxgrid,nygrid)

        # Now normalize the raw data using the fiterpolated model

        # Get useful things and set up
        
        startcol = xranges[i,0]
        stopcol  = xranges[i,1]

        order_ncols = stopcol-startcol+1
        order_cols = np.arange(order_ncols,dtype=int)+startcol
        
        botedge = np.polynomial.polynomial.polyval(order_cols,edgecoeffs[i,0,:])
        topedge = np.polynomial.polynomial.polyval(order_cols,edgecoeffs[i,1,:])
        dif = topedge-botedge

        # Loop over each column, interpolate the model onto the data sampling, and
        # divide

        mask = np.full((nrows,ncols),0)
        for j in range(order_ncols):

            # get linear conversion from pixels to arcseconds

            b = slith_arc/dif[j]
            m = -1*b*botedge[j]

            # get range over which the slit falls and create y values

            yrange = [np.ceil(botedge[j]).astype('int')+ybuffer,\
                      np.floor(topedge[j]).astype('int')-ybuffer]
            
            ypix_slit = np.arange(yrange[0],yrange[1]+1)
        
            # Do the linterpolation

            f = interpolate.interp1d(spatmap[:,0],model[:,j])
            tmp = np.polynomial.polynomial.polyval(ypix_slit,[m,b])

            
            slit_model = f(np.polynomial.polynomial.polyval(ypix_slit,[m,b]))

            # divide the data by the interpolated model
            
            nimg[yrange[0]:yrange[1]+1,j+startcol] = \
              img[yrange[0]:yrange[1]+1,j+startcol]/slit_model

            # divide the variance by the interpolated model if need be

            if var is not None:

                nvar[yrange[0]:yrange[1]+1,j+startcol] = \
                var[yrange[0]:yrange[1]+1,j+startcol]/slit_model                
            
            # fill the mask

            mask[yrange[0]:yrange[1]+1,j+startcol] = 1

        # Compute the RMS of the result

        z = np.where(mask == 1)
        m = moments(nimg[z],robust=4.5)
        rms[i] = m['stddev']
        
        if clupdate: loopprogress(i,0,norders,message='Normalizing the flat...')


    return (nimg,nvar,rms)


def rectifyorder1d(img,edgecoeffs,xranges,slith_arc,oversamp=1,ybuffer=0):

    """
    To rectify a spectral order 

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
        (norders,2,ncoeffs) array giving the polynomial coefficients 
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

    if ybuffer > 0:
    
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
        
        f = interpolate.interp1d(slit_arc,img[slitbot_pix:slittop_pix+1,\
                                 i+startcol])
        order[:,i] = f(rslit_arc)

    # Create the spatial map using rslit_arc

    spatmap = np.rot90(np.tile(rslit_arc,(order_ncols,1)),k=3)
        
    return (spatmap,order)
    
