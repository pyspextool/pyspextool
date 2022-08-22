"""Functions for spectroscopy utlitiy management."""


import os
import numpy as np
from numpy import zeros as npzeros
from numpy import sum as npsum
from numpy import int8 as npint8
import matplotlib.pyplot as pl
np.set_printoptions(threshold=np.inf)
import scipy
from scipy import ndimage, signal

from pyspextool.utils.ioutils import bitset
from pyspextool.utils.math import *
from pyspextool.utils.text import forprint
from pyspextool.fit.polyfit import robustpolyfit1d,polyfit1d
from pyspextool.utils.fitsutils import getimgrange
from pyspextool.plotting.plot import getyrange
from pyspextool.fit.fitpeak import fitpeak1d
from pyspextool.fit.robust import robustsavgol


def combflagstack(stack,nbits=8):

    '''
    To combine bit-set flag arrays together.

    Input Parameters
    ----------------
    stack : numpy.ndarray 
        The stack of bit-set flag arrays to combine.  The stack
        can either be a stack of spectra [nspec,ndat] or a stack
        of images [nimg,nx,ny].   

    nbits : int, optional
        The number of bits that can potentially be set.  This
        routine assumes the bits are set sequentially, starting
        with the zeroth bit.  So if nbits is 2, then it will
        check the 0th and 1st bit.  The default is to check all
        eight bits

    Output Parameters
    ------------------
    numpy.ndarray
        A bit-set flag array that reflects the bit-set flags from all of
        the spectra or images.

    Procedure
    ---------
    Just some basic math.

    Example
    -------
    Consider a two spectra masks
    
    > spec1 = np.array([0,4,4,2])
    > spec2 = np.array([0,0,3,1])
    > stack = np.stack((spec1,spec2))
    > combflagstack(stack)

    [0 4 7 3]

    Consider two image masks

    > img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
    > img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
    > stack = np.stack((img1,img2))     
    > combflagstack(stack)

    [[1 2 0]
     [3 0 4]
     [0 0 0]]

    Modification History
    --------------------
    2022-03-09 - Written by M. Cushing, University of Toledo.  
    Based on the mc_combflagstack.pro IDL program.

    '''

    # Determine whether this is a spectral stack or image stack
    
    ndim = stack.ndim
    shape = stack.shape

    # Set up the output array
    # Replace with match case statement when you upgrade to 3.10.

    if ndim == 2: comb = npzeros(shape[1],dtype=npint8)
    if ndim == 3: comb = npzeros(shape[1:2],dtype=npint8)        
        
    # Now just loop over each bit requested.

    for i in range(0,nbits):

        #  Identify the pixels with the particular bit set

        set = bitset(stack,i)

        #  Collapse everything down one dimension

        sum = npsum(set,axis=0)

        #  Identify which pixels are set

        mask = sum > 0

        #  Set them to the proper bit value and add to the comb

        comb = comb+mask*2**i 

    return(comb)
    

def extr_xs1dxd(img,var,ordermask,orders,wavecal,spatcal,appos,apradii,\
                    linmax_bitmask=None,badpixel_mask=None,bginfo=None,
                    clupdate=True):

    '''
    To extract and extended source.


    Input Parameters
    ----------------
    img : numpy.ndarray of float
        An (nrows,ncols) image with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    var : The variance image for `img`.

    ordermask : numpy.ndarray of int
        (nrows,cols) order mask image where each pixel is set to its 
        order number.

    wavecal : numpy.ndarray of float
        (nrows,ncols) array where each pixel is set to its wavelength.

    spatcal : numpy.ndarray of float
        (nrows,ncols) array where each pixel is set to its angular position 
        on the sky (arcseconds).

    appos : numpy.ndarray of float or float
        (naps,) array of aperture positions (arcseconds).

    apradius : numpy.ndarray of float or float
        (naps,) array of aperture radii (arcseconds).

    linmax_bitmask : numpy.ndarray of int, optional
        (nrows,ncols) array where a pixel is set to unity if its value is 
        beyond the linearity maximum.  

    badpixel_mask : numpy.ndarray of int, optional
        (nrows,ncols) array where good pixels are set to unity and bad pixels
        are set to zero.  

    bginfo : later

    clupdate : {True, False}, optional
        Set to True for command line updates during execution. 

    Returns
    -------
    list
        A list of integers giving the individual  numbers


    Notes
    -----


    Examples
    --------

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mx_extxdspec.pro.

    '''    

    ##### Do basic things
    
    nrows,ncols = img.shape
    
    norders = len(orders)

    # Convert to numpy arrays to deal with float/int versus list/array problem

    appos = np.array(appos,dtype='float',ndmin=1)
    apradii = np.array(apradii,dtype='float',ndmin=1)    

    naps = len(appos)

    # Deal with bginfo

    if bginfo is not None:

        bgregions = bginfo['regions']
        bgdeg = bginfo['bgdeg']
        
    else:

        bgregions = None
        bgdeg = None
    
    # Create pixel coordinate arrays
    
    xx,yy = mkimgidxs(nrows,ncols)
    
    #  Create linearity maximum mask

    if linmax_bitmask is None:

        linmax_bitmask = np.zeros((nrows,ncols)).astype(int)

    #  Create linearity maximum mask

    if badpixel_mask is None:

        badpixel_mask = np.ones((nrows,ncols)).astype(int)        
        
    ##### Start the order loop

    output_dict = {}
    for i in range(norders):

        if clupdate is not None and i == 0:

            print('Extracting spectra...')
        
        zordr = np.where(ordermask == orders[i])
        xmin = np.min(xx[zordr])
        xmax = np.max(xx[zordr])        
        nwaves = xmax-xmin+1

        # Create the output arrays

        owave = np.full((nwaves),np.nan)
        oflux = np.full((naps,nwaves),np.nan)
        ounc = np.full((naps,nwaves),np.nan)
        omask = np.zeros((naps,nwaves),dtype=int)                

        ##### Start the wavelength loop

        for j in range(nwaves):

            # get the slit mask
            
            colmask = ordermask[:,xmin+j]
            zslit = colmask == orders[i]

            # Store the wavelength

            owave[j] = wavecal[zslit,xmin+j][0]

            # Carve out the slit values            

            slit_pix = yy[zslit,xmin+j]
            slit_arc = spatcal[zslit,xmin+j]
            slit_img = img[zslit,xmin+j]
            slit_var = var[zslit,xmin+j]
            slit_bpm = badpixel_mask[zslit,xmin+j]
            slit_lmm = linmax_bitmask[zslit,xmin+j]            
            
            # Generate the aperture mask

            slitmask = mkslitmask(slit_arc,appos,apradii,\
                                xsbginfo=bgregions)

            # Do the background subtraction

            if bgregions is not None:

                print('do subtraction later.')

            # Do the sum extraction

            for k in range(naps):

                z = (slitmask > float(k)) & (slitmask <= float(k+1))
                oflux[k,j] = np.sum(slit_img[z]*(slitmask[z]-float(k)))
                varval = np.sum(slit_var[z]*(slitmask[z]-float(k+1)))
                ounc[k,j] = np.sqrt(np.abs(varval))

            # Check the bitmask for linearity

            z = slit_lmm == 1
            if sum(z) is True: omask[k,j] = 1

        ##### Store the results

        # Trim the NaNs if need be

        nonan = nantrim(owave,flag=2)
            
        # Generate the key

        for k in range(naps):
    
            key = 'OR'+str(orders[i]).zfill(3)+'_AP'+str(k+1).zfill(2)
            arr=np.stack((owave[nonan],oflux[k,nonan],ounc[k,nonan],\
                    omask[k,nonan]))

            output_dict[key] = arr
            
                
        if clupdate is not None:

            loopprogress(i,0,norders)

    return output_dict


def findorders(img,guesspos,sranges,step,slith_range,deg,bufpix,frac,comwidth,
               qafig=None):
    '''
    Locats orders in a (cross-dispersed) spectral image


    Input Parameters
    ----------------
    img : array_like of float
        An (nrows,ncols) image with (cross-dispersed) spectral orders.  It is 
        assumed that the dispersion direction is roughly aligned with the 
        rows of `img` and the spatial axis is roughly aligned with the 
        columns of `img.  That is, orders go left-right and not up-down.  

    guesspos : array_like of int
        An (norders,2) array giving the positions to start the search.  
        guesspos[0,0] gives the column number for the order closest to the 
        bottom of `img` and guesspos[0,1] gives the row number for said 
        order.  Typically the positions are near the center of the image 
        and center of the slit.

    sranges : array_like of float
        An (norders,2) array giving the column numbers over which to search.
        sranges[0,0] gives the starting column number for the order closest 
        to the bottom of `img` and sranges[0,1] gives the end column number 
        for said order.

    step : int
        The step size in the dispersion or column direction

    slith_range: array_like of float
        (2,) array giving the minimum and maxmimum possible slit heights in 
        units of pixels.  

    deg : int
        Polynomial fit degree for the edges of the orders

    bufpix : int
        Number of pixels to buffer from the top and bottom of `img` during 
        the search

    frac : float
        The fraction of the of the flux of the center of the slit used to 
        identify the location of the edge of the order

    comwidth : int
        The window in units of pixels used to compute the center-of-mass (COM) 
        
    qaflag : {None, str}
        The fullpath of a filename for the quality assurance figure

        
    Output Parameters
    -----------------
    edgecoeffs : array_like of float
        (norders,2,ncoeffs) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of `img` and edgecoeffs[0,1,:] gives the coefficients for 
        the top of said order.  

    Procedure
    ---------
        The IDL sobel function is used to enhance the edges of the orders.  
        Within an order, and at a particular column, the positions above 
        and below the guess position at which the flux falls to `frac` the 
        value at the guess position is determined.  Centroids are computed 
        at these positions on the Sobeled image to determine the location 
        of the edge of the order.  As the program steps away from the guess 
        position, the new guess position is estimated using the previously 
        determined edges.  After the entire array has been checked, the 
        results are fit with a robust least-squares polynomial of degree 
        `deg`.
    
    Examples
    --------
    later?

    Modification History
    --------------------
        2022-03-09 - Written by M. Cushing, University of Toledo.  
                     Based on the Spextool IDL program mc_findorders.pro
    '''

    # Check for qaplotting and set up empty lists

    if qafig is not None:

        plcols = []        
        pledges = []
        pledgecoeffs = []
        plogoodbad = []
    
    # Get basic info and do basic things

    nrows,ncols = img.shape
    norders = len(guesspos)

    rownum = np.arange(nrows)

    halfwin = int(comwidth/2.)

    edgecoeffs = np.empty((norders,2,deg+1))    

    # Sobel the image

    scl = np.max(img)
    simg1 = ndimage.sobel(img/scl,axis=0)
    simg2 = ndimage.sobel(img/scl,axis=1)
    simg = np.sqrt(simg1**2 + simg2**2)

    # Start looping over each order
    
    for i in range(0,norders):

        # Determine the start and stop column
            
        start = sranges[i,0]+step-1
        stop = sranges[i,1]-step+1

        # Generate an array of columns where you will find the slit

        fcols = np.arange(start,stop,step) # "f"ind columns
        nfcols = len(fcols)

        # Create some empty arrays to fill
        
        edges = np.full((2,nfcols),np.nan) # edges of slit
        cens = np.full(nfcols,np.nan)      # center point of slit

        # Fill the first few points in the cens array with the guess position
    
        gpidx = np.asscalar(np.where(np.absolute(fcols-guesspos[i,0])==\
                            np.absolute(fcols-guesspos[i,0]).min())[0])

        cens[(gpidx-deg):(gpidx+deg)] = guesspos[i,1]

        #
        # --------------Now move left from the guess position-------------------
        #

        for j in range(gpidx,0,-1):

            imgcol = img[:,fcols[j]].ravel()       # the image column
            sobcol = simg[:,fcols[j]].ravel()      # the sobel image column

            # Fit the centers so you can project away from the guess position
            
            r = polyfit1d(fcols,cens,max(1,deg-2),justfit=True,silent=True)

            # Find the new guess position yguess 

            yguess = np.polynomial.polynomial.polyval(fcols[j],r['coeffs'])

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess,bufpix,(nrows-bufpix-1)))
            iguess = imgcol[yguess]

            bot,top = find_top_bot(rownum,imgcol,sobcol,yguess,iguess,frac,\
                                 halfwin,debug=True)

            # Confirm that both COMs were computed

            if np.isnan(bot+top) == False:

                # Now check to make sure calculated slit height falls within the
                # parameter slith_pix
                
                if (((top-bot) > slith_range[0]) and
                    ((top-bot) < slith_range[1])):

                    # Store the results, update the cens array

                    edges[:,j] = np.array([bot,top])
                    cens[j] = (bot+top)/2

        #
        # --------------Now move right from the guess position-------------------
        #
                        
        for j in range(gpidx,nfcols,1):

            imgcol = img[:,fcols[j]].ravel()       # the image column
            sobcol = simg[:,fcols[j]].ravel()      # the sobel image column

            # Fit the centers so you can project away from the guess position
            
            r = polyfit1d(fcols,cens,max(1,deg-2),silent=True)

            # Find the new guess position yguess

            yguess = np.polynomial.polynomial.polyval(fcols[j],r['coeffs'])

            # Clip it to avoid the edges

            yguess = int(np.clip(yguess,bufpix,(nrows-bufpix-1)))
            iguess = imgcol[yguess]

            bot,top = find_top_bot(rownum,imgcol,sobcol,yguess,iguess,frac,\
                                 halfwin)                        
            # Confirm that both COMs were computed

            if np.isnan(bot+top) == False:

                # Now check to make sure calculated slit height falls within the
                # parameter slith_pix
                
                if (((top-bot) > slith_range[0]) and
                   ((top-bot) < slith_range[1])):

                    # Store the results, update the cens array

                    edges[:,j] = np.array([bot,top])
                    cens[j] = (bot+top)/2

        # Now fit the results

        tmp = np.empty([2,deg+1])            
        for j in range(0,2):

            fit = robustpolyfit1d(fcols,edges[j,:],deg,4,0.1,\
                                  justfit=True,silent=True)
            tmp[j,:] = fit['coeffs']

            # Store the results for possible plotting         

            if qafig is not None:

                plcols.append(fcols)                
                pledges.append(edges[j,:])
                pledgecoeffs.append(fit['coeffs'])
                plogoodbad.append(fit['ogoodbad'])
                    
        edgecoeffs[i,:,:] = tmp


    if qafig is not None:

        minmax = getimgrange(img,'zscale')        
        fig = pl.figure(figsize=(7,7))
        pl.imshow(img,vmin=minmax[0],vmax=minmax[1],cmap='gray',origin='lower')
        pl.xlabel('Columns (pixels)')
        pl.ylabel('Rows (pixels)')

        # Overplot everything
        
        for i in range(len(pledges)):

            pl.plot(plcols[i],pledges[i],'go',markersize=1)

            z = (np.where(plogoodbad[i] == 0))[0]
            if z.size != 0:

                pl.plot(plcols[i][z],pledges[i][z],'bo',markersize=1)
                
            pl.plot(plcols[i],\
                np.polynomial.polynomial.polyval(plcols[i],pledgecoeffs[i]),\
                'r-',linewidth=0.3)
            
        pl.savefig(qafig)
                
    return edgecoeffs


def find_top_bot(rownum,imgcol,sobcol,yguess,imgguess,frac,halfwin,\
                 debug=False):

    # Moving up from the guess position, find the index of the first pixel in
    # the column that falls below the flux threshold, i.e. the edge of the slit

    ztop = np.where((imgcol < frac*imgguess) & (rownum > yguess))
    if (len(ztop[0]) > 0):
                
        # Compute the index of the bottom and top of the center-of-mass window
        # being careful not to fall off of the array

        bidx = max(0,ztop[0][0]-halfwin)
        tidx = min(ztop[0][0]+halfwin,max(rownum))

        # Yank out the pixels to compute the center of mass

        yvals = rownum[bidx:tidx]
        zvals = sobcol[bidx:tidx]

        # Compute the center of mass

        com_top = np.sum(yvals*zvals)/np.sum(zvals)
        
    else: com_top = np.nan

    # Moving down from the guess position, find the index of the first pixel in
    # the column that falls below the flux threshold, i.e. the edge of the slit

    zbot = np.where((imgcol < frac*imgguess) & (rownum < yguess))
    if (len(zbot[0]) > 0):
                
        # Compute the index of the bottom and top of the center-of-mass window
        # being careful not to fall off of the array

        bidx = max(0,zbot[0][-1]-halfwin)
        tidx = min(zbot[0][-1]+halfwin,max(rownum))
                
        # Yank out the pixels to compute the center of mass

        yvals = rownum[bidx:tidx]
        zvals = sobcol[bidx:tidx]
        
        # Compute the center of mass

        com_bot = np.sum(yvals*zvals)/np.sum(zvals)

    else: com_bot = np.nan

#    if debug is not False:
#
#        plotfig = pl.figure(2)                
#        pl.plot(rownum,sobcol)
#        pl.xlim([com_bot-10,com_top+10])
#        pl.axvline(com_top,color='r')
#        pl.axvline(com_bot,color='r')
#        pl.axvline(yguess,color='g')
#        pl.title(fcols[j])
#        val = input("-->")
#        plotfig.clear()


    return(com_bot,com_top)


def getspecpixshift(xanchor,yanchor,xsource,ysource,savitzky_golay=True,\
                    qafileinfo=None):

    '''
    To determine the pixel shift between two spectra.

    Input Parameters
    ----------------
    xanchor : array-like of int
        (nanchor,) array of pixel values.

    yanchor : array-like
        (nanchor,) array of intensity values.

    xsource : array-like of int
        (nsource,) array of pixel values.

    ysource : array-like
        (nsource,) array of intensity values.

    savitzky_golay : {False, True}, optional
        Set to smooth the `yanchor` and `ysource` with a savitzky-golay 
        filter before cross correlation to protect against bad pixels.

    qafileinfo : dict, optional
        `"figsize"` : tuple
            (2,) tuple of the figure size (inches).

        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.  

    Returns
    -------
    int
        The shift of `ysource` relative to `yanchor` in pixels.  
        if `ysource` is shifted to the right of `yanchor` the returned 
        value is positive.  

    Notes
    -----
    Uses scipy.signal.correlate for the cross correlation.  

    The savitzky-golay parameters are window_length=5 and polydeg=2.

    Examples
    --------
    later

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program xmc_corspec.pro.

    '''

    # Convert to numpy arrays and do basic things

    xanchor = np.array(xanchor,dtype=int)
    yanchor = np.array(yanchor)   
    xsource = np.array(xsource,dtype=int)
    ysource = np.array(ysource)

    ndat = len(xanchor)    
    
    # Find the intersection

    junk, zanchor, zsource = np.intersect1d(xanchor,xsource,\
                                            return_indices=True)

    # clip the results

    xanchor = xanchor[zanchor]
    yanchor = yanchor[zanchor]

    xsource = xsource[zsource]
    ysource = ysource[zsource]

    # Savitzky-Golay the results to protect against bad pixels

    if savitzky_golay is not False:
    
        sgyanchor = robustsavgol(xanchor,yanchor,5)['fit']
        sgysource = robustsavgol(xsource,ysource,5)['fit']

    else:
        
        sgyanchor = yanchor
        sgysource = ysource
        
    # Run the cross correlation
    
    xcor = scipy.signal.correlate(sgysource,sgyanchor,\
                                mode='same',method='fft')
    xcor = xcor/np.nanmax(xcor)

    lag = scipy.signal.correlation_lags(ndat,ndat,mode='same')

    ##### Now lets find the fit window

    # Determine the pixel when the derivative becomes positive again

    maxidx = np.argmax(xcor)

    dif = -1
    npix = 0
    while dif < 0 and maxidx+npix+1 < ndat-1:
     
        dif = xcor[maxidx+npix+1]-xcor[maxidx+npix]
        npix += 1

    halfwin = np.around(npix*1.5).astype('int')

    # Clip out the fit zone
        
    fitlag = lag[maxidx-halfwin:maxidx+halfwin]
    fitxcor = xcor[maxidx-halfwin:maxidx+halfwin]
    
    # Do the fit

    r = fitpeak1d(fitlag,fitxcor,nparms=4,positive=True)
    offset = r['parms'][1]
    
    if qafileinfo is not None:

        # Normalize the spectra for plotting purposes

        np.divide(ysource,np.median(ysource),out=ysource)
        np.divide(yanchor,np.median(yanchor),out=yanchor)        

        if 'figsize' in qafileinfo.keys():

            figsize = qafileinfo['figsize']

        else: figsize=(8.5,11)

        # A 3 panel figure

        pl.rcParams['font.size'] = '8'
        pl.rcParams['lines.linewidth'] = 0.75
            
        fig, (axes1, axes2, axes3) = pl.subplots(3,figsize=figsize,\
                                        constrained_layout=True)
        # Plot the two spectra

        yrange = getyrange([yanchor,ysource],frac=0.1)
        
        axes1.margins(x=0)
        axes1.step(xanchor,yanchor,'#1f77b4')
        axes1.step(xsource,ysource,'r')
        axes1.set_ylim(ymin=yrange[0],ymax=yrange[1])
        axes1.set(xlabel='Column Number',ylabel='Relative Intensity')

        axes1.text(0.95,0.8,'anchor',color='#1f77b4',ha='right',\
                       transform=axes1.transAxes)

        axes1.text(0.95,0.7,'source',color='r',ha='right',\
                       transform=axes1.transAxes)        
        
        # Plot the entire cross correlation results

        yrange = getyrange(xcor,frac=0.1)

        axes2.margins(x=0)        
        axes2.tick_params(axis='x')
        axes2.tick_params(axis='y')
        axes2.set_title('Cross Correlation')
        axes2.step(lag,xcor)
        axes2.set_ylim(ymin=yrange[0],ymax=yrange[1])
        axes2.set(xlabel='Lag (pixels)',ylabel='Relative Intensity')        

        # Plot a zoom in of the cross correlation and the fit        

        yrange = getyrange([fitxcor,r['fit']],frac=0.1)

        axes3.margins(x=0)        
        
        axes3.step(fitlag,fitxcor)
        axes3.step(fitlag,r['fit'],'r')
        axes3.set_title('Fit of Cross Correlation')
        axes3.set_ylim(ymin=yrange[0],ymax=yrange[1])

        axes3.set(xlabel='Offset (pixels)',ylabel='Relative Intensity')
        axes3.axvline(x=offset,linestyle='dotted',color='r')

        axes3.text(0.95,0.8,'offset='+"{:.1f}".format(offset)+' pixels',\
                       ha='right',c='r',transform=axes3.transAxes)
        
        # Save the figure
        
        pl.savefig(os.path.join(qafileinfo['filepath'],\
                                qafileinfo['filename'])+\
                                qafileinfo['extension'])
        
    return offset


def mkapmask(slit_arc,appos,apradius,psbginfo=None,xsbginfo=None):

    '''
    To create an aperture mask for use in extraction.  

    Input Parameters
    ----------------
    slit_arc : array_like
        An array of angular positions along the slit (typically in arcseconds).

    appos : array_like or float or int
        The aperture positions (same units as `slit_arc`).

    apradius : array_like or float or int
        The aperture radius (same units as `slit_arc`).

    psbginfo : list, optional
        (2,) list giving the angular distance from the aperture to begin the 
        background region and the angular width of the background region. 

    xsbginfo : list of list, optional
        (n_bg,) list where each element is a 2-element list of angular 
        positions that indicate a background region.

    Returns
    -------
    numpy.ndarray
         An array of the aperture mask.  The pixel values denote their status:

         val = 0 - nothing
         val = -1 - a background pixel
         0 < val <= 1 - a pixel belonging to aperture 1
         1 < val <= 2 - a pixel belonging to aperture 2
         etc.

    Notes
    -----
    The function generates an array where each pixel value denotes its status.
    To construct the mask, the background regions are identified without 
    regard to potential overlap.  If `psbginfo` is passed, then overlap is 
    dealt with by clearing out from `appos`-`bgstart` to `appos`+`bgstart`.  
    Finally, the apertures are labeled by their aperture number, e.g. 1, 2, etc.
    Pixels at the edge of an aperture are adjust to reflect their fractional 
    pixel value, e.g. 1 -> 0.5 or 2 -> 1.6.

    Examples
    --------
    > slit_arc = np.arange(100)
    > mkapmask(slit_arc,[50,60],[4,4.1],psbginfo=[11,20])
      [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.
        0.   0.   0.   0.   0.5  1.   1.   1.   1.   1.   1.   1.   0.5  0.
        1.6  2.   2.   2.   2.   2.   2.   2.   1.6  0.   0.   0.   0.   0.
        0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.   0.   0.   0.   0.
        0.   0. ]

    > slit_arc = np.arange(100)
    > mkapmask(slit_arc,[50,60],[4,4.1],xsbginfo=[[7,15],[80,90]])
      [ 0.   0.   0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.5  1.   1.   1.   1.   1.   1.   1.   0.5  0.
        1.6  2.   2.   2.   2.   2.   2.   2.   1.6  0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.   0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0. ]

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_mkapmask.pro.

    '''

    # Convert to numpy arrays

    appos = np.array(appos,dtype='float',ndmin=1)
    apradius = np.array(apradius,dtype='float',ndmin=1)    
    
    # Now do basic things
        
    npix = len(slit_arc)
    naps = len(appos)
    
    # Now check to see if the apertures overlap

    mask = np.zeros([npix])
    for i in range(naps):

        ap = [appos[i]-apradius[i],appos[i]+apradius[i]]
        idx = findidx(slit_arc,ap)
        mask[int(idx[0]):int(idx[1])] += 1

    test = mask > 1
    if sum(test) >=1:

        exception = 'Apertures overlap.'
        raise ValueError(exception)

    # Create the mask

    mask = np.zeros([npix])

    # Let's deal with the background regions, point source?

    if psbginfo is not None:

        bgstart = psbginfo[0]
        bgwidth = psbginfo[1]

        # Set all pixels that fall in the background region around each aperture
        # to -1 irregardless of whether there is overlap.

        for i in range(naps):

            if bgstart <= apradius[i]:

                error = 'psbginfo[0] (bgstart) must be greater than '+\
                  'the aperture radius.'
                raise ValueError(error)
                    
            leftbg = [appos[i]-bgstart-bgwidth,appos[i]-bgstart]
            rghtbg = [appos[i]+bgstart,appos[i]+bgstart+bgwidth]
            
            leftidx = findidx(slit_arc,leftbg)
            rghtidx = findidx(slit_arc,rghtbg)        
            
            mask[int(leftidx[0]):int(leftidx[1])] = -1
            mask[int(rghtidx[0]):int(rghtidx[1])] = -1

        # Now clear out from -bgstart to +bgstart around each aperture position to
        # deal with any overlap

        for i in range(naps):

            clr = [appos[i]-bgstart,appos[i]+bgstart]
            clridx = findidx(slit_arc,clr)
            mask[int(clridx[0]):int(clridx[1])] = 0
                    
    # extended source?

    if xsbginfo is not None:

        # Let's check against a single background region

        if any(isinstance(el,list) for el in xsbginfo) is False:

            xsbginfo = [xsbginfo]
            
        for i in range(len(xsbginfo)):

            idx = findidx(slit_arc,xsbginfo[i])
            mask[int(idx[0]):int(idx[1])] = -1

    # Now fill in the apertures

    for i in range(naps):

        ap = [appos[i]-apradius[i],appos[i]+apradius[i]]
        idx = findidx(slit_arc,ap)
        mask[int(idx[0]):int(idx[1])]+=i+1

        # Fix end points to reflect fractional pixels

        if idx[0]-np.floor(idx[0]) >= 0.5:

            
            mask[int(idx[0])] = 0
            mask[int(idx[0])+1] = 0.5 + round_tntafz(idx[0]) - idx[0] + i
            
        else:
            
            mask[int(idx[0])] = 0.5 - (idx[0] - np.floor(idx[0])) + i

        if idx[1]-np.floor(idx[1]) >= 0.5:

            mask[int(idx[0])+1] = 0.5 - (round_tntafz(idx[1]) - idx[1])  + i
            
        else:
            
            mask[int(idx[1])] = 0.5 + (idx[1] - round_tntafz(idx[1])) + i            
    

    return mask


def simwavecal1dxd(ncols,nrows,edgecoeffs,xranges,slith_arc):

    '''
    To simulate Spextool wavecal and spatcal arrays.

    Will generate wavecal and spatcal files in the 1DXD case with the
    wavelengths replaced with the column numbers.


    Input Parameters
    ----------------
    ncols : int
        The number of columns of the image.

    nrows : int
        The number of rows of the image.

    edgecoeffs : array_like of float
        (norders,`edgedeg`+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  

    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.

    slith_arc : float
        The nominal slit height (arcseconds).

    Returns
    -------
    wavecal, spatcal : numpy.ndarray, numpy.ndarray
        - wavecal (nrows,ncols) array where each pixel is set to its 
          wavelength which in this case is the column number.
        - spatcal (nrows,ncols) array where each pixel is set to its 
          angular position on the sky (in arcseconds).

    Notes
    -----
    None

    Examples
    --------
    later

    Modification History
    --------------------
    2022-06-28 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_simwavecal2d.pro.

    '''
    
    # Get basic info
    
    ndimen = edgecoeffs.ndim
    if ndimen == 2: norders = 1
    if ndimen == 3: norders = edgecoeffs.shape[0]
    
    # Make some NaN arrays
    
    wavecal = np.full([nrows,ncols],np.nan)
    spatcal = np.full_like(wavecal,np.nan)

    # start the loop over order and column number

    y = np.arange(nrows)

    for i in range(norders):

        start = xranges[i,0]
        stop = xranges[i,1]

        x = np.arange(stop-start+1)+start

        # Get the top and bottom positions of the slit
        
        botedge = np.polynomial.polynomial.polyval(x,edgecoeffs[i,0,:])
        topedge = np.polynomial.polynomial.polyval(x,edgecoeffs[i,1,:])

        # Creat the pixel to arcsecond transformation
        
        pixtoarc = np.empty([2,stop-start+1])
        pixtoarc[1,:] = slith_arc/(topedge-botedge)
        pixtoarc[0,:] = -1*pixtoarc[1,:]*botedge

        # Fill things in
        
        for j in range(stop-start+1):

            wavecal[np.floor(botedge[j]).astype('int'):\
                    np.ceil(topedge[j]).astype('int'),x[j]] = x[j]

            # Create ysub to make things readable...
            
            ysub = y[np.floor(botedge[j]).astype('int'):\
                    np.ceil(topedge[j]).astype('int')]
            
            spatcal[np.floor(botedge[j]).astype('int'):\
                    np.ceil(topedge[j]).astype('int'),x[j]] = \
                    np.polynomial.polynomial.polyval(ysub,pixtoarc[:,j])

            
    return wavecal,spatcal