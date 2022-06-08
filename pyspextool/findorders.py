import numpy as np
from scipy import ndimage
import matplotlib.pyplot as pl
from polyfit1d import polyfit1d
from robustpolyfit1d import robustpolyfit1d
from getimgrange import getimgrange

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
        guesspos[0,0] gives the column number for the first order and 
        guesspos[0,1] gives the row number for the first order.  Typically
        the positions are near the center of the image and center of the slit.

    sranges : array_like of float
        An (norders,2) array giving the column numbers over which to search.
        sranges[0,0] gives the starting column number for the first order and 
        sranges[0,1] gives the end column number for the first order.

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
        (norders,ncoeffs+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,:,0]
        gives the coefficients for the bottom of the order closests to the 
        bottom of `img` and edgecoeffs[0,:,1] gives the coefficients for 
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

#
#============================================================================
#
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
