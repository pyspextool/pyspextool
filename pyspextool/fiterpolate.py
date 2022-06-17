import numpy as np
from imgquadfit import imgquadfit
from mkimgidxs import mkimgidxs
from bicuval import bicuval

def fiterpolate(img,ncg,nrg):

    '''
    Not sure how to explain!


    Input Parameters
    ----------------
    img : array_like
        The image to be fiterpolate-d.  

    ncg : int
        The number of columnn grid cells
       
    nrg : int
        The number of row grid cells

    Returns
    --------
    numpy.ndarray
        The fitted image

    Notes
    -----
    Based on John Tonry's fiterpolate program.  Going to require
    me going into this again and remembering exactly what it does.
    

    Examples
    --------
    later?

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_fiterpolate.pro.

    '''

# Get basic infor and create basic things
    
    nrows,ncols = img.shape

    fitimg = np.empty((nrows,ncols))

# Determine the number of grid points
    
    nxgrid = ncg+1
    nygrid = nrg+1

# Compute the actually grid points

    gx = np.rint(np.arange(nxgrid)*(ncols-1)/(nxgrid-1)).astype(int)
    gy = np.rint(np.arange(nygrid)*(nrows-1)/(nygrid-1)).astype(int)

# Set up the grid information for ease of use later
    
    gridinfo = []

    for i in range(ncg):

        for j in range(nrg):

            dict = {}        
            dict.update({'xrng':[gx[i],gx[i+1]]})
            dict.update({'yrng':[gy[j],gy[j+1]]})
            dict.update({'xsize':gx[i+1]-gx[i]+1})
            dict.update({'ysize':gy[j+1]-gy[j]+1})
            gridinfo.append(dict)

# Now determine the values of z, dy1, dy2, dy12 at the grid points.
# They will be stored in the val array

    idx = 0
    vals = np.empty((nxgrid*nygrid,4))
    
    for i in range(nxgrid):

        x1 = np.rint((gx[np.max([i-1,0])]+gx[i])/2).astype(int)
        x2 = np.rint((gx[i]+gx[np.min([nxgrid-1,i+1])])/2).astype(int)

        for j in range(nygrid):

            y1 = np.rint((gy[np.max([j-1,0])]+gy[j])/2).astype(int)
            y2 = np.rint((gy[j]+gy[np.min([nygrid-1,j+1])])/2).\
                         astype(int)

            c = imgquadfit(img[y1:y2+1,x1:x2+1])

            x = gx[i]-x1 # offset to account for zero-based fit
            y = gy[j]-y1


            vals[idx,0]= c[0]+c[1]*x+c[2]*y+c[3]*x**2+c[4]*y**2+c[5]*x*y
            vals[idx,1]= c[1]+2.*c[3]*x+c[5]*y
            vals[idx,2]= c[2]+2.*c[4]*y+c[5]*x
            vals[idx,3]= c[5]

            idx +=1 

# Now perform the bicubic interpolation and reconstruct the image


            
    k = 0
    for i in range(ncg):

        for j in range(nrg):

            idx = [i*nygrid+j,(i+1)*nygrid+j,(i+1)*nygrid+j+1,\
                   (1+i*nygrid)+j]

# Create coordinate images

            ximg, yimg = mkimgidxs(gridinfo[k]['ysize'],\
                                   gridinfo[k]['xsize'])

            nimg = bicuval(vals[idx,0],vals[idx,1],vals[idx,2],\
                           vals[idx,3],0.0,gridinfo[k]['xsize'],0.0,\
                           gridinfo[k]['ysize'],ximg,yimg)

            fitimg[gridinfo[k]['yrng'][0]:gridinfo[k]['yrng'][1]+1,\
                   gridinfo[k]['xrng'][0]:gridinfo[k]['xrng'][1]+1] = \
                   nimg

            k +=1 

    return fitimg
