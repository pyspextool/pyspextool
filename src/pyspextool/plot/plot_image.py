import numpy as np
import matplotlib.pyplot as pl

from pyspextool.fit.polyfit import poly_1d
from pyspextool.io.check_parameter import check_parameter
from pyspextool.plot.limits import get_image_range

def plot_image(image, edgecoeffs, xranges, orders):

    '''
    To plot a spectral image along with the edges and order numbers

    Parameters
    ----------
    img : numpy.ndarray 
        An (nrows, ncols) float array with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    edgecoeffs : numpy.ndarray 
        (norders,`edgedeg`+1,2) float array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  

    xranges : array_like
        An (norders, 2) float array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.

    orders : list of int
	    (norders,) int array of the order numbers.  By Spextool convention, 
        orders[0] is the order closest to the bottom of the array.

    Returns
    -------
        None

    Notes
    -----
        Just basic plotting stuff.
    
    '''
    
    #
    # Check parameters
    #
    
    check_parameter('plot_image', 'image', image, 'ndarray', 2)

    check_parameter('plot_image', 'edgecoeffs', edgecoeffs, 'ndarray', [2, 3])

    check_parameter('plot_image', 'xranges', xranges, 'ndarray', [1, 2])

    check_parameter('plot_image', 'orders', orders, ['int', 'ndarray'], [0,1])

    #
    # Just plot it up
    #
    
    norders = len(orders)
        
    minmax = get_image_range(image, 'zscale')
    fig = pl.figure(figsize=(7, 7))
    pl.imshow(image, vmin=minmax[0], vmax=minmax[1], cmap='gray',
                  origin='lower')
    pl.xlabel('Columns (pixels)')
    pl.ylabel('Rows (pixels)')
    
    for i in range(norders):

        x = np.arange(xranges[i,0],xranges[i,1]+1)
        bot = poly_1d(x,edgecoeffs[i,0,:])
        top = poly_1d(x,edgecoeffs[i,1,:])
        
        pl.plot(x,bot,color='purple')
        pl.plot(x,top,color='purple')                
        pl.fill_between(x,bot,y2=top,color='purple',alpha=0.15)
        
        delta = xranges[i,1] - xranges[i,0]
        idx = np.fix(delta*0.02).astype(int)
        
        pl.text(x[idx],(top[idx]+bot[idx])/2., str(orders[i]),
                    color='yellow', verticalalignment='center')
        
                
    pl.show()

    
    
