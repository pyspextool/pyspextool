import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator

from pyspextool.fit.polyfit import poly_1d
from pyspextool.io.check import check_parameter
from pyspextool.plot.limits import get_image_range

def plot_image(image:npt.ArrayLike,
               zrange:str | list | float=None,
               mask:npt.ArrayLike=None,
               orders_plotinfo:dict=None,
               trace_plotinfo:dict=None,
               locateorders_plotinfo:dict=None,
               figure_size:tuple=(7,7),
               font_size:int=12,               
               output_fullpath:str=None,
               showblock:bool=False,
               showscale:float | int=1.0,
               plot_number:int=None):
    
    """
    To plot a spectral image along with the edges and order numbers

    Parameters
    ----------
    image : ndarray
        An (nrows, ncols) array with (cross-dispersed) spectral orders.  
        It is assumed that the dispersion direction is roughly aligned 
        with the rows of `img` and the spatial axis is roughly aligned 
        with the columns of `img.  That is, orders go left-right and 
        not up-down. 

    zrange : str, list, float, None


    mask : ndarray, optional
        An (nrows, ncols) array where "bad" pixels are set to unity and 
        "good" pixels are zero. If passed, "bad" pixels are colored red.

    orders_plotinfo : dict, optional
        'edgecoeffs' : ndarray 
            (norders,`edgedeg`+1,2) float array giving the polynomial 
            coefficients delineating the top and bottom of each order.  
            edgecoeffs[0,0,:] gives the coefficients for the bottom of the 
            order closest to the bottom of the image and edgecoeffs[0,1,:] 
            gives the coefficients for the top of said order.  

        'xranges' : ndarray
            An (norders, 2) float array giving the column numbers over which to 
            operate.  xranges[0,0] gives the starting column number for the 
            order nearest the bottom of the image and xranges[0,1] gives 
            the end column number for said order.

        'orders' : list of int
            (norders,) int array of the order numbers.  By Spextool convention, 
            orders[0] is the order closest to the bottom of the array.

    trace_plotinfo : dict, optional
        'x' : list
            An (npoints,) list of x positions for the `y` fit positions.

        'y' : list
            An (npoints,) list of fitted positions.

        'goodbad' : list
            An (npoints,) list of goodbad values.  0=bad, 1=good.

        'fits' : list
            An (norders*napertures,) list of fitted values to the `x`s and `y`s.
            fits[i] is a (2,ndata) array where fits[i][0,:] = x values, and
            fits[i][1,:] = y values.  
       
    locateorders_plotinfo : dict, optional
        'guess_positions': list
            An (norders,) list where each element is a two-element list 
            giving the (x,y) position of the guess position

        'x' : list
            An (2*norders,) list where each element is the x positions of 
            either the top or the bottom of an order.

        'y' : list
            An (2*norders,) list where each element is the x positions of 
            either the top or the bottom of an order.

        'goodbad' : list
            An (2*norders,) list where each element is the x positions of 
            either the top or the bottom of an order.

        'coefficients' : list
            An (2*norders,) list where each element is (ncoeffs,) ndarray of 
            the polynomial coefficients of the top or bottom of an order.

    figure_size : tuple of int or float, default (7,7)
         The figure size in inches.

    font_size : int, default 12
         The font size to use in the plot.
    
    output_fullpath : str, optional
        A string giving the fullpath filename to write the QA to disk.
    
    showblock : {False, True}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.

    showscale : int or float, default 1
        A scale factor by which to scale the size of the plot displayed to the
        screen and the font size.
    
    plot_number : int, default None
        A plot number to pass to matplotlib.  

    Returns
    -------
        None
        Displays an image to the screen or writes the image to disk.
    
    """

    #
    # Check parameters
    #

    check_parameter('plot_image', 'image', image, 'ndarray', 2)

    check_parameter('plot_image', 'mask', mask, ['NoneType', 'ndarray'], 2)

    check_parameter('plot_image', 'locateorders_plotinfo',
                    locateorders_plotinfo, ['NoneType', 'dict'])

    check_parameter('plot_image', 'orders_plotinfo', orders_plotinfo,
                    ['NoneType', 'dict'])

    check_parameter('plot_image', 'trace_plotinfo', trace_plotinfo,
                    ['NoneType', 'dict'])

    check_parameter('plot_image', 'output_fullpath', output_fullpath,
                    ['NoneType', 'str'])

    check_parameter('plot_image', 'plot_number', plot_number,
                    ['NoneType', 'int'])

    check_parameter('plot_image', 'showblock', showblock, 'bool')

    check_parameter('plot_image', 'showscale', showscale, ['int','float'])

    
    #
    # Make the plot
    #
 
    if output_fullpath is not None:

        # Write the plot to disk.

        doplot(None,
               figure_size, 
               font_size,
               image,   
               zrange=zrange,
               mask=mask,
               locateorders_plotinfo=locateorders_plotinfo,
               orders_plotinfo=orders_plotinfo,
               trace_plotinfo=trace_plotinfo)

        pl.savefig(output_fullpath)
        pl.close()

    if output_fullpath is None:
    
        # Display the image to the screen.

        doplot(plot_number,
               (figure_size[0]*showscale,figure_size[1]*showscale),
               font_size*showscale,
               image,
               zrange=zrange,
               mask=mask,
               locateorders_plotinfo=locateorders_plotinfo,
               orders_plotinfo=orders_plotinfo,
               trace_plotinfo=trace_plotinfo)

        pl.show(block=showblock)
        if showblock is False: pl.pause(1)

        return plot_number

        
    
def doplot(plot_number:int,
           figure_size:tuple,
           font_size:float | int,
           image:npt.ArrayLike,
           zrange=None,
           mask:npt.ArrayLike=None,
           locateorders_plotinfo:dict=None,
           orders_plotinfo:dict=None,
           trace_plotinfo:dict=None):

    """
    To plot the image "independent of the device"


    """

    # Set the fonts

    # removed helvetica - problem for windows OS
    font = {'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)
    
    # Set the color map

    cmap = pl.cm.gray

    # Get the plot range

    if isinstance(zrange, list):
        
        minmax = np.array(zrange)
        
    else:

        type = 'zscale' if zrange == None else zrange

        minmax = get_image_range(image, type)
        if minmax[0] > minmax[1]:
            
            minmax = (np.nanmin(image), np.nanmax(image))
    
    # Now check to see if the mask is passed.

    pimage = np.copy(image)
    if mask is not None:
        cmap.set_bad((1, 0, 0, 1))
        z = np.where(mask == 1)
        pimage[z] = np.nan

    # Now draw the figure

    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.clf()
    axes1 = fig.add_subplot(111)

    axes1.imshow(pimage, vmin=minmax[0], vmax=minmax[1], cmap=cmap,
              origin='lower')
    
    pl.xlabel('Columns (pixels)')
    pl.ylabel('Rows (pixels)')

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())

    # Disable autoscaling
    axes1.autoscale(False)

       
    #
    # Overplot orders if requested
    #

    if orders_plotinfo is not None:

        xranges = orders_plotinfo['xranges']
        edgecoeffs = orders_plotinfo['edgecoeffs']
        orders = orders_plotinfo['orders']
        norders = len(orders)

        for i in range(norders):
            x = np.arange(xranges[i, 0], xranges[i, 1] + 1)
            bot = poly_1d(x, edgecoeffs[i, 0, :])
            top = poly_1d(x, edgecoeffs[i, 1, :])

            pl.plot(x, bot, color='purple')
            pl.plot(x, top, color='purple')
            pl.fill_between(x, bot, y2=top, color='purple', alpha=0.15)

            delta = xranges[i, 1] - xranges[i, 0]
            idx = np.fix(delta * 0.02).astype(int)

            pl.text(x[idx], (top[idx] + bot[idx]) / 2., str(orders[i]),
                    color='yellow', verticalalignment='center')

    #
    # Overplot traces if requested
    #

    if trace_plotinfo is not None:

        if 'x' in trace_plotinfo and \
                'y' in trace_plotinfo and \
                'goodbad' in trace_plotinfo:
            pl.plot(trace_plotinfo['x'], trace_plotinfo['y'], 'go',
                    markersize=2)
            bad = trace_plotinfo['goodbad'] == 0
            pl.plot(trace_plotinfo['x'][bad], trace_plotinfo['y'][bad], 'bo',
                    markersize=2)

        if 'fits' in trace_plotinfo:

            for i in range(len(trace_plotinfo['fits'])):
                pl.plot(trace_plotinfo['fits'][i][0, :],
                        trace_plotinfo['fits'][i][1, :], color='cyan',
                        linewidth=0.5)

    #
    # Overplot order locating data if requested.
    #

    if locateorders_plotinfo is not None:

        guess = locateorders_plotinfo['guess_positions']
        x = locateorders_plotinfo['x']
        y = locateorders_plotinfo['y']
        goodbad = locateorders_plotinfo['goodbad']
        coeffs = locateorders_plotinfo['edgecoeffs']

        for i in range(len(guess)):
            pl.plot(guess[i][0], guess[i][1], 'go')

        for i in range(len(x)):

            pl.plot(x[i], y[i], 'go', markersize=1)

            bad = goodbad[i] == 0
            pl.plot(x[i][bad], y[i][bad], 'bo', markersize=1)

            pl.plot(x[i], 
                    np.polynomial.polynomial.polyval(x[i], coeffs[i]),
                    'r-',linewidth=0.5)
                
    # fix layout overlap
    
    pl.tight_layout()

