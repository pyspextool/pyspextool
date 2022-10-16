import numpy as np
from scipy import interpolate

from pyspextool.fit.polyfit import poly_1d
from pyspextool.io.check import check_parameter
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.arrays import find_index
from pyspextool.utils.loop_progress import loop_progress

def trace_to_xy(order_mask, wavecal, spatcal, xranges, orders, doorders,
                naps, tracecoeffs, clupdate=True):

    '''
    To convert a trace in wavelength/arcseconds to x/y in an image.

    Parameters
    ----------
    order_mask : numpy.ndarray
        An int (nrows, cols) order mask.  A pixel is set to its 
        order number.  Interorder pixels are set to zero.

    wavecal : numpy.ndarray
        A float (nrows, cols) wavelength image.  A pixel is set to its 
        wavelength.  Interorder pixels are set to NaN.  

    spatcal : numpy.ndarray
        A float (nrows, cols) spatial image.  A pixel is set to its 
        spatial position.  Interorder pixels are set to NaN.  

    xranges : numpy.ndarray
        A float (norders, 2) array giving the column numbers over which 
        to operate.  xranges[0,0] gives the starting column number for 
        the order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.
    
    orders : numpy.ndarray 
        An int (norders,) array giving the order numbers.  By Spextool 
        convention, `orders[0]` is the order closest to the bottom 
        of the array.

    doorders : numpy.ndarray 
        An int (norders,) array indicating whether an order should be
        operated on.  1=yes, 0=no.

    naps : int
        The number of apertures.

    tracecoeffs : numpy.ndarray
        A float (norders*naps, ncoeffs) array giving the polynomial 
        coefficents to convert from wavelength to spatial angle in the slit.

    clupdate : {True, False}, optional
        Set to True for command line updates during execution. 

    Returns
    -------
    list
    A (norders*naps) list where each element is a (2, ndat) numpy.ndarray
    giving the trace points in (x,y) image coordinates.
    x values = array[0,:], y values = array[0,:].  

    Notes
    -----
    
    Example
    -------

    '''

    #
    # Check parameters
    #
    check_parameter('trace_to_xy', 'order_mask', order_mask, 'ndarray')

    check_parameter('trace_to_xy', 'wavecal', wavecal, 'ndarray')

    check_parameter('trace_to_xy', 'spatcal', spatcal, 'ndarray')

    check_parameter('trace_to_xy', 'xranges', xranges, 'ndarray')

    #
    # Get set up
    #

    nrows, ncols = np.shape(order_mask)
    
    norders = len(orders)
    donorders = np.sum(doorders)
    
    xx, yy = make_image_indices(nrows, ncols)

    plot_fits = []
    l = 0

    #
    # Start the loop over each 
    #
    
    for i in range(norders):

        # Check the doorders array

        if doorders[i] == 0:
            continue

        # Define the start and end points
        
        starts = xranges[i,0]
        stops = xranges[i,1]
        ncols = stops-starts+1

        # Create the x array for this array
        
        x = np.arange(ncols)+starts

        # Get the wavelengths

        z = order_mask == orders[i]
        waves = np.unique(wavecal[z])

        # Now determine the trace position in angular units at each wavelength
        
        trace_ang = np.full((naps, ncols), np.nan)

        for j in range(naps):

            trace_ang[j,:] = poly_1d(waves, tracecoeffs[l,:])            

            loop_progress(l, 0, donorders*naps,
                          message='Collecting plotting data...')
            l += 1

        # Convert these angular positons to pixel coordinates
            
        trace_pix = np.full((naps, ncols),np.nan)            
        for j in range(ncols):

            omaskz = order_mask[:, x[j]]
            z = omaskz == orders[i]

            slits = spatcal[z, x[j]]
            slity = yy[z, x[j]]            

            f = interpolate.interp1d(slits,slity)
            trace_pix[:,j] =f(trace_ang[:,j])

        # Store them in a list
            
        for j in range(naps):

            plot_fits.append(np.stack([x, trace_pix[j,:]]))
            
    return plot_fits
