import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)

from pyspextool.io.check import check_parameter
from pyspextool.plot.limits import get_spectra_range

def plot_shifts(plot_number:int,
                subplot_size:tuple,
                subplot_stackmax:int,
                font_size:int,
                scale:int | float,
                spectrum_linewidth:int | float,
                spine_linewidth:int | float,
                xlabel:str,
                orders:npt.ArrayLike,
                object_spectra:npt.ArrayLike,
                rawtc_spectra:npt.ArrayLike,
                shiftedtc_spectra:npt.ArrayLike,                
                shift_ranges:npt.ArrayLike,
                shifts:npt.ArrayLike,
                reverse_order:bool=True):

    """
    To create a QA plot for the standard shifts

    Parmaeters
    ----------
    plot_number : int or None
        The plot number to pass to matplotlib

    subplot_size : tuple
        A (2,) tuple giving the plot size for a single order.

    subplot_stackmax : int
        The maximum number of orders to plot in the vertical direction

    font_size : int
        An int giving the font size to pass to matplotlib

    scale : int or float
        A scale factor to multiple the font and final figure size by.
    
    spectrum_linewidth : int or float
        An int or float giving the spectrum line width to pass to matplotlib

    spine_linewidth : int or float
        An int or float giving the spine line width to pass to matplotlib

    xlabel : str
        A string giving the x axis label.

    orders : ndarray
        An (norders,) array of order numbers.

    object_spectra : ndarray
        An (norders*napertures, 4, nwavelength) array of object spectra.
    
    rawtc_spectra : ndarray
        An (norders*napertures, 4, nwavelength) array of raw telluric
        correction spectra.

    shiftedtc_spectra : ndarray
        An (norders*napertures, 4, nwavelength) array of shifted
        telluric correction spectra.

    shift_ranges : ndarray
        An (norders,2) array of shift ranges.  shift_ranges[0,:] gives the
        lower and upper wavelength limit over which the shift was determined.
        If no shift is requsted, the values are np.nan.

    shifts : ndarray
        An (norders,napertures) array of shift ranges.  shift_ranges[0,:]
        gives the lower and upper wavelength limit over which the shift was
        determined.  If no shift is requsted, the values are np.nan.

    reverse_order : {True, False}
        Set to True to plot the orders in decreasing order, e.g. 8,7,6,5
        Set to False to plot the orders in increasing order, e.g. 5,6,7,8
      
    Returns
    -------
    None
    
    """

    #
    # Check parameters
    #

    check_parameter('plot_shifts', 'plot_number', plot_number,
                    ['int', 'NoneType'])

    check_parameter('plot_shifts', 'subplot_size', subplot_size, 'tuple')

    check_parameter('plot_shifts', 'subplot_stackmax', subplot_stackmax, 'int')

    check_parameter('plot_shifts', 'font_size', font_size, 'int')

    check_parameter('plot_shifts', 'scale', scale, ['int','float'])

    check_parameter('plot_shifts', 'spectrum_linewidth', spectrum_linewidth,
                    ['int','float'])        

    check_parameter('plot_shifts', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        

    check_parameter('plot_shifts', 'xlabel', xlabel, 'str')

    check_parameter('plot_shifts', 'orders', orders, 'ndarray')
    
    check_parameter('plot_shifts', 'object_spectra', object_spectra, 'ndarray')

    check_parameter('plot_shifts', 'rawtc_spectra', rawtc_spectra, 'ndarray')

    check_parameter('plot_shifts', 'shiftedtc_spectra', shiftedtc_spectra,
                    'ndarray')    
    
    check_parameter('plot_shifts', 'shift_ranges', shift_ranges, 'ndarray')

    check_parameter('plot_shifts', 'shifts', shifts, 'ndarray')        

    check_parameter('plot_shifts', 'reverse_order', reverse_order, 'bool') 
       
    #
    # Get important information
    #

    # Get norders and napertures for the object spectra.
    
    norders = len(orders)

    napertures = int(np.shape(object_spectra)[0]/norders)
           
    # Determine which orders were shifted

    zshifted = ~np.isnan(shift_ranges[:,0])    
    shifted_norders = len(orders[zshifted])
    
    #
    # Now start the plotting
    #
    
    # Set the fonts

    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size*scale}

    rc('font', **font)
    
    # Determine the plot size
    
    ncols = np.ceil(shifted_norders*napertures / subplot_stackmax).astype(int)

    nrows = np.min([shifted_norders*napertures, subplot_stackmax]).astype(int)

    plot_index = np.arange(1,nrows*ncols+1)
    
    plot_index = np.reshape(np.reshape(plot_index,(nrows,ncols)),
                            ncols*nrows,order='F')
    
    figure_size = (subplot_size[0]*ncols*scale, subplot_size[1]*nrows*scale)
    
    #
    # Make the figure
    #
    
    pl.figure(num=plot_number,
              figsize=figure_size)
    pl.clf()    
    pl.subplots_adjust(hspace=0.5,
                       wspace=0.2,
                       left=0.1,
                       right=0.95,
                       bottom=0.075,
                       top=0.95)

    m = 0
    for i in range(norders):

        if reverse_order is True:

            # We are going to plot things in reverse order.
            
            order_idx = norders-i-1

        else:

            order_idx = i

        # Did this order get shifted?


        if zshifted[order_idx].item() is False:

            continue

        for j in range(napertures):

            # Find the correct order index

            if reverse_order is True:

                object_idx = norders*napertures-i*napertures-napertures


            else:

                object_idx = i*napertures

            # Now add the aperture values

            object_idx = object_idx + j

            # Create the telluric corrected spectra
                
            wavelength = object_spectra[object_idx,0,:]
            object_flux = object_spectra[object_idx,1,:]
            raw_telluric = rawtc_spectra[object_idx,1,:]
            shifted_telluric = shiftedtc_spectra[object_idx,1,:]
            
            raw_ratio = object_flux *raw_telluric
            shifted_ratio = object_flux*shifted_telluric            
            
            # Clip out the shift range.

            zshift = np.where((wavelength >= shift_ranges[order_idx,0]) &
                              (wavelength <= shift_ranges[order_idx,1]))[0]

            wavelength = wavelength[zshift]
            raw_ratio = raw_ratio[zshift]
            shifted_ratio = shifted_ratio[zshift]                        

            # Normalize intensity

            raw_ratio /= np.nanmedian(raw_ratio)
            shifted_ratio /= np.nanmedian(shifted_ratio)            
            
            # Get the plot range

            xrange = get_spectra_range(wavelength)
            yrange = get_spectra_range(shifted_ratio, raw_ratio, frac=0.1)

            # Do the plot

            axe = pl.subplot(nrows, ncols, plot_index[m])

            axe.step(wavelength,
                     raw_ratio,
                     color='grey',
                     where='mid',
                     label='Raw')
            
            axe.step(wavelength,
                     shifted_ratio,
                     color='green',
                     where='mid',
                     label='Shifted')
           
            axe.set_xlim(xrange)
            axe.set_ylim(yrange)

            axe.set_title('Order ' + str(orders[order_idx])+', aperture '+\
                          str(j+1)+r', $\Delta x$='+'$'+\
                          str(shifts[order_idx,j])+'$ pixels')
            axe.set_ylabel('Relative Intensity')
            axe.set_xlabel(xlabel)            

            axe.xaxis.set_minor_locator(AutoMinorLocator())    
            axe.tick_params(right=True, left=True, top=True, bottom=True,
                            which='both', direction='in', width=1.5)
            axe.tick_params(which='minor', length=3)
            axe.tick_params(which='major', length=5)
            axe.yaxis.set_minor_locator(AutoMinorLocator())
            
            if m == 0:

                axe.legend()

            m += 1
