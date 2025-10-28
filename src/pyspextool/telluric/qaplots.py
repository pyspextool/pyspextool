import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)
import matplotlib

from pyspextool.io.check import check_parameter
from pyspextool.plot.limits import get_spectra_range

def plot_deconvolve_line(
    plot_number:int,
    figure_size:tuple,
    font_size:float,
    spectrum_linewidth:float,
    spine_linewidth:float,
    line_wavelength:npt.ArrayLike,
    line_data_fluxdensity:npt.ArrayLike,
    line_model_fluxdensity:npt.ArrayLike,
    line_model_convolved_fluxdensity:npt.ArrayLike,
    line_fluxdensity_ratio:npt.ArrayLike,
    rms_deviation:float,
    maximum_deviation:float,
    plot_xlabel:str=None,
    plot_title:str=None):
    

    """
    To create the QA plot for ps.telluric.core.deconvolve_line.

    Creates a vertical two-panel plot.  The upper panel shows the normalized
    spectrum and the region over which the deconvolution was done and the
    lower plot shows the results of the deconvolution process.

    Parameters
    ----------
    plot_number : int or None
        The plot number to pass to matplotlib

    figure_size : tuple
        A (2,) tuple giving the figure size to pass to matplotlib

    font_size : int
        An int giving the font size to pass to matplotlib

    spectrum_linewidth : int or float
        An int or float giving the spectrum line width to pass to matplotlib

    spine_linewidth : int or float
        An int or float giving the spine line width to pass to matplotlib
    
    line_wavelength : ndarray
        A (nwave2,) array of wavelengths used in the deconvolution
        
    line_data_fluxdensity : ndarray
        A (nwave2,) array of data continuum-normalized, continuum-subtracted
        flux densities.

    line_model_fluxdensity : ndarray
        A (nwave2,) array of model continuum-normalized, continuum-subtracted
        flux densities.

    line_model_convolved_fluxdensity : ndarray
        A (nwave2,) array of model continuum-normalized, continuum-subtracted,
        convolved flux densities.

    line_fluxdensity_ratio : ndarray
        A (nwave2,) array of the ratio of line_data_fluxdensity and
        line_model_convolved_fluxdensity.
      
    rms_deviation: float
        A float giving the rms deviation of line_fluxdensity_ratio.
    
    maximum_deviation : float
        A float giving giving the maximum deviation of line_fluxdensity_ratio.

    plot_scale : float or int, default is 1
        A scale factor by which to increase the QA plots size.

    plot_number : int, default is None
        The plot number to pass to pl.figure.
    
    plot_xtitle : string, default is None
        A string giving the xtitle to pass directly to matplotlib.

    plot_title : string, default is None
        A string giving the title to pass directly to matplotlib.

       
    Returns
    -------
    None

    """

    #
    # Check the parameters
    #

    check_parameter('plot_deconvolve_line', 'plot_number', plot_number,
                    ['int','NoneType'])

    check_parameter('plot_deconvolve_line', 'figure_size', figure_size, 'tuple')

    check_parameter('plot_deconvolve_line', 'font_size', font_size,
                    ['int','float'])

    check_parameter('plot_deconvolve_line', 'spectrum_linewidth',
                    spectrum_linewidth, ['int','float'])        

    check_parameter('plot_deconvolve_line', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        
    
    check_parameter('plot_deconvolve_line', 'line_wavelength',
                    line_wavelength, 'ndarray',1)

    check_parameter('plot_deconvolve_line', 'line_data_fluxdensity',
                    line_data_fluxdensity, 'ndarray',1)

    check_parameter('plot_deconvolve_line', 'line_model_fluxdensity',
                    line_model_fluxdensity, 'ndarray',1)
    
    check_parameter('plot_deconvolve_line', 'line_model_convolved_fluxdensity',
                    line_model_convolved_fluxdensity, 'ndarray',1)

    check_parameter('plot_deconvolve_line', 'line_fluxdensity_ratio',
                    line_fluxdensity_ratio, 'ndarray',1)
        
    check_parameter('plot_deconvolve_line', 'plot_xlabel', plot_xlabel,
                    ['str','NoneType'])

    check_parameter('plot_deconvolve_line', 'plot_title', plot_title,
                    ['str','NoneType'])

    
    #
    # Set up the plot
    #    

    # removed helvetica - problem for windows OS
    font = {'weight' : 'normal',
            'size'   : font_size}

    matplotlib.rc('font', **font)
    
    fig = pl.figure(num=plot_number,
                    figsize=figure_size)

    pl.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.2)

    #
    # Do the top plot
    #
    
    axes1 = fig.add_subplot(111)

    # Get the yrange based on the data in that wavelength range
    
    yrange = get_spectra_range(line_data_fluxdensity,
                            line_model_fluxdensity,
                            line_model_convolved_fluxdensity,
                            line_fluxdensity_ratio-1, frac=0.1)
    
    axes1.set_ylim(yrange)
 
    # Plot the spectrum

    axes1.step(line_wavelength, line_data_fluxdensity,color='black',
               where='mid', lw=spectrum_linewidth)
    axes1.step(line_wavelength, line_model_fluxdensity, color='green',
               where='mid', lw=spectrum_linewidth)

    axes1.step(line_wavelength, line_model_convolved_fluxdensity,color='red',
               where='mid', lw=spectrum_linewidth)
    
    axes1.step(line_wavelength,line_fluxdensity_ratio-1,color='blue',
               where='mid', lw=spectrum_linewidth)

    # Plot the rms limit lines of 0.01.
    
    axes1.axhline(y=0.01,color='magenta',linestyle='dotted')
    axes1.axhline(y=-0.01,color='magenta',linestyle='dotted')

    # Deal with the tickmarks

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=spine_linewidth)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Label the axes
    
    axes1.set(xlabel=plot_xlabel, ylabel='Residual Normalized Flux Density')

    # Add all the information texts
        
    axes1.text(0.02, 0.2, 'Max Deviation = '+"{:.4f}".format(maximum_deviation),
               ha='left', va='bottom', transform=axes1.transAxes, color='black')

    axes1.text(0.02, 0.15, 'RMS Deviation = '+"{:.4f}".format(rms_deviation),
               ha='left', va='bottom', transform=axes1.transAxes, color='black')


    axes1.text(0.95, 0.3, 'A0 V', color='black', ha='right', va='bottom',
               transform=axes1.transAxes)

    axes1.text(0.95, 0.25, 'Scaled Vega', color='green', ha='right',
               va='bottom', transform=axes1.transAxes)

    axes1.text(0.95, 0.2, 'Scaled & Convolved Vega', color='red', ha='right',
               va='bottom', transform=axes1.transAxes)

    axes1.text(0.95, 0.15, 'Ratio', color='blue', ha='right',
               va='bottom', transform=axes1.transAxes)

        # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(spine_linewidth)

def plot_estimate_ewscales(
    plot_number:int,
    figure_size:tuple,
    font_size:int,
    spectrum_linewidth:int | float,
    spine_linewidth:int | float,
    xlabel:str,
    title:str,
    wavelength:npt.ArrayLike,
    default_fluxdensity:npt.ArrayLike,
    new_fluxdensity:npt.ArrayLike,
    atmospheric_transmission:npt.ArrayLike,
    line_wavelengths:npt.ArrayLike | float,
    line_scales:npt.ArrayLike | float,
    default_scale:float):

    #
    # Now start the plotting
    #
    
    # Set the fonts

    font = {'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)

    fig = pl.figure(num=plot_number,
                    figsize=figure_size)

    pl.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.2)

    xrange = [np.nanmin(wavelength),np.nanmax(wavelength)]

    #
    # Do the top plot
    #
    
    axes1 = fig.add_subplot(211)

    axes1.step(wavelength,default_fluxdensity,where='mid',color='black',
               label='Raw Correction')
    axes1.step(wavelength,new_fluxdensity,where='mid',color='red',
               label='Adjusted Correction')

    
    axes1.set_xlim(xrange)

    for i in range(len(line_scales)):

        axes1.axvline(x=line_wavelengths[i],linestyle='dotted')

    axes1.set_xticklabels([]) 

    axes1.set_ylabel('Relative Intensity')
    axes1.set_title(title)



    axes1b = axes1.twinx() 
    axes1b.step(wavelength, atmospheric_transmission, 
               where='mid',color='purple',label='Atmospheric Transmission')
    axes1b.set_ylim(ymin=0, 
                    ymax=1)
    axes1b.set_yticklabels([]) 
    axes1b.tick_params(right = False) 
    
    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                    which='both', direction='in', width=1.5)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())
    
    lines, labels = axes1.get_legend_handles_labels()
    lines2, labels2 = axes1b.get_legend_handles_labels()
    axes1.legend(lines + lines2, labels + labels2, loc=0)


    axes2 = fig.add_subplot(212)

    axes2.plot(line_wavelengths,line_scales,'or')

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                    which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())

    axes2.set_ylabel('EW Scale Factor')
    axes2.set_xlabel(xlabel)
    axes2.set_xlim(xrange)


    axes2.axhline(y=default_scale,linestyle='dashed',color='black')

    yrange = get_spectra_range(line_scales, [default_scale], frac=0.1)
    axes2.set_ylim(yrange)


def plot_measure_linerv(
    plot_number:int,
    figure_size:tuple,
    font_size:int,
    spectrum_linewidth:int | float,
    spine_linewidth:int | float,
    wavelength:npt.ArrayLike,
    object_flux:npt.ArrayLike,
    model_flux:npt.ArrayLike,
    lag:npt.ArrayLike,
    xcorrelation:npt.ArrayLike,
    offset:float,
    velocity:float,
    redshift:float,
    fit:npt.ArrayLike=None,
    plot_xlabel:str=None,
    plot_title:str=None):
    
    """
    To create a plot for the cross correlation in a device-independent way

    Parameters
    ----------
    plot_number : int or None
        The plot number to pass to matplotlib

    figure_size : tuple
        A (2,) tuple giving the figure size to pass to matplotlib

    font_size : int
        An int giving the font size to pass to matplotlib

    spectrum_linewidth : int or float
        An int or float giving the spectrum line width to pass to matplotlib

    spine_linewidth : int or float
        An int or float giving the spine line width to pass to matplotlib
    
    wavelength : ndarray
        A (nwave,) array of wavelengths.

    object_flux : ndarray
        A (nwave,) array of flux density values for the object.

    model_flux : ndarray
        A (nwave,) array of flux density values for the model.
    
    lag : ndarray
        A (nlag,) array of lag values used to perform the cross correlation

    xcorrelation : ndarray
        A (nlag,) array of cross correlation values.

    offset : float
        The number of pixels the xcorrelation peak is offset from 0.

    velocity : float
        The velocity in km s-1 corresponding to 'offset'.

    redshift : float
        (1 + 'velocity'/c)

    fit : ndarray, default=None
        A (nwave,) array of fitted values to the xcorrelation array.

    plot_xlabel : str, default=None
        A string given an optional x label.  Useful because the wavelength
        units may change.

    plot_title : str, default=None
        A string giving an optional title.

    
    Returns
    -------
    None
        May write a file to disk.  
    
    """
        
    #
    # Check the parameters
    #

    check_parameter('plot_measure_linerv', 'plot_number', plot_number,
                    ['int','NoneType'])

    check_parameter('plot_measure_linerv', 'figure_size', figure_size, 'tuple')

    check_parameter('plot_measure_linerv', 'font_size', font_size,
                    ['int','float'])

    check_parameter('plot_measure_linerv', 'spectrum_linewidth',
                    spectrum_linewidth, ['int','float'])        

    check_parameter('plot_measure_linerv', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        
    
    check_parameter('plot_measure_linerv', 'wavelength', wavelength,
                    'ndarray')

    check_parameter('plot_measure_linerv', 'object_flux', object_flux,
                    'ndarray')
    
    check_parameter('plot_measure_linerv', 'model_flux', model_flux, 'ndarray')

    check_parameter('plot_measure_linerv', 'lag', lag, 'ndarray')

    check_parameter('plot_measure_linerv', 'xcorrelation', xcorrelation,
                    'ndarray')

    check_parameter('plot_measure_linerv', 'offset', offset,
                    ['float64','float','int'])

    check_parameter('plot_measure_linerv', 'velocity', velocity,
                    ['float','float64'])

    check_parameter('plot_measure_linerv', 'fit', fit, ['NoneType','ndarray'])
    
    check_parameter('plot_measure_linerv', 'plot_xlabel', plot_xlabel,
                    ['Nonetype','str'])    

    check_parameter('plot_measure_linerv', 'plot_title', plot_title,
                    ['NoneType','str'])    

    
    #
    # Set the fonts
    #
    
    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    matplotlib.rc('font', **font)

    #
    # Create the figure
    #
    
    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.2)

    #
    # Create the cross correlation plot
    #

    # Get the y range
    
    if fit is not None:

        yrange = get_spectra_range(xcorrelation, fit, frac=0.1)

    else:

        yrange = get_spectra_range(xcorrelation, frac=0.1)        

    # Build the figure
        
    axes1 = fig.add_subplot(211)    
    axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes1.step(lag, xcorrelation, 'black',lw=spectrum_linewidth)
    axes1.set(xlabel='lag (pixels)')
    axes1.set_title(plot_title)    

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in',width=spine_linewidth)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())

    axes1.axvline(x=offset,color='black',linestyle='dotted')
    
    
    axes1.text(0.05, 0.9, 'X Correlation', color='black', ha='left',
                transform=axes1.transAxes)

    axes1.text(0.95, 0.9, 'Offset='+'%+.2f' % offset+' pixels',
               color='black', ha='right', transform=axes1.transAxes)

    axes1.text(0.95, 0.8, 'Velocity='+'%+.2f' % velocity+' km s$^{-1}$',
               color='black', ha='right', transform=axes1.transAxes)    

    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(spine_linewidth)

    if fit is not None:
        
        axes1.step(lag, fit, 'r')
        axes1.text(0.05, 0.8, 'Fit', color='r', ha='left',
                   transform=axes1.transAxes)

    #     
    # Show the data, spectrum, and shifted spectrum
    #

    # Get the y range
    
    yrange = get_spectra_range(object_flux, model_flux, frac=0.1)

    # Build the figure
    
    axes2 = fig.add_subplot(212)
    axes2.margins(x=0)
    axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes2.step(wavelength, object_flux, 'black', label='spectrum',
               lw=spectrum_linewidth)
    axes2.step(wavelength, model_flux, 'r',linestyle='dashed',
               label='Model (0 km s$^{-1}$)',
               lw=spectrum_linewidth)
    axes2.step(wavelength*(1+redshift), model_flux, 'r',
               label='Model (%+.2f' % velocity+' km s$^{-1}$)',
               lw=spectrum_linewidth)  
    axes2.set(xlabel=plot_xlabel, ylabel='Relative Flux Density')

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=spine_linewidth)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())    


    axes2.legend(bbox_to_anchor=(0.5,0.1),
                 ncols=1,
                 frameon=False,
                 loc='lower left',
                 handlelength=1)


    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(spine_linewidth)


def plot_shifts(
    plot_number:int,
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
    # removed helvetica

    font = {
    #'family' : 'helvetica',
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



        


