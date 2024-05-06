import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.ticker import (AutoMinorLocator)
import typing
import os
import scipy

from pyspextool.setup_utils import pySpextoolError
from pyspextool.io.check import check_parameter
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.plot.limits import get_spec_range
from pyspextool.utils.arrays import find_index

def model_xcorrelate(object_wavelength:npt.ArrayLike,
                     object_fluxdensity:npt.ArrayLike,
                     model_wavelength:npt.ArrayLike,
                     model_fluxdensity:npt.ArrayLike,
                     minimum_wavelength:typing.Optional[float]=None,
                     maximum_wavelength:typing.Optional[float]=None,
                     resolving_power:typing.Optional[float]=None,
                     qashow_info:typing.Optional[dict]=None,
                     qafile_info:typing.Optional[dict]=None):

    """
    To determine the velocity shift between a model and data using a line

    Parameters
    ----------

    object_wavelength : ndarray
        A (nwave1,) array of wavelengths for the object in microns.

    object_fluxdensity : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the object.

    model_wavelength : ndarray
        A (nwave2,) array of wavelengths for a model in microns.

    model_fluxdensity : ndarray
        A (nwave2,) array of continuum-normalized flux densities for a model.

    minimum_wavelength : int or float, default=None
        The minimum wavelength value over which to perform the xcorrelation.

    maximum_wavelength : int or float, default=None
        The maximum wavelength value over which to perform the xcorrelation.

    resolving_power : int or float, default=None
        The resolving power of the observations.  A not-perfect attempt to 
        convolve the model to the same resolving power is then made.

    qashow_info : None or dict

        `'plot_number'` : int
            The plot number to be passed back in to update the same plot.
            Useful if you are doing to-screen plotting.
    
        `'plot_scale'` : float or int
            A scale factor to increase the size of the plot.  

        `'block'`: {False, True}, optional
            Set to make the plot block access to the command line, e.g.
            pl.ioff().

        `'plot_xlabel'` : str, optional
            A latex string giving the xlabel.

        `'plot_title'` : str, optional
            A latex string giving the title of the plot.

        
    qafile_info : dict, optional
        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

        `'plot_xlabel'` : str, optional
            A latex string giving the xlabel.

        `'plot_title'` : str, optional
            A latex string giving the title of the plot.

    
    
    Returns
    -------
    ndarray
        The velocity shift of the data relative to the model in km s-1.  

    Any resulting QA plot is written to disk as,
    os.path.join(qa_fileinfo['filepath'],qa_fileinfo['filename'] + \
                 qa_fileinfo['extension'])

           
    """

    #
    # Check Parameters
    #

    check_parameter('model_xcorrelate', 'object_wavelength', object_wavelength,
                    'ndarray', 1)

    check_parameter('model_xcorrelate', 'object_fluxdensity',
                    object_fluxdensity, 'ndarray', 1)

    check_parameter('model_xcorrelate', 'model_wavelength', model_wavelength,
                    'ndarray', 1)

    check_parameter('model_xcorrelate', 'model_nflux', model_fluxdensity,
                    'ndarray', 1)

    check_parameter('model_xcorrelate', 'minimum_wavelength',
                    minimum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('model_xcorrelate', 'maximum_wavelength',
                    maximum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('model_xcorrelate', 'resolving_power',
                    resolving_power, ['NoneType', 'int', 'float'])        
    
    check_parameter('model_xcorrelate', 'qashow_info', qashow_info,
                    ['NoneType', 'dict'])
    
    check_parameter('model_xcorrelate', 'qafile_info', qafile_info,
                    ['NoneType', 'dict'])
    
    #
    # Get set up
    #

    cspeed = 2.99792458E5 # km/s

    # zero the continuum

    model_zflux = model_fluxdensity - 1
    object_zflux = object_fluxdensity-1

    #
    # Now get the wavelength range over which to do the x-correlation
    #
    
    # Get wavelength range from the object spectrum

    if minimum_wavelength is None:

        minimum_wavelength = np.min(object_wavelength)

    if maximum_wavelength is None:

        maximum_wavelength = np.max(object_wavelength)        

    # Now confirm the range is valid.

    if minimum_wavelength >= maximum_wavelength:

        message = '`minimum_wavelenth` >= `maximum_wavelength`.'
        raise pySpextoolError(message)

    #
    # Now check to make sure the object has fewer pixels than the model
    #

    nobject = np.size(object_wavelength)

    z = np.logical_and((model_wavelength > minimum_wavelength),
                       (model_wavelength < maximum_wavelength))

    nmodel = np.sum(z)

    if nobject > nmodel:

        message = 'Data has a higher resolution than the model.'
        raise pySpextoolError(message)

    npixels_lnlambda = np.max([nobject, nmodel])

    #
    # Smooth the model to get close
    #

    if resolving_power is not None:

        # Determine the wavelength mid point

        wavelength = (minimum_wavelength+maximum_wavelength)/2.

        # Find the dispersion of the Vega model at this wavelength
        
        idx = int(find_index(model_wavelength, wavelength))

        model_dispersion = model_wavelength[idx]-model_wavelength[idx-1] 
        
        # Determine the number of pixels for the resolving power.

        fwhm_kernel = int(wavelength/resolving_power/model_dispersion)

        # Create a Gaussian kernel

        npixels_kernel = fwhm_kernel*5
        if npixels_kernel % 2 == 0:

            npixels_kernel +=1

        x = np.arange(npixels_kernel)-npixels_kernel//2
        
        sigma = fwhm_kernel/2.354

        gaussian = np.exp(-(x/sigma)**2 / 2)
        gaussian /= np.sum(gaussian)

        model_zflux = np.convolve(model_zflux,gaussian, mode='same')
        
    #
    # Resampling to a constant spacing in ln lambda: v/c = d(ln lambda)
    #

    # Create the wavelength array
    
    acoeff = (npixels_lnlambda - 1)/(np.log(maximum_wavelength) -
                                     np.log(minimum_wavelength))
    bcoeff  = npixels_lnlambda - (acoeff * (np.log(maximum_wavelength)))
  
    xpon = np.arange(npixels_lnlambda)+1.0
    lnlambda_wavelengths = np.exp((xpon-bcoeff)/acoeff)

    # Do the resampling

    f = scipy.interpolate.interp1d(object_wavelength, object_zflux)
    object_resampled_zflux = f(lnlambda_wavelengths)

    f = scipy.interpolate.interp1d(model_wavelength, model_zflux)
    model_resampled_zflux = f(lnlambda_wavelengths)    

    #
    # Do the cross correlation
    #

    xcor = scipy.signal.correlate(object_resampled_zflux, model_resampled_zflux,
                                  mode='same', method='fft')
    xcor = xcor / np.nanmax(xcor)

    lag = scipy.signal.correlation_lags(npixels_lnlambda, npixels_lnlambda,
                                        mode='same')

    #
    # Fit the cross correlation
    #

    fit = fit_peak1d(lag, xcor, nparms=4, positive=True)
    offset_pixels = fit['parms'][1]

    velocity_shift = (cspeed * offset_pixels)/acoeff
    redshift = velocity_shift/cspeed

    #
    # Make the QA plot
    #

    if qashow_info is not None:

        if qashow_info['block'] is True:

            pl.ioff()

        else:

            pl.ion()

        plotnum = plot_model_xcorrelate(lnlambda_wavelengths,
                                        object_resampled_zflux,
                                        model_resampled_zflux,
                                        lag, xcor, fit['fit'],
                                        offset_pixels, velocity_shift,
                                        redshift,
                                        plot_scale=qashow_info['plot_scale'],
                                        plot_number=qashow_info['plot_number'],
                                        plot_xlabel=qashow_info['plot_xlabel'],
                                        plot_title=qashow_info['plot_title'])
                               
        pl.show()
        if qashow_info['block'] is False:
            pl.pause(1)
        
    if qafile_info is not None:

        plot_model_xcorrelate(lnlambda_wavelengths,
                              object_resampled_zflux,
                              model_resampled_zflux,
                              lag, xcor, fit['fit'], offset_pixels,
                              velocity_shift, redshift,
                              plot_xlabel=qashow_info['plot_xlabel'],
                              plot_title=qashow_info['plot_title'])
                
        pl.savefig(os.path.join(qafile_info['filepath'],
                                qafile_info['filename'] + \
                                qafile_info['extension']))
        pl.close()

        plotnum = None

        
    return velocity_shift, redshift, plotnum
    

        
def normalize_continuum(wavelength:npt.ArrayLike,
                        fluxdensity:npt.ArrayLike,
                        ranges:npt.ArrayLike,
                        degree:int,
                        robust:typing.Optional[dict]=None,
                        qashow_info:typing.Optional[str]=None,
                        qafile_info:typing.Optional[dict]=None):

    
    """
    To normalize the continuum of a spectrum using a robust polynomial

    Parameters
    ----------
    wavelength : ndarray
        An (nwave,) array of wavelengths.

    fluxdensity : ndarray
        An (nwave,) array of flux densities.

    ranges : ndarray
        A (2*nrange,) list of wavelength ranges over which to fit the continuum.

    degree : int
        The polynomial degree of the fit.

    robust : dict, optional
        `"threshold"`: int or float
            The sigma threshold over which a pixel is identified as bad in the 
            sigma clipping loop, e.g. 

                     |x_i - model|
            bad if  ---------------   >  threshold
                         sigma 

        `"epsilon"` : int or float
            The fractional change in the standard deviation below which the 
            sigma clipping loop stops.  That is,

                      old_stddev - new_stddev
            stop if  -------------------------   < epsilon
                           new_stddev

    qashow_info : None or dict
        `'plot_number'` : int
            The plot number.  Useful if you are doing to-screen plotting
            because you can plot to the same window multiple times.
    
        `'plot_scale'` : float or int
            A scale factor to increase the size of the plot over the default.

        `'block'`: {False, True}, optional
            Set to make the plot block access to the command line, e.g.
            pl.ioff().

        `'plot_xlabel'` : str, optional
            A latex string giving the xlabel.

        `'plot_title'` : str, optional
            A latex string giving the title of the plot.
        
    qafile_info : dict, optional    
        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

        `'plot_xlabel'` : str, optional
            A latex string giving the xlabel.

        `'plot_title'` : str, optional
            A latex string giving the title of the plot.

        
    Returns
    -------
    ndarray, int
    

    Notes
    -----
    TODO:  Must include proper error propagation.  


    """

    #
    # Check the parameters
    #

    check_parameter('normalize_continuum', 'wavelength', wavelength,
                    'ndarray', 1)

    check_parameter('normalize_continuum', 'fluxdensity', fluxdensity,
                    'ndarray', 1)

    check_parameter('normalize_continuum', 'ranges', ranges, 'ndarray', 1)

    check_parameter('normalize_continuum', 'degree', degree, 'int')

    check_parameter('normalize_continuum', 'robust', robust,
                    ['NoneType', 'dict'])

    check_parameter('normalize_continuum', 'qashow_info', qashow_info,
                    ['NoneType', 'dict'])
    
    check_parameter('normalize_continuum', 'qafile_info', qafile_info,
                    ['NoneType', 'dict'])

    #
    # Determine if the wavelengths are monotonically increasing
    #

    if np.all(np.diff(ranges) > 0) == False:
        
        message = "The normalization ranges, "+str(ranges)+\
            ", are not monotonic."
        raise pySpextoolError(message)

    #
    # Make sure all wavelengths fall within the input range of `wavelength`.
    #

    min_wavelength_spectrum = np.nanmin(wavelength)

    max_wavelength_spectrum = np.nanmax(wavelength)

    if np.all(ranges > min_wavelength_spectrum) == False:

        message = "The noramlization ranges, "+str(ranges)+\
            ", do not fall within the wavelength range of the spectrum, "+\
        str([min_wavelength_spectrum,max_wavelength_spectrum])+"."
        raise pySpextoolError(message)

    if np.all(ranges < max_wavelength_spectrum) == False:

        message = "The noramlization ranges, "+str(ranges)+\
            ", do not fall within the wavelength range of the spectrum, "+\
        str([min_wavelength_spectrum,max_wavelength_spectrum])+"."
        raise pySpextoolError(message)
            
    #
    # Identify the wavelength to be fit.
    #

    nranges = np.size(ranges)//2

    # Generate a list where each region gets its element which is a
    # True/False mask
    
    masks = []
    for i in range(nranges):

        zleft = (wavelength > ranges[i*2])
        zright = (wavelength < ranges[i*2+1])
        
        zselection = np.logical_and(zleft,zright)
        masks.append(zselection)

    #
    # Now combine True/False mask to create a single True/False mask
    #
    
    zregions = np.array(masks[0])
    
    for i in range(len(masks)-1):

        np.logical_or(zregions, masks[i+1], out=zregions)
        
    #
    # Now fit the continuum and normalize the results
    # 

    result = poly_fit_1d(wavelength[zregions], fluxdensity[zregions], degree)

    continuum = poly_1d(wavelength,result['coeffs'])
    
    normalized_fluxdensity = fluxdensity/continuum

    #
    # Make the QA plot
    #

    figure_size = (9,6)
    
    if qashow_info is not None:
        
        if qashow_info['block'] is True:

            pl.ioff()

        else:

            pl.ion()

        plotnum = plot_normalization(wavelength,
                                     fluxdensity,
                                     continuum,
                                     zregions,
                                     plot_scale=qashow_info['plot_scale'],
                                     plot_number=qashow_info['plot_number'],
                                     plot_xlabel=qashow_info['plot_xlabel'],
                                     plot_title=qashow_info['plot_title'])
            
        pl.show()

        if qashow_info['block'] is False:
            pl.pause(1)
        
    if qafile_info is not None:

        plot_normalization(wavelength,
                           fluxdensity,
                           continuum,
                           zregions,
                           plot_xlabel=qafile_info['plot_xlabel'],
                           plot_title=qafile_info['plot_title'])
        

        pl.savefig(os.path.join(qafile_info['filepath'],
                                qafile_info['filename'] + \
                                qafile_info['extension']))
        pl.close()

        plotnum = None

    return normalized_fluxdensity, plotnum



def plot_normalization(wavelength:npt.ArrayLike,
                       intensity:npt.ArrayLike,
                       continuum:npt.ArrayLike,
                       normalization_mask:npt.ArrayLike,
                       plot_scale:float=1.0,
                       plot_number:typing.Optional[int]=None,
                       plot_xlabel:typing.Optional[str]=None,
                       plot_title:typing.Optional[str]=None):

    
    """
    To plot the results of normalize_continuum in a device independent way

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths.

    intensity : ndarray
        A (nwave,) array of "intensities".

    continuum : ndarray
        A (nwave,) array of the fitted continuum values.

    normalization_mask : ndarray
        A (nwave,) True/False array where pixels used in the fit are True.
    
    figure_size : tuple
        A (2,) tuple giving the figure size.

    plot_number : int or None
        The plot number to be passed back in to update the same plot.
        Useful if you are doing to-screen plotting.

    plot_xlabel : str, optional, default=None
        A latex string giving the xlabel.

    plot_title : str, optional, default=None
        A latex string giving the title of the plot.
    
    Returns
    -------
    int
    The plot number of the window.

    """

    #
    # Check the parameters
    #
    
    check_parameter('plot_normalization', 'wavelength', wavelength, 'ndarray')

    check_parameter('plot_normalization', 'intensity', intensity, 'ndarray')
    
    check_parameter('plot_normalization', 'continuum', continuum, 'ndarray')

    check_parameter('plot_normalization', 'normalization_mask',
                    normalization_mask, 'ndarray')

    check_parameter('plot_normalization', 'plot_scale', plot_scale,
                    ['int','float'])

    check_parameter('plot_normalization', 'plot_number', plot_number,
                    ['int','NoneType'])

    check_parameter('plot_normalization', 'plot_xlabel', plot_xlabel, 'str')

    check_parameter('plot_normalization', 'plot_title', plot_title, 'str')    
        
    #
    # Make the two-panel figure
    #

    figure_size = (9,6)
    scaled_figure_size = [figure_size[0]*plot_scale,figure_size[1]*plot_scale]

    font_size = 12
    scaled_font_size = font_size*plot_scale
    
    # Set the fonts

    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : scaled_font_size}

    matplotlib.rc('font', **font)

    # Start the figure, and set the spacing
    
    fig = pl.figure(num=plot_number, figsize=scaled_figure_size)
    pl.clf()
    pl.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.05)
    
    # Get the plot range for x axis
    
    xrange = get_spec_range(wavelength[normalization_mask], frac=0.1)

    # Determine which pixels fall entirely within the plot range
    
    zleft = (wavelength > xrange[0])
    zright = (wavelength < xrange[1])

    zplotrange = np.logical_and(zleft,zright)    

    #
    # Create the spectral plot with the fit
    #

    # Determine the yrange for spectral plot
    
    yrange = get_spec_range([intensity[zplotrange], continuum[zplotrange]],
                            frac=0.1)

    axes1 = fig.add_subplot(211)    
    axes1.step(wavelength, intensity, 'black')
    axes1.set_title(plot_title)
    axes1.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes1.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes1.set_ylabel('Flux Density')

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5,
                      labelbottom=False)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())


    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(1.5)
    
    # Now plot the fitted pixels in red
    
    tmp = np.copy(intensity)
    tmp[~normalization_mask] = np.nan
    
    axes1.step(wavelength, tmp, 'red')

    # Plot the continuum
    
    axes1.step(wavelength, continuum, 'green')    

    #
    # Now plot the normalized spectrum
    #

    # Normalize and get the plot range
    
    normalized = intensity/continuum
    yrange = get_spec_range(normalized[zplotrange], frac=0.1)


    axes2 = fig.add_subplot(212)    
    axes2.step(wavelength, normalized, 'black')
    axes2.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes2.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes2.set_xlabel(plot_xlabel)
    axes2.set_ylabel('Normalized Flux Density')    

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())    
    
    # Plot the fitted pixels in red
    
    tmp = np.copy(normalized)
    tmp[~normalization_mask] = np.nan
    
    axes2.step(wavelength, tmp, 'red')
    axes2.axhline(y=1, linestyle='--', color='green')

    # change all spines
    for axis in ['top','bottom','left','right']:
        axes2.spines[axis].set_linewidth(1.5)

    
    #
    # Get the plot number and return the results
    #
    
    plot_number = pl.gcf().number
    return plot_number

    

def plot_model_xcorrelate(wavelength:npt.ArrayLike,
                          object_flux:npt.ArrayLike,
                          model_flux:npt.ArrayLike,
                          lag:npt.ArrayLike,
                          xcorrelation:npt.ArrayLike,
                          fit:npt.ArrayLike,
                          offset:float,
                          velocity:float,
                          redshift:float,
                          plot_scale:float=1.0,
                          plot_number:typing.Optional[int]=None,
                          plot_xlabel:typing.Optional[str]=None,
                          plot_title:typing.Optional[str]=None):
    """
    To create a plot for the cross correlation in device independent way

    Parameters
    ----------
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

    fit : ndarray
        A (nlag,) array of fitted values to the xcorrelation array.

    offset : float
        The number of pixels the xcorrelation peak is offset from 0.

    velocity : float
        The velocity in km s-1 corresponding to 'offset'.

    redshift : float
        (1 + 'velocity'/c)

    figure_size : tuple
        A (2,) tuple giving the figure size.

    Returns
    -------
    int
    The plot number of the window.

    
    """
        
    #
    # Check the parameters
    #

    check_parameter('plot_model_xcorrelate', 'wavelength', wavelength,
                    'ndarray')

    check_parameter('plot_model_xcorrelate', 'object_flux', object_flux,
                    'ndarray')

    check_parameter('plot_model_xcorrelate', 'model_flux', model_flux,
                    'ndarray')

    check_parameter('plot_model_xcorrelate', 'lag', lag, 'ndarray')

    check_parameter('plot_model_xcorrelate', 'xcorrelation', xcorrelation,
                    'ndarray')

    check_parameter('plot_model_xcorrelate', 'fit', fit, 'ndarray')

    check_parameter('plot_model_xcorrelate', 'offset', offset,
                    ['float64','float','int'])

    check_parameter('plot_model_xcorrelate', 'velocity', velocity,
                    ['float','float64'])

#    check_parameter('plot_model_xcorrelate', 'figure_size', figure_size,
#                    'tuple')

    
    #
    # Make the figure
    #

    figure_size = (6,9)
    scaled_figure_size = (figure_size[0]*plot_scale,
                          figure_size[1]*plot_scale)
    
    font_size = 12
    scaled_font_size = font_size*plot_scale
    
    # Set the fonts

    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : scaled_font_size}

    matplotlib.rc('font', **font)
    
    fig = pl.figure(num=plot_number, figsize=scaled_figure_size)
    pl.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                       hspace=0.2)


    # Create the spectral plot

    # Create the cross correlation plot

    yrange = get_spec_range(xcorrelation, fit, frac=0.1)
    
    axes1 = fig.add_subplot(211)    
    axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes1.step(lag, xcorrelation, 'black')
    axes1.step(lag, fit, 'r')
    axes1.set(xlabel='lag (pixels)')
    axes1.set_title(plot_title)    

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())

    axes1.axvline(x=0,color='black',linestyle='dotted')
    
    
    axes1.text(0.05, 0.9, 'X Correlation', color='black', ha='left',
                transform=axes1.transAxes)

    axes1.text(0.05, 0.8, 'Fit', color='r', ha='left',
              transform=axes1.transAxes)

    axes1.text(0.95, 0.9, 'Offset='+'%+.2f' % offset+' pixels',
               color='black', ha='right', transform=axes1.transAxes)

    axes1.text(0.95, 0.8, 'Velocity='+'%+.2f' % velocity+' km s$^{-1}$',
               color='black', ha='right', transform=axes1.transAxes)    

    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(1.5)


    # Show the data, spectrum, and shifted spectrum

    
    yrange = get_spec_range(object_flux, model_flux, frac=0.1)
    
    axes2 = fig.add_subplot(212)
    axes2.margins(x=0)
    axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes2.step(wavelength, object_flux, 'black', label='spectrum')
    axes2.step(wavelength, model_flux, 'r',linestyle='dashed',
               label='Model (0 km s$^{-1}$)')
    axes2.step(wavelength*(1+redshift), model_flux, 'r',
               label='Model (%+.2f' % velocity+' km s$^{-1}$)')  
    axes2.set(xlabel=plot_xlabel, ylabel='Relative Flux Density')

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())    



    axes2.legend(bbox_to_anchor=(0.5,0.1),
                 ncols=1,
#                 mode='expand',
                 frameon=False,
                 loc='lower left',
                 handlelength=1)


    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(1.5)

    
        
    
    #
    # Get the plot number and return the results
    #
    
    plot_number = pl.gcf().number
    return plot_number
