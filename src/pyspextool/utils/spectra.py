import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import typing
import os
import scipy

from pyspextool.io.check import check_parameter
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.plot.limits import get_spec_range
from pyspextool.utils.arrays import find_index

def model_xcorrelate(object_wavelength:npt.ArrayLike,
                     object_nflux:npt.ArrayLike,
                     model_wavelength:npt.ArrayLike,
                     model_nflux:npt.ArrayLike,
                     minimum_wavelength:typing.Optional[float]=None,
                     maximum_wavelength:typing.Optional[float]=None,
                     resolving_power:typing.Optional[float]=None,
                     qa_show:bool=False, qa_show_plotsize:tuple=(6,10),
                     qa_fileinfo:typing.Optional[dict]=None,
                     block:bool=False):                    

    """
    To determine the velocity shift between a model and data using a line

    Parameters
    ----------

    object_wavelength : ndarray
        A (nwave1,) array of wavelengths for the object in microns.

    object_nflux : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the object.

    model_wavelength : ndarray
        A (nwave2,) array of wavelengths for a model in microns.

    model_nflux : ndarray
        A (nwave2,) array of continuum-normalized flux densities for a model.

    minimum_wavelength : int or float, default=None
        The minimum wavelength value over which to perform the xcorrelation.

    maximum_wavelength : int or float, default=None
        The maximum wavelength value over which to perform the xcorrelation.

    resolving_power : int or float, default=None
        The resolving power of the observations.  A not-perfect attempt to 
        convolve the model to the same resolving power is then made.

    qa_show : {False, True}
        If set to True, a quality assurance plot will be shown to the user.

    qa_show_plotsize : tuple, default=(6,10)
        if qa_show is True, then the quality assurance plot will be created
        with this size.

    qa_fileinfo : dict, optional

        `'figsize'` : tuple
            A (2,) tuple giving the figure size.

        `'filepath'` : str
            A string giving the write path.
    
        `'filename'` : str
            The root of the file.

        `'extension'` : str
            The file extension.

        The result is os.path.join(qa_fileinfo['filepath'],
                                   qa_fileinfo['filename'] + \
                                   qa_fileinfo['extension'])
    
    block : {False, True}, optional
        Set to make the plot block access to the command line, e.g. pl.ioff().
    
    Returns
    -------
    ndarray
        The velocity shift of the model relative to the data in km s-1.  

        
    """

    #
    # Check Parameters
    #

    check_parameter('model_xcorrelate', 'object_wavelength', object_wavelength,
                    'ndarray', 1)

    check_parameter('model_xcorrelate', 'object_nflux', object_nflux, 'ndarray',
                    1)

    check_parameter('model_xcorrelate', 'model_wavelength', model_wavelength,
                    'ndarray', 1)

    check_parameter('model_xcorrelate', 'model_nflux', model_nflux, 'ndarray',
                    1)

    check_parameter('model_xcorrelate', 'minimum_wavelength',
                    minimum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('model_xcorrelate', 'maximum_wavelength',
                    maximum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('model_xcorrelate', 'resolving_power',
                    resolving_power, ['NoneType', 'int', 'float'])        

    check_parameter('model_xcorrelate', 'qa_show', qa_show, 'bool')

    check_parameter('model_xcorrelate', 'qa_fileinfo', qa_fileinfo,
                    ['NoneType','dict'])        

    check_parameter('model_xcorrelate', 'block', block, 'bool')
    
    
    #
    # Get set up
    #

    cspeed = 2.99792458E5 # km/s

    # subtract the continuum

    model_zflux = model_nflux - 1
    object_zflux = object_nflux-1

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
        raise ValueError(message)

    #
    # Now check to make sure the object has fewer pixels than the model
    #

    nobject = np.size(object_wavelength)

    z = np.logical_and((model_wavelength > minimum_wavelength),
                       (model_wavelength < maximum_wavelength))

    nmodel = np.sum(z)

    if nobject > nmodel:

        message = 'Data has a higher resolution than the model.'
        raise ValueError(message)

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


    if qa_show is True:


        if block is True:

            pl.ioff()

        else:

            pl.ion()
    
        number = plot_model_xcorrelate(lnlambda_wavelengths,
                                       object_resampled_zflux,
                                       model_resampled_zflux,
                                       lag, xcor, fit, offset_pixels,
                                       velocity_shift, redshift,
                                       qa_show_plotsize)


#        telluric.state['xcorrelate_plotnum'] = number
        pl.show()


    if qa_fileinfo is not None:

        pl.ioff()
        plot_model_xcorrelate(lnlambda_wavelengths,
                              object_resampled_zflux,
                              model_resampled_zflux,
                              lag, xcor, fit, offset_pixels,
                              velocity_shift, redshift,
                              qa_fileinfo['figsize'])

        pl.savefig(os.path.join(qa_fileinfo['filepath'],
                                qa_fileinfo['filename'] + \
                                qa_fileinfo['extension']))
        pl.close()
    
    return velocity_shift, redshift
    

    
    
def normalize_continuum(wavelength:npt.ArrayLike, flux:npt.ArrayLike,
                        ranges:npt.ArrayLike, degree:int,
                        robust:typing.Optional[dict]=None,
                        qa_show:bool=False, qa_show_scale:float=1.0,
                        qa_fileinfo:typing.Optional[dict]=None,
                        latex_xlabel:typing.Optional[str]=None,
                        qa_block:bool=True):

    
    """
    To normalize the continuum of a spectrum using a robust polynomial

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths.

    flux : ndarray
        A (nwave,) array of "intensities".

    ranges : ndarray
        A (2*nrange,) list of wavelength ranges over which to fit the contimuum.

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

    qa_show : {False, True}
        Set to True to show a QA plot on the screen.

    qa_show_scale : float or int, default=1.0
        The scale factor for the QA plot size.  The baseline size is (10, 6).

    qa_fileinfo : dict, optional
        `"figsize"` : tuple
            (2,) tuple of the figure size (inches).

        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

    latex_xlabel : str, default='Wavelength'
        A latex string giving the xlabel.

    qa_block : {False, True}, optional
        Set to make the plot block access to the command line, e.g. pl.ioff().
    
    Returns
    -------

    Notes
    -----
    TODO:  Must include proper error propagation.  


    """

    #
    # Check the parameters
    #

    check_parameter('normalize_continuum', 'wavelength', wavelength,
                    'ndarray', 1)

    check_parameter('normalize_continuum', 'flux', flux, 'ndarray', 1)

    check_parameter('normalize_continuum', 'ranges', ranges, 'ndarray', 1)

    check_parameter('normalize_continuum', 'degree', degree, 'int')

    check_parameter('normalize_continuum', 'robust', robust,
                    ['NoneType', 'dict'])

    check_parameter('normalize_continuum', 'qa_show', qa_show, 'bool')

    check_parameter('normalize_continuum', 'qa_show_scale', qa_show_scale,
                    ['int','float'])    

    check_parameter('normalize_continuum', 'qa_fileinfo', qa_fileinfo,
                    ['NoneType', 'dict'])

    check_parameter('normalize_continuum', 'latex_xlabel', latex_xlabel, 'str')

    check_parameter('normalize_continuum', 'block', block, 'bool')    

    
    #
    # Determine if the wavelengths are monotonically increasing
    #

    if np.all(np.diff(ranges) > 0) == False:

        message = "`ranges` is not monotonic."
        raise ValueError(message)

    #
    # Make sure all wavelengths fall within the input range of `wavelength`.
    #

    min_wavelength_spectrum = np.nanmin(wavelength)

    max_wavelength_spectrum = np.nanmax(wavelength)

    if np.all(ranges > min_wavelength_spectrum) == False:

        message = "`ranges` is out of range."
        raise ValueError(message)

    if np.all(ranges < max_wavelength_spectrum) == False:

        message = "`ranges` is out of range."
        raise ValueError(message)    
            
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

    result = poly_fit_1d(wavelength[zregions], flux[zregions], degree)

    continuum = poly_1d(wavelength,result['coeffs'])
    
    flux_normalized = flux/continuum

    #
    # Make the QA plot
    #
    
    plot_number = None
    if qa_show is True:

        # This is to the screen

        if block is True:

            pl.ioff()

        else:

            pl.ion()

        plot_number = plot_normalization(wavelength,
                                         flux,
                                         continuum,
                                         zregions,
                                         (10*qa_show_scale, 6*qa_show_scale),
                                         latex_xlabel=latex_xlabel)
        pl.show()
        pl.pause(1)

    if qa_fileinfo is not None:

        pl.ioff()
        plot_normalization(wavelength,
                           flux,
                           continuum,
                           zregions,
                           qa_fileinfo['figsize'],
                           latex_xlabel=latex_xlabel)

        pl.savefig(os.path.join(qa_fileinfo['filepath'],
                                qa_fileinfo['filename'] + \
                                qa_fileinfo['extension']))
        pl.close()

            
    return flux_normalized, plot_number



def plot_normalization(wavelength:npt.ArrayLike,
                       intensity:npt.ArrayLike,
                       continuum:npt.ArrayLike,
                       zregions:npt.ArrayLike,
                       figure_size:tuple,
                       plot_number:typing.Optional[float]=None,
                       latex_xlabel:str='Wavelength',
                       file_info:typing.Optional[str]=None,
                       block:bool=False):

    
    """
    To plot the results of normalize_continuum

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths.

    intensity : ndarray
        A (nwave,) array of "intensities".

    continuum : ndarray
        A (nwave,) array of the fitted continuum values.

    zranges : ndarray
        A (nwave,) True/False array where pixels used in the fit are True.

    Returns
    -------
    Nothing

    """

    #
    # Make the two-panel figure
    #
    
    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.clf()
    
    # Get the plot range for x axis
    
    xrange = get_spec_range(wavelength[zregions], frac=0.1)

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
    axes1.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes1.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes1.set_ylabel('Intensity')    

    # Now plot the fitted pixels in red
    
    tmp = np.copy(intensity)
    tmp[~zregions] = np.nan
    
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
    axes2.set_xlabel(latex_xlabel)
    axes2.set_ylabel('Normalized Intensity')    
    
    # Plot the fitted pixels in red
    
    tmp = np.copy(normalized)
    tmp[~zregions] = np.nan
    
    axes2.step(wavelength, tmp, 'red')
    axes2.axhline(y=1, linestyle='--', color='green')




def plot_model_xcorrelate(wavelength, object_flux, model_flux, lag,
                          xcorrelation, fit, offset, velocity, redshift,
                          plot_size):


    """
    To create a plot for the cross correlation device independently

    Parameters
    ----------
    
    

    """

    #
    # Make the figure
    #
    
    fig = pl.figure(figsize=plot_size)

    # Create the spectral plot

    yrange = get_spec_range(object_flux, model_flux, frac=0.1)
    
    axes1 = fig.add_subplot(211)
    axes1.margins(x=0)
    axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes1.step(wavelength, object_flux, '#1f77b4')
    axes1.step(wavelength, model_flux, 'r')
    axes1.set(xlabel='Wavelength ($\mu$m)', ylabel='Relative Intensity')

    axes1.text(0.95, 0.2, 'data', color='#1f77b4', ha='right',
                transform=axes1.transAxes)

    axes1.text(0.95, 0.1, 'Model', color='r', ha='right',
              transform=axes1.transAxes)
    
    # Create the cross correlation plot

    yrange = get_spec_range(xcorrelation, fit['fit'], frac=0.1)
    
    axes2 = fig.add_subplot(212)    
    axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes2.step(lag, xcorrelation, '#1f77b4')
    axes2.step(lag, fit['fit'], 'r')
    axes2.set(xlabel='lag (pixels)', ylabel='Relative Intensity')

    axes2.text(0.05, 0.9, 'X Correlation', color='#1f77b4', ha='left',
                transform=axes2.transAxes)

    axes2.text(0.05, 0.85, 'Fit', color='r', ha='left',
              transform=axes2.transAxes)

    axes2.text(0.95, 0.9, 'Offset='+'%.2f' % offset+' pixels',
               color='black', ha='right', transform=axes2.transAxes)

    axes2.text(0.95, 0.85, 'Velocity='+'%.2f' % velocity+' km s$^{-1}$',
               color='black', ha='right', transform=axes2.transAxes)    
    
    #
    # Get the plot number and return the results
    #
    
    plot_number = pl.gcf().number
    return plot_number
