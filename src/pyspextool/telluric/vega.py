import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import scipy
import os

from pyspextool.io.check import check_parameter
from pyspextool.utils.arrays import find_index
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.telluric import config as telluric
from pyspextool.plot.limits import get_spec_range

import typing

def vega_xcorrelate(object_wavelength:npt.ArrayLike, object_flux:npt.ArrayLike,
                    vega_wavelength:npt.ArrayLike, vega_flux:npt.ArrayLike,
                    vega_continuum:npt.ArrayLike,
                    minimum_wavelength:typing.Optional[float]=None,
                    maximum_wavelength:typing.Optional[float]=None,
                    resolving_power:typing.Optional[float]=None,
                    qa_show:bool=False, qa_show_plotsize:tuple=(6,10),
                    qa_fileinfo:typing.Optional[dict]=None):

    """
    To determine the velocity shift between an A0 V star and a Vega model

    Parameters
    ----------

    object_wavelength : ndarray
        A (nwave1,) array of wavelengths for the object in microns.

    object_flux : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the object.

    vega_wavelength : ndarray
        A (nwave2,) array of wavelengths for a Vega model in microns.

    vega_flux : ndarray
        A (nwave2,) array of flux densities for a Vega model.

    vega_continuum : ndarray
        A (nwave2,) array of continuum flux densities for a Vega model.

    minimum_wavelength : int or float, default=None
        The minimum wavelength value over which to perform the xcorrelation.

    maximum_wavelength : int or float, default=None
        The maximum wavelength value over which to perform the xcorrelation.

    resolving_power : int or float, default=None
        The resolving power of the observations.  A not-perfect attempt to 
        convolve the Vega model to the same resolving power is then made.

    qa_show : {False, True}
        Set to True to show a QA plot on the screen.

    qa_show_plotsize : tuple of int, default=(6,10)
        The plot size in inches for a screen QA plot.

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


    Returns
    -------
    ndarray
        The velocity shift of the Vega model relative to the data in km s-1.  

    Notes
    -----
    While in principle any Vega model can be passed, the code was designed 
    to work with a Kurucz Vega model.
        

    """

    #
    # Check Parameters
    #

    check_parameter('vega_xcorrelate', 'object_wavelength', object_wavelength,
                    'ndarray', 1)

    check_parameter('vega_xcorrelate', 'object_flux', object_flux, 'ndarray',
                    1)

    check_parameter('vega_xcorrelate', 'vega_wavelength', vega_wavelength,
                    'ndarray', 1)

    check_parameter('vega_xcorrelate', 'vega_flux', vega_flux, 'ndarray',
                    1)

    check_parameter('vega_xcorrelate', 'vega_continuum', vega_continuum,
                    'ndarray', 1)

    check_parameter('vega_xcorrelate', 'minimum_wavelength',
                    minimum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('vega_xcorrelate', 'maximum_wavelength',
                    maximum_wavelength, ['NoneType', 'int', 'float'])

    check_parameter('vega_xcorrelate', 'resolving_power',
                    resolving_power, ['NoneType', 'int', 'float'])        

    check_parameter('vega_xcorrelate', 'qa_show', qa_show, 'bool')

    check_parameter('vega_xcorrelate', 'qa_show_plotsize', qa_show_plotsize,
                    'tuple')

    check_parameter('vega_xcorrelate', 'qa_fileinfo', qa_fileinfo,
                    ['NoneType', 'dict'])    

    
    #
    # Get set up
    #

    cspeed = 2.99792458E5 # km/s

    # Normalize the Vega model

    vega_flux_normalized = vega_flux/vega_continuum - 1

    # Subtract the continuum from the data

    object_flux -= 1

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

    z = np.logical_and((vega_wavelength > minimum_wavelength),
                       (vega_wavelength < maximum_wavelength))

    nvega = np.sum(z)

    if nobject > nvega:

        message = 'Data has a higher resolution than the Vega model.'
        raise ValueError(message)

    npixels_lnlambda = np.max([nobject, nvega])

    #
    # Smooth the Vega model to get close
    #

    if resolving_power is not None:

        # Determine the wavelength mid point

        wavelength = (minimum_wavelength+maximum_wavelength)/2.

        # Find the dispersion of the Vega model at this wavelength
        
        idx = int(find_index(vega_wavelength, wavelength))

        vega_dispersion = vega_wavelength[idx]-vega_wavelength[idx-1] 
        
        # Determine the number of pixels for the resolving power.

        fwhm_kernel = int(wavelength/resolving_power/vega_dispersion)

        # Create a Gaussian kernel

        npixels_kernel = fwhm_kernel*5
        if npixels_kernel % 2 == 0:

            npixels_kernel +=1

        x = np.arange(npixels_kernel)-npixels_kernel//2
        
        sigma = fwhm_kernel/2.354

        gaussian = np.exp(-(x/sigma)**2 / 2)
        gaussian /= np.sum(gaussian)

        vega_flux_normalized = np.convolve(vega_flux_normalized,gaussian,
                                           mode='same')
        
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

    f = scipy.interpolate.interp1d(object_wavelength, object_flux)
    object_resampled_flux = f(lnlambda_wavelengths)

    f = scipy.interpolate.interp1d(vega_wavelength, vega_flux_normalized)
    vega_resampled_flux = f(lnlambda_wavelengths)    

    #
    # Do the cross correlation
    #

    xcor = scipy.signal.correlate(object_resampled_flux, vega_resampled_flux,
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
    # Do the QA plotting
    #

#    if qa_show is True:
#
#        pl.ion()
#        number = make_xcorrelate_plot(lnlambda_wavelengths,
#                                      object_resampled_flux,
#                                      vega_resampled_flux,
#                                      xcor, lag, fit, offset_pixels,
#                                      velocity_shift, redshift,
#                                      qa_show_plotsize,
#                            plot_number=telluric.state['xcorrelate_plotnum'])
#
#
#        telluric.state['xcorrelate_plotnum'] = number
#        pl.show()
#        pl.pause(1)
#
#    if qa_fileinfo is not None:
# 
#        pl.ioff()
#
#        make_xcorrelate_plot(lnlambda_wavelengths, object_resampled_flux,
#                             vega_resampled_flux, xcor, lag, fit, offset_pixels,
#                             velocity_shift, redshift, qa_show_plotsize)
#
#        pl.savefig(os.path.join(qa_fileinfo['filepath'],
#                                qa_fileinfo['filename']) + \
#                   '_xcorrelate' + qa_fileinfo['extension'])
#        pl.close()

    #
    # Return the velocity shift
    #
            
    return velocity_shift

#def make_xcorrelate_plot(wavelengths, object_flux, vega_flux, xcor, lag,
#                         fit, offset, velocity, redshift, plot_size,
#                         plot_number=None):
#
#    """
#    To create a plot for the cross correlation device independently
#
#
#    """
#
#    #
#    # Make the figure
#    #
#    
#    fig = pl.figure(num=plot_number, figsize=plot_size)
#
#    # Create the spectral plot
#
#    yrange = get_spec_range(object_flux, vega_flux, frac=0.1)
#    
#    axes1 = fig.add_subplot(211)
#    axes1.margins(x=0)
#    axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
#    axes1.step(wavelengths, object_flux, '#1f77b4')
#    axes1.step(wavelengths, vega_flux, 'r')
#    axes1.set(xlabel='Wavelength ($\mu$m)', ylabel='Relative Intensity')
#
#    axes1.text(0.95, 0.2, 'data', color='#1f77b4', ha='right',
#                transform=axes1.transAxes)
#
#    axes1.text(0.95, 0.1, 'Vega', color='r', ha='right',
#              transform=axes1.transAxes)
#    
#    # Create the cross correlation plot
#
#    yrange = get_spec_range(xcor, fit['fit'], frac=0.1)
#    
#    axes2 = fig.add_subplot(212)    
#    axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
#    axes2.step(lag, xcor, '#1f77b4')
#    axes2.step(lag, fit['fit'], 'r')
#    axes2.set(xlabel='lag (pixels)', ylabel='Relative Intensity')
#
#    axes2.text(0.05, 0.9, 'X Correlation', color='#1f77b4', ha='left',
#                transform=axes2.transAxes)
#
#    axes2.text(0.05, 0.85, 'Fit', color='r', ha='left',
#              transform=axes2.transAxes)
#
#    axes2.text(0.95, 0.9, 'Offset='+'%.2f' % offset+' pixels',
#               color='black', ha='right', transform=axes2.transAxes)
#
#    axes2.text(0.95, 0.85, 'Velocity='+'%.2f' % velocity+' km s$^{-1}$',
#               color='black', ha='right', transform=axes2.transAxes)    
#    
#    #
#    # Get the plot number and return the results
#    #
#    
#    plot_number = pl.gcf().number
#    return plot_number
    
