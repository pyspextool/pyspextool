import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.ticker import (AutoMinorLocator)

import logging
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.special import erf
import typing
import os

from pyspextool.utils.for_print import for_print
from pyspextool.io.check import check_parameter
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.utils.math import moments
from pyspextool.plot.limits import get_spec_range


def deconvolve_line(data_wavelength:npt.ArrayLike,
                    data_fluxdensity:npt.ArrayLike,
                    model_wavelength:npt.ArrayLike,
                    model_fluxdensity:npt.ArrayLike,
                    line_wavelength_range:npt.ArrayLike,
                    tapered_window_factor:int|float=10,
                    verbose:bool=False,
                    qashow_info:typing.Optional[dict]=None,
                    qafile_info:typing.Optional[dict]=None):

    """
    To deconvolve an absorption line using a model spetrum.

    Parameters
    ----------

    data_wavelength : ndarray
        A (nwave1,) array of wavelengths for the object.

    data_fluxdensity : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the object

    model_wavelength : ndarray
        A (nwave2,) array of wavelengths for the model. 

    model_fluxdensity : ndarray
        A (nwave12,) array of continuum-normalized flux densities for the model

    line_wavelength_range : ndarray
        A (2,) array of of wavelengths giving the range over which to do the 
        deconvolution.
    

    Returns
    -------

    Notes
    ----- 
    The units of `data_wavelength`, `model_wavelength`, and 
    `line_wavelength_range` must be the same.

    """
    

    #
    # Check Parameters
    #

    check_parameter('deconvolve_line', 'data_wavelength', data_wavelength,
                    'ndarray', 1)

    check_parameter('deconvolve_line', 'data_fluxdensity', data_fluxdensity,
                    'ndarray', 1)    

    check_parameter('deconvolve_line', 'model_wavelength', model_wavelength,
                    'ndarray', 1)

    check_parameter('deconvolve_line', 'model_fluxdensity', model_fluxdensity,
                    'ndarray', 1)    

    check_parameter('deconvolve_line', 'line_wavelength_range',
                    line_wavelength_range, 'ndarray', 1)    

    #
    # Check verbose
    #

    if verbose is True:
        logging.getLogger().setLevel(logging.INFO)
        
    elif verbose is False:
        logging.getLogger().setLevel(logging.ERROR)


    #
    # Do basic set up stuff
    #

    # Get the wavelength range of the data, which we will assume is always
    # smaller than the model


    data_wavelength_range = [np.min(data_wavelength),np.max(data_wavelength)]
    
    # Locate the deconvolution region for the object and Vega

    zdata = np.where((data_wavelength >= line_wavelength_range[0]) &
                     (data_wavelength <= line_wavelength_range[1]))[0]
    
    zmodel = np.where((model_wavelength >= line_wavelength_range[0]) & 
                      (model_wavelength <= line_wavelength_range[1]))[0]
    
    # Require the number of points to be even

    if np.size(zdata) % 2 == 1:

        zdata = np.append(zdata, zdata[-1]+1)

    if np.size(zmodel) % 2 == 1:

        zmodel = np.append(zmodel, zmodel[-1]+1)

        
    # Clip the line out of the data

    data_line_wavelength = data_wavelength[zdata]
    data_line_fluxdensity = data_fluxdensity[zdata]
    model_line_wavelength = model_wavelength[zmodel]
    model_line_fluxdensity = model_fluxdensity[zmodel]    
    
    ndata = np.size(data_line_wavelength)
    nmodel = np.size(model_line_wavelength)

    #
    # Fit absorption line
    #
    
    r = fit_peak1d(model_line_wavelength, model_line_fluxdensity, nparms=4,
                   negative=True)
    line_wavelength = r['parms'][1]

    logging.info(' Absorption line at '+str(line_wavelength))
    
    #
    # Determine the sampling frequency
    #

    data_step = (np.max(data_line_wavelength)-np.min(data_line_wavelength))/\
                  (ndata-1)
    model_step = (np.max(model_line_wavelength)-np.min(model_line_wavelength))/\
                  (nmodel-1)

    logging.info(' Model sampling = '+str(model_step))
    logging.info(' Data sampling = '+str(data_step))

    #
    # Determine EWs of the absorption features and set the scale factor
    #

    data_ew = data_step*(np.sum(1.0-data_line_fluxdensity))
    model_ew = model_step*(np.sum(1.0-model_line_fluxdensity))    
    scale_ew = data_ew/model_ew

    logging.info(' Model Line EW = '+str(model_ew))
    logging.info(' Data Line EW = '+str(data_ew))
    logging.info(' EW Scale Factor= '+str(scale_ew))            
    
    #
    # Zero the spectra and scale the model
    #

    data_zeroed_line_fluxdensity = data_line_fluxdensity-1
    model_zeroed_scaled_line_fluxdensity = (model_line_fluxdensity-1)*\
        scale_ew
    
    #
    # Interpolate the data to the model wavelength grid.  I set fill_value
    # To match the default behavior or interpol in IDL.
    #

    f = interp1d(data_line_wavelength, data_zeroed_line_fluxdensity,
                 bounds_error=False, fill_value="extrapolate")
    rdata_zeroed_line_fluxdensity = f(model_line_wavelength)

    #
    # Perform the deconvolution
    #

    fft_rdata = fft(rdata_zeroed_line_fluxdensity)
    fft_model = fft(model_zeroed_scaled_line_fluxdensity)
    fft_kernel = fft_rdata/fft_model

    #
    # Build the frequency array.  I did not use fftfreq because it doesn't
    # quite give the same value as below (doesn't give the zero at the end).
    #

    dw_model = (np.max(model_line_wavelength)-\
                np.min(model_line_wavelength))/(nmodel-1)
    dw_data = (np.max(model_line_wavelength)-\
               np.min(model_line_wavelength))/(ndata-1)    

    freq = np.arange(nmodel/2)/(nmodel*dw_model)
    frequency = np.append(np.append(freq, -1*np.flip(freq[1:nmodel//2])),0)
    
    #
    # Apply a tapered window function
    #

    r = fit_peak1d(frequency, np.square(fft_rdata.real), nparms=4,
                   positive=True)
    sigma_fft_rdata = r['parms'][2]
    
    r = fit_peak1d(frequency, np.square(fft_model.real), nparms=4,
                   positive=True)
    sigma_fft_model = r['parms'][2]

    sigma_fft = np.min((sigma_fft_rdata, sigma_fft_model))

    FWHM_fft = sigma_fft*2*np.sqrt(2*np.log(2))

    max_frequency = tapered_window_factor * FWHM_fft
        
    tapered_window = 1/(1+np.abs(frequency/max_frequency)**10)

    np.multiply(fft_kernel, tapered_window, out=fft_kernel)

    #
    # Inverse FFT to get kernel and then reorder and normalize
    #

    kernel = ifft(fft_kernel).real    
    kernel = np.roll(kernel, nmodel//2)
    kernel /= np.sum(kernel)

#    f = open("kernel_python.dat","w+")
#    for i in range(len(kernel)):
#        f.write(str(kernel[i])+'\n')

#    f.close()
        
#    return
    
    #
    # Convolve the model over the data_wavelength range
    #
    
    zconvolve = np.where((model_wavelength >= data_wavelength_range[0]) &
                     (model_wavelength <= data_wavelength_range[1]))[0]
    

    input = (model_fluxdensity[zconvolve]-1)*scale_ew
    model_zeroed_scaled_convolved_fluxdensity = np.convolve(input, kernel,
                                                            mode='same')

    # Interpolate the convolved model onto data_line_wavelength
    
    f = interp1d(model_wavelength[zconvolve],
                 model_zeroed_scaled_convolved_fluxdensity,
                 bounds_error=False, fill_value="extrapolate")
    rmodel_zeroed_scaled_convolved_line_fluxdensity = f(data_line_wavelength)
    
    # Compute the statistics

    line_fluxdensity_ratio = (data_zeroed_line_fluxdensity+1)/\
        (rmodel_zeroed_scaled_convolved_line_fluxdensity+1)
    maximum_deviation = np.max(np.abs(line_fluxdensity_ratio-1))

    results = moments(line_fluxdensity_ratio)
    rms_deviation = results['stddev']

    logging.info(' RMS deviation '+str(rms_deviation))        
    logging.info(' Maximum deviation '+str(maximum_deviation))    

    #
    # Get things ready to return in the dictionary
    #

    # Interpolate the model onto data_line_wavelength
    
    f = interp1d(model_wavelength, model_fluxdensity-1, bounds_error=False,
                 fill_value="extrapolate")
    rmodel_zeroed_line_fluxdensity = f(data_line_wavelength)

    
    # Resample kernel onto data wavelengths

    f = interp1d(model_line_wavelength ,kernel, bounds_error=False,
                 fill_value=0.0)
    rkernel = f(data_line_wavelength)

    dict = {'kernel':rkernel, 'wavelengths':data_line_wavelength,
            'data':data_zeroed_line_fluxdensity,
            'model':rmodel_zeroed_line_fluxdensity,
            'convolved_model':rmodel_zeroed_scaled_convolved_line_fluxdensity,
            'residuals':line_fluxdensity_ratio-1,
            'maximum_deviation':maximum_deviation,
            'rms_deviation':rms_deviation}

    #
    # Make the QA plot
    #

    figure_size = (6,10)
    
    if qashow_info is not None:


        
        if qashow_info['block'] is True:

            pl.ioff()

        else:

            pl.ion()

        plotnum = plot_deconvolve_line(data_wavelength,
                                       data_fluxdensity,
                                qashow_info['normalization_wavelength_range'],
                                       line_wavelength_range,
                                       data_line_wavelength,
                                       data_zeroed_line_fluxdensity,
                                       rmodel_zeroed_line_fluxdensity,
                                rmodel_zeroed_scaled_convolved_line_fluxdensity,
                                       line_fluxdensity_ratio,
                                       rms_deviation,
                                       maximum_deviation,
                                       plot_scale=qashow_info['figure_scale'],
                                       plot_xlabel=qashow_info['plot_xlabel'],
                                       plot_title=qashow_info['plot_title'])
                                

            
        pl.show()
        if qashow_info['block'] is False:
            pl.pause(1)

        
    if qafile_info is not None:

        plot_deconvolve_line(data_wavelength,
                             data_fluxdensity,
                             qashow_info['normalization_wavelength_range'],
                             line_wavelength_range,
                             data_line_wavelength,
                             data_zeroed_line_fluxdensity,
                             rmodel_zeroed_line_fluxdensity,
                             rmodel_zeroed_scaled_convolved_line_fluxdensity,
                             line_fluxdensity_ratio,
                             rms_deviation,
                             maximum_deviation,
                             plot_number=qashow_info['plot_number'],
                             plot_xlabel=qashow_info['plot_xlabel'],
                             plot_title=qashow_info['plot_title'])
        

        pl.savefig(os.path.join(qafile_info['filepath'],
                                qafile_info['filename'] + \
                                qafile_info['extension']))
        pl.close()


        plotnum = None

    return plotnum


    
def make_instrument_profile(x:npt.ArrayLike, parameters:npt.ArrayLike):

    """
    To create a pySpextool instrument profile.

    Parameters
    ----------
    x : ndarray
        An (ndata,) array of x

    parameters : ndarray
        An (3,) array of parameters.

    Returns
    -------
    ndarray
        An (ndata,) array of values that give the instrument profile at x.

    Notes
    -----

    The instrument profile P is given by,

    P(x) = C * {  erf[ (x + parameter[1] - parameter[0])/parameter[2] -
                  erf[ (x - parameter[1] - parameter[0])/parameter[2] }

    where parameter[0] is typically zero for the profile to be centered in x
    and C is determined by normalizing the profile by its sum.  
    
    """

    #
    # Check Parameters
    #

    check_parameter('make_instrument_profile', 'x', x, 'ndarray', 1)

    check_parameter('make_instrument_profile', 'parameters', parameters,
                    'ndarray', 1)    

    #
    # Just do it
    #

    ip = erf( (x+parameters[1]-parameters[0])/parameters[2]) - \
         erf( (x-parameters[1]-parameters[0])/parameters[2])

    ip /= np.sum(ip)

    return ip



def plot_deconvolve_line(order_wavelength:npt.ArrayLike,
                         order_fluxdensity:npt.ArrayLike,
                         normalization_wavelength_range:npt.ArrayLike,
                         line_wavelength_range:npt.ArrayLike,
                         line_wavelength:npt.ArrayLike,
                         line_data_fluxdensity:npt.ArrayLike,
                         line_model_fluxdensity:npt.ArrayLike,
                         line_model_convolved_fluxdensity:npt.ArrayLike,
                         line_fluxdensity_ratio:npt.ArrayLike,
                         rms_deviation:float,
                         maximum_deviation:float,
                         plot_scale:typing.Optional[int]=1.0,
                         plot_number:typing.Optional[int]=None,
                         plot_xlabel:typing.Optional[str]=None,
                         plot_title:typing.Optional[str]=None):


    """
    To plot the results of the deconvolution.

    Creates a vertical two-panel plot.  The upper panel shows the normalized
    spectrum and the region over which the deconvolution was done and the
    lower plot shows the results of the deconvolution process.

    Parameters
    ----------
    order_wavelength : ndarray
        A (nwave1,) array of wavelengths.

    order_fluxdensity : ndarray
        A (nwave1,) array of data continuum-normalized flux densities.
    
    normalization_wavelength_range : ndarray
        A (2,) array giving the wavelength range over which the normalization
        of the continuum was done.

    line_wavelength_range : ndarray
        A (2,) array giving the wavelength range over which the line was
    deconvolved.

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
    int
    The plot number

    """

    #
    # Check the parameters
    #
    
    check_parameter('plot_deconvolve_line', 'order_wavelength',
                    order_wavelength, 'ndarray', 1)

    check_parameter('plot_deconvolve_line', 'order_fluxdensity',
                    order_fluxdensity, 'ndarray', 1)

    check_parameter('plot_deconvolve_line',
                    'normalization_wavelength_range',
                    normalization_wavelength_range, ['ndarray','list'])

    check_parameter('plot_deconvolve_line', 'line_wavelength_range',
                    line_wavelength_range, 'ndarray', 1)
    
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

    check_parameter('plot_deconvolve_line', 'plot_scale', plot_scale,
                    ['float','int'])

    check_parameter('plot_deconvolve_line', 'plot_number', plot_number,
                    ['int','NoneType'])
        
    check_parameter('plot_deconvolve_line', 'plot_xlabel', plot_xlabel,
                    ['str','NoneType'])

    check_parameter('plot_deconvolve_line', 'plot_title', plot_title,
                    ['str','NoneType'])

    
    #
    # Set up the plot
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
    
    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.2)

    #
    # Do the top plot
    #
    
    axes1 = fig.add_subplot(211)

    # Get the xrange based on 'line_wavelength_range'
    
    delta = normalization_wavelength_range[1]-normalization_wavelength_range[0]
    xrange = [normalization_wavelength_range[0]-delta*0.1,
              normalization_wavelength_range[1]+delta*0.1]

    pl.xlim(xrange)

    # Get the yrange based on te data in that wavelength range
    
    zleft = (order_wavelength > normalization_wavelength_range[0])
    zright = (order_wavelength < normalization_wavelength_range[1])
        
    zselection = np.logical_and(zleft,zright)

    yrange = get_spec_range(order_fluxdensity[zselection], frac=0.1)    
    
    pl.ylim(yrange)
 
    # Plot the spectrum
    
    axes1.step(order_wavelength, order_fluxdensity,color='black', where='mid')

    # Plot the normalization level and line region
    
    axes1.axvline(x=line_wavelength_range[0],color='red',linestyle='dotted')
    axes1.axvline(x=line_wavelength_range[1],color='red',linestyle='dotted')
    axes1.hlines(y=1,color='green',xmin=normalization_wavelength_range[0],
                  xmax=normalization_wavelength_range[1],linestyle='dashed')

    # Deal with the tickmarks

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())

    # Label axes
    
    axes1.set_title(plot_title)   
    axes1.set_ylabel('Normalized Flux Density')
    
    #
    # Do the second figure showing the convolution results
    #
    
    axes2 = fig.add_subplot(212)

    # Get the xrange 
    
    xrange = [np.min(line_wavelength),np.max(line_wavelength)]
    pl.xlim(xrange)

    # Get the yrange based on all three spectra to be plotted.
    
    yrange = get_spec_range(line_data_fluxdensity, line_model_fluxdensity,
                            line_fluxdensity_ratio-1, frac=0.1)    
    
    pl.ylim(yrange)    

    # Plot the data 
    
    axes2.step(line_wavelength, line_data_fluxdensity,color='black',
               where='mid')
    axes2.step(line_wavelength, line_model_fluxdensity, color='green',
               where='mid')
    axes2.step(line_wavelength,line_fluxdensity_ratio-1,color='blue',
               where='mid')
    axes2.step(line_wavelength, line_model_convolved_fluxdensity,color='red',
               where='mid')

    
    # Plot the rms limit lines of 0.01.
    
    axes2.axhline(y=0.01,color='magenta',linestyle='dotted')
    axes2.axhline(y=-0.01,color='magenta',linestyle='dotted')

    # Deal with the tickmarks

    # Deal with the tickmarks

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Label the axes
    
    axes2.set(xlabel=plot_xlabel, ylabel='Residual Normalized Flux Density')

    # Add all the information texts
        
    axes2.text(0.02, 0.2, 'Max Deviation = '+"{:.4f}".format(maximum_deviation),
               ha='left', va='bottom', transform=axes2.transAxes, color='black')

    axes2.text(0.02, 0.15, 'RMS Deviation = '+"{:.4f}".format(rms_deviation),
               ha='left', va='bottom', transform=axes2.transAxes, color='black')


    axes2.text(0.95, 0.3, 'A0 V', color='black', ha='right', va='bottom',
               transform=axes2.transAxes)

    axes2.text(0.95, 0.25, 'Vega', color='green', ha='right', va='bottom',
               transform=axes2.transAxes)

    axes2.text(0.95, 0.2, 'Scaled & Convolved Vega', color='red', ha='right',
               va='bottom', transform=axes2.transAxes)

    axes2.text(0.95, 0.15, 'Residuals', color='blue', ha='right',
               va='bottom', transform=axes2.transAxes)
    

    #
    # Get the plot number and return the results
    #
    
    plot_number = pl.gcf().number

    return plot_number





    
    
