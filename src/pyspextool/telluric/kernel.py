import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl

import logging
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.special import erf
import typing

from pyspextool.utils.for_print import for_print
from pyspextool.io.check import check_parameter
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.utils.math import moments


def deconvolve_line(data_wavelength:npt.ArrayLike,
                    data_normalized_flux:npt.ArrayLike,
                    model_wavelength:npt.ArrayLike,
                    model_normalized_flux:npt.ArrayLike,
                    line_wavelength_range:npt.ArrayLike,
                    tapered_window_factor:int|float=10,
                    verbose:bool=False, qa_show:bool=False,
                    qa_show_plotsize:tuple=(6,4),
                    qa_fileinfo:typing.Optional[dict]=None,
                    qa_wavelength_units:str='$\mu$m'):

    """
    To deconvolve an absorption line using a model spetrum.

    Parameters
    ----------

    data_wavelength : ndarray
        A (nwave1,) array of wavelengths for the object.

    data_normalized_flux : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the object

    model_wavelength : ndarray
        A (nwave2,) array of wavelengths for the model. 

    model_normalized_flux : ndarray
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

    check_parameter('deconvolve_line', 'data_wavelength',
                    data_wavelength, 'ndarray', 1)

    check_parameter('deconvolve_line', 'data_normalized_flux',
                    data_normalized_flux, 'ndarray', 1)    

    check_parameter('deconvolve_line', 'model_wavelength',
                    model_wavelength, 'ndarray', 1)

    check_parameter('deconvolve_line', 'model_normalized_flux',
                    model_normalized_flux, 'ndarray', 1)    

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
    data_normalized_line_flux = data_normalized_flux[zdata]        
    model_line_wavelength = model_wavelength[zmodel]
    model_normalized_line_flux = model_normalized_flux[zmodel]    
    
    ndata = np.size(data_line_wavelength)
    nmodel = np.size(model_line_wavelength)

    #
    # Fit absorption line
    #
    
    r = fit_peak1d(model_line_wavelength, model_normalized_line_flux, nparms=4,
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

    data_ew = data_step*(np.sum(1.0-data_normalized_line_flux))
    model_ew = model_step*(np.sum(1.0-model_normalized_line_flux))    
    scale_ew = data_ew/model_ew

    logging.info(' Model Line EW = '+str(model_ew))
    logging.info(' Data Line EW = '+str(data_ew))
    logging.info(' EW Scale Factor= '+str(scale_ew))            
    
    #
    # Zero the spectra and scale the model
    #

    data_zeroed_line_flux = data_normalized_line_flux-1
    model_zeroed_scaled_line_flux = (model_normalized_line_flux-1)*scale_ew
    
    #
    # Interpolate the data to the model wavelength grid.  I set fill_value
    # To match the default behavior or interpol in IDL.
    #

    f = interp1d(data_line_wavelength, data_zeroed_line_flux,
                 bounds_error=False, fill_value="extrapolate")
    rdata_zeroed_line_flux = f(model_line_wavelength)

    #
    # Perform the deconvolution
    #

    fft_rdata = fft(rdata_zeroed_line_flux)
    fft_model = fft(model_zeroed_scaled_line_flux)
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
    

    input = (model_normalized_flux[zconvolve]-1)*scale_ew
    model_zeroed_scaled_convolved_flux = np.convolve(input, kernel,mode='same')

    # Interpolate the convolved model onto data_line_wavelength
    
    f = interp1d(model_wavelength[zconvolve],
                 model_zeroed_scaled_convolved_flux,
                 bounds_error=False, fill_value="extrapolate")
    rmodel_zeroed_scaled_convolved_line_flux = f(data_line_wavelength)
    
    # Compute the statistics

    ratio = (data_zeroed_line_flux+1)/\
        (rmodel_zeroed_scaled_convolved_line_flux+1)
    maximum_deviation = np.max(np.abs(ratio-1))

    results = moments(ratio)
    rms_deviation = results['stddev']

    #
    # Get things ready to return in the dictionary
    #

    # Interpolate the model onto data_line_wavelength
    
    f = interp1d(model_wavelength, model_normalized_flux-1,
                 bounds_error=False, fill_value="extrapolate")
    rmodel_zeroed_line_flux = f(data_line_wavelength)

    
    # Resample kernel onto data wavelengths

    f = interp1d(model_line_wavelength ,kernel, bounds_error=False,
                 fill_value=0.0)
    rkernel = f(data_line_wavelength)

    dict = {'kernel':rkernel, 'wavelengths':data_line_wavelength,
            'data':data_zeroed_line_flux, 'model':rmodel_zeroed_line_flux,
            'convolved_model':rmodel_zeroed_scaled_convolved_line_flux,
            'residuals':ratio-1,'maximum_deviation':maximum_deviation,
            'rms_deviation':rms_deviation}

    #
    # Create the QA plot
    #

    # Create the resampled Vega model just for show

    
    
    
    fig = pl.figure()
    axes1 = fig.add_subplot(111)

    xrange = [np.min(data_line_wavelength),np.max(data_line_wavelength)]
    pl.xlim(xrange)
    pl.ylim([-0.4,0.05])
    
    axes1.step(data_line_wavelength, data_zeroed_line_flux,color='black',
               where='mid')
    axes1.step(data_line_wavelength, rmodel_zeroed_scaled_convolved_line_flux,
               color='green', where='mid')
    axes1.step(data_line_wavelength,ratio-1,color='blue', where='mid')
    axes1.axhline(y=0.01,color='magenta',linestyle='dotted')
    axes1.axhline(y=-0.01,color='magenta',linestyle='dotted')   

    axes1.set(xlabel='Wavelength ('+qa_wavelength_units+')',
              ylabel='Residual Intensity')
    
    axes1.step(data_line_wavelength, rmodel_zeroed_line_flux,color='red',
               where='mid')

#    axes1.step(model_wavelength, model_normalized_flux-1,color='red',
#               where='mid')


    
    axes1.text(0.02, 0.2, 'Max Deviation = '+"{:.4f}".format(maximum_deviation),
               ha='left', va='bottom', transform=axes1.transAxes, color='black')

    axes1.text(0.02, 0.15, 'RMS Deviation = '+"{:.4f}".format(rms_deviation),
               ha='left', va='bottom', transform=axes1.transAxes, color='black')


    axes1.text(0.95, 0.3, 'A0 V', color='black', ha='right', va='bottom',
               transform=axes1.transAxes)

    axes1.text(0.95, 0.25, 'Vega', color='red', ha='right', va='bottom',
               transform=axes1.transAxes)

    axes1.text(0.95, 0.2, 'Scaled & Convolved Vega', color='green', ha='right',
               va='bottom', transform=axes1.transAxes)

    axes1.text(0.95, 0.15, 'Residuals', color='blue', ha='right',
               va='bottom', transform=axes1.transAxes)
    
    pl.show()
    
    
    logging.info(' Maximum deviation '+str(maximum_deviation))


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
    
