import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import scipy

import logging
from scipy import signal
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline
from scipy.special import erf
from dust_extinction.parameter_averages import G23
import astropy.units as u

from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.fit.polyfit import polyfit_1d, poly_1d
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.math import moments
from pyspextool.utils.units import convert_fluxdensity
from pyspextool.telluric.qaplots import plot_measure_linerv
from pyspextool.telluric.qaplots import plot_deconvolve_line
from pyspextool.telluric.qaplots import plot_estimate_ewscales
import pyspextool.telluric.qaplots as qaplots


def deconvolve_line(
    data_wavelength:npt.ArrayLike,
    data_fluxdensity:npt.ArrayLike,
    model_wavelength:npt.ArrayLike,
    model_fluxdensity:npt.ArrayLike,
    line_wavelength_range:npt.ArrayLike,
    tapered_window_factor:int|float=10,
    verbose:bool=False,
    qashow_info:dict=None,
    qafile_info:dict=None):

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

    tapered_window_factor : float, default=10

    verbose : {None, True, False}
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file.

    qashow_info : None or dict
        `'number'` : int
            The plot number.  Useful if you are doing to-screen plotting
            because you can plot to the same window multiple times.
    
        `'scale'` : float or int
            A scale factor to increase the size of the plot over the default.

        `'block'`: {False, True}, optional
            Set to make the plot block access to the command line, e.g.
            pl.ioff().

        `'xlabel'` : str, optional
            A latex string giving the xlabel.

        `'title'` : str, optional
            A latex string giving the title of the plot.
        
    qafile_info : dict, optional    
        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

        `'xlabel'` : str, optional
            A latex string giving the xlabel.

        `'title'` : str, optional
            A (latex) string giving the title of the plot.
      
    Returns
    -------
    dict
        `"object_pixels"` : ndarray
            The pixels associated with the kernel in data space

        `"kernel"` : ndarray
            The kernel

        `"plot_number"` : int or NoneType
            The plot number if qashow_info is passed.
            
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

    check_parameter('deconvolve_line', 'tapered_window_factor',
                    tapered_window_factor, ['float','int'])

    check_parameter('deconvolve_line', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('deconvolve_line', 'qashow_info', qashow_info,
                    ['NoneType', 'dict'])
    
    check_parameter('deconvolve_line', 'qafile_info', qafile_info,
                    ['NoneType', 'dict'])

    check_qakeywords(verbose=verbose)
    
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

    if np.size(zdata) % 2 == 0:

        zdata = np.append(zdata, zdata[-1]+1)

    if np.size(zmodel) % 2 == 0:

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
    # Determine the sampling frequency of the data and the model
    #

    data_step = (np.max(data_line_wavelength)-np.min(data_line_wavelength))/\
                  (ndata-1)
    model_step = (np.max(model_line_wavelength)-np.min(model_line_wavelength))/\
                  (nmodel-1)

    logging.info(' Model sampling = '+str(model_step))
    logging.info(' Data sampling = '+str(data_step))

    #
    # Determine EWs of the absorption features in the data and the model and
    # set the EW scale factor
    #

    data_ew = data_step*(np.sum(1.0-data_line_fluxdensity))
    model_ew = model_step*(np.sum(1.0-model_line_fluxdensity))    
    scale_ew = data_ew/model_ew

    logging.info(' Model Line EW = '+str(model_ew))
    logging.info(' Data Line EW = '+str(data_ew))
    logging.info(' EW Scale Factor = '+str(scale_ew))
    
    #
    # Zero the continuum of the data and model and scale the model to adjust
    # its equivalent width to match that of the object
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
    # Build a tapered window function
    #
    
    dw_model = (np.max(model_line_wavelength)-\
                np.min(model_line_wavelength))/(nmodel-1)

    freq = np.arange(nmodel/2)/(nmodel*dw_model)
    frequency = np.append(np.append(freq, -1*np.flip(freq[1:nmodel//2])),0)
    
    
    # Apply a tapered window function to the FFT of the data and model
    

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

       
    # Inverse FFT to get kernel and then reorder and normalize
    
    kernel = ifft(fft_kernel).real    
    kernel = np.roll(kernel, nmodel//2)
    kernel /= np.sum(kernel)

    
    #
    # Build the dictionary for return
    #

    # Get the kernel, which is in integer model pixels into data pixels
        
    model_pixels = np.arange(len(kernel))

    r = fit_peak1d(model_pixels,kernel, nparms=3, positive=True)

    object_relative_pixels = (model_pixels-r['parms'][1])*model_step/data_step

    # Create the dictionary for later return
    
    dict = {'data_pixels':object_relative_pixels,
            'kernel': kernel,
            'ew_scale':scale_ew}
    
    #
    # Now start the plotting process
    #    
    
    # Convolve the model over the data_wavelength range
        
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

    dict['rms_deviation'] = rms_deviation
    dict['max_deviation'] = maximum_deviation    

    logging.info(' RMS deviation = '+str(rms_deviation))        
    logging.info(' Maximum deviation = '+str(maximum_deviation))    

    #
    # Get things ready to return in the dictionary
    #

    # Interpolate the model onto data_line_wavelength
    
    f = interp1d(model_wavelength, model_fluxdensity-1, bounds_error=False,
                 fill_value="extrapolate")
    rmodel_zeroed_line_fluxdensity = f(data_line_wavelength)

    #
    # Make the QA plot
    #

    if qashow_info is not None:

        plot_deconvolve_line(
            qashow_info['plot_number'],
            qashow_info['figure_size'],
            qashow_info['font_size'],
            qashow_info['spectrum_linewidth'],
            qashow_info['spine_linewidth'],
            data_line_wavelength,
            data_zeroed_line_fluxdensity,
            rmodel_zeroed_line_fluxdensity,
            rmodel_zeroed_scaled_convolved_line_fluxdensity,
            line_fluxdensity_ratio,
            rms_deviation,
            maximum_deviation,
            plot_xlabel=qashow_info['xlabel'],
            plot_title=qashow_info['title'])

        pl.show(block=qashow_info['block'])
        if qashow_info['block'] is False:

            pl.pause(1)

    if qafile_info is not None:

        plot_deconvolve_line(
            None,
            qafile_info['figure_size'],
            qafile_info['font_size'],
            qafile_info['spectrum_linewidth'],
            qafile_info['spine_linewidth'],
            data_line_wavelength,
            data_zeroed_line_fluxdensity,
            rmodel_zeroed_line_fluxdensity,
            rmodel_zeroed_scaled_convolved_line_fluxdensity,
            line_fluxdensity_ratio,
            rms_deviation,
            maximum_deviation,
            plot_xlabel=qafile_info['xlabel'],
            plot_title=qafile_info['title'])
        
        
        pl.savefig(qafile_info['file_fullpath'])
        pl.close()

    result = linear_interp1d(
        model_line_wavelength,
        kernel,
        data_line_wavelength)
    
    result = np.nan_to_num(result)
    
    dict['kernel'] = result/np.sum(result)

    return dict




#def estimate_ewscales(
#    standard_wavelength:npt.ArrayLike,
#    standard_fluxdensity:npt.ArrayLike,
#    standard_rv:int |float,
#    standard_vmag:int | float,
#    standard_bmag:int | float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    ew_scale:float | int,
#    line_wavelengths:npt.ArrayLike,
#    atmospheric_transmission:npt.ArrayLike,
#    poly_degree:int,
#    tolerance:float=0.01,
#    pixelshift_model:bool=False,
#    include_edges:bool=False,
#    qashow_info:dict=None,
#    qafile_info:dict=None):
#
#
#    """
#    To estimate the H line EW scale factors .
#
#    Parameters
#    ----------
#    standard_wavelength : ndarray
#        An (nstd,) array of standard star wavelengths.  
#
#    standard_fluxdensity : ndarray
#        An (nstd,) array of standard star flux densities.
#
#    standard_rv : float
#        A float giving the standard star radial velocities in km s-1.
#
#    standard_bmag : float
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float
#        The V-band magnitude of the standard star.
#
#    vega_wavelength : ndarray
#        An (nvega,) array of Vega model wavelengths.  
#
#    vega_fluxdensity : ndarray
#        An (nvega,) array of Vega model flux densities.
#        
#    vega_continuum : ndarray
#        An (nvega,) array of Vega model continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        An (nvega,) array of Vega model fitted continuum flux densities.
#
#    kernel : ndarray
#        An (nkernel,) kernel used to smooth the Vega model.
#
#    ew_scale : float
#        The EW scale factor.
#
#    poly_degree : int
#        The polynomial degree used to model the instrument throughput.
#
#    include_edges : {False, True}
#        Set to True to include the min() and max() values of 
#        `standard_wavelength` as control points.
#        Set to False to not include the min() and max() values of 
#        `standard_wavelength` as control points.
#
#    qashow_info : None or dict
#        `'number'` : int
#            The plot number.  Useful if you are doing to-screen plotting
#            because you can plot to the same window multiple times.
#    
#        `'scale'` : float or int
#            A scale factor to increase the size of the plot over the default.
#
#        `'block'`: {False, True}, optional
#            Set to make the plot block access to the command line, e.g.
#            pl.ioff().
#
#        `'xlabel'` : str, optional
#            A latex string giving the xlabel.
#
#        `'title'` : str, optional
#            A latex string giving the title of the plot.
#        
#    qafile_info : dict, optional    
#        `"filepath"` : str
#            The directory to write the QA figure.
#
#        `"filename"` : str
#            The name of the file, sans suffix/extension.
#
#        `"extension"` : str
#            The file extension.  Must be compatible with the savefig
#            function of matplotlib.
#
#        `'xlabel'` : str, optional
#            A latex string giving the xlabel.
#
#        `'title'` : str, optional
#            A (latex) string giving the title of the plot.
#
#    Returns
#    -------
#
#    """
#
#    #
#    # Check parameters
#    #
#
#    check_parameter('estimate_ewscales', 'standard_wavelength', 
#                    standard_wavelength, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'standard_fluxdensity', 
#                    standard_fluxdensity, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'standard_rv', standard_rv, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'standard_bmag', standard_bmag, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'standard_vmag', standard_vmag, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'vega_wavelength', 
#                    vega_wavelength, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_fluxdensity', 
#                    vega_fluxdensity, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_continuum', 
#                    vega_continuum, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_fitted_continuum', 
#                    vega_fitted_continuum, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'kernel', kernel, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'line_wavlenegths', line_wavelengths, 
#                    ['float', 'ndarray'])
#
#    check_parameter('estimate_ewscales', 'ew_scale', ew_scale, 'float')
#
#    check_parameter('estimate_ewscales', 'poly_degree', poly_degree, 
#                    'int')
#
#    check_parameter('estimate_ewscales', 'tolerance', tolerance, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'include_edges', include_edges, 
#                    'bool')
#
#    check_parameter('adjust_ews', 'qashow_info', qashow_info, 
#                    ['NoneType','dict'])
#
#    check_parameter('adjust_ews', 'qafile_info', qafile_info, 
#                    ['NoneType','dict'])
#
#
#    #
#    # Normalize the data
#    #
#
#    standard_fluxdensity /= np.nanmedian(standard_fluxdensity)
#    
#    #
#    # Create a default scales array
#    #
#
#    if include_edges is True:
#        
#        points = np.append(np.insert(line_wavelengths, 0, 
#                                     np.min(standard_wavelength)),
#                           np.max(standard_wavelength))
#        
#    else:
#
#        points = np.array(line_wavelengths)
#
#    default_scales = np.full_like(points, ew_scale) 
#
#    #
#    # Create the scales array on the Vega wavelength grid
#    #
#
#    f = scipy.interpolate.interp1d(
#        points, default_scales, 
#        fill_value=(default_scales[0],default_scales[-1]),
#        bounds_error=False)
#
#    vega_scales = f(vega_wavelength)
#
#    #
#    # Create default Vega model
#    #
#
#    result = modify_kuruczvega(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        vega_scales,
#        units='erg s-1 cm-2 A-1',
#        new=False)
#
#    # Normalize the vega model
#
#    default_vegamodel = result[0]/np.median(result[0])
#
#    #
#    # Estimate model parameters p0.
#    # p0[0] = atmospheric scale factor
#    # p0[1:poly_degree+2] = coeffs
#    # p0[poly_degree+2:] = EW scale factors
#    #
#
#    default_ratio = standard_fluxdensity/default_vegamodel
#
#
#    coeffs = polyfit_1d(
#        standard_wavelength, 
#        default_ratio, 
#        poly_degree)['coeffs']
#
#    p0 = np.insert(coeffs, 0, 1)
#
#    p0 = np.append(p0, default_scales)
#
#    bounds = ((None,None),)*np.size(p0)
#
#    pixelshift_model = True
#
#    if pixelshift_model is True:
#
#        p0 = np.append(p0,0.0)
#        bounds = (*bounds, (0,1))
#
##    print(bounds)
#
#    #
#    # Run the fit
#    #
#
#    # Build the arg array to pass to the optimizer
#    
#    args = (standard_wavelength,
#            standard_fluxdensity,
#            standard_rv,
#            standard_vmag,
#            standard_bmag,
#            vega_wavelength,
#            vega_fluxdensity,
#            vega_continuum,
#            vega_fitted_continuum,
#            kernel,
#            ew_scale,
#            line_wavelengths,
#            atmospheric_transmission,
#            poly_degree,
#            include_edges,
#            pixelshift_model)
#
#    # Do the fit
#
#    result = scipy.optimize.minimize(
#        ewscale_objfunction, 
#        p0, 
##        method='SLSQP',
#        method='BFGS',
#        args=args, 
#        bounds=bounds,
#        tol=1e-5)
#
#    # upack the results
#
#    p = result.x
#
##    print('p0',p0)
##    print('p',p)
#
#    idx = poly_degree+2
#    atm_scale = p[0]
#    coeffs = p[1:idx]
#
#    if pixelshift_model is True:
#
#        optimized_scales = p[idx:-1]
#        offset = p[-1]
#
#    else:
#
#        optimized_scales = p[idx:]
#
#    #
#    # plot the results
#    #
#
#    # Create a new vega model and ratio
#
#    if include_edges is True:
#                
#        scales = np.append(np.insert(optimized_scales, 0, ew_scale),
#                              ew_scale)
#    else:
#
#        scales = optimized_scales
#
#    # Create the new vega_scales array
#
#    f = scipy.interpolate.interp1d(
#        points, 
#        scales,
#        fill_value=(scales[0],scales[-1]),
#        bounds_error=False)
#
#    vega_scales = f(vega_wavelength)
#
#    # Create the new Vega model
#
#    result = modify_kuruczvega(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        vega_scales,
#        units='erg s-1 cm-2 A-1',
#        new=False)
#
#    # Create the new data model
#
#    new_vegamodel = result[0]/np.median(result[0])
#
##    atmosphere = (atm_scale*(atmospheric_transmission-1))+1
#
##    model = new_vegamodel*atmosphere*poly_1d(standard_wavelength,coeffs)
#
#
#    if pixelshift_model is True:
#
#        x = np.arange(len(standard_wavelength))
#        new_vegamodel = linear_interp1d(x,new_vegamodel,x-offset)
#
#
#    new_ratio = standard_fluxdensity/new_vegamodel
#
#    #
#    # Do the QA plotting
#    #
#
#    if qashow_info is not None:
#
#        plot_estimate_ewscales(
#            qashow_info['plot_number'],
#            qashow_info['figure_size'],
#            qashow_info['font_size'],
#            qashow_info['spectrum_linewidth'],
#            qashow_info['spine_linewidth'],
#            qashow_info['xlabel'],
#            qashow_info['title'],
#            standard_wavelength,
#            default_ratio,
#            new_ratio,
#            atmospheric_transmission,
#            line_wavelengths,
#            optimized_scales,
#            ew_scale)
#        
#        pl.show(block=qashow_info['block'])
#        if qashow_info['block'] is False:
#
#            pl.pause(1)
#
#    if qafile_info is not None:
#
#        plot_estimate_ewscales(
#            None,
#            qafile_info['figure_size'],
#            qafile_info['font_size'],
#            qafile_info['spectrum_linewidth'],
#            qafile_info['spine_linewidth'],
#            qafile_info['xlabel'],
#            qafile_info['title'],
#            standard_wavelength,
#            default_ratio,
#            new_ratio,
#            atmospheric_transmission,
#            line_wavelengths,
#            optimized_scales,
#            ew_scale)
#        
#        pl.savefig(qafile_info['file_fullpath'])
#        pl.close()
#    
#    return atm_scale, coeffs, optimized_scales
#
#
#
#def estimate_ewscales2(
#    standard_wavelength:npt.ArrayLike,
#    standard_fluxdensity:npt.ArrayLike,
#    standard_rv:int |float,
#    standard_vmag:int | float,
#    standard_bmag:int | float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    ew_scale:float | int,
#    line_wavelengths:npt.ArrayLike,
#    atmospheric_transmission:npt.ArrayLike,
#    poly_degree:int,
#    tolerance:float=0.01,
#    pixelshift_model:bool=False,
#    include_edges:bool=False,
#    qashow_info:dict=None,
#    qafile_info:dict=None):
#
#
#    """
#    To estimate the H line EW scale factors .
#
#    Parameters
#    ----------
#    standard_wavelength : ndarray
#        An (nstd,) array of standard star wavelengths.  
#
#    standard_fluxdensity : ndarray
#        An (nstd,) array of standard star flux densities.
#
#    standard_rv : float
#        A float giving the standard star radial velocities in km s-1.
#
#    standard_bmag : float
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float
#        The V-band magnitude of the standard star.
#
#    vega_wavelength : ndarray
#        An (nvega,) array of Vega model wavelengths.  
#
#    vega_fluxdensity : ndarray
#        An (nvega,) array of Vega model flux densities.
#        
#    vega_continuum : ndarray
#        An (nvega,) array of Vega model continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        An (nvega,) array of Vega model fitted continuum flux densities.
#
#    kernel : ndarray
#        An (nkernel,) kernel used to smooth the Vega model.
#
#    ew_scale : float
#        The EW scale factor.
#
#    poly_degree : int
#        The polynomial degree used to model the instrument throughput.
#
#    include_edges : {False, True}
#        Set to True to include the min() and max() values of 
#        `standard_wavelength` as control points.
#        Set to False to not include the min() and max() values of 
#        `standard_wavelength` as control points.
#
#    qashow_info : None or dict
#        `'number'` : int
#            The plot number.  Useful if you are doing to-screen plotting
#            because you can plot to the same window multiple times.
#    
#        `'scale'` : float or int
#            A scale factor to increase the size of the plot over the default.
#
#        `'block'`: {False, True}, optional
#            Set to make the plot block access to the command line, e.g.
#            pl.ioff().
#
#        `'xlabel'` : str, optional
#            A latex string giving the xlabel.
#
#        `'title'` : str, optional
#            A latex string giving the title of the plot.
#        
#    qafile_info : dict, optional    
#        `"filepath"` : str
#            The directory to write the QA figure.
#
#        `"filename"` : str
#            The name of the file, sans suffix/extension.
#
#        `"extension"` : str
#            The file extension.  Must be compatible with the savefig
#            function of matplotlib.
#
#        `'xlabel'` : str, optional
#            A latex string giving the xlabel.
#
#        `'title'` : str, optional
#            A (latex) string giving the title of the plot.
#
#    Returns
#    -------
#
#    """
#
#    #
#    # Check parameters
#    #
#
#    check_parameter('estimate_ewscales', 'standard_wavelength', 
#                    standard_wavelength, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'standard_fluxdensity', 
#                    standard_fluxdensity, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'standard_rv', standard_rv, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'standard_bmag', standard_bmag, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'standard_vmag', standard_vmag, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'vega_wavelength', 
#                    vega_wavelength, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_fluxdensity', 
#                    vega_fluxdensity, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_continuum', 
#                    vega_continuum, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'vega_fitted_continuum', 
#                    vega_fitted_continuum, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'kernel', kernel, 'ndarray')
#
#    check_parameter('estimate_ewscales', 'line_wavlenegths', line_wavelengths, 
#                    ['float', 'ndarray'])
#
#    check_parameter('estimate_ewscales', 'ew_scale', ew_scale, 'float')
#
#    check_parameter('estimate_ewscales', 'poly_degree', poly_degree, 
#                    'int')
#
#    check_parameter('estimate_ewscales', 'tolerance', tolerance, 
#                    'float')
#
#    check_parameter('estimate_ewscales', 'include_edges', include_edges, 
#                    'bool')
#
#    check_parameter('adjust_ews', 'qashow_info', qashow_info, 
#                    ['NoneType','dict'])
#
#    check_parameter('adjust_ews', 'qafile_info', qafile_info, 
#                    ['NoneType','dict'])
#
#
#    #
#    # Normalize the data
#    #
#
#    standard_fluxdensity /= np.nanmedian(standard_fluxdensity)
#    
#    #
#    # Create a default scales array
#    #
#
#    if include_edges is True:
#        
#        points = np.append(np.insert(line_wavelengths, 0, 
#                                     np.min(standard_wavelength)),
#                           np.max(standard_wavelength))
#        
#    else:
#
#        points = np.array(line_wavelengths)
#
#        
#    npoints = np.size(points)
#    default_scales = np.full_like(points, ew_scale) 
#
#    #
#    # Create the scales array on the Vega wavelength grid
#    #
#
#    f = scipy.interpolate.interp1d(
#        points, default_scales, 
#        fill_value=(default_scales[0],default_scales[-1]),
#        bounds_error=False)
#
#    vega_scales = f(vega_wavelength)
#
#    #
#    # Create default Vega model
#    #
#
#    result = modify_kuruczvega2(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        vega_scales,
#        units='erg s-1 cm-2 A-1')
#
#
#    # Normalize the vega model
#
#    default_vegamodel = result[0]/np.median(result[0])
#
#    #
#    # Estimate model parameters p0.
#    # p0[0] = atmospheric scale factor
#    # p0[1:poly_degree+2] = coeffs
#    # p0[poly_degree+2:] = EW scale factors
#    #
#
#    default_ratio = standard_fluxdensity/default_vegamodel
#
#
 #    results = polyfit_1d(
#        standard_wavelength, 
#        default_ratio, 
#        poly_degree,
#        robust={'thresh':4,'eps':0.1})
#
#    continuum=results['yfit']
#
#    coeffs = results['coeffs']
#
#    p0 = np.array([1])
#    bounds = (0.8,1)
#
#    p0 = np.append(p0, default_scales)
#
#    bounds = (bounds, *((0.8,1.2),)*npoints)
#
#    pixelshift_model = True
#
#    if pixelshift_model is True:
#
#        p0 = np.append(p0,0.5)
#        bounds = (*bounds, (0.0,1))
#
#
#    args = (standard_wavelength,
#            standard_fluxdensity,
#            standard_rv,
#            standard_vmag,
#            standard_bmag,
#            vega_wavelength,
#            vega_fluxdensity,
#            vega_continuum,
#            vega_fitted_continuum,
#            kernel,
#            ew_scale,
#            line_wavelengths,
#            atmospheric_transmission,
#            coeffs,
#            include_edges,
#            pixelshift_model)
#
##    # Do the fit
##
##
##    result = scipy.optimize.brute(
##        ewscale_objfunction2, 
##        bounds,
##        Ns=10,
##        args=args)
#
#
#
#
#    result = scipy.optimize.minimize(
#        ewscale_objfunction2, 
#        p0, 
#        method='Powell',
##        method='SLSQP',
##        method='Nelder-Mead',
##        method='L-BFGS-B',
##        method='TNC',
##        method='COBYLA',
##        method='COBYQA',
#        args=args, 
#        bounds=bounds,
#        tol=1e-3)
#
#
#
#
##    result = scipy.optimize.differential_evolution(
##        ewscale_objfunction2, 
##        bounds,
##        x0=p0, 
##        args=args, 
##        tol=tolerance)
##
##
##
##
##    # upack the results
##
#    p = result.x
##    print('fun',result.fun)
#
#
#    atm_scale = p[0]
#
#    if pixelshift_model is True:
#
#        optimized_scales = p[1:-1]
#        offset = p[-1]
#
#    else:
#
#        optimized_scales = p[1:]
#
##    print('scales',optimized_scales)
#
#    #
#    # plot the results
#    #
#
#    # Create a new vega model and ratio
#
#    if include_edges is True:
#                
#        scales = np.append(np.insert(optimized_scales, 0, ew_scale),
#                              ew_scale)
#    else:
#
#        scales = optimized_scales
#
#    # Create the new vega_scales array
#
#    f = scipy.interpolate.interp1d(
#        points, 
#        scales,
#        fill_value=(scales[0],scales[-1]),
#        bounds_error=False)
#
#    vega_scales = f(vega_wavelength)
#
#    # Create the new Vega model
#
#    result = modify_kuruczvega2(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        vega_scales,
#        units='erg s-1 cm-2 A-1')
#
#    # Create the new data model
#
#    new_vegamodel = result[0]/np.median(result[0])
#
##    atmosphere = (atm_scale*(atmospheric_transmission-1))+1
#
##    model = new_vegamodel*atmosphere*poly_1d(standard_wavelength,coeffs)
#
#
#    if pixelshift_model is True:
#
#        x = np.arange(len(standard_wavelength))
#        new_vegamodel = linear_interp1d(x,new_vegamodel,x-offset)
#
#
#    new_ratio = standard_fluxdensity/new_vegamodel
#
#    #
#    # Do the QA plotting
#    #
#
#    if qashow_info is not None:
#
#        plot_estimate_ewscales(
#            qashow_info['plot_number'],
#            qashow_info['figure_size'],
#            qashow_info['font_size'],
#            qashow_info['spectrum_linewidth'],
#            qashow_info['spine_linewidth'],
#            qashow_info['xlabel'],
#            qashow_info['title'],
#            standard_wavelength,
#            default_ratio,
#            new_ratio,
#            atmospheric_transmission,
#            line_wavelengths,
#            optimized_scales,
#            ew_scale)
#        
#        pl.show(block=qashow_info['block'])
#        if qashow_info['block'] is False:
#
#            pl.pause(1)
#
#    if qafile_info is not None:
#
#        plot_estimate_ewscales(
#            None,
#            qafile_info['figure_size'],
#            qafile_info['font_size'],
#            qafile_info['spectrum_linewidth'],
#            qafile_info['spine_linewidth'],
#            qafile_info['xlabel'],
#            qafile_info['title'],
#            standard_wavelength,
#            default_ratio,
#            new_ratio,
#            atmospheric_transmission,
#            line_wavelengths,
#            optimized_scales,
#            ew_scale,
#            continuum)
#        
#        pl.savefig(qafile_info['file_fullpath'])
#        pl.close()
#    
#    return atm_scale, coeffs, optimized_scales
#
#
#
def estimate_ewscales(
    standard_wavelength:npt.ArrayLike,
    standard_fluxdensity:npt.ArrayLike,
    standard_rv:int |float,
    standard_vmag:int | float,
    standard_bmag:int | float,
    vega_wavelength:npt.ArrayLike,
    vega_fluxdensity:npt.ArrayLike,
    vega_continuum:npt.ArrayLike,
    vega_fitted_continuum:npt.ArrayLike,
    kernel:npt.ArrayLike,
    vega_pixelshift:float,
    default_ewscale:float | int,
    ew_wavelengths:npt.ArrayLike,
    atmospheric_transmission:npt.ArrayLike,
    poly_degree:int,
    tolerance:float=0.01,
    include_edges:bool=True,
    qashow_info:dict=None,
    qafile_info:dict=None):


    """
    To estimate the H line EW scale factors.

    Parameters
    ----------
    standard_wavelength : ndarray
        An (nstd,) array of standard star wavelengths.  

    standard_fluxdensity : ndarray
        An (nstd,) array of standard star flux densities.

    standard_rv : float
        A float giving the standard star radial velocities in km s-1.

    standard_bmag : float
        The B-band magnitude of the standard star.

    standard_vmag : float
        The V-band magnitude of the standard star.

    vega_wavelength : ndarray
        An (nvega,) array of Vega model wavelengths.  

    vega_fluxdensity : ndarray
        An (nvega,) array of Vega model flux densities.
        
    vega_continuum : ndarray
        An (nvega,) array of Vega model continuum flux densities.

    vega_fitted_continuum : ndarray
        An (nvega,) array of Vega model fitted continuum flux densities.

    kernel : ndarray
        An (nkernel,) kernel used to smooth the Vega model.

    vega_pixelshift : float
        The number of pixels to shift the Vega model post convolution. 

    ew_scale : float
        An estimate of the EW scale factor.

    ew_wavelengths : ndarray
        An (nlines,) array of wavelengths of lines to adjust EWs.

    atmospheric_transmission: ndarray
        An (nstd,) array giving the atmospheric transmission.  

    poly_degree : int
        The polynomial degree used to model the instrument throughput.

    include_edges : {False, True}
        Set to True to include the min() and max() values of 
        `standard_wavelength` as control points.
        Set to False to not include the min() and max() values of 
        `standard_wavelength` as control points.

    qashow_info : None or dict
        `'number'` : int
            The plot number.  Useful if you are doing to-screen plotting
            because you can plot to the same window multiple times.
    
        `'scale'` : float or int
            A scale factor to increase the size of the plot over the default.

        `'block'`: {False, True}, optional
            Set to make the plot block access to the command line, e.g.
            pl.ioff().

        `'xlabel'` : str, optional
            A latex string giving the xlabel.

        `'title'` : str, optional
            A latex string giving the title of the plot.
        
    qafile_info : dict, optional    
        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

        `'xlabel'` : str, optional
            A latex string giving the xlabel.

        `'title'` : str, optional
            A (latex) string giving the title of the plot.

    Returns
    -------

    """

    #
    # Check parameters
    #

    check_parameter('estimate_ewscales', 'standard_wavelength', 
                    standard_wavelength, 'ndarray')

    check_parameter('estimate_ewscales', 'standard_fluxdensity', 
                    standard_fluxdensity, 'ndarray')

    check_parameter('estimate_ewscales', 'standard_rv', 
                    standard_rv, 'float')

    check_parameter('estimate_ewscales', 'standard_bmag', 
                    standard_bmag, 'float')

    check_parameter('estimate_ewscales', 'standard_vmag', 
                    standard_vmag, 'float')

    check_parameter('estimate_ewscales', 'vega_wavelength', 
                    vega_wavelength, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_fluxdensity', 
                    vega_fluxdensity, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_continuum', 
                    vega_continuum, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_fitted_continuum', 
                    vega_fitted_continuum, 'ndarray')

    check_parameter('estimate_ewscales', 'kernel', 
                    kernel, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_pixelshift', 
                    vega_pixelshift, 'float')

    check_parameter('estimate_ewscales', 'ew_wavlengths', 
                    ew_wavelengths, ['float', 'ndarray'])

    check_parameter('estimate_ewscales', 'default_ewscale', 
                    default_ewscale, 'float')

    check_parameter('estimate_ewscales', 'atmospheric transmission', 
                    atmospheric_transmission, 'ndarray')

    check_parameter('estimate_ewscales', 'poly_degree', 
                    poly_degree, 'int')

    check_parameter('estimate_ewscales', 'tolerance', 
                    tolerance, 'float')

    check_parameter('estimate_ewscales', 'include_edges', 
                    include_edges, 'bool')

    check_parameter('adjust_ews', 'qashow_info', 
                    qashow_info, ['NoneType','dict'])

    check_parameter('adjust_ews', 'qafile_info', 
                    qafile_info, ['NoneType','dict'])

    #
    # Create default Vega model
    #

    vega_setup = _modify_kuruczvega_setup(
        standard_wavelength,
        standard_rv,
        standard_vmag,
        standard_bmag,
        vega_wavelength,
        vega_fluxdensity,
        vega_continuum,
        vega_fitted_continuum,
        kernel,
        vega_pixelshift)

    result = _modify_kuruczvega_execute(vega_setup)

    # parse results and normalize the model

    defaultvega_fluxdensity = result[0]
    defaultvega_fluxdensity /= np.nanmedian(defaultvega_fluxdensity) 

    #
    # Estimate model parameters p0.
    # p0[0] = atmospheric scale factor
    # p0[1:poly_degree+2] = coeffs
    # p0[poly_degree+2:-2] = EW scale factors
    # p0[-1] = Pixel offset
    #

    normalized_standard_fluxdensity = standard_fluxdensity/np.nanmedian(standard_fluxdensity)


    default_ratio = normalized_standard_fluxdensity/defaultvega_fluxdensity

    results = polyfit_1d(
        standard_wavelength, 
        default_ratio, 
        poly_degree,
        robust={'thresh':4,'eps':0.1})
    
    p0 = []

    # Atmosphere parameter


    p0 = np.array([1])
    p0 = np.append(p0, results['coeffs'])
    p0 = np.append(p0, np.full(np.size(ew_wavelengths),default_ewscale))
    p0 = np.append(p0,0.0)

    args = (vega_setup,
            normalized_standard_fluxdensity,
            default_ewscale,
            ew_wavelengths,
            atmospheric_transmission,
            poly_degree)

    result = scipy.optimize.minimize(
        ewscale_objfunction, 
        p0, 
        method='Powell',
#        method='SLSQP',
#        method='Nelder-Mead',
#        method='L-BFGS-B',
#        method='TNC',
#        method='COBYLA',
#        method='COBYQA',
        args=args, 
#        bounds=bounds,
        tol=tolerance)

#
    # upack the results
#

    p = result.x



    idx = poly_degree+2
    atm_scale = p[0]
    coeffs = p[1:idx]
    optimized_scales = p[idx:-1]


    continuum = poly_1d(vega_setup['standard_wavelength'],coeffs)

    

    #  Add the order edges to the control points

    points = np.insert(ew_wavelengths, 0, vega_setup['min_vega_wavelength'])
    points = np.append(points, vega_setup['max_vega_wavelength'])
        
    values = np.insert(optimized_scales, 0, default_ewscale)
    values = np.append(values, default_ewscale)

    ew_info = {'wavelengths':points, 'scales':values}

    result = _modify_kuruczvega_execute(
        vega_setup,
        ew_info=ew_info)
        
    #
    # Build the model of the data
    #

    vega_model = result[0]
    vega_model /= np.median(vega_model)
    new_ratio = normalized_standard_fluxdensity/vega_model

    #
    # Do the QA plotting
    #

    if qashow_info is not None:

        plot_estimate_ewscales(
            qashow_info['plot_number'],
            qashow_info['figure_size'],
            qashow_info['font_size'],
            qashow_info['spectrum_linewidth'],
            qashow_info['spine_linewidth'],
            qashow_info['xlabel'],
            qashow_info['title'],
            standard_wavelength,
            default_ratio,
            new_ratio,
            atmospheric_transmission,
            ew_wavelengths,
            optimized_scales,
            default_ewscale,
            continuum)
        
        pl.show(block=qashow_info['block'])
        if qashow_info['block'] is False:

            pl.pause(1)

    if qafile_info is not None:

        plot_estimate_ewscales(
            None,
            qafile_info['figure_size'],
            qafile_info['font_size'],
            qafile_info['spectrum_linewidth'],
            qafile_info['spine_linewidth'],
            qafile_info['xlabel'],
            qafile_info['title'],
            standard_wavelength,
            default_ratio,
            new_ratio,
            atmospheric_transmission,
            ew_wavelengths,
            optimized_scales,
            default_ewscale,
            continuum)
        
        pl.savefig(qafile_info['file_fullpath'])
        pl.close()
    
    return atm_scale, coeffs, optimized_scales





#def ewscale_objfunctionold(
#    p:npt.ArrayLike, 
#    standard_wavelength:npt.ArrayLike,
#    standard_fluxdensity:npt.ArrayLike,
#    standard_rv:int | float,
#    standard_vmag:int | float,
#    standard_bmag:int | float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    startew_scale:int | float,
#    line_wavelengths:npt.ArrayLike,
#    atmospheric_transmission:npt.ArrayLike,
#    poly_degree:int,
#    include_edges:bool,
#    pixelshift_model:bool):
#    
#    """
#    The objective function used to determine the H line EW scale factors.
#
#    The function is used along with scipy.optmize.minimize in order to 
#    determine optimal model parameters of a standard star by minimizing 
#    this function.
#
#    Parameters
#    ----------
#    p : ndarray
#        A array of model parameters.  The length depends on the number of 
#        hydrogen lines whose EWs need adjustment and the polynomial degree 
#        used to model the instrument throughput.  
#
#       atm_scale = p[0]
#       coeffs = p[1:`poly_degree`+2]
#       ew_scales = p[`poly_degree`+2:]
#
#    standard_wavelength : ndarray
#        An (nstd,) array of standard star wavelengths.  
#
#    standard_fluxdensity : ndarray
#        An (nstd,) array of standard star flux densities.
#
#    standard_rv : float
#        A float giving the standard star radial velocities in km s-1.
#
#    standard_bmag : float
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float
#        The V-band magnitude of the standard star.
#
#    vega_wavelength : ndarray
#        An (nvega,) array of Vega model wavelengths.  
#
#    vega_fluxdensity : ndarray
#        An (nvega,) array of Vega model flux densities.
#        
#    vega_continuum : ndarray
#        An (nvega,) array of Vega model continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        An (nvega,) array of Vega model fitted continuum flux densities.
#
#    kernel : ndarray
#        An (nkernel,) kernel used to smooth the Vega model.
#
#    start_ewscale : float
#        The EW scale factor.
#
#    poly_degree : int
#        The polynomial degree used to model the instrument throughput.
#
#    include_edges : {False, True}
#        Set to True to include the min() and max() values of 
#        `standard_wavelength` as control points.
#        Set to False to not include the min() and max() values of 
#        `standard_wavelength` as control points.
#
#    Returns
#    -------
#    float 
#    
#        The summed squared residuals between the data and model.  
#
#    Notes
#    -----
#    We model the A0 V standard star spectrum over some limited wavelength
#    range as the product of a Vega model, the atmospheric transmission, and 
#    a polynomial to account for the smooth throughput of the instrument.  The 
#    EWs of the hydrogen lines in the Vega model are allowed to vary and 
#    the scale factors used to adjust the EWs are parameters of the model.
#
#    The length of the parameter array is a function of the number of hydrogen 
#    list whose EWs need adjustment, and the degree of the polynomial used to 
#    model the instrument throughput.  `poly_degree` is the polynomial degree
#    and so the model parameters are given by:
#
#    atm_scale = p[0]
#    coeffs = p[1:poly_degree+2]
#    ew_scales = p[poly_degree+2:]
#
#    """
#
#    #
#    # Check parameters
#    #
#
#    check_parameter('ewscale_objfunction', 'p', p, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_wavelength', 
#                    standard_wavelength, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_fluxdensity', 
#                    standard_fluxdensity, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_rv', standard_rv, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'standard_bmag', standard_bmag, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'standard_vmag', standard_vmag, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'vega_wavelength', 
#                    vega_wavelength, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_fluxdensity', 
#                    vega_fluxdensity, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_continuum', 
#                    vega_continuum, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_fitted_continuum', 
#                    vega_fitted_continuum, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'kernel', kernel, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'line_wavlenegths', line_wavelengths, 
#                    ['float', 'ndarray'])
#
#    check_parameter('ewscale_objfunction', 'startew_scale', startew_scale, 'float')
#
#    check_parameter('ewscale_objfunction', 'atmospheric_transmission', 
#                    atmospheric_transmission, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'poly_degree', poly_degree, 
#                    'int')
#
#    check_parameter('ewscale_objfunction', 'include_edges', include_edges, 
#                    'bool')
#
#    #
#    # Upack the parameters
#    #
#
#    idx = poly_degree+2
#    atm_scale = p[0]
#    coeffs = p[1:idx]
#
#    if pixelshift_model is True:
#
#        ew_scales = p[idx:-1]
#        offset = p[-1]
#
#    else:
#
#        ew_scales = p[idx:]
#
#    #
#    # Does the user request adding the order edges to the control points?
#    #
#    
#    if include_edges is True:
#        
#        points = np.append(np.insert(line_wavelengths, 0, 
#                                     np.min(standard_wavelength)),
#                           np.max(standard_wavelength))
#        
#        
#        
#        scales = np.append(np.insert(ew_scales, 0, startew_scale),
#                           startew_scale)
#
#    else:
#
#        points = line_wavelengths
#        scales = ew_scales
#
##    cs = CubicSpline(ew_points, ew_values)
##    ewscales = cs(vega_wavelength)
#
#    #
#    # Interpolate the control points onto the vega model wavelength grid
#    #
#
#    f = scipy.interpolate.interp1d(
#        points, ew_scales, 
#        fill_value=(ew_scales[0], ew_scales[-1]),
#        bounds_error=False)
#    scales = f(vega_wavelength)
#
#    #
#    #  Generate the Vega model
#    #
#
#    result = modify_kuruczvega2(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        scales,
#        units='erg s-1 cm-2 A-1')
#
#    #
#    # Build the model of the data
#    #
#
#    vega_model = result[0]
#    vega_model /= np.median(vega_model)
#
#    if pixelshift_model is True:
#
#        print(offset)
#        x = np.arange(len(standard_wavelength))
#        rvega_model = linear_interp1d(x,vega_model,x-offset)
#
#    else:
#
#        rvega_model = vega_model
#
#    atmosphere = (atm_scale*(atmospheric_transmission-1))+1
#
#    model = rvega_model*atmosphere*poly_1d(standard_wavelength,coeffs)
#
#    #
#    # Compute the sum of the squared residuals
#    #
#
#    objective = np.nansum((standard_fluxdensity-model)**2)
#
#    #
#    # Return the value
#    #
#
#    return objective
#
#
#def ewscale_objfunction2(
#    p:npt.ArrayLike, 
#    standard_wavelength:npt.ArrayLike,
#    standard_fluxdensity:npt.ArrayLike,
#    standard_rv:int | float,
#    standard_vmag:int | float,
#    standard_bmag:int | float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    startew_scale:int | float,
#    line_wavelengths:npt.ArrayLike,
#    atmospheric_transmission:npt.ArrayLike,
#    ratio_coeffs:npt.ArrayLike,
#    include_edges:bool,
#    pixelshift_model:bool):
#    
#    """
#    The objective function used to determine the H line EW scale factors.
#
#    The function is used along with scipy.optmize.minimize in order to 
#    determine optimal model parameters of a standard star by minimizing 
#    this function.
#
#    Parameters
#    ----------
#    p : ndarray
#        A array of model parameters.  The length depends on the number of 
#        hydrogen lines whose EWs need adjustment and the polynomial degree 
#        used to model the instrument throughput.  
#
#       atm_scale = p[0]
#       coeffs = p[1:`poly_degree`+2]
#       ew_scales = p[`poly_degree`+2:]
#
#    standard_wavelength : ndarray
#        An (nstd,) array of standard star wavelengths.  
#
#    standard_fluxdensity : ndarray
#        An (nstd,) array of standard star flux densities.
#
#    standard_rv : float
#        A float giving the standard star radial velocities in km s-1.
#
#    standard_bmag : float
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float
#        The V-band magnitude of the standard star.
#
#    vega_wavelength : ndarray
#        An (nvega,) array of Vega model wavelengths.  
#
#    vega_fluxdensity : ndarray
#        An (nvega,) array of Vega model flux densities.
#        
#    vega_continuum : ndarray
#        An (nvega,) array of Vega model continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        An (nvega,) array of Vega model fitted continuum flux densities.
#
#    kernel : ndarray
#        An (nkernel,) kernel used to smooth the Vega model.
#
#    start_ewscale : float
#        The EW scale factor.
#
#    poly_degree : int
#        The polynomial degree used to model the instrument throughput.
#
#    include_edges : {False, True}
#        Set to True to include the min() and max() values of 
#        `standard_wavelength` as control points.
#        Set to False to not include the min() and max() values of 
#        `standard_wavelength` as control points.
#
#    Returns
#    -------
#    float 
#    
#        The summed squared residuals between the data and model.  
#
#    Notes
#    -----
#    We model the A0 V standard star spectrum over some limited wavelength
#    range as the product of a Vega model, the atmospheric transmission, and 
#    a polynomial to account for the smooth throughput of the instrument.  The 
#    EWs of the hydrogen lines in the Vega model are allowed to vary and 
#    the scale factors used to adjust the EWs are parameters of the model.
#
#    The length of the parameter array is a function of the number of hydrogen 
#    list whose EWs need adjustment, and the degree of the polynomial used to 
#    model the instrument throughput.  `poly_degree` is the polynomial degree
#    and so the model parameters are given by:
#
#    atm_scale = p[0]
#    coeffs = p[1:poly_degree+2]
#    ew_scales = p[poly_degree+2:]
#
#    """
#
#    #
#    # Check parameters
#    #
#
#    check_parameter('ewscale_objfunction', 'p', p, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_wavelength', 
#                    standard_wavelength, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_fluxdensity', 
#                    standard_fluxdensity, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'standard_rv', standard_rv, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'standard_bmag', standard_bmag, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'standard_vmag', standard_vmag, 
#                    'float')
#
#    check_parameter('ewscale_objfunction', 'vega_wavelength', 
#                    vega_wavelength, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_fluxdensity', 
#                    vega_fluxdensity, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_continuum', 
#                    vega_continuum, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'vega_fitted_continuum', 
#                    vega_fitted_continuum, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'kernel', kernel, 'ndarray')
#
#    check_parameter('ewscale_objfunction', 'line_wavlenegths', line_wavelengths, 
#                    ['float', 'ndarray'])
#
#    check_parameter('ewscale_objfunction', 'startew_scale', startew_scale, 'float')
#
#    check_parameter('ewscale_objfunction', 'atmospheric_transmission', 
#                    atmospheric_transmission, 'ndarray')
#
##    check_parameter('ewscale_objfunction', 'poly_degree', poly_degree, 
##                    'int')
#
#    check_parameter('ewscale_objfunction', 'include_edges', include_edges, 
#                    'bool')
#
#    #
#    # Upack the parameters
#    #
#
#    atm_scale = p[0]
#    coeffs = ratio_coeffs
#
#    if pixelshift_model is True:
#
#        ew_scales = p[1:-1]
#        offset = p[-1]
#
#    else:
#
#        ew_scales = p[1:]
#
#    print(ew_scales)
#
#    #
#    # Does the user request adding the order edges to the control points?
#    #
#    
#    if include_edges is True:
#        
#        points = np.append(np.insert(line_wavelengths, 0, 
#                                     np.min(standard_wavelength)),
#                           np.max(standard_wavelength))
#        
#        
#        
#        scales = np.append(np.insert(ew_scales, 0, startew_scale),
#                           startew_scale)
#
#    else:
#
#        points = line_wavelengths
#        scales = ew_scales
#
##    cs = CubicSpline(ew_points, ew_values)
##    ewscales = cs(vega_wavelength)
#
#    #
#    # Interpolate the control points onto the vega model wavelength grid
#    #
#
#    f = scipy.interpolate.interp1d(
#        points, ew_scales, 
#        fill_value=(ew_scales[0], ew_scales[-1]),
#        bounds_error=False)
#    scales = f(vega_wavelength)
#
#    #
#    #  Generate the Vega model
#    #
#
#    result = modify_kuruczvega2(
#        standard_wavelength,
#        standard_rv,
#        standard_vmag,
#        standard_bmag,
#        vega_wavelength,
#        vega_fluxdensity,
#        vega_continuum,
#        vega_fitted_continuum,
#        kernel,
#        scales,
#        units='erg s-1 cm-2 A-1')
#
#    #
#    # Build the model of the data
#    #
#
#    vega_model = result[0]
#    vega_model /= np.median(vega_model)
#
#    if pixelshift_model is True:
#
#        print(offset)
#        x = np.arange(len(standard_wavelength))
#        rvega_model = linear_interp1d(x,vega_model,x-offset)
#
#    else:
#
#        rvega_model = vega_model
#
#    atmosphere = (atm_scale*(atmospheric_transmission-1))+1
#
#    model = rvega_model*atmosphere*poly_1d(standard_wavelength,coeffs)
#
#    #
#    # Compute the sum of the squared residuals
#    #
#
#    objective = np.nansum((standard_fluxdensity-model)**2)
#    print('objective', objective)
#
#    #
#    # Return the value
#    #
#
#    return objective


def ewscale_objfunction(
    p:npt.ArrayLike,
    setup:dict,
    standard_fluxdensity:npt.ArrayLike,
    default_ewscale:float,
    ew_wavelengths:npt.ArrayLike,
    atmospheric_transmission,
    poly_degree:int):

#    args = (vega_setup,
#            standard_fluxdensity,
#            default_ewscale,
#            ew_wavelengths,
#            atmospheric_transmission,
#            poly_degree)


    """
    The objective function used to optimize the EWs of standard star lines

    Parameters
    ----------
    p : ndarray
        
    setup : dict
        The set up dictionary generated by `_modify_kuruczvega_setup`

    ew_scale : float
        The default EW scale factor

    ew_wavelengths : ndarray


    Notes
    -----
    We model the A0 V standard star spectrum over some limited wavelength
    range as the product of a Vega model, the atmospheric transmission, and 
    a polynomial to account for the smooth throughput of the instrument.  The 
    EWs of the hydrogen lines in the Vega model are allowed to vary and 
    the scale factors used to adjust the EWs are parameters of the model.

    The length of the parameter array is a function of the number of hydrogen 
    list whose EWs need adjustment, and the degree of the polynomial used to 
    model the instrument throughput.  `poly_degree` is the polynomial degree
    and so the model parameters are given by:

    atm_scale = p[0]
    coeffs = p[1:poly_degree+2]
    ew_scales = p[poly_degree+2:]

    """

    #
    # Upack the parameters
    #


    idx = poly_degree+2
    atm_scale = p[0]
    coeffs = p[1:idx]
    ew_scales = p[idx:-1]

    #  Add the order edges to the control points

    points = np.insert(ew_wavelengths, 0, setup['min_vega_wavelength'])
    points = np.append(points, setup['max_vega_wavelength'])
        
    values = np.insert(ew_scales, 0, default_ewscale)
    values = np.append(values, default_ewscale)

    ew_info = {'wavelengths':points, 'scales':values}

    result = _modify_kuruczvega_execute(
        setup,
        ew_info=ew_info)
        
    #
    # Build the model of the data
    #

    vega_model = result[0]
    vega_model /= np.median(vega_model)

    atmosphere = (atm_scale*(atmospheric_transmission-1))+1

    model = vega_model*atmosphere*poly_1d(setup['standard_wavelength'],coeffs)

    #
    # Compute the sum of the squared residuals and return the results
    #

    objective = np.nansum((standard_fluxdensity-model)**2)

    return objective



def find_modelshift(
    standard_wavelength:npt.ArrayLike,
    standard_fluxdensity:npt.ArrayLike,
    vega_wavelength:npt.ArrayLike,
    vega_fluxdensity:npt.ArrayLike,
    wavelength_range:npt.ArrayLike,
    qafile_info:dict=None):

    """
    To shift the telluric standard model 

    Parameters
    ----------
    standard_wavelength : ndarray
        An (nstd,) array of standard star wavelengths.  

    standard_fluxdensity : ndarray
        An (nstd,) array of standard star flux densities.

    vega_wavelength : ndarray
        An (nvega,) array of Vega wavelengths.  

    vega_fluxdensity : ndarray
        An (nvega,) array of Vega flux densities.

    qafile_info : dict, optional    
        
        `"figure_size"` : tuple
            A (2,) tuple of the figure size in inches.

        `"font_size"` : int
            The font size.

        `"spectrum_linewidth"` : float
            The line thickness of the spectrum.

        `"spine_linewidth"` : float
            The spine thickness of the spectrum.

        `"file_fullpath"`: str
            The fullpath to the file to be written to.

        `'xlabel'` : str, optional
            A latex string giving the xlabel.

        `'title'` : str, optional
            A (latex) string giving the title of the plot.
    
    Returns
    -------
    float
        The shift of the model relative to the data.  
    
    """

    #
    # Check Parameters
    #

    check_parameter('find_modelshift', 'standard_wavelength', 
                    standard_wavelength, 'ndarray', 1)

    check_parameter('find_modelshift', 'standard_flulxdensity', 
                    standard_fluxdensity, 'ndarray', 1)

    check_parameter('find_modelshift', 'vega_wavelength', 
                    vega_wavelength, 'ndarray', 1)

    check_parameter('find_modelshift', 'vega_flulxdensity', 
                    vega_fluxdensity, 'ndarray', 1)

    check_parameter('find_modelshift', 'wavelength_range', 
                    wavelength_range, 'list')

    check_parameter('find_modelshift', 'qafile_info', 
                    qafile_info, ['dict','NoneType'])


    #
    # Start the process
    #

    # Clip the standard to 'wavelength_range' and normalize

    z = np.where((standard_wavelength > wavelength_range[0]) & 
                 (standard_wavelength < wavelength_range[1]))[0]
    
    standard_wavelength = standard_wavelength[z]
    standard_fluxdensity = standard_fluxdensity[z]
    standard_fluxdensity /= np.median(standard_fluxdensity)
    
    # Clip the model to 'wavelength_range' and normalize

    z = np.where((vega_wavelength > wavelength_range[0]) & 
                 (vega_wavelength < wavelength_range[1]))[0]

    vega_wavelength = vega_wavelength[z]
    vega_fluxdensity = vega_fluxdensity[z]
    
    vega_fluxdensity /= np.median(vega_fluxdensity)

    # Resample the model onto the standard

    rvega_fluxdensity = linear_interp1d(
        vega_wavelength, 
        vega_fluxdensity, 
        standard_wavelength)
    
    # Find the shift

    result = find_shift(
        standard_wavelength,
        standard_fluxdensity,
        rvega_fluxdensity,
        wavelength_range,
        pixel_range=[-3,3])

    #
    # Plot the results
    #

    if qafile_info is not None:

        x = np.arange(len(standard_wavelength))
        xshift = x+result

        rsvega_fluxdensity = linear_interp1d(xshift,rvega_fluxdensity,x)

        qaplots.plot_find_modelshift(
            None,
            qafile_info['figure_size'],
            qafile_info['font_size'],
            qafile_info['spectrum_linewidth'],
            qafile_info['spine_linewidth'],
            standard_wavelength,
            standard_fluxdensity,
            rvega_fluxdensity,
            result,
            rsvega_fluxdensity,
            qafile_info['xlabel'])

        pl.savefig(qafile_info['file_fullpath'])
        pl.close()
            
    return result
    


def find_shift(
    wavelength_object:npt.ArrayLike,
    intensity_object:npt.ArrayLike,
    intensity_telluric:npt.ArrayLike,
    wavelength_range:list | tuple,
    pixel_range:list=[-1.5,1.5],
    nsteps:int=301,
    qa_show=False,
    qa_showblock=True):

    """
    To determine the shift between spectra that minimizes telluric noise.


    Parameters
    ----------
    wavelength_object: ndarray
        An (nwave,) array of wavelength values of the object spectrum.

    intensity_object: ndarray
        An (nwave,) array of intensity values  of the object spectrum.

    intensity_telluric: ndarray
        An (nwave,) array of intensity values  of the telluric spectrum.
    
    wavelength_range : list
        An (2,) list of wavelengths over which the noise minimization is
        desired.

    pixel_range : list, default [-1.5,1.5]
        A (2,) list giving the number of pixels over which to do the
        noise minimization.

    nsteps : int
        The number of steps between `pixel_range` to evaluate the rms at.

    Returns
    -------
    float
        The shift in pixels that minimizes the noise in `wavelength_range`


    
    """

    #
    # Check parameters
    #

    check_parameter('find_shift', 'wavelength_object', 
                    wavelength_object, 'ndarray')

    check_parameter('find_shift', 'intensity_object', 
                    intensity_object, 'ndarray')

    check_parameter('find_shift', 'intensity_telluric', 
                    intensity_telluric, 'ndarray')
    
    check_parameter('find_shift', 'wavelength_range', 
                    wavelength_range, ['list','tuple'])

    check_parameter('find_shift', 'nsteps', 
                    nsteps, 'int')

    #
    # Get set up
    #

    wrange = np.sort(wavelength_range)
    
    zdata = np.where(
        (wavelength_object >= wrange[0]) &
        (wavelength_object <= wrange[1]))[0]
    
    shifts = np.linspace(pixel_range[0], pixel_range[1], num=nsteps)
    nshifts = len(shifts)

    rms = np.zeros(nshifts)
    
    x = np.arange(len(wavelength_object))

    #
    # Start the loop
    #
    

#    pl.figure()
    for i in range(nshifts):

        xshift = x+shifts[i]

        telluric_shifted = np.interp(
            x,
            xshift,
            intensity_telluric,
            left=np.nan,
            right=np.nan)

        ratio = np.multiply(intensity_object,telluric_shifted)

        rms[i] = np.nanstd(ratio[zdata])

#        pl.step(wavelength_object[zdata],ratio[zdata],label='Ratio')
#        pl.step(wavelength_object,intensity_object,color='red',
#                label='object')
#        pl.step(wavelength_object,telluric_shifted,color='green',
#                label='telluric')
#        pl.step(wavelength_object,intensity_telluric,color='green',
#                label='telluric')
#        pl.title(i)
#        pl.legend()
#        pl.pause(0.01)
#        pl.clf()
#                
#    pl.plot(shifts,rms)
#    pl.show()


    #
    # Find the shift
    #

    # Find the minimum rms value
    
    minimum_idx = np.argmin(rms)

    # Fit a 2nd order polynomial around this point

    # temporary fix for end points
    if minimum_idx<5 or minimum_idx>len(shifts)-6: 
        logging.info(' find shift found minimum near end point {:.0f}, assuming zero shift '.format(minimum_idx))
        return 0.    

    result = polyfit_1d(shifts[(minimum_idx-5):(minimum_idx+6)],
                        rms[(minimum_idx-5):(minimum_idx+6)],2)

    # Set the derivative equal to zero to find the minimum shift
    
    shift = -result['coeffs'][1]/2./result['coeffs'][2]

    return shift


def get_kernelwavelengths(
    standard_wavelengths:npt.ArrayLike,
    nkernel:int):

    """
    To construct an array containing the wavelengths associated with the kernel.

    Parameters
    ----------
    standard_wavelength : ndarray
        An (nstd,) array of data wavelengths.

    nkernel: int
        The length of the kernel

    Returns
    -------
    ndarray
        An (nkernel,nstd) array.


    """

    #
    # Check parameters
    #

    check_parameter('get_kernelwavelengths', 'standard_wavelengths', 
                    standard_wavelengths, 'ndarray')

    check_parameter('get_kernelwavelengths', 'nkernel', 
                    nkernel, 'int')

    #  Is nkernel odd?

    if nkernel % 2 == 0:

        message = 'Parameter `nkernel` must be odd.'
        raise pySpextoolError(message)


    half_nkernel = nkernel // 2

    nwavelengths = np.size(standard_wavelengths)

    output = np.full((nkernel,nwavelengths),np.nan)

    for i in range(nkernel):

        shift = half_nkernel-i
        
        test = np.roll(standard_wavelengths,shift)
        if shift < 0:  
            
            test[shift:] = np.nan
            
        if shift > 0:
            
            test[:shift] = np.nan
    
        output[i,:] = test

    # create mask 

    mask = np.tile(standard_wavelengths*0.0+1.0, (nkernel, 1))
    output *= mask

    # Loop over each wavelength to search for incomplete columns.

    for i in range(nwavelengths):

        # If there is no NaN, then the column is complete.

        if np.sum(np.isnan(output[:,i])) == 0:

            continue

        # Found an incomplete column.

        for j in range(half_nkernel):

            # First we move up from the middle


            idx = half_nkernel+j+1
            test = output[idx,i]
            if np.isnan(test):

                delta = output[idx-1,i]-output[idx-2,i]
                output[idx,i] = output[idx-1,i]+delta

            # Now we move down from the middle

            idx = half_nkernel-j-1
            test = output[idx,i]

            if np.isnan(test):

                delta = output[idx+2,i]-output[idx+1,i]
                output[idx,i] = output[idx+1,i]-delta
                

    return output

    
def make_instrument_profile(
        x:npt.ArrayLike,
        parameters:npt.ArrayLike):

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


def make_telluric_spectrum(
    standard_wavelength:npt.ArrayLike,
    standard_fluxdensity:npt.ArrayLike,
    standard_uncertainty:npt.ArrayLike,
    standard_rv:float,
    standard_vmag:float,
    standard_bmag:float,
    vega_wavelength:npt.ArrayLike,
    vega_fluxdensity:npt.ArrayLike,
    vega_continuum:npt.ArrayLike,
    vega_fitted_continuum:npt.ArrayLike,
    kernel:npt.ArrayLike,
    vega_pixelshift:float,
    control_points:npt.ArrayLike,
    control_values:npt.ArrayLike,
    order,
    intensity_unit:str='erg s-1 cm-2 A-1'):

    """
    To create a telluric correction spectrum


    Parameters
    ----------
    standard_wavelength : ndarray
        A (nwave1,) array of data wavelengths.

    standard_fluxdensity : ndarray
        A (nwave1,) array of data flux densities.

    standard_uncertainty : ndarray
        A (nwave1,) array of data uncertainties.

    standard_rv : float, int
        The radial velocities of the standard star in km s-1.

    standard_bmag : float, int
        The B-band magnitude of the standard star.

    standard_vmag : float, int
        The V-band magnitude of the standard star.
    
    vega_wavelength : ndarray
        A (nwave2,) array of Vega-model wavelengths.

    vega_fluxdensity : ndarray
        A (nwave2,) array of Vega-model flux densities.

    vega_continuum : ndarray
        A (nwave2,) array of Vega-continuum flux densities.

    vega_fitted_continuum : ndarray
        A (nwave2,) array of Vega-fitted-continuum flux densities.

    kernel : ndarray
        An array of of kernel values.
    
    Returns
    -------
    ndarray

        The telluric correction spectrum.

    ndarray

        The uncertainty on the telluric correction spectrum.

    ndarray

        The modified Vega model spectrum.

    
    """

    #
    # Check the parameters
    #
    
    check_parameter('make_telluric_spectrum', 'standard_wavelength',
                    standard_wavelength, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'standard_fluxdensity',
                    standard_fluxdensity, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'standard_uncertainty',
                    standard_uncertainty, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'standard_bmag',
                    standard_bmag, ['float','int'])

    check_parameter('make_telluric_spectrum', 'standard_vmag',
                    standard_vmag, ['float','int'])

    check_parameter('make_telluric_spectrum', 'standard_rv',
                    standard_rv, ['float64','float','int'])
    
    check_parameter('make_telluric_spectrum', 'vega_wavelength',
                    vega_wavelength, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'vega_fluxdensity',
                    vega_fluxdensity, 'ndarray', 1)
    
    check_parameter('make_telluric_spectrum', 'vega_continuum',
                    vega_continuum, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'vega_fitted_continuum',
                    vega_fitted_continuum, 'ndarray', 1)

    check_parameter('make_telluric_spectrum', 'kernel', 
                    kernel, 'ndarray', 1)


    #
    # Determine the EW scale factor array using the control points
    #

#    f = scipy.interpolate.interp1d(
#        control_points,
#        control_values,
#        fill_value=(control_values[0],control_values[-1]),
#        bounds_error=False)
#    ewscales = f(vega_wavelength)

    #
    # Modify the Kurucz Vega model 
    #

#    standard_wavelength:npt.ArrayLike,
#    standard_rv:float,
#    standard_vmag:float,
#    standard_bmag:float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    vega_pixelshift:float,
#    ew_wavelengths:npt.ArrayLike,
#    ew_scales:npt.ArrayLike,
#    edge_scale:float,
#    Rv:float=3.1,
#    vega_vmag:float=0.03,
#    vega_bminv:float=0.0,
#    units:str='erg s-1 cm-2 A-1',
#    return_raw_model:bool=False):





    vega_fd, vega_cont = modify_kuruczvega(
        standard_wavelength,
        standard_rv,
        standard_vmag,
        standard_bmag,
        vega_wavelength,
        vega_fluxdensity,
        vega_continuum,
        vega_fitted_continuum,
        kernel,
        vega_pixelshift,
        ew_info={'wavelengths':control_points, 'scales':control_values},
        units=intensity_unit)

#    standard_wavelength:npt.ArrayLike,
#    standard_rv:float,
#    standard_vmag:float,
#    standard_bmag:float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    vega_pixelshift:float,
#    ew_info:dict,
#    Rv:float=3.1,
#    vega_vmag:float=0.03,
#    vega_bminv:float=0.0,
#    units:str='erg s-1 cm-2 A-1',
#    return_raw_model:bool=False):
#



    #
    # Create the telluric correction spectrum and return results
    #

    telluric = vega_fd/standard_fluxdensity
    
    telluric_var = (vega_fd/standard_fluxdensity**2)**2 * \
        standard_uncertainty**2
    
    telluric_unc = np.sqrt(telluric_var)
    
    return telluric, telluric_unc, vega_fd, vega_cont


def measure_linerv(
    data_wavelength:npt.ArrayLike,
    data_fluxdensity:npt.ArrayLike,
    model_wavelength:npt.ArrayLike,
    model_fluxdensity:npt.ArrayLike,
    peak_method:str='max',
    resolving_power:float=None,
    qashow_info:dict=None,
    qafile_info:dict=None):
    
    """
    To determine the velocity shift of a line using a model

    Parameters
    ----------

    data_wavelength : ndarray
        A (nwave1,) array of wavelengths for the data.

    data_fluxdensity : ndarray
        A (nwave1,) array of continuum-normalized flux densities for the data.

    model_wavelength : ndarray
        A (nwave2,) array of wavelengths for a model.  Must be the same units
        as `data_wavelength`.  

    model_fluxdensity : ndarray
        A (nwave2,) array of continuum-normalized flux densities for a model.

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
            A (latex) string giving the xlabel.

        `'plot_title'` : str, optional
            A (latex) string giving the title of the plot.

        
    qafile_info : dict, default=None
        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

        `'plot_xlabel'` : str, optional
            A (latex) string giving the xlabel.

        `'plot_title'` : str, optional
            A (latex) string giving the title of the plot.

    
    
    Returns
    -------
    ndarray
        The velocity shift of the data relative to the model in km s-1.  

           
    """

    #
    # Check Parameters
    #

    check_parameter('measure_linerv', 'data_wavelength', data_wavelength,
                    'ndarray', 1)

    check_parameter('measure_linerv', 'data_fluxdensity', data_fluxdensity,
                    'ndarray', 1)

    check_parameter('measure_linerv', 'model_wavelength', model_wavelength,
                    'ndarray', 1)

    check_parameter('measure_linerv', 'model_nflux', model_fluxdensity,
                    'ndarray', 1)
    
    check_parameter('measure_linerv', 'peak_method', peak_method,
                    'str', possible_values=['max','fit'])

    check_parameter('measure_linerv', 'resolving_power',
                    resolving_power, ['NoneType', 'int', 'float'])        
    
    check_parameter('measure_linerv', 'qashow_info', qashow_info,
                    ['NoneType', 'dict'])
    
    check_parameter('measure_linerv', 'qafile_info', qafile_info,
                    ['NoneType', 'dict'])

    #
    # Get set up
    #

    cspeed = 2.99792458E5 # km/s
    
    #
    # Check that the data has fewer pixels than the model
    #

    ndata = len(data_wavelength)
    nmodel = len(model_wavelength)

    if ndata > nmodel:

        message = 'Data has a higher sampling frequency than the model.'
        raise pySpextoolError(message)

    if ndata == 0:

        message = 'The wavelength array is empty.'
        raise pySpextoolError(message)

    #
    # Compute the EW scale factor
    #

    data_step = (np.max(data_wavelength)-np.min(data_wavelength))/(ndata-1)
    model_step = (np.max(model_wavelength)-np.min(model_wavelength))/\
                  (nmodel-1)

    #
    # Determine EWs of the absorption features in the data and the model and
    # set the EW scale factor
    #

    data_ew = (np.sum(data_step*(1.0-data_fluxdensity)))
    model_ew = (np.sum(model_step*(1.0-model_fluxdensity)))    
    scale_ew = data_ew/model_ew
    
    #
    # Now we can start the process
    #
    
    model_zflux = model_fluxdensity - 1
    data_zflux = data_fluxdensity - 1
            
    #
    # Resampling to a constant spacing in ln lambda: v/c = d(ln lambda)
    #

    # Create the wavelength array

    min_wavelength = np.min(np.concatenate((data_wavelength, model_wavelength)))
    max_wavelength = np.max(np.concatenate((data_wavelength, model_wavelength)))
        
    acoeff = (nmodel - 1)/(np.log(max_wavelength) -
                                     np.log(min_wavelength))
    bcoeff  = nmodel - (acoeff * (np.log(max_wavelength)))
  
    xpon = np.arange(nmodel)+1.0
    lnlambda_wavelengths = np.exp((xpon-bcoeff)/acoeff)

    # Do the resampling

    f = interp1d(data_wavelength, data_zflux, bounds_error=False, fill_value=0)
    data_resampled_zflux = f(lnlambda_wavelengths)

    f = interp1d(model_wavelength, model_zflux, bounds_error=False,
                 fill_value=0)
    model_resampled_zflux = f(lnlambda_wavelengths)
    
    #
    # Do the cross correlation
    #

    xcor = signal.correlate(data_resampled_zflux, model_resampled_zflux,
                                  mode='same', method='fft')
    xcor = xcor / np.nanmax(xcor)

    lag = signal.correlation_lags(nmodel, nmodel, mode='same')

    #
    # Find the peak of the cross correlation
    #

    if peak_method == 'max':

        max_idx = np.argmax(xcor)
        offset_pixels = int(lag[max_idx])
        fit = None
                
    if peak_method == 'fit':

    
        # Fit the cross correlation

        fit = fit_peak1d(lag, xcor, nparms=3, positive=True)
        offset_pixels = fit['parms'][1]
        fit = fit['fit']

    #
    # Compute the velocity shift
    #
        
    velocity_shift = (cspeed * offset_pixels)/acoeff
    redshift = velocity_shift/cspeed
    
    #
    # Make the QA plot
    #

    plotnum = None

    if qashow_info is not None:

        plot_measure_linerv(
            qashow_info['plot_number'],
            qashow_info['figure_size'],
            qashow_info['font_size'],
            qashow_info['spectrum_linewidth'],
            qashow_info['spine_linewidth'],
            lnlambda_wavelengths,
            data_resampled_zflux,
            model_resampled_zflux,
            lag,
            xcor,
            offset_pixels,
            velocity_shift,
            redshift,
            fit=fit,
            plot_xlabel=qashow_info['xlabel'],
            plot_title=qashow_info['title'])

        pl.show(block=qashow_info['block'])
        if qashow_info['block'] is False:

            pl.pause(1)

    if qafile_info is not None:

        plot_measure_linerv(
            None,
            qafile_info['figure_size'],
            qafile_info['font_size'],
            qafile_info['spectrum_linewidth'],
            qafile_info['spine_linewidth'],
            lnlambda_wavelengths,
            data_resampled_zflux,
            model_resampled_zflux,
            lag,
            xcor,
            offset_pixels,
            velocity_shift,
            redshift,
            fit=fit,
            plot_xlabel=qafile_info['xlabel'],
            plot_title=qafile_info['title'])

        pl.savefig(qafile_info['file_fullpath'])
        pl.close()
               
    return {'ew_scale':scale_ew, 'rv':velocity_shift, 'z':redshift,
            'plotnum':plotnum}


def modify_kuruczvega(
    standard_wavelength:npt.ArrayLike,
    standard_rv:float,
    standard_vmag:float,
    standard_bmag:float,
    vega_wavelength:npt.ArrayLike,
    vega_fluxdensity:npt.ArrayLike,
    vega_continuum:npt.ArrayLike,
    vega_fitted_continuum:npt.ArrayLike,
    kernel:npt.ArrayLike,
    vega_pixelshift:float,
    ew_info:dict,
    Rv:float=3.1,
    vega_vmag:float=0.03,
    vega_bminv:float=0.0,
    units:str='erg s-1 cm-2 A-1',
    return_raw_model:bool=False):

    """
    To modify the Kurucz Vega model for use in telluric correction.

    Parameters
    ----------
    standard_wavelength : ndarray
        An (nstd,) array of standard star wavelengths.

    standard_rv : float, int
        The radial velocities of the standard star in km s-1.

    standard_bmag : float, int
        The B-band magnitude of the standard star.

    standard_vmag : float, int
        The V-band magnitude of the standard star.
    
    vega_wavelength : ndarray
        An (nvega,) array of Vega-model wavelengths.

    vega_fluxdensity : ndarray
        An (nvega,) array of Vega-model flux densities.

    vega_continuum : ndarray
        An (nvega,) array of Vega-continuum flux densities.

    vega_fitted_continuum : ndarray
        An (nvega,) array of Vega-fitted-continuum flux densities.

    kernel : ndarray
        An (nkernel,) array of the kernel (in standard_wavelength pixels).

    ew_info : dict, default=None
        A dictionary giving the wavelengths and scale factors used to adjust the EWs
        of H lines in the standard star.  The information is used to compute an 
        (nvega,) array of scale factors.  Linear interpolation is used.


        'wavelengths' : ndarray
            An (nwaves,) array of wavelengths at which scale factors are known.

        'scales' : ndarray
                An (nwaves,) array of scale factors.
    
    Rv : float, default=3.1
        The R_V value used for reddening the Vega model.

    vega_vmag : float, default=0.03
        The Vega V-band magnitude

    vega_bminv : float, default=0.0
        The Vega B-V color.

    units : str

    return_raw_model : {False, True}
        Set to True to return only RV and convolution modifications.
        Set to False to also return modifications for color and brightness.

    Returns
    -------
    ndarray : An (nstd,) array with modified flux densities of the Vega model.

    ndarray : An (nstd,) array with modified continuum flux densities of the Vega model.

    """

    check_parameter('modify_kuruczvega', 'standard_wavelength',
                    standard_wavelength, 'ndarray', 1)

    check_parameter('modify_kuruczvega', 'standard_bmag',
                    standard_bmag, ['float','int'])

    check_parameter('modify_kuruczvega', 'standard_vmag',
                    standard_vmag, ['float','int'])

    check_parameter('modify_kuruczvega', 'standard_rv',
                    standard_rv, ['float64','float','int'])
    
    check_parameter('modify_kuruczvega', 'vega_wavelength',
                    vega_wavelength, 'ndarray', 1)

    check_parameter('modify_kuruczvega', 'vega_fluxdensity',
                    vega_fluxdensity, 'ndarray', 1)
    
    check_parameter('modify_kuruczvega', 'vega_continuum',
                    vega_continuum, 'ndarray', 1)

    check_parameter('modify_kuruczvega', 'vega_fitted_continuum',
                    vega_fitted_continuum, 'ndarray', 1)

    check_parameter('modify_kuruczvega', 'kernel', kernel, 'ndarray', 1)

    check_parameter('modify_kuruczvega', 'Rv', Rv, 'float')

    check_parameter('modify_kuruczvega', 'vega_vmag', vega_vmag, 'float')

    check_parameter('modify_kuruczvega', 'vega_bminv', vega_bminv, 'float')

    #
    # Make the two calls
    #

    setup = _modify_kuruczvega_setup(
        standard_wavelength,
        standard_rv,
        standard_vmag,
        standard_bmag,
        vega_wavelength,
        vega_fluxdensity,
        vega_continuum,
        vega_fitted_continuum,
        kernel,
        vega_pixelshift,
        Rv=Rv,
        vega_vmag=vega_vmag,
        vega_bminv=vega_bminv)


    f, c = _modify_kuruczvega_execute(
        setup,
        ew_info=ew_info,
        units=units,
        return_raw_model=return_raw_model)

    return f, c




def _modify_kuruczvega_setup(
    standard_wavelength:npt.ArrayLike,
    standard_rv:float,
    standard_vmag:float,
    standard_bmag:float,
    vega_wavelength:npt.ArrayLike,
    vega_fluxdensity:npt.ArrayLike,
    vega_continuum:npt.ArrayLike,
    vega_fitted_continuum:npt.ArrayLike,
    kernel:npt.ArrayLike,
    vega_pixelshift:float,
    Rv:float=3.1,
    vega_vmag:float=0.03,
    vega_bminv:float=0.0):

    """
    To obtain non-changing parts to pass to modify_kuruczvega_execute.

    Parameters
    ----------
    standard_wavelength : ndarray
        An (nstd,) array of data wavelengths.

    standard_rv : float, int
        The radial velocities of the standard star in km s-1.

    standard_bmag : float, int
        The B-band magnitude of the standard star.

    standard_vmag : float, int
        The V-band magnitude of the standard star.
    
    vega_wavelength : ndarray
        An (nvega,) array of Vega-model wavelengths.

    vega_fluxdensity : ndarray
        An (nvega,) array of Vega-model flux densities.

    vega_continuum : ndarray
        An (nvega,) array of Vega-continuum flux densities.

    vega_fitted_continuum : ndarray
        An (nvega,) array of Vega-fitted-continuum flux densities.

    kernel : ndarray
        An (nkernel,) array matched to the standard_wavelength array.

    vega_pixelshift : float
        The number of pixels to shift the Vega model after convolution.  

    Rv : float, default=3.1
        The R_V value used for reddening the Vega model.

    vega_vmag : float, default=0.03
        The Vega V-band magnitude

    vega_bminv : float, default=0.0
        The Vega B-V color.

    Returns
    -------
    dict : A dictionary of non-changing things.
    
        'standard_wavelength' : ndarray
            `standard_wavelength`

        'kernel' : ndarray
            `kernel`

        'ndata' : int
            The number of data points in `standard_wavelength`.

        'nkernel' : int
            The number of data points in `kernel`.

        'news' : int
            The number lines to adjust EWs.

        'standard_extwavelength' : ndarray
            `standard_wavelength` with extra wavelengths on either side to 
            accomodate the kernel when convolving.

        'rvcorrected_vega_wavelength' : ndarray
            `vega_wavelength` corrected for the RV of the standard.

        'vega_fluxdensity' : ndarray
            The Vega flux density array covering min/max of `standard_extwavelength`.

        'vega_continuum' : ndarray
            The Vega continuum array covering min/max of `standard_extwavelength`.

        'vega_fitted_continuum' : ndarray
            The Vega fitted continuum array covering min/max of 
            `standard_extwavelength`.

        'min_vega_wavelength' : ndarray
            The minimum of 'rvcorrected_vega_wavelength'
    
        'max_vega_wavelength' : ndarray
            The maximum of 'rvcorrected_vega_wavelength'

        'n_vega_wavelength' : ndarray
            The number of elements in 'rvcorrected_vega_wavelength'


        'vega_extinction_correction' : ndarray
             The extinction array on `standard_extwavelength`.

        'scale_fluxcalibration' : float
             A scale factor for flux calibration.
    
    """

    #
    # Check the parameters
    #
    
    check_parameter('_modify_kuruczvega_setup', 'standard_wavelength',
                    standard_wavelength, 'ndarray', 1)

    check_parameter('_modify_kuruczvega_setup', 'standard_bmag',
                    standard_bmag, ['float','int'])

    check_parameter('_modify_kuruczvega_setup', 'standard_vmag',
                    standard_vmag, ['float','int'])

    check_parameter('_modify_kuruczvega_setup', 'standard_rv',
                    standard_rv, ['float64','float','int'])
    
    check_parameter('_modify_kuruczvega_setup', 'vega_wavelength',
                    vega_wavelength, 'ndarray', 1)

    check_parameter('_modify_kuruczvega_setup', 'vega_fluxdensity',
                    vega_fluxdensity, 'ndarray', 1)
    
    check_parameter('_modify_kuruczvega_setup', 'vega_continuum',
                    vega_continuum, 'ndarray', 1)

    check_parameter('modify_kuruczvega_setup', 'vega_fitted_continuum',
                    vega_fitted_continuum, 'ndarray', 1)

    check_parameter('modify_kuruczvega_setup', 'vega_pixelshift',
                    vega_pixelshift, 'float')

    check_parameter('_modify_kuruczvega_setup', 'kernel', 
                    kernel, 'ndarray', 1)

    check_parameter('_modify_kuruczvega', 'Rv', 
                    Rv, 'float')

    check_parameter('_modify_kuruczvega', 'vega_vmag', 
                    vega_vmag, 'float')

    check_parameter('_modify_kuruczvega', 'vega_bminv', 
                    vega_bminv, 'float')
    
    #
    # Fill in the setup dictionary that gets passed to _modify_kuruczvega_execute
    #

    setup = {}

    setup['standard_wavelength'] = standard_wavelength
    setup['kernel'] = kernel
    setup['ndata'] = np.size(standard_wavelength)
    setup['nkernel'] = np.size(kernel)
    setup['vega_pixelshift'] = vega_pixelshift
        
    #
    # Extend the `setup['standard_wavelength']` to account for the width of the kernel.
    #

    setup['kernel_wavelengths'] = get_kernelwavelengths(
        standard_wavelength,
        setup['nkernel'])

    #
    # Shift to the radial velocity of the standard
    #

    rvcorrected_vega_wavelength = vega_wavelength*(1+standard_rv/2.99792458e5)

    #
    # Clip the Vega model to the wavelengths of `setup['standard_extwavelength']`
    #

    min = np.nanmin(setup['kernel_wavelengths'])
    max = np.nanmax(setup['kernel_wavelengths'])
    
    zleft = np.where(rvcorrected_vega_wavelength > min)
    zleft_idx = zleft[0][0]
    
    zright = np.where(rvcorrected_vega_wavelength < max)
    zright_idx = zright[0][-1]

    z_idx = np.arange(zleft_idx,zright_idx,dtype=int)

    setup['rvcorrected_vega_wavelength'] = rvcorrected_vega_wavelength[z_idx]
    setup['vega_fluxdensity'] = vega_fluxdensity[z_idx]
    setup['vega_continuum'] = vega_continuum[z_idx]
    setup['vega_fitted_continuum'] = vega_fitted_continuum[z_idx]
    setup['min_vega_wavelength'] = np.min(setup['rvcorrected_vega_wavelength'])
    setup['max_vega_wavelength'] = np.max(setup['rvcorrected_vega_wavelength'])
    setup['n_vega_wavelength'] = np.size(setup['rvcorrected_vega_wavelength'])

    #
    # Redden the model to match that of the standard
    #
    
    ebminv = (standard_bmag-standard_vmag - vega_bminv)
    ebminv = np.max([ebminv,0])  # ensure >= 0
    Av = Rv*ebminv
    ext = G23(Rv=Rv)
    
    setup['vega_extinction_correction'] = \
        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
    
    #
    # Scale the model to the observed magnitude of the standard
    #

    setup['scale_fluxcalibration'] = 10.**(-0.4*(standard_vmag-vega_vmag))*10**(0.4*Av)

    return setup
    



def _modify_kuruczvega_execute(
    s,
    ew_info:dict=None,
    units:str='erg s-1 cm-2 A-1',
    return_raw_model:bool=False):

    """
    To modify the Kurucz Vega model...

    Parameters
    ----------
    s : dict    
        'standard_wavelength' : ndarray
            `standard_wavelength`

        'kernel' : ndarray
            `kernel`

        'ndata' : int
            The number of data points in `standard_wavelength`.

        'nkernel' : int
            The number of data points in `kernel`.

        'news' : int
            The number lines to adjust EWs.

        'standard_extwavelength' : ndarray
            `standard_wavelength` with extra wavelengths on either side to 
            accomodate the kernel when convolving.

        'rvcorrected_vega_wavelength' : ndarray
            `vega_wavelength` corrected for the RV of the standard.

        'vega_fluxdensity' : ndarray
            The Vega flux density array covering min/max of `standard_extwavelength`.

        'vega_continuum' : ndarray
            The Vega continuum array covering min/max of `standard_extwavelength`.

        'vega_fitted_continuum' : ndarray
            The Vega fitted continuum array covering min/max of 
            `standard_extwavelength`.

        'min_vega_wavelength' : ndarray
            The minimum of 'rvcorrected_vega_wavelength'
    
        'max_vega_wavelength' : ndarray
            The maximum of 'rvcorrected_vega_wavelength'

        'n_vega_wavelength' : ndarray
            The number of elements in 'rvcorrected_vega_wavelength'

        'vega_extinction_correction' : ndarray
             The extinction array on `standard_extwavelength`.

        'scale_fluxcalibration' : float
             A scale factor for flux calibration.

        A dictionary of non-changing things that gets generated by 
        _modify_kuruczvega_setup

    ew_info : dict, default=None
        'wavelengths' : ndarray
            An (nwaves,) array of wavelengths at which scale factors are known.

        'scales' : ndarray
                An (nwaves,) array of scale factors.

    units : str

    raw_model : {False, True}
        Set to True to return only RV, convolution, and pixel shift modifications.
        Set to False to also return modifications for color and brightness.

    Returns
    -------
    ndarray : s[`standard_wavelength`]
    
    ndarray : The modified flux density of the Vega model.

    ndarray : The modified continuum density of the Vega model.

    """

    #
    # Generate the EW scale factor array
    #

#    start_time = time.perf_counter()
    if ew_info is not None:

        # User passed values.  Interpolate onto the vega wavelength grid

        f = scipy.interpolate.interp1d(
            ew_info['wavelengths'],
            ew_info['scales'], 
            fill_value=(ew_info['scales'][0], ew_info['scales'][-1]),
            bounds_error=False)
        vega_scales = f(s['rvcorrected_vega_wavelength'])

    else:

        vega_scales = np.full(s['n_vega_wavelength'],1)

#    end_time = time.perf_counter()
#    execution_time = end_time - start_time
#    print(f"Execution time: {execution_time:.6f} seconds")

    # 
    # Do the convolution
    #

#    start_time = time.perf_counter()

    flipped_kernel = np.flip(s['kernel'])

    c_n_z_vega_fd = np.full(s['ndata'],np.nan) # convolved, normalized, zeroed
    c_n_z_vega_cont = np.full(s['ndata'],np.nan)

    for i in range(s['ndata']):

        if np.isnan(s['standard_wavelength'][i]):

            continue
    
        kernel_wavelength = s['kernel_wavelengths'][:,i]

        rkernel = np.interp(
            s['rvcorrected_vega_wavelength'],
            kernel_wavelength,
            flipped_kernel,
            left=0.0,
            right=0.0)

        rkernel /= np.nansum(rkernel)

        z_ratio = s['vega_fluxdensity']/s['vega_continuum']-1
        c_n_z_vega_fd[i] = np.nansum(rkernel*z_ratio*vega_scales)

        z_ratio = s['vega_continuum']/s['vega_fitted_continuum']-1
        c_n_z_vega_cont[i] = np.nansum(rkernel*z_ratio*vega_scales)

#    end_time = time.perf_counter()
#    execution_time = end_time - start_time
#    print(f"Execution time: {execution_time:.6f} seconds")


#    start_time = time.perf_counter()

    # Reconstruct the flux density and continuum
    
    vega_fitted_rcontinuum = np.interp(
        s['standard_wavelength'],
        s['rvcorrected_vega_wavelength'],
        s['vega_fitted_continuum'],
        left=0.0,
        right=0.0)

    c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_rcontinuum
    
    c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont

    if s['vega_pixelshift'] != 0:

        # We use the spline (k=1 -> linear interp) because it allows us to 
        # extrapolate when we do the pixel shift which is not accounted for 
        # when we determine the ranges used to cut the Vega model

        x = np.arange(np.size(s['standard_wavelength']))

        func = make_interp_spline(x+s['vega_pixelshift'],c_vega_fd,k=1)
        c_vega_fd = func(x)

        func = make_interp_spline(x+s['vega_pixelshift'],c_vega_cont,k=1)
        c_vega_cont = func(x)

    if return_raw_model:

        return c_vega_fd, c_vega_cont

    #
    # Redden the model to match that of the standard
    #

    e_c_vega_fd = c_vega_fd*s['vega_extinction_correction']

    e_c_vega_cont = c_vega_cont*s['vega_extinction_correction']

    #
    # Scale the matched observed brightness
    #

    s_e_c_vega_fd = e_c_vega_fd*s['scale_fluxcalibration']
    s_e_c_vega_cont = e_c_vega_cont*s['scale_fluxcalibration']

    #
    # Convert the model to the requested flux density units
    #

    vega_fd = convert_fluxdensity(
        s['standard_wavelength'],
        s_e_c_vega_fd,
        'um','erg s-1 cm-2 A-1',
        units)

    vega_cont = convert_fluxdensity(
        s['standard_wavelength'],
        s_e_c_vega_cont,
        'um','erg s-1 cm-2 A-1',
        units)

    #
    # Returm the results
    #
#    end_time = time.perf_counter()
#    execution_time = end_time - start_time
#    print(f"Execution time: {execution_time:.6f} seconds")

    return vega_fd, vega_cont




#def modify_kuruczvega2(
#    standard_wavelength:npt.ArrayLike,
#    standard_rv:float,
#    standard_vmag:float,
#    standard_bmag:float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    vega_ewscalefactors:float | npt.ArrayLike,
#    units:str='erg s-1 cm-2 A-1'):
#
#    """
#    To modify the Kurucz Vega model to match the standard star spectrum
#
#    Parameters
#    ----------
#    standard_wavelength : ndarray
#        A (nwave1,) array of data wavelengths.
#
#    standard_rv : float, int
#        The radial velocities of the standard star in km s-1.
#
#    standard_bmag : float, int
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float, int
#        The V-band magnitude of the standard star.
#    
#    vega_wavelength : ndarray
#        A (nwave2,) array of Vega-model wavelengths.
#
#    vega_fluxdensity : ndarray
#        A (nwave2,) array of Vega-model flux densities.
#
#    vega_continuum : ndarray
#        A (nwave2,) array of Vega-continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        A (nwave2,) array of Vega-fitted-continuum flux densities.
#
#    vega_ewscalefactors : ndarray or float
#        A (nwave2,) array of EW scale factors or a float giving a 
#        single scale factor.
#
#    units : str
#        The requested flux density units.
#    
#    Returns
#    -------
#    ndarray
#
#        The modified Vega flux density
#
#    ndarray
#
#        The modified Vega continiuum
#
#    """
#
#    #
#    # Check the parameters
#    #
#    
#    check_parameter('modify_kuruczvega', 'standard_wavelength',
#                    standard_wavelength, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'standard_bmag',
#                    standard_bmag, ['float','int'])
#
#    check_parameter('modify_kuruczvega', 'standard_vmag',
#                    standard_vmag, ['float','int'])
#
#    check_parameter('modify_kuruczvega', 'standard_rv',
#                    standard_rv, ['float64','float','int'])
#    
#    check_parameter('modify_kuruczvega', 'vega_wavelength',
#                    vega_wavelength, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_fluxdensity',
#                    vega_fluxdensity, 'ndarray', 1)
#    
#    check_parameter('modify_kuruczvega', 'vega_continuum',
#                    vega_continuum, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_fitted_continuum',
#                    vega_fitted_continuum, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'kernel', kernel, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_ewscalefactors',
#                    vega_ewscalefactors, ['float','ndarray'])
#
#    check_parameter('modify_kuruczvega', 'units', units, 'str')
#
#    
#    #
#    # Get set up 
#    #
#
#    ndata = np.size(standard_wavelength)
#    nkernel = np.size(kernel)
#
#    # Extend the wavelength coverage of the data to account for the kernel
#
#    standard_extwavelength = extend_wavelengths(
#        standard_wavelength,
#        int(np.floor(nkernel/2)),
#        int(np.floor(nkernel/2)))
#
#    #
#    # Clip the Vega model to avoid interpolating over the entire thing.
#    #
#
#    min = np.nanmin(standard_extwavelength)
#    max = np.nanmax(standard_extwavelength)
#    
#    zleft = np.where(vega_wavelength > min)
#    zleft_idx = zleft[0][0]
#    
#    zright = np.where(vega_wavelength < max)
#    zright_idx = zright[0][-1]
#
#    z_idx = np.arange(zleft_idx,zright_idx,dtype=int)
#   
#    svega_wavelength = vega_wavelength[z_idx]
#    svega_fluxdensity = vega_fluxdensity[z_idx]
#    svega_continuum = vega_continuum[z_idx]
#    svega_fitted_continuum = vega_fitted_continuum[z_idx]
#    svega_ewscalefactors = vega_ewscalefactors[z_idx]
#    
#    #
#    # Do the convolutions
#    #
#
#    # Convolve the flux density array
#
#    c_n_z_vega_fd = np.empty(ndata) # convolved, normalized, zeroed
#    c_n_z_vega_cont = np.empty(ndata)
#    for i in range(ndata):
#
#        kernel_wavelength = standard_extwavelength[0+i:nkernel+i]
#
#        rkernel = linear_interp1d(
#            kernel_wavelength,
#            np.flip(kernel),
#            svega_wavelength)
#
#        rkernel /= np.nansum(rkernel)
#
#
#        c_n_z_vega_fd[i] = np.nansum(
#            rkernel*(svega_fluxdensity/svega_continuum-1)*svega_ewscalefactors)
#
#        c_n_z_vega_cont[i] = np.nansum(
#            rkernel*(svega_continuum/svega_fitted_continuum-1)*svega_ewscalefactors)
#
#
#    # Reconstruct the flux density and continuum
#    
#    vega_fitted_rcontinuum = linear_interp1d(
#        svega_wavelength,
#        svega_fitted_continuum,
#        standard_wavelength)
#
#    c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_rcontinuum
#    
#    c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont
#
#
#
#        
##    pl.figure()
##    pl.step(standard_wavelength, c_n_z_vega_fd,where='mid')
##    pl.step(vega_wavelength, vega_fluxdensity/vega_continuum-1,where='mid')
##    pl.xlim(1,2)
##    pl.ylim(-0.5,0.1)
##
##    pl.show()
##    return 1, 1
#
##    # Reconstruct the flux density and continuum
##    
##    vega_fitted_rcontinuum = linear_interp1d(
##        vega_wavelength,
##        vega_fitted_continuum,
##        standard_wavelength)
##
##    c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_rcontinuum
##    
##    c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont
##
##    return c_vega_fd, c_vega_cont
##
#
#
#
##    c_n_z_vega_fd = np.empty(ndata) # convolved, normalized, zeroed
##    for i in range(ndata):
##
##        kernel_wavelength = standard_extwavelength[0+i:nkernel+i]
##
##        model_rintensity = linear_interp1d(
##            vega_wavelength,
##            (vega_fluxdensity/vega_continuum-1)*vega_ewscalefactors,
##            kernel_wavelength)
##
##        c_n_z_vega_fd[i] = np.sum(model_rintensity*np.flip(kernel))
##
##
##    pl.figure()
##    pl.step(standard_wavelength, c_n_z_vega_fd,where='mid')
##    pl.step(vega_wavelength, vega_fluxdensity/vega_continuum-1,where='mid')
###    pl.step(vega_wavelength, vega_fluxdensity/vega_continuum,where='mid')
##    pl.xlim(1,2)
##    pl.ylim(-0.5,0.1)
##    pl.show()
##    return 1, 1
##
##
##
##
##
##    # Convolve the continuum array
##
##    c_n_z_vega_cont = np.empty(ndata)
##    for i in range(ndata):
##
##        kernel_wavelength = standard_extwavelength[0+i:nkernel+i]
##
##        model_rintensity = linear_interp1d(
##            vega_wavelength,
##            (vega_continuum/vega_fitted_continuum-1)*vega_ewscalefactors,
##            kernel_wavelength)
##
##        c_n_z_vega_cont[i] = np.sum(model_rintensity*np.flip(kernel))
##
##    # Reconstruct the flux density and continuum
##    
##    vega_fitted_rcontinuum = linear_interp1d(
##        vega_wavelength,
##        vega_fitted_continuum,
##        standard_wavelength)
##
##    c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_rcontinuum
##    
##    c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont
#
##    return c_vega_fd, c_vega_cont
#        
#    #
#    # Shift to the radial velocity of the standard
#    #
#
#    shifted_vega_wavelength = standard_wavelength*(1+standard_rv/2.99792458e5)
#
#    #
#    # Redden the model to match that of the standard
#    #
#
#    Rv=3.1
#    vega_vmag = 0.03
#    vega_bminv = 0.0
#    
#    ebminv = (standard_bmag-standard_vmag - vega_bminv)
#    ebminv = np.max([ebminv,0])
#    Av = Rv*ebminv
#    
#    #
#    # Extinguish the model to match that of the standard
#    #
#
#    ext = G23(Rv=Rv)
#    
#    e_c_vega_fd = c_vega_fd*\
#        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
#
#    e_c_vega_cont = c_vega_cont*\
#        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
#    
#    #
#    # Scale the model to the observed magnitude of the standard
#    #
#
#    scale = 10.**(-0.4*(standard_vmag-vega_vmag))*10**(0.4*Av)
#
#    s_e_c_vega_fd = e_c_vega_fd*scale
#    s_e_c_vega_cont = e_c_vega_cont*scale
#
#    #
#    # Convert the model to the requested flux density units
#    #
#
#    vega_fd = convert_fluxdensity(
#        standard_wavelength,
#        s_e_c_vega_fd,
#        'um','erg s-1 cm-2 A-1',
#        units)
#
#    vega_cont = convert_fluxdensity(
#        standard_wavelength,
#        s_e_c_vega_cont,
#        'um','erg s-1 cm-2 A-1',
#        units)
#
#    #
#    # Returm the results
#    #
#
#    return vega_fd, vega_cont
#    
#
#
#
#def modify_kuruczvega(
#    standard_wavelength:npt.ArrayLike,
#    standard_rv:float,
#    standard_vmag:float,
#    standard_bmag:float,
#    vega_wavelength:npt.ArrayLike,
#    vega_fluxdensity:npt.ArrayLike,
#    vega_continuum:npt.ArrayLike,
#    vega_fitted_continuum:npt.ArrayLike,
#    kernel:npt.ArrayLike,
#    vega_ewscalefactors:float | npt.ArrayLike,
#    units:str='erg s-1 cm-2 A-1',
#    new=False):
#
#    """
#    To modify the Kurucz Vega model to match the standard star spectrum
#
#    Parameters
#    ----------
#    standard_wavelength : ndarray
#        A (nwave1,) array of data wavelengths.
#
#    standard_rv : float, int
#        The radial velocities of the standard star in km s-1.
#
#    standard_bmag : float, int
#        The B-band magnitude of the standard star.
#
#    standard_vmag : float, int
#        The V-band magnitude of the standard star.
#    
#    vega_wavelength : ndarray
#        A (nwave2,) array of Vega-model wavelengths.
#
#    vega_fluxdensity : ndarray
#        A (nwave2,) array of Vega-model flux densities.
#
#    vega_continuum : ndarray
#        A (nwave2,) array of Vega-continuum flux densities.
#
#    vega_fitted_continuum : ndarray
#        A (nwave2,) array of Vega-fitted-continuum flux densities.
#
#    vega_ewscalefactors : ndarray or float
#        A (nwave2,) array of EW scale factors or a float giving a 
#        single scale factor.
#
#    units : str
#        The requested flux density units.
#    
#    Returns
#    -------
#    ndarray
#
#        The modified Vega flux density
#
#    ndarray
#
#        The modified Vega continiuum
#
#    """
#
#    #
#    # Check the parameters
#    #
#    
#    check_parameter('modify_kuruczvega', 'standard_wavelength',
#                    standard_wavelength, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'standard_bmag',
#                    standard_bmag, ['float','int'])
#
#    check_parameter('modify_kuruczvega', 'standard_vmag',
#                    standard_vmag, ['float','int'])
#
#    check_parameter('modify_kuruczvega', 'standard_rv',
#                    standard_rv, ['float64','float','int'])
#    
#    check_parameter('modify_kuruczvega', 'vega_wavelength',
#                    vega_wavelength, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_fluxdensity',
#                    vega_fluxdensity, 'ndarray', 1)
#    
#    check_parameter('modify_kuruczvega', 'vega_continuum',
#                    vega_continuum, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_fitted_continuum',
#                    vega_fitted_continuum, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'kernel', kernel, 'ndarray', 1)
#
#    check_parameter('modify_kuruczvega', 'vega_ewscalefactors',
#                    vega_ewscalefactors, ['float','ndarray'])
#
#    check_parameter('modify_kuruczvega', 'units', units, 'str')
#
#    #
#    # Get the Vega pixels that exactly cover the order
#    #
#
#    min = np.nanmin(standard_wavelength)
#    max = np.nanmax(standard_wavelength)
#    
#    zleft = np.where(vega_wavelength > min)
#    zleft_idx = zleft[0][0]
#    
#    zright = np.where(vega_wavelength < max)
#    zright_idx = zright[0][-1]
#
#    # Now figure out how many pixels you need to add to both sides to 
#    # account for the width of the kernel during the convolution
#
#    nadd = np.ceil(len(kernel)/2.)
#
#    z_idx = np.arange(zleft_idx-nadd,zright_idx+nadd+1,dtype=int)
#    
#    #
#    # Do the convolution on the normalized, zeroed, and EW scaled spectrum 
#    # and continuum
#    #
#
#    scale = vega_ewscalefactors
#    if new is True:
#
#        input = (vega_fluxdensity[z_idx]/vega_fitted_continuum[z_idx]-1)*scale[z_idx]
#
#        input = (input+1)*vega_fitted_continuum[z_idx]
#
#        c_vega_fd = np.convolve(input, kernel, mode='same')
#
#        input = vega_fitted_continuum[z_idx]
#        c_vega_cont = np.convolve(input, kernel, mode='same')
##        
##        # Reconstruct the flux density and continuum
##        
##        c_vega_fd = (c_n_z_vega_fd+1)*vega_fitted_continuum[z_idx]
##        c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont
#
#
##        raw_vega_conv = np.convolve(vega_fluxdensity[z_idx], kernel, mode='same')
#
#
#
##        pl.figure()
##        pl.plot(vega_wavelength, vega_fluxdensity, color='red', 
##                label='Raw Vega')
##        pl.plot(vega_wavelength[z_idx], raw_vega_conv, color='black', 
##                label='Convolved Vega')
##
##        pl.plot(vega_wavelength[z_idx], c_vega_fd,label='Scaled, Convolved Vega')
##        pl.xlabel('Wavelength ($\mu$m)')
##        pl.xlim(np.min(vega_wavelength[z_idx]), 
##                  np.max(vega_wavelength[z_idx]))
##
##        pl.ylim(np.min(c_vega_fd), np.max(c_vega_fd))
##
##        pl.title('Order'+str(order))
##        pl.legend()
##        pl.savefig('Order '+str(order)+'.pdf')
#
#
#    
#    if new is False:
#
#        input = (vega_fluxdensity[z_idx]/vega_continuum[z_idx]-1)*scale[z_idx]
##\
##            vega_ewscalefactors[z_idx]
#        c_n_z_vega_fd = np.convolve(input, kernel, mode='same')
#        
#        input = (vega_continuum[z_idx]/vega_fitted_continuum[z_idx]-1)*scale[z_idx]
##\
##            vega_ewscalefactors[z_idx]
#        c_n_z_vega_cont = np.convolve(input, kernel, mode='same')
#        
#        # Reconstruct the flux density and continuum
#        
#        c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_continuum[z_idx]
#        
#        c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont
#
#
#
#
##        raw_vega_conv = np.convolve(vega_fluxdensity[z_idx], kernel, mode='same')
#
#
#
##        pl.figure()
##
##
##        pl.plot(vega_wavelength, vega_fluxdensity, color='red', 
##                label='Raw Vega')
##        pl.plot(vega_wavelength[z_idx], raw_vega_conv, color='black', 
##                label='Convolved Vega')
##
##        pl.plot(vega_wavelength[z_idx], c_vega_fd,label='Scaled, Convolved Vega')
##
##        pl.xlim(np.min(vega_wavelength[z_idx]), 
##                  np.max(vega_wavelength[z_idx]))
##        pl.ylim(np.min(c_vega_fd), np.max(c_vega_fd))
##        pl.title('Order'+str(order))
##        pl.legend()
##        pl.savefig('Order '+str(order)+'.pdf')
#
#
#        
#    #
#    # Shift to the radial velocity of the standard
#    #
#
#    shifted_vega_wavelength = vega_wavelength*(1+standard_rv/2.99792458e5)
#
#    #
#    # Interpolate onto the wavelength grid of the standard
#    #
#
#    r_c_vega_fd = linear_interp1d(shifted_vega_wavelength[z_idx],
#                                  c_vega_fd, standard_wavelength)
#
#    r_c_vega_cont = linear_interp1d(shifted_vega_wavelength[z_idx],
#                                    c_vega_cont, standard_wavelength)
#
#    return r_c_vega_fd, r_c_vega_cont
#    
#    #
#    # Redden the model to match that of the standard
#    #
#
#    Rv=3.1
#    vega_vmag = 0.03
#    vega_bminv = 0.0
#    
#    ebminv = (standard_bmag-standard_vmag - vega_bminv)
#    ebminv = np.max([ebminv,0])
#    Av = Rv*ebminv
#    
#    #
#    # Extinguish the model to match that of the standard
#    #
#
#    ext = G23(Rv=Rv)
#    
#    e_r_c_vega_fd = r_c_vega_fd*\
#        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
#
#    e_r_c_vega_cont = r_c_vega_cont*\
#        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
#    
#    #
#    # Scale the model to the observed magnitude of the standard
#    #
#
#    scale = 10.**(-0.4*(standard_vmag-vega_vmag))*10**(0.4*Av)
#
#    s_e_r_c_vega_fd = e_r_c_vega_fd*scale
#    s_e_r_c_vega_cont = e_r_c_vega_cont*scale
#
#    #
#    # Convert the model to the requested flux density units
#    #
#
#    vega_fd = convert_fluxdensity(
#        standard_wavelength,
#        s_e_r_c_vega_fd,
#        'um','erg s-1 cm-2 A-1',
#        units)
#
#    vega_cont = convert_fluxdensity(
#        standard_wavelength,
#        s_e_r_c_vega_cont,
#        'um','erg s-1 cm-2 A-1',
#        units)
#
#    #
#    # Returm the results
#    #
#
#    return vega_fd, vega_cont
    

