import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)
import scipy


import logging
from scipy import signal
from scipy.interpolate import CubicSpline
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.special import erf
from dust_extinction.parameter_averages import G23
import astropy.units as u

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.fit.polyfit import polyfit_1d, poly_1d
from pyspextool.utils.math import moments
from pyspextool.plot.limits import get_spectra_range
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.units import convert_fluxdensity
from pyspextool.utils.arrays import find_index



def deconvolve_line(data_wavelength:npt.ArrayLike,
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
    logging.info(' EW Scale Factor= '+str(scale_ew))
    
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

    logging.info(' RMS deviation '+str(rms_deviation))        
    logging.info(' Maximum deviation '+str(maximum_deviation))    

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

        plot_deconvolve_line(qashow_info['plot_number'],
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

        plot_deconvolve_line(None,
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


    return dict



def estimate_ewscales(standard_wavelength:npt.ArrayLike,
                      standard_fluxdensity:npt.ArrayLike,
                      standard_rv:int |float,
                      standard_vmag:int | float,
                      standard_bmag:int | float,
                      vega_wavelength:npt.ArrayLike,
                      vega_fluxdensity:npt.ArrayLike,
                      vega_continuum:npt.ArrayLike,
                      vega_fitted_continuum:npt.ArrayLike,
                      kernel:npt.ArrayLike,
                      ew_scale:float | int,
                      line_wavelengths:npt.ArrayLike,
                      atmospheric_transmission:npt.ArrayLike,
                      poly_degree:int,
                      tolerance:float=0.01,
                      include_edges:bool=False,
                      qashow_info:dict=None,
                      qafile_info:dict=None):


    """
    To estimate the H line EW scale factors .

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

    ew_scale : float
        The EW scale factor.

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

    check_parameter('estimate_ewscales', 'standard_rv', standard_rv, 
                    'float')

    check_parameter('estimate_ewscales', 'standard_bmag', standard_bmag, 
                    'float')

    check_parameter('estimate_ewscales', 'standard_vmag', standard_vmag, 
                    'float')

    check_parameter('estimate_ewscales', 'vega_wavelength', 
                    vega_wavelength, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_fluxdensity', 
                    vega_fluxdensity, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_continuum', 
                    vega_continuum, 'ndarray')

    check_parameter('estimate_ewscales', 'vega_fitted_continuum', 
                    vega_fitted_continuum, 'ndarray')

    check_parameter('estimate_ewscales', 'kernel', kernel, 'ndarray')

    check_parameter('estimate_ewscales', 'line_wavlenegths', line_wavelengths, 
                    ['float', 'ndarray'])

    check_parameter('estimate_ewscales', 'ew_scale', ew_scale, 'float')

    check_parameter('estimate_ewscales', 'poly_degree', poly_degree, 
                    'int')

    check_parameter('estimate_ewscales', 'tolerance', tolerance, 
                    'float')

    check_parameter('estimate_ewscales', 'include_edges', include_edges, 
                    'bool')

    check_parameter('adjust_ews', 'qashow_info', qashow_info, 
                    ['NoneType','dict'])

    check_parameter('adjust_ews', 'qafile_info', qafile_info, 
                    ['NoneType','dict'])


    #
    # Normalize the data
    #

    standard_fluxdensity /= np.nanmedian(standard_fluxdensity)
    
    #
    # Create a default scales array
    #

    if include_edges is True:
        
        points = np.append(np.insert(line_wavelengths, 0, 
                                     np.min(standard_wavelength)),
                           np.max(standard_wavelength))
        
    else:

        points = np.array(line_wavelengths)

    default_scales = np.full_like(points, ew_scale) 

#    print(default_scales)
#    print(default_scales[0])

    #
    # Create the scales array on the Vega wavelength grid
    #

    f = scipy.interpolate.interp1d(points, default_scales, 
                                   fill_value=(default_scales[0], 
                                               default_scales[-1]),
                                   bounds_error=False)
    vega_scales = f(vega_wavelength)

    #
    # Create default Vega model
    #

    result = modify_kuruczvega(standard_wavelength,
                               standard_rv,
                               standard_vmag,
                               standard_bmag,
                               vega_wavelength,
                               vega_fluxdensity,
                               vega_continuum,
                               vega_fitted_continuum,
                               kernel,
                               vega_scales,
                               units='erg s-1 cm-2 A-1',
                               new=False)

    # Normalize the vega model

    fd = scipy.interpolate.interp1d(standard_wavelength, standard_fluxdensity)
    fr = scipy.interpolate.interp1d(standard_wavelength, result[0])
    fv = scipy.interpolate.interp1d(vega_wavelength, vega_fluxdensity)

    default_vegamodel = result[0]/np.median(result[0])

    #
    # Estimate model parameters p0.
    # p0[0] = atmospheric scale factor
    # p0[1:poly_degree+2] = coeffs
    # p0[poly_degree+2:] = EW scale factors
    #

    default_ratio = standard_fluxdensity/default_vegamodel


    coeffs = polyfit_1d(standard_wavelength, default_ratio, 
                        poly_degree)['coeffs']

    p0 = np.insert(coeffs, 0, 1)

    p0 = np.append(p0, default_scales)

    #
    # Run the fit
    #

    # Build the arg array to pass to the optimizer
    
    args = (standard_wavelength,
            standard_fluxdensity,
            standard_rv,
            standard_vmag,
            standard_bmag,
            vega_wavelength,
            vega_fluxdensity,
            vega_continuum,
            vega_fitted_continuum,
            kernel,
            ew_scale,
            line_wavelengths,
            atmospheric_transmission,
            poly_degree,
            include_edges)


    # Do the fit

    result = scipy.optimize.minimize(ewscale_objfunction, 
                                     p0, 
                                     args=args, 
                                     tol=tolerance)

    # upack the results

    p = result.x

    idx = poly_degree+2
    atm_scale = p[0]
    coeffs = p[1:idx]
    optimized_scales = p[idx:]

    #
    # plot the results
    #

    # Create a new vega model and ratio

    if include_edges is True:
                
        scales = np.append(np.insert(optimized_scales, 0, ew_scale),
                              ew_scale)
    else:

        scales = optimized_scales

    # Create the new vega_scales array

    f = scipy.interpolate.interp1d(points, scales,
                                   fill_value=(scales[0], 
                                               scales[-1]),
                                   bounds_error=False)

    vega_scales = f(vega_wavelength)

    # Create the new Vega model

    result = modify_kuruczvega(standard_wavelength,
                               standard_rv,
                               standard_vmag,
                               standard_bmag,
                               vega_wavelength,
                               vega_fluxdensity,
                               vega_continuum,
                               vega_fitted_continuum,
                               kernel,
                               vega_scales,
                               units='erg s-1 cm-2 A-1',
                               new=False)

    # Create the new data model

    new_vegamodel = result[0]/np.median(result[0])

    atmosphere = (atm_scale*(atmospheric_transmission-1))+1

    model = new_vegamodel*atmosphere*poly_1d(standard_wavelength,coeffs)

    new_ratio = standard_fluxdensity/new_vegamodel

    #
    # Do the QA plotting
    #

    if qashow_info is not None:

        plot_estimate_ewscales(qashow_info['plot_number'],
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
                               line_wavelengths,
                               optimized_scales,
                               ew_scale)
        
        pl.show(block=qashow_info['block'])
        if qashow_info['block'] is False:

            pl.pause(1)

    if qafile_info is not None:

        plot_estimate_ewscales(None,
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
                               line_wavelengths,
                               optimized_scales,
                               ew_scale)
               
        pl.savefig(qafile_info['file_fullpath'])
        pl.close()
    
    return atm_scale, coeffs, optimized_scales




def ewscale_objfunction(p:npt.ArrayLike, 
                        standard_wavelength:npt.ArrayLike,
                        standard_fluxdensity:npt.ArrayLike,
                        standard_rv:int | float,
                        standard_vmag:int | float,
                        standard_bmag:int | float,
                        vega_wavelength:npt.ArrayLike,
                        vega_fluxdensity:npt.ArrayLike,
                        vega_continuum:npt.ArrayLike,
                        vega_fitted_continuum:npt.ArrayLike,
                        kernel:npt.ArrayLike,
                        ew_scale:int | float,
                        line_wavelengths:npt.ArrayLike,
                        atmospheric_transmission:npt.ArrayLike,
                        poly_degree:int,
                        include_edges:bool):

    """
    The objective function used to determine the H line EW scale factors.

    The function is used along with scipy.optmize.minimize in order to 
    determine optimal model parameters of a standard star by minimizing 
    this function.

    Parameters
    ----------
    p : ndarray
        A array of model parameters.  The length depends on the number of 
        hydrogen lines whose EWs need adjustment and the polynomial degree 
        used to model the instrument throughput.  

       atm_scale = p[0]
       coeffs = p[1:`poly_degree`+2]
       ew_scales = p[`poly_degree`+2:]

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

    start_ewscale : float
        The EW scale factor.

    poly_degree : int
        The polynomial degree used to model the instrument throughput.

    include_edges : {False, True}
        Set to True to include the min() and max() values of 
        `standard_wavelength` as control points.
        Set to False to not include the min() and max() values of 
        `standard_wavelength` as control points.

    Returns
    -------
    float 
    
        The summed squared residuals between the data and model.  

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
    # Check parameters
    #

    check_parameter('ewscale_objfunction', 'p', p, 'ndarray')

    check_parameter('ewscale_objfunction', 'standard_wavelength', 
                    standard_wavelength, 'ndarray')

    check_parameter('ewscale_objfunction', 'standard_fluxdensity', 
                    standard_fluxdensity, 'ndarray')

    check_parameter('ewscale_objfunction', 'standard_rv', standard_rv, 
                    'float')

    check_parameter('ewscale_objfunction', 'standard_bmag', standard_bmag, 
                    'float')

    check_parameter('ewscale_objfunction', 'standard_vmag', standard_vmag, 
                    'float')

    check_parameter('ewscale_objfunction', 'vega_wavelength', 
                    vega_wavelength, 'ndarray')

    check_parameter('ewscale_objfunction', 'vega_fluxdensity', 
                    vega_fluxdensity, 'ndarray')

    check_parameter('ewscale_objfunction', 'vega_continuum', 
                    vega_continuum, 'ndarray')

    check_parameter('ewscale_objfunction', 'vega_fitted_continuum', 
                    vega_fitted_continuum, 'ndarray')

    check_parameter('ewscale_objfunction', 'kernel', kernel, 'ndarray')

    check_parameter('ewscale_objfunction', 'line_wavlenegths', line_wavelengths, 
                    ['float', 'ndarray'])

    check_parameter('ewscale_objfunction', 'ew_scale', ew_scale, 'float')

    check_parameter('ewscale_objfunction', 'atmospheric_transmission', 
                    atmospheric_transmission, 'ndarray')

    check_parameter('ewscale_objfunction', 'poly_degree', poly_degree, 
                    'int')

    check_parameter('ewscale_objfunction', 'include_edges', include_edges, 
                    'bool')

    #
    # Upack the parameters
    #

    idx = poly_degree+2
    atm_scale = p[0]
    coeffs = p[1:idx]
    ew_scales = p[idx:]

    #
    # Does the user request adding the order edges to the control points?
    #
    
    if include_edges is True:
        
        points = np.append(np.insert(line_wavelengths, 0, 
                                     np.min(standard_wavelength)),
                           np.max(standard_wavelength))
        
        
        
        scales = np.append(np.insert(ew_scales, 0, ew_scale),
                           ew_scale)

    else:

        points = line_wavelengths
        scales = ew_scales

#    cs = CubicSpline(ew_points, ew_values)
#    ewscales = cs(vega_wavelength)

    #
    # Interpolate the control points onto the vega model wavelength grid
    #

    f = scipy.interpolate.interp1d(points, ew_scales, 
                                   fill_value=(ew_scales[0], 
                                               ew_scales[-1]),
                                   bounds_error=False)
    scales = f(vega_wavelength)

    #
    #  Generate the Vega model
    #

    result = modify_kuruczvega(standard_wavelength,
                               standard_rv,
                               standard_vmag,
                               standard_bmag,
                               vega_wavelength,
                               vega_fluxdensity,
                               vega_continuum,
                               vega_fitted_continuum,
                               kernel,
                               scales,
                               units='erg s-1 cm-2 A-1',
                               new=False)

    #
    # Build the model of the data
    #

    vega_model = result[0]
    vega_model /= np.median(vega_model)

    atmosphere = (atm_scale*(atmospheric_transmission-1))+1

    model = vega_model*atmosphere*poly_1d(standard_wavelength,coeffs)

    #
    # Compute the sum of the squared residuals
    #

    difference = (standard_fluxdensity-model)**2

    objective = np.sum((standard_fluxdensity-model)**2)

    #
    # Return the value
    #

    return objective





def find_shift(wavelength_object:npt.ArrayLike,
               intensity_object:npt.ArrayLike,
               intensity_telluric:npt.ArrayLike,
               wavelength_range:list,
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

    check_parameter('find_shift', 'wavelength_object', wavelength_object,
                    'ndarray')

    check_parameter('find_shift', 'intensity_object', intensity_object,
                    'ndarray')

    check_parameter('find_shift', 'intensity_telluric', intensity_telluric,
                    'ndarray')
    
    check_parameter('find_shift', 'wavelength_range', wavelength_range,
                    'list')

    check_parameter('find_shift', 'nsteps', nsteps, 'int')


    #
    # Get set up
    #

    wrange = np.sort(wavelength_range)
    
    zdata = np.where((wavelength_object >= wrange[0]) &
                     (wavelength_object <= wrange[1]))[0]
    
    shifts = np.linspace(pixel_range[0], pixel_range[1], num=nsteps)
    nshifts = len(shifts)

    rms = np.zeros(nshifts)
    
    x = np.arange(len(wavelength_object))

    #
    # Start the loop
    #
    
    for i in range(nshifts):

        xshift = x+shifts[i]

        telluric_shifted = linear_interp1d(xshift,intensity_telluric,x)

        ratio = np.multiply(intensity_object,telluric_shifted)

        rms[i] = np.std(ratio[zdata])

    #
    # Find the shift
        

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


    
def make_instrument_profile(x:npt.ArrayLike,
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


def make_telluric_spectrum(standard_wavelength:npt.ArrayLike,
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
                           control_points:npt.ArrayLike,
                           control_values:npt.ArrayLike,
                           order,
                           intensity_unit:str='erg s-1 cm-2 A-1',
                           new=False):

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

    check_parameter('make_telluric_spectrum', 'kernel', kernel, 'ndarray', 1)


    #
    # Determine the EW scale factor array using the control points
    #

#    cs = CubicSpline(control_points, control_values)
#    ewscales = cs(vega_wavelength)

    f = scipy.interpolate.interp1d(control_points,
                                   control_values,
                                   fill_value=(control_values[0], 
                                               control_values[-1]),
                                   bounds_error=False)
    ewscales = f(vega_wavelength)



    #
    # Modify the Kurucz Vega model 
    #


    vega_fd, vega_cont = modify_kuruczvega(standard_wavelength,
                                           standard_rv,
                                           standard_vmag,
                                           standard_bmag,
                                           vega_wavelength,
                                           vega_fluxdensity,
                                           vega_continuum,
                                           vega_fitted_continuum,
                                           kernel,
                                           ewscales,
                                           units=intensity_unit,
                                           new=new)

    #
    # Create the telluric correction spectrum and return results
    #

    telluric = vega_fd/standard_fluxdensity
    

    telluric_var = (vega_fd/standard_fluxdensity**2)**2 * \
        standard_uncertainty**2

    telluric_unc = np.sqrt(telluric_var)
    
    return telluric, telluric_unc, vega_fd, vega_cont


def measure_linerv(data_wavelength:npt.ArrayLike,
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

        plot_measure_linerv(qashow_info['plot_number'],
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

        plot_measure_linerv(None,
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




def modify_kuruczvega(standard_wavelength:npt.ArrayLike,
                      standard_rv:float,
                      standard_vmag:float,
                      standard_bmag:float,
                      vega_wavelength:npt.ArrayLike,
                      vega_fluxdensity:npt.ArrayLike,
                      vega_continuum:npt.ArrayLike,
                      vega_fitted_continuum:npt.ArrayLike,
                      kernel:npt.ArrayLike,
                      vega_ewscalefactors:float | npt.ArrayLike,
                      units:str='erg s-1 cm-2 A-1',
                      new=False):

    """
    To modify the Kurucz Vega model to match the standard star spectrum

    Parameters
    ----------
    standard_wavelength : ndarray
        A (nwave1,) array of data wavelengths.

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

    vega_ewscalefactors : ndarray or float
        A (nwave2,) array of EW scale factors or a float giving a 
        single scale factor.

    units : str
        The requested flux density units.
    
    Returns
    -------
    ndarray

        The modified Vega flux density

    ndarray

        The modified Vega continiuum

    """

    #
    # Check the parameters
    #
    
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

    check_parameter('modify_kuruczvega', 'vega_ewscalefactors',
                    vega_ewscalefactors, ['float','ndarray'])

    check_parameter('modify_kuruczvega', 'units', units, 'str')

    #
    # Get the Vega pixels that exactly cover the order
    #

    min = np.nanmin(standard_wavelength)
    max = np.nanmax(standard_wavelength)
    
    zleft = np.where(vega_wavelength > min)
    zleft_idx = zleft[0][0]
    
    zright = np.where(vega_wavelength < max)
    zright_idx = zright[0][-1]

    # Now figure out how many pixels you need to add to both sides to 
    # account for the width of the kernel during the convolution

    nadd = np.ceil(len(kernel)/2.)

    z_idx = np.arange(zleft_idx-nadd,zright_idx+nadd+1,dtype=int)
    
    #
    # Do the convolution on the normalized, zeroed, and EW scaled spectrum 
    # and continuum
    #

#    scale = 1
    scale = vega_ewscalefactors
    if new is True:

        input = (vega_fluxdensity[z_idx]/vega_fitted_continuum[z_idx]-1)*scale[z_idx]

        input = (input+1)*vega_fitted_continuum[z_idx]

        c_vega_fd = np.convolve(input, kernel, mode='same')



#\
#            vega_ewscalefactors[z_idx]
#        c_n_z_vega_fd = np.convolve(input, kernel, mode='same')
#        
        input = vega_fitted_continuum[z_idx]
        c_vega_cont = np.convolve(input, kernel, mode='same')
#        
#        # Reconstruct the flux density and continuum
#        
#        c_vega_fd = (c_n_z_vega_fd+1)*vega_fitted_continuum[z_idx]
#        c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont


        raw_vega_conv = np.convolve(vega_fluxdensity[z_idx], kernel, mode='same')



#        pl.figure()
#        pl.plot(vega_wavelength, vega_fluxdensity, color='red', 
#                label='Raw Vega')
#        pl.plot(vega_wavelength[z_idx], raw_vega_conv, color='black', 
#                label='Convolved Vega')
#
#        pl.plot(vega_wavelength[z_idx], c_vega_fd,label='Scaled, Convolved Vega')
#        pl.xlabel('Wavelength ($\mu$m)')
#        pl.xlim(np.min(vega_wavelength[z_idx]), 
#                  np.max(vega_wavelength[z_idx]))
#
#        pl.ylim(np.min(c_vega_fd), np.max(c_vega_fd))
#
#        pl.title('Order'+str(order))
#        pl.legend()
#        pl.savefig('Order '+str(order)+'.pdf')


    
    if new is False:

        input = (vega_fluxdensity[z_idx]/vega_continuum[z_idx]-1)*scale[z_idx]
#\
#            vega_ewscalefactors[z_idx]
        c_n_z_vega_fd = np.convolve(input, kernel, mode='same')
        
        input = (vega_continuum[z_idx]/vega_fitted_continuum[z_idx]-1)*scale[z_idx]
#\
#            vega_ewscalefactors[z_idx]
        c_n_z_vega_cont = np.convolve(input, kernel, mode='same')
        
        # Reconstruct the flux density and continuum
        
        c_vega_cont = (c_n_z_vega_cont+1)*vega_fitted_continuum[z_idx]
        
        c_vega_fd = (c_n_z_vega_fd+1)*c_vega_cont


        raw_vega_conv = np.convolve(vega_fluxdensity[z_idx], kernel, mode='same')



#        pl.figure()
#
#
#        pl.plot(vega_wavelength, vega_fluxdensity, color='red', 
#                label='Raw Vega')
#        pl.plot(vega_wavelength[z_idx], raw_vega_conv, color='black', 
#                label='Convolved Vega')
#
#        pl.plot(vega_wavelength[z_idx], c_vega_fd,label='Scaled, Convolved Vega')
#
#        pl.xlim(np.min(vega_wavelength[z_idx]), 
#                  np.max(vega_wavelength[z_idx]))
#        pl.ylim(np.min(c_vega_fd), np.max(c_vega_fd))
#        pl.title('Order'+str(order))
#        pl.legend()
#        pl.savefig('Order '+str(order)+'.pdf')


        
    #
    # Shift to the radial velocity of the standard
    #

    shifted_vega_wavelength = vega_wavelength*(1+standard_rv/2.99792458e5)

    #
    # Interpolate onto the wavelength grid of the standard
    #

    r_c_vega_fd = linear_interp1d(shifted_vega_wavelength[z_idx],
                                  c_vega_fd, standard_wavelength)

    r_c_vega_cont = linear_interp1d(shifted_vega_wavelength[z_idx],
                                    c_vega_cont, standard_wavelength)
    
    #
    # Redden the model to match that of the standard
    #

    Rv=3.1
    vega_vmag = 0.03
    vega_bminv = 0.0
    
    ebminv = (standard_bmag-standard_vmag - vega_bminv)
    ebminv = np.max([ebminv,0])
    Av = Rv*ebminv
    
    #
    # Extinguish the model to match that of the standard
    #

    ext = G23(Rv=Rv)
    
    e_r_c_vega_fd = r_c_vega_fd*\
        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)

    e_r_c_vega_cont = r_c_vega_cont*\
        ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
    
    #
    # Scale the model to the observed magnitude of the standard
    #

    scale = 10.**(-0.4*(standard_vmag-vega_vmag))*10**(0.4*Av)

    s_e_r_c_vega_fd = e_r_c_vega_fd*scale
    s_e_r_c_vega_cont = e_r_c_vega_cont*scale

    #
    # Convert the model to the requested flux density units
    #

    vega_fd = convert_fluxdensity(standard_wavelength,
                                  s_e_r_c_vega_fd,
                                  'um','erg s-1 cm-2 A-1',
                                  units)

    vega_cont = convert_fluxdensity(standard_wavelength,
                                    s_e_r_c_vega_cont,
                                    'um','erg s-1 cm-2 A-1',
                                    units)

    #
    # Returm the results
    #

    return vega_fd, vega_cont
    



def plot_measure_linerv(plot_number:int,
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


        
def plot_deconvolve_line(plot_number,
                         figure_size,
                         font_size,
                         spectrum_linewidth,
                         spine_linewidth,
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
    To plot the results of the deconvolution.

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
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
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


def plot_estimate_ewscales(plot_number:int,
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

    font = {'family' : 'helvetica',
            'weight' : 'normal',
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
                shifts:npt.ArrayLike):

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
    
    ncols = np.ceil(shifted_norders / subplot_stackmax).astype(int)

    nrows = np.min([shifted_norders,subplot_stackmax]).astype(int)

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

        # Did this order get shifted?

        if not zshifted[norders-i-1]:

            continue
        
        for j in range(napertures):

            k = i+j*napertures

            wavelength = object_spectra[norders-k-1,0,:]
            object_flux = object_spectra[norders-k-1,1,:]
            raw_telluric = rawtc_spectra[norders-k-1,1,:]
            shifted_telluric = shiftedtc_spectra[norders-k-1,1,:]
            
            raw_ratio = object_flux*raw_telluric
            shifted_ratio = object_flux*shifted_telluric            
            
            zshift = np.where((wavelength >= shift_ranges[norders-k-1,0]) &
                              (wavelength <= shift_ranges[norders-k-1,1]))[0]

            wavelength = wavelength[zshift]
            raw_ratio = raw_ratio[zshift]
            shifted_ratio = shifted_ratio[zshift]                        

            raw_ratio /= np.nanmedian(shifted_ratio)
            shifted_ratio /= np.nanmedian(shifted_ratio)            
            
            
            # Get the plot range

            xrange = get_spectra_range(wavelength)
            yrange = get_spectra_range(shifted_ratio, raw_ratio, frac=0.1)

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



            axe.set_title('Order ' + str(orders[norders-k-1])+', aperture '+\
                          str(j+1)+r', $\Delta x$='+'$'+str(shifts[norders-k-1,j])+\
                          '$ pixels')
            axe.set_ylabel('Relative Intensity')
            axe.set_xlabel(xlabel)            

            axe.xaxis.set_minor_locator(AutoMinorLocator())    
            axe.tick_params(right=True, left=True, top=True, bottom=True,
                            which='both', direction='in', width=1.5)
            axe.tick_params(which='minor', length=3)
            axe.tick_params(which='major', length=5)
            axe.yaxis.set_minor_locator(AutoMinorLocator())

#            axe2 = axe.twinx()
#            axe2.plot(atmosphere_wavelength,
#                      atmosphere_transmission,
#                      color='grey',
#                      label='Atmospheric Transmission')
#            axe2.set_ylim((-1,1))
            
            
            if m == 0:

                axe.legend()

            m += 1
    





    
    
