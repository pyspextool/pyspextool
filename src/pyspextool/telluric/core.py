import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.ticker import (AutoMinorLocator)

import logging
from scipy import signal
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.special import erf
from os.path import join as osjoin
from dust_extinction.parameter_averages import G23
import astropy.units as u

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.utils.math import moments
from pyspextool.utils.arrays import find_index
from pyspextool.plot.limits import get_spectra_range
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.interpolate import linear_bitmask_interp1d
from pyspextool.utils.math import combine_flag_stack
from pyspextool.pyspextoolerror import pySpextoolError

def correct_spectrum(object_wavelength:npt.ArrayLike,
                     object_intensity:npt.ArrayLike,
                     correction_wavelength:npt.ArrayLike,
                     correction_intensity:npt.ArrayLike,
                     uncertainties:list=None,
                     masks:list=None):

    """
    To correct a spectrum for telluric absorption and flux calibration.

    The function also optionally corrects uncertainty arrays and will merge
    bit masks.

    Parameters
    ----------

    object_wavelength : ndarray
        A (nwave1,) ndarray of wavelength values.
    
    object_intensity : ndarray
        A (nwave1,) ndarray of intensity values, typically in units of ADU s-1
        or DN s-1.

    correction_wavelength : ndarray
        A (nwave2,) ndarray of wavelength values.
    
    correction_intensity : ndarray
        A (nwave2,) ndarray of intensity values, typically in units of ADU s-1
        or DN s-1 per "flux density".

    uncertainties : list, default=None
        A (2,) list where uncertainties[0] is an (nwave1,) ndarray of
        uncertainties associated with the object_intensity ndarray and
        uncertainties[1] is an (nwave2,) ndarray of
        uncertainties associated with the correction_intensity ndarray

    masks : list, default=None
        A (2,) list where masks[0] is an (nwave1,) ndarray mask associated
        with the object_intensity ndarray and mask[1] is an (nwave2,) ndarray
        mask associated with the correction_intensity ndarray
       
    Returns
    -------
    ndarray, ndarray, ndarray

    
   
    """
    
    #
    # Check the parameters
    #

    check_parameter('correct_spectrum', 'object_wavelength', object_wavelength,
                    'ndarray')

    check_parameter('correct_spectrum', 'object_intensity', object_intensity,
                    'ndarray')

    check_parameter('correct_spectrum', 'correction_wavelength',
                    correction_wavelength, 'ndarray')

    check_parameter('correct_spectrum', 'correction_intensity',
                    correction_intensity, 'ndarray')

    check_parameter('correct_spectrum', 'uncertainties', uncertainties,
                    ['NoneType','list'], list_types=['ndarray', 'ndarray'])

    check_parameter('correct_spectrum', 'masks', masks, ['ndarray', 'list'],
                    list_types=['ndarray', 'ndarray'])


    #
    # Unpack the uncertainties
    #

    object_uncertainty = None
    correction_uncertainty = None
    
    if uncertainties is not None:

        object_uncertainty = uncertainties[0]
        correction_uncertainty = uncertainties[1]

    #
    # Do the correction
    #
    
    # Interpolate the correction spectrum onto to wavelength grid of 
    
    
    rtc_i, rtc_u = linear_interp1d(correction_wavelength,
                                   correction_intensity,
                                   object_wavelength,
                                   input_u=correction_uncertainty)

    # Correct and propagate errors
        
    fluxdensity = object_intensity*rtc_i

    if uncertainties is not None:
    
        fluxdensity_variance = object_intensity**2 * rtc_u**2 + \
            rtc_i**2 * object_uncertainty**2
        fluxdensity_uncertainty = np.sqrt(fluxdensity_variance)
        
    else:

        fluxdensity_uncertainty = None
        
    #
    # Do the masks
    #

    if masks is not None:

        object_mask = masks[0]
        correction_mask = masks[1]

        # Interpolate the masks and combine
    
        rmask = linear_bitmask_interp1d(correction_wavelength,
                                        correction_mask.astype(np.uint8),
                                        object_wavelength)
    
        stack = np.stack((object_mask.astype(np.uint8),rmask))
        fluxdensity_mask = combine_flag_stack(stack)


    else:

        fluxdensity_mask = None


    return fluxdensity, fluxdensity_uncertainty, fluxdensity_mask



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
        if qashow_info['block'] is False: pl.pause(1)

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
                           scale:float):

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

    scale : float

    
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

    check_parameter('make_telluric_spectrum', 'scale', scale, ['int','float'])


    #
    # Determine the range over which to convolved the Vega model
    #

    # Get the pixels that exactly cover the order
    
    min = np.nanmin(standard_wavelength)
    max = np.nanmax(standard_wavelength)
    
    zleft = np.where(vega_wavelength > min)
    zleft_idx = zleft[0][0]
    
    zright = np.where(vega_wavelength < max)
    zright_idx = zright[0][-1]

    # Now figure out how may you have to add to both sides to account for
    # the width of the kernel

    nadd = np.ceil(len(kernel)/2.)

    z_idx = np.arange(zleft_idx-nadd,zright_idx+nadd+1,dtype=int)

    # Do the convolution
    
    con_norm_vega_fd = np.convolve((vega_fluxdensity[z_idx]/\
                                   vega_continuum[z_idx]-1)*scale,
                                   kernel, mode='same')+1

    con_norm_vega_cont = np.convolve((vega_continuum[z_idx]/\
                                     vega_fitted_continuum[z_idx]-1)*scale,
                                     kernel, mode='same')+1
    
    con_vega_cont = con_norm_vega_cont*vega_fitted_continuum[z_idx]

    con_vega_fd = con_norm_vega_fd*con_vega_cont
    
    #
    # Shift to the radial velocity of the standard
    #

    shifted_vega_wavelength = vega_wavelength*(1+standard_rv/2.99792458e5)

    #
    # Interpolate onto the wavelength grid of the standard
    #

    r_con_vega_fd = linear_interp1d(shifted_vega_wavelength[z_idx],
                                    con_vega_fd, standard_wavelength)
    
    #
    # Redden the model to match that of the standard
    #

    Rv=3.1
    vega_vmag = 0.03
    vega_bminv = 0.0
    
    ebminv = (standard_bmag-standard_vmag - vega_bminv)
    ebminv = np.max([ebminv,0])
    Av = Rv*ebminv
    
    # Load the extinction package and extinguish
    
    ext = G23(Rv=Rv)
    
    r_con_vega_fd *= ext.extinguish(standard_wavelength*1e4*u.AA,Ebv=ebminv)
    
    #
    # Scale the Vega model to the observed magnitude of the standard
    #

    scale = 10.**(-0.4*(standard_vmag-vega_vmag))*10**(0.4*Av)
    r_con_vega_fd *= scale

    #
    # Calculate the ratio:  model/standard
    #

    telluric_spectrum = r_con_vega_fd/standard_fluxdensity
    
    telluric_spectrum_var = (r_con_vega_fd/standard_fluxdensity**2)**2 * \
                                    standard_uncertainty**2
    telluric_spectrum_unc = np.sqrt(telluric_spectrum_var)
    
    
    return telluric_spectrum, telluric_spectrum_unc, r_con_vega_fd




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
        if qashow_info['block'] is False: pl.pause(1)

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
    
    font = {'family' : 'helvetica',
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
    int
    The plot number

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

    font = {'family' : 'helvetica',
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


    





    
    
