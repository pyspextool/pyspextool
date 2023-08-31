import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import typing
import os

from pyspextool.io.check import check_parameter
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.plot.limits import get_spec_range


def normalize_continuum(wavelength:npt.ArrayLike, flux:npt.ArrayLike,
                        ranges:npt.ArrayLike, degree:int,
                        robust:typing.Optional[dict]=None,
                        qa_show:bool=False, qa_show_plotsize:tuple=(6,10),
                        qa_fileinfo:typing.Optional[dict]=None):

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

    if qa_show is True:

        pl.ion()
        plot_normalize_continuum(qa_show_plotsize, wavelength, flux,
                                 continuum, zregions)

        pl.show()
        pl.pause(1)

    if qa_fileinfo is not None:
        
        pl.ioff()

        plot_normalize_continuum((8.5,11), wavelength, flux, continuum,
                                 zregions)
    
        pl.savefig(os.path.join(qa_fileinfo['filepath'],
                                qa_fileinfo['filename']) + \
                                qa_fileinfo['extension'])
        pl.close()
    
    return flux_normalized


def plot_normalize_continuum(figsize:tuple, wavelength:npt.ArrayLike,
                             flux:npt.ArrayLike, continuum:npt.ArrayLike,
                             zregions:npt.ArrayLike):

    """
    To plot the results of normalize_continuum

    Parameters
    ----------
    figsize : tuple
        A (2,) tuple giving the figure size in inches.

    wavelength : ndarray
        A (nwave,) array of wavelengths.

    flux : ndarray
        A (nwave,) array of "intensities".

    flux : ndarray
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
    
    fig = pl.figure(figsize=figsize)

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
    
    yrange = get_spec_range([flux[zplotrange], continuum[zplotrange]], frac=0.1)

    axes1 = fig.add_subplot(211)    
    axes1.step(wavelength, flux, '#1f77b4')
    axes1.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes1.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes1.set_xlabel('Wavlength')
    axes1.set_ylabel('Intensity')    

    # Now plot the fitted pixels in red
    
    tmp = np.copy(flux)
    tmp[~zregions] = np.nan
    
    axes1.plot(wavelength, tmp, 'red')

    # Plot the continuum
    
    axes1.step(wavelength, continuum, 'green')    

    #
    # Now plot the normalized spectrum
    #

    # Normalize and get the plot range
    
    normalized = flux/continuum
    yrange = get_spec_range(normalized[zplotrange], frac=0.1)


    axes2 = fig.add_subplot(212)    
    axes2.step(wavelength, normalized, '#1f77b4')
    axes2.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes2.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes2.set_xlabel('Wavlength')
    axes2.set_ylabel('Normalized Intensity')    
    
    # Plot the fitted pixels in red
    
    tmp = np.copy(normalized)
    tmp[~zregions] = np.nan
    
    axes2.plot(wavelength, tmp, 'red')
    axes2.axhline(y=1, linestyle='--', color='green')
