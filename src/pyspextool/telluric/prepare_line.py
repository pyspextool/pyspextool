import numpy as np
import numpy.typing as npt
import logging
import matplotlib.pyplot as pl
import matplotlib
from matplotlib.ticker import (AutoMinorLocator)
import os

from pyspextool import config as setup
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_range, check_qakeywords
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.plot.limits import get_spectra_range
from pyspextool.utils.arrays import find_index

def prepare_line(order:int,
                 line_wavelength:float,
                 resolving_power:int | float,
                 wavelength_range:npt.ArrayLike,
                 fit_type:str,
                 poly_degree:int,
                 verbose:bool=None,
                 qa_show:bool=None,
                 qa_showscale:float=None,
                 qa_showblock:bool=None,
                 qa_write:bool=None):

    """
    To normalize a spectrum around a line for RV/deconvolution work

    Parameters
    ----------
    order : int
        The order number with the line to be fit.

    line_wavelength : float
        The wavelength of the line to be fit.

    resolving_power : int or float
        The resolving power of the spectrum.
    
    wavelength_range : ndarray
        An (2,) array of wavelenths that encompose the region to be fit.

    fit_type : {'lorentzian', 'gaussian'}
        The type of function to fit to the line.
    
    poly_degree : int
        The polynomial degree to fit the continuum with (deg=1 -> line).

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
           
    Returns
    -------
    None
    Load data into memory and writes QA plots to disk.

        tc.state['normalized_line_wavelength']
        tc.state['normalized_line_flux']
        tc.state['line_center']
        tc.state['line_fwhm']
        tc.state['prepare_done']
    
    """

    #
    # Check the load_done variable
    #
    
    if tc.state['load_done'] is False:

        message = "Spectra have not been loaded.  Please run load_spectra.py."
        raise pySpextoolError(message)

    #
    # Check the parameters and keywords
    #

    check_parameter('prepare_line', 'order', order, 'int')

    check_parameter('prepare_line', 'line_wavelength', line_wavelength,
                    'float')

    check_parameter('prepare_line', 'resolving_power', resolving_power,
                    ['int', 'float'])
    
    check_parameter('prepare_line', 'wavelength_range', wavelength_range,
                    ['ndarray','list'])
    
    check_parameter('prepare_line', 'fit_type', fit_type, 'str',
                    possible_values=['gaussian','lorentzian'])    
    
    check_parameter('prepare_line', 'poly_degree', poly_degree, 'int')

    check_parameter('prepare_line', 'verbose', verbose, ['NoneType','bool'])

    check_parameter('prepare_line', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('prepare_line', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])

    check_parameter('prepare_line', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
    
    check_parameter('prepare_line', 'qa_write', qa_write, ['NoneType','bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # log the operation
    #
    
    logging.info(" Normalizing continuum in order "+str(order)+".")

    #
    # Get set up for the normalization
    #
    
    # Find the order given the modeinfo file

    z_order = np.where(tc.state['standard_orders'] == order)


    # Store values in shorter variable names for ease
    
    wavelength = np.squeeze(tc.state['standard_spectra'][z_order,0,:])
    flux = np.squeeze(tc.state['standard_spectra'][z_order,1,:])
    xlabel = tc.state['standard_hdrinfo']['LXLABEL'][0]
    title = tc.state['standard_name']+', '+\
        tc.state['mode']+' Order '+str(order)+', degree='+str(poly_degree)
    
    #
    # Determine if wavelength_range is monotonically increasing
    #

    if wavelength_range[1] <= wavelength_range[0]:
        
        message = "The fit range, "+str(wavelength_range)+\
            ", is not monotonic."
        raise pySpextoolError(message)

    #
    # Check to make sure the two values fall within the wavelength range of
    # the order
    #

    min = float(np.nanmin(wavelength))

    max = float(np.nanmax(wavelength))
    
    check_range(wavelength_range[0],[min,max],'gtlt','wavelength_range[0]')

    check_range(wavelength_range[1],[min,max],'gtlt','wavelength_range[1]')    
       
    # Cut the region out of the order

    zleft = (wavelength > wavelength_range[0])
    zright = (wavelength < wavelength_range[1])
        
    zselection = np.logical_and(zleft,zright)

    wavelength = wavelength[zselection]
    flux = flux[zselection]
    
    # Fit the line and continuum

    idx = find_index(wavelength,line_wavelength)
    p0 = np.array([flux[int(idx)],line_wavelength,1/float(resolving_power)])
    p0 = np.pad(p0,(0,poly_degree+1),constant_values=0.0)

    result = fit_peak1d(wavelength,
                        flux,
                        type=fit_type,
                        negative=True,
                        p0=p0,
                        nparms=3+poly_degree+1,
                        robust={'thresh':4, 'eps':0.1})

    continuum_coefficients = result['parms'][3:]
    line_center = float(result['parms'][1])
    line_halfwidth = float(result['parms'][2])

    continuum = poly_1d(wavelength, continuum_coefficients)    
    normalized_flux = flux/continuum

    # added error message
#    if len(normalized_flux[np.isfinite(normalized_flux)==False]==0):
#        message = "Failure to fit a H I line in the fit range "+str(wavelength_range)+\
#            "; verify there is a line in this spectrum."
#        raise pySpextoolError(message)


    #
    # Store the reults
    #

    tc.state['normalized_line_wavelength'] = wavelength
    tc.state['normalized_line_flux'] = normalized_flux
    tc.state['line_center'] = line_center

    if fit_type == 'gaussian':

        tc.state['line_fwhm'] = 2.354*line_halfwidth

    else:

        tc.state['line_fwhm'] = 2*line_halfwidth        
        
    #
    # Make the QA plot
    #

    if qa['show'] is True:

        figure_size = (setup.plots['portrait_size'][0]*qa['showscale'],
                       setup.plots['portrait_size'][1]*qa['showscale'])

        font_size = setup.plots['font_size']*qa['showscale']
        
        plot_normalization(setup.plots['normalize_order'],
                           figure_size,
                           font_size,
                           setup.plots['zoomspectrum_linewidth'],
                           setup.plots['spine_linewidth'],                     
                           wavelength,
                           flux,
                           result['fit'],
                           result['goodbad'],
                           line_center,
                           line_halfwidth,
                           continuum,
                           plot_xlabel=xlabel,
                           plot_title=title)

        pl.show(block=qa['showblock'])
        if qa['showblock'] is False:

            pl.pause(1)
        
    if qa['write'] is True:

        plot_normalization(None,
                           setup.plots['portrait_size'],
                           setup.plots['font_size'],
                           setup.plots['zoomspectrum_linewidth'],
                           setup.plots['spine_linewidth'],                     
                           wavelength,
                           flux,
                           result['fit'],
                           result['goodbad'],
                           line_center,
                           line_halfwidth,
                           continuum,
                           plot_xlabel=xlabel,
                           plot_title=title)

        pl.savefig(os.path.join(setup.state['qa_path'],
                                tc.state['output_filename']+ \
                                '_normalization' + \
                                setup.state['qa_extension']))
        pl.close()

    #
    # Set the done variable
    #
        
    tc.state['prepare_done'] = True

        
def plot_normalization(plot_number:int,
                       figure_size:tuple,
                       font_size:int,
                       spectrum_linewidth:int | float,
                       spine_linewidth:int | float,        
                       wavelength:npt.ArrayLike,
                       intensity:npt.ArrayLike,
                       fit:npt.ArrayLike,
                       goodbad:npt.ArrayLike,
                       line_center:float,
                       line_halfwidth:float,
                       continuum:npt.ArrayLike,
                       plot_xlabel:str=None,
                       plot_title:str=None):

    
    """
    To plot the results of prepare_line in a device independent way

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths.

    intensity : ndarray
        A (nwave,) array of "intensities".

    fit : ndarray
        A (nwave,) array of the fitted values of the line+continuum.
        
    line_center : float
        The line center in units of `wavelength`.

    line_halfwidth : float
        The line "half width" either the HWHM of the Lorentzian fit or standard
        deviation of the gaussisn fit.  In units of `wavelength`.

    continuum : ndarray
        A (nwave,) array of the fitted continuum values.

    plot_scale : float
        A value by which to scale the default size of the figure.  
    
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
    
    check_parameter('plot_normalization', 'fit', fit, 'ndarray')

    check_parameter('plot_normalization', 'line_center', line_center, 'float')
    
    check_parameter('plot_normalization', 'line_halfwidth', line_halfwidth,
                    'float')
    
    check_parameter('plot_normalization', 'continuum', continuum, 'ndarray')    

    check_parameter('plot_normalization', 'plot_xlabel', plot_xlabel, 'str')

    check_parameter('plot_normalization', 'plot_title', plot_title, 'str')    
        
    #
    # Make the two-panel figure
    #
    
    # Set the fonts

    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    matplotlib.rc('font', **font)

    # Start the figure, and set the spacing
    
    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.clf()
    pl.subplots_adjust(left=0.1,
                       bottom=0.1, 
                       right=0.95, 
                       top=0.9, 
                       hspace=0.05)
    
    # Get the plot range for x axis

    xrange = [np.min(wavelength), np.max(wavelength)]
    
    #
    # Create the spectral plot with the fit
    #

    # Determine the yrange for spectral plot
    
    yrange = get_spectra_range([intensity,continuum],frac=0.1)


    
    axes1 = fig.add_subplot(211)    
    axes1.step(wavelength, intensity, 'black',lw=spectrum_linewidth)
    axes1.set_title(plot_title)
    axes1.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes1.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes1.set_ylabel('Normalized Intensity')

    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=spine_linewidth,
                      labelbottom=False)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())
    axes1.step(wavelength, fit, 'red')

    z = goodbad == 0
    axes1.plot(wavelength[z], intensity[z], 'or')

    
    # change all spines
    for axis in ['top','bottom','left','right']:
        axes1.spines[axis].set_linewidth(spine_linewidth)
    
    # Plot the continuum
    
    axes1.step(wavelength, continuum, 'green')    

    #
    # Now plot the normalized spectrum
    #

    # Normalize and get the plot range
    
    normalized = intensity/continuum
    yrange = get_spectra_range(normalized, frac=0.1)

    
    axes2 = fig.add_subplot(212)    
    axes2.step(wavelength, normalized, 'black',lw=spectrum_linewidth)
    axes2.set_ylim(ymin = yrange[0], ymax=yrange[1])
    axes2.set_xlim(xmin = xrange[0], xmax=xrange[1])    
    axes2.set_xlabel(plot_xlabel)
    axes2.set_ylabel('Normalized Intensity')    
    axes2.axvline(x=line_center, linestyle='--', color='red')
    axes2.axvline(x=line_center-line_halfwidth, linestyle='--', color='red')
    axes2.axvline(x=line_center+line_halfwidth, linestyle='--', color='red')
        
    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=spine_linewidth)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())    
    axes2.axhline(y=1, linestyle='--', color='green')


    # change all spines
    for axis in ['top','bottom','left','right']:
        axes2.spines[axis].set_linewidth(spine_linewidth)

    
    #
    # Get the plot number and return the results
    #
    

