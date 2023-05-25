import numpy as np
import matplotlib.pyplot as pl
import os

from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.plot.limits import get_spec_range
from pyspextool.io.check import check_parameter, check_file
from pyspextool.io.read_spectra_fits import read_spectra_fits


def plot_spectra(file, plot_type='continuous', plot_size=(10, 6), y='flux',
                 aperture=0, xlabel=None, ylabel=None, title=None,
                 colors='green', line_width=1, yrange_buffer=0.05,
                 order_numbers=True, file_info=None, display=True):
    """
    To plot a pyspextool FITS file.

    Parameters
    ----------
        file : str
            The full path to a pySpextool spectral fits files.

        plot_type : {'continuous', 'ladder'}
            Currently only continuous is allowed.

        plot_size : tuple of (float, float), optional
            A (2,) tuple giving the page size.

        y : {'flux', 'uncertainty', 'snr', 'flux and uncertainty'}, optional
            Which spectrum to plot.

        aperture : int, default=0
            The aperture to plot.  It is zero-based.

        xlabel : str, optional
            The label for the x axis.

        ylabel : str, optional
            The label for the y axis.

        title : str, optional
            The title.

        colors : str or list of (str, str) or list of (str, str, str), 
                 default='green'
            The spectra will be plotted with alternative colors depending on 
            how many colors are given.
            
        line_width : float or int, default=1
            The line width value passed to maplotlib.

        yrange_buffer : float, default=0.05
            The fraction by which to expand the y range if desired.

        order_numbers : bool, default=True
            Set to True for the order numbers to be plotted on the plot.

        file_info : dict, optional

            `'figsize'` : tuple
                A (2,) tuple giving the figure size.

            `'filepath'` : str, optional, default=''

            `'filename'` : str, optional, default=root of spectra.

            `'extension'` : str, optional, default='.pdf'

        display : bool, default=True
            Set to True to display figure with call.
            
    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter('plot_spectra', 'file', file, 'str')

    check_parameter('plot_spectra', 'plot_type', plot_type, 'str',
                    ['continuous', 'ladder'])

    check_parameter('plot_spectra', 'plot_size', plot_size, 'tuple')

    check_parameter('plot_spectra', 'y', y, 'str',
                    ['flux', 'uncertainty', 'snr'])

    check_parameter('plot_spectra', 'aperture', aperture, 'int')

    check_parameter('plot_spectra', 'xlabel', xlabel, ['NoneType', 'str'])

    check_parameter('plot_spectra', 'ylabel', ylabel, ['NoneType', 'str'])

    check_parameter('plot_spectra', 'title', title, ['NoneType', 'str'])

    check_parameter('plot_spectra', 'colors', colors,
                    ['NoneType', 'str', 'list'])

    check_parameter('plot_spectra', 'line_width', line_width, ['float', 'int'])
    
    check_parameter('plot_spectra', 'yrange_buffer', yrange_buffer, 'float')

    check_parameter('plot_spectra', 'order_numbers', order_numbers, 'bool')

    check_parameter('plot_spectra', 'file_info', file_info,
                    ['NoneType', 'dict'])

    #
    # Check various other things
    #

    # Does the file exist?

    check_file(file)

    #
    # Read the file
    #

    spectra, info = read_spectra_fits(file)
    norders = info['norders']
    napertures = info['napertures']
    orders = info['orders']
    latex_yunits = info['lyunits']
    latex_xlabel = info['lxlabel']
    latex_ylabel = info['lylabel']

    #
    # Get the plot ranges
    #

    wranges, franges, uranges, furanges, sranges = get_ranges(spectra,
                                                      info['norders'],
                                                      info['napertures'],
                                                      aperture,
                                                      fraction=yrange_buffer)

    xrange = [np.min(wranges), np.max(wranges)]

    #
    # Get the index for the spectrum
    #

    if y == 'flux':

        yrange = [np.nanmin(franges), np.nanmax(franges)]
        ylabel = latex_ylabel

    elif y == 'uncertainty':

        ylabel = 'Uncertainty (' + latex_yunits + ')'
        yrange = [np.nanmin(uranges), np.nanmax(uranges)]

    elif y == 'snr':

        ylabel = 'Signal to Noise'
        yrange = [np.nanmin(sranges), np.nanmax(sranges)]


    elif y == 'flux and uncertainty':

        yrange = [np.nanmin(furanges), np.nanmax(furanges)]

    else:

        raise ValueError('Do not recognize y axis name {}'.format(y))


    if plot_type == 'continuous':

        if file_info is None:

            figure_size = plot_size

        else:

            figure_size = file_info['figsize']

        #
        # Make the plot
        #

        figure = pl.figure(figsize=figure_size)
        ax = figure.add_axes([0.125, 0.11, 0.775, 0.77])
        ax.plot([np.nan], [np.nan])
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_xlabel(latex_xlabel)
        ax.set_ylabel(ylabel)

        # Get the coordinate transform ready for order number labels

        axis_to_data = ax.transAxes + ax.transData.inverted()
        data_to_axis = axis_to_data.inverted()

        for i in range(norders):

            # Get the plotting values

            xvalues = spectra[i * napertures + aperture, 0, :]
            y2values = [np.nan]*len(xvalues)

            if y == 'flux':
                
                yvalues = spectra[i*napertures+aperture,1,:]

            elif y == 'uncertainty':
            
                yvalues = spectra[i*napertures+aperture,2,:]
                
            elif y == 'snr':
                
                yvalues = spectra[i*napertures+aperture,1,:]/\
                spectra[i*napertures+aperture,2,:]

            elif y == 'flux and uncertainty':

                yvalues = spectra[i*napertures+aperture,1,:]
                y2values = spectra[i*napertures+aperture,2,:]


            # Get the colors

            if isinstance(colors, list):

                ncolors = len(colors)

                if ncolors == 2:

                    color = colors[0] if i % 2 == 0 else colors[1]

                elif ncolors == 3:

                    mod = i % 3

                    if mod == 0: color = colors[0]

                    if mod == 1: color = colors[1]

                    if mod == 2:  color = colors[2]

            else:

                color = colors

            # Plot the spectrum
                
            ax.plot(xvalues, yvalues, color=color, ls='-', lw=line_width)
            ax.plot(xvalues, y2values, color='grey', lw=line_width)

            # Now label the order numbers

            if order_numbers is True:

                # Get the order number position

                ave_wave = np.mean(wranges[i, :])
                xnorm = data_to_axis.transform((ave_wave, 0))[0]

                # Print the text

                ax.text(xnorm, 1.01 + 0.01 * (i % 2), str(orders[i]),
                        color=color, transform=ax.transAxes)

                # Add the title

                ax.set_title(title, pad=20.0)

            else:

                ax.set_title(title)

    if file_info is not None:

        keys = file_info.keys()

        if 'filepath' in keys:

            filepath = file_info['filepath']

        else:

            filepath = ''

        if 'filename' in keys:

            filename = file_info['filename']

        else:

            filename = os.path.basename(filename).removesuffix('.fits')

        if 'extension' in keys:

            extension = file_info['extension']

        else:

            extension = '.pdf'

        pl.savefig(os.path.join(filepath, filename + extension))

    if display==True: pl.show()

    pl.close()
    return ax



def get_ranges(spectra, norders, napertures, aperture, fraction=0.05):

    """
    To determine the x and y ranges of each order

    Parameters 
    ----------
        spectra : nd.array

        norders : int
            The number of orders.

        napertures : int
            The number of apertures.

        aperture : int
            The aperture requested.  The value is zero-based.

        fraction : float, default=0.05
            The fraction by which to expand the range if desired.

    Returns
    -------
        tuple : (ndarray, ndarray, ndarray, ndarray)

    """

    wave_ranges = np.empty((norders, 2))
    flux_ranges = np.empty((norders, 2))
    unc_ranges = np.empty((norders, 2))
    fluxunc_ranges = np.empty((norders, 2))    
    snr_ranges = np.empty((norders, 2))

    for i in range(norders):
        # Do the wavelengths

        wave = spectra[i * napertures + aperture, 0, :]
        wave_ranges[i, :] = [np.nanmin(wave), np.nanmax(wave)]

        # Do the various intensities.  Smooth with a Savitzky-Golay to avoid
        # bad pixels.

        ndat = len(wave)
        x_values = np.arange(ndat)

        flux = spectra[i * napertures + aperture, 1, :]
        sg_flux = robust_savgol(x_values, flux, 11)['fit']

        flux_ranges[i, :] = get_spec_range(sg_flux, frac=fraction)

        unc = spectra[i * napertures + aperture, 2, :]
        sg_unc = robust_savgol(x_values, unc, 11)['fit']

        unc_ranges[i, :] = get_spec_range(sg_unc, frac=fraction)

        fluxunc_ranges[i,0] = min(flux_ranges[i,0], unc_ranges[i,0])
        fluxunc_ranges[i,1] = max(flux_ranges[i,1], unc_ranges[i,1])        
        
        snr = flux / unc
        sg_snr = robust_savgol(x_values, snr, 11)['fit']

        snr_ranges[i, :] = get_spec_range(sg_snr, frac=fraction)


        
    return wave_ranges, flux_ranges, unc_ranges, fluxunc_ranges, snr_ranges
