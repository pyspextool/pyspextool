import numpy as np
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
import os

from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.plot.limits import get_spec_range
from pyspextool.io.check import check_parameter, check_file
from pyspextool.io.read_spectra_fits import read_spectra_fits


def plot_spectra(file, plot_size=(10, 6), ytype='flux', aperture=None,
                 title=None, colors='green', line_width=0.5,
                 yrange_buffer=0.05, order_numbers=True, file_info=None,
                 plot_number=None):

    """
    To plot a pyspextool FITS file.

    Parameters
    ----------
        file : str
            The full path to a pySpextool spectral fits files.

        plot_size : tuple of (float, float), optional
            A (2,) tuple giving the page size.

        ytype : {'flux', 'uncertainty', 'snr', 'flux and uncertainty'}
            Which spectrum to plot.

        aperture : int or NoneType, default=None
            The aperture to plot.  It is one-indexed.  if None, all apertures
            are plotted.

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

        order_numbers : {True, False}
            Set to True for the order numbers to be plotted on the plot.

        file_info : dict, optional

            `'figsize'` : tuple
                A (2,) tuple giving the figure size.

            `'filepath'` : str, optional, default=''

            `'filename'` : str, optional, default=root of spectra.

            `'extension'` : str, optional, default='.pdf'
            
    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter('plot_spectra', 'file', file, 'str')

    check_parameter('plot_spectra', 'plot_size', plot_size, 'tuple')

    check_parameter('plot_spectra', 'ytype', ytype, 'str',
                    ['flux', 'uncertainty', 'snr'])

    check_parameter('plot_spectra', 'aperture', aperture, ['NoneType', 'int'])

    check_parameter('plot_spectra', 'title', title, ['NoneType', 'str'])

    check_parameter('plot_spectra', 'colors', colors,
                    ['NoneType', 'str', 'list'])

    check_parameter('plot_spectra', 'line_width', line_width, ['float', 'int'])
    
    check_parameter('plot_spectra', 'yrange_buffer', yrange_buffer, 'float')

    check_parameter('plot_spectra', 'order_numbers', order_numbers, 'bool')

    check_parameter('plot_spectra', 'file_info', file_info,
                    ['NoneType', 'dict'])

    # Does the file exist?
    
    check_file(file)

    #
    # Read the file
    #

    spectra, info = read_spectra_fits(file)
    norders = info['norders']
    napertures = info['napertures']
    latex_yunits = info['lyunits']
    latex_xlabel = info['lxlabel']
    latex_ylabel = info['lylabel']

    if order_numbers is True:

        orders = info['orders']

    else:

        orders = None

    if plot_number is None:

        plot_numbers = [None]*napertures

    else:

        plot_numbers = plot_number

    #
    # Deal with file or window
    #

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

        if 'figsize' in keys:

            figsize = file_info['figsize']

        else:

            figsize = (10,6)


        if extension == '.pdf':

            # Create a multipage pdf file
            
            pdf = PdfPages(os.path.join(filepath,filename+extension))

        #
        # Loop over each aperture
        #

        for i in range(napertures):

            if napertures != 1:

                label = 'aperture_'+str(i+1).zfill(2)

            else:

                label = ''

            # Only plot a single aperture if requested.
            
            if aperture is not None:

                if i != aperture-1:

                    continue
        
            # create the index array to get the spectra for each aperture
        
            aperture_index = np.arange(norders)*napertures+i
            
            sub_spectra = spectra[aperture_index, :, :]

            fig = plot_aperture(sub_spectra, i, ytype, colors, line_width,
                                    latex_xlabel, latex_ylabel, latex_yunits,
                                    yrange_buffer, figsize, plot_numbers[i],
                                    orders=orders, title=title)

            if extension == '.pdf':

                pdf.savefig(fig)
                pl.close(fig)            

            else:
                
                pl.savefig(os.path.join(filepath,filename+label+extension))
                pl.close()    
                    
        if extension == '.pdf':
                    
            pdf.close()
            
    else:
        
        #
        # Loop over each aperture
        #

        for i in range(napertures):

            # Only plot a single aperture if requested.
        
            if aperture is not None:

                if i != aperture-1:

                    continue
        
            # create the index array to get the spectra for each aperture
        
            aperture_index = np.arange(norders)*napertures+i
            
            sub_spectra = spectra[aperture_index, :, :]

            pl.ion()
            number = plot_aperture(sub_spectra, i, ytype, colors, line_width,
                                   latex_xlabel, latex_ylabel, latex_yunits,
                                   yrange_buffer, plot_size, plot_numbers[i],
                                   orders=orders, title=title)
            pl.show()
            pl.pause(1)
            plot_numbers[i] = number

        return plot_numbers


        
def plot_aperture(spectra, aperture, ytype, colors, line_width, latex_xlabel,
                  latex_ylabel, latex_yunits, yrange_buffer, figure_size,
                  plot_number, orders=None, title=None):

    """
    Actually does the plot independent of device.

    


    """

    #
    # Get the plot ranges
    #
    
    wranges, yranges = get_ranges(spectra, ytype, fraction=yrange_buffer)

    xrange = [np.nanmin(wranges), np.nanmax(wranges)]
    yrange = [np.nanmin(yranges), np.nanmax(yranges)]

    #
    # Deal with the labels
    #
    
    if ytype == 'flux':

        ylabel = latex_ylabel

    elif ytype == 'uncertainty':

        ylabel = 'Uncertainty (' + latex_yunits + ')'

    elif ytype == 'snr':

        ylabel = 'Signal to Noise'


    elif ytype == 'flux and uncertainty':

        ylabel = latex_ylabel       
        
    else:

        raise ValueError('Do not recognize y axis name {}'.format(ytype))

    #
    # Make the plot
    #

    figure = pl.figure(num=plot_number, figsize=figure_size)
    pl.clf()
    ax = figure.add_axes([0.125, 0.11, 0.775, 0.77])
    ax.plot([np.nan], [np.nan])
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.set_xlabel(latex_xlabel)
    ax.set_ylabel(ylabel)
    
    # Get the coordinate transform ready for order number labels
    
    axis_to_data = ax.transAxes + ax.transData.inverted()
    data_to_axis = axis_to_data.inverted()

    norders = np.shape(spectra)[0]
    
    for i in range(norders):
            
        # Get the plotting values
            
        xvalues = spectra[i, 0, :]
            
        if ytype == 'flux':
                
            yvalues = spectra[i,1,:]
            y2values = [np.nan]*len(yvalues)

        elif ytype == 'uncertainty':
            
            yvalues = spectra[i,2,:]
            y2values = [np.nan]*len(yvalues)
                
        elif ytype == 'snr':
                
            yvalues = spectra[i,1,:]/spectra[i,2,:]
            y2values = [np.nan]*len(yvalues)

        elif ytype == 'flux and uncertainty':

            yvalues = spectra[i,1,:]
            y2values = spectra[i,2,:]

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

        if orders is not None:

            # Get the order number position

            ave_wave = np.mean(wranges[i, :])
            xnorm = data_to_axis.transform((ave_wave, 0))[0]

            # Print the text

            ax.text(xnorm, 1.01 + 0.01 * (i % 2), str(orders[i]),
                    color=color, transform=ax.transAxes)

            # Add the title

            if title is not None:

                ax.set_title(title+' - Aperture '+str(aperture+1), pad=20.0)

            # else:

            #     ax.set_title('Aperture '+str(aperture+1), pad=20.0)                

        else:

            if title is not None:

                ax.set_title('Aperture '+str(aperture+1))

    #
    # Get plot number
    #

    plot_number = pl.gcf().number
    return plot_number



def get_ranges(spectra, ytype, fraction=0.05):

    """


    """
    
    norders = np.shape(spectra)[0]

    wranges = np.empty((norders, 2))
    yranges = np.empty((norders, 2))

    for i in range(norders):

        # Do the wavelengths

        wave = spectra[i, 0, :]
        wranges[i, :] = [np.nanmin(wave), np.nanmax(wave)]

        ndat = len(wave)
        x_values = np.arange(ndat)
        
        if ytype == 'flux':

            array = spectra[i, 1, :]
            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spec_range(sg_array, frac=fraction)

        if ytype == 'uncertainty':

            array = spectra[i, 2, :]
            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spec_range(sg_array, frac=fraction)

        if ytype == 'snr':

            array = spectra[i, 1, :]/spectra[i, 2, :]
            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spec_range(sg_array, frac=fraction)

        if ytype == 'flux and uncertainty':

            array1 = spectra[i, 1, :]
            sg_array1 = robust_savgol(x_values, array1, 11)['fit']
            
            array2 = spectra[i, 2, :]
            sg_array2 = robust_savgol(x_values, array2, 11)['fit']            
                    
            yranges[i, :] = get_spec_range(sg_array1, sg_array2, frac=fraction)                                                
    return wranges, yranges
    
    
