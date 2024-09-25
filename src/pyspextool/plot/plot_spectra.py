import numpy as np
import matplotlib.pyplot as pl
import numpy.typing as npt
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator


from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.plot.limits import get_spectra_range
from pyspextool.io.check import check_parameter, check_file
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.math import bit_set


def plot_spectra(file:str,
                 ytype:str='flux',
                 apertures:int=None,
                 title:str=None,
                 colors:str | list=['green','black'],
                 line_width:int | float=0.5,
                 yrange_buffer:int | float=0.05,
                 order_numbers:bool=True,
                 figure_size:tuple=(6,4),
                 font_size:int=12,               
                 output_fullpath:str=None,
                 showblock:bool=False,
                 showscale:float | int=1.0,
                 plot_number:int=None,
                 flag_linearity:bool=False,
                 flag_optimal:bool=False,
                 flag_replace:bool=False,
                 flag_fix:bool=False):

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

    block : {False, True}, optional
        Set to make the plot block access to the command line, e.g. pl.ioff().
            
    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter('plot_spectra', 'file', file, 'str')

    check_parameter('plot_spectra', 'ytype', ytype, 'str',
                    ['flux', 'uncertainty', 'snr', 'flux and uncertainty'])

    check_parameter('plot_spectra', 'apertures', apertures,
                    ['NoneType', 'int', 'list'])
    
    check_parameter('plot_spectra', 'title', title, ['NoneType', 'str'])

    check_parameter('plot_spectra', 'colors', colors,
                    ['NoneType', 'str', 'list'])

    check_parameter('plot_spectra', 'line_width', line_width, ['float', 'int'])
    
    check_parameter('plot_spectra', 'yrange_buffer', yrange_buffer, 'float')

    check_parameter('plot_spectra', 'order_numbers', order_numbers, 'bool')

    check_parameter('plot_spectra', 'figure_size', figure_size, 'tuple')

    check_parameter('plot_spectra', 'font_size', font_size, ['int','float'])

    check_parameter('plot_spectra', 'output_fullpath', output_fullpath,
                    ['NoneType','str'])            

    check_parameter('plot_spectra', 'showblock', showblock, 'bool')

    check_parameter('plot_spectra', 'plot_number', plot_number,
                    ['NoneType', 'int'])    

    check_parameter('plot_spectra', 'flag_linearity', flag_linearity, 'bool')
    
    check_parameter('plot_spectra', 'flag_optimal', flag_optimal, 'bool')

    check_parameter('plot_spectra', 'flag_replace', flag_replace, 'bool')

    check_parameter('plot_spectra', 'flag_fix', flag_fix, 'bool')    
    
    #
    # Read the file
    #

    check_file(file)
    
    spectra, info = read_spectra_fits(file)
    norders = info['norders']
    napertures = info['napertures']
    latex_yunits = info['lyunits']
    latex_xlabel = info['lxlabel']
    latex_ylabel = info['lylabel']

    orders = info['orders'] if order_numbers is True else None
    
    #
    # Get aperture parameter ready for use
    #

    if apertures is None:

        # Plot them all
        
        plot_apertures = list(range(1,1+napertures))

    else:

        plot_apertures = [apertures] if isinstance(apertures,int) else apertures

    #
    # Check to make the request is is doable
    #

    napertures_plot = len(plot_apertures)
    if napertures_plot > 4:

        message = 'Cannot plot more than four apertures.'
        raise pySpextoolError(message)
        
    #
    # Do the plot
    #

    if output_fullpath is not None:

        # Write the plot to disk.
        
        doplot(spectra,
               norders,
               napertures,
               plot_apertures,
               ytype,
               colors,
               line_width,
               latex_xlabel,
               latex_ylabel,
               latex_yunits,
               yrange_buffer,
               figure_size,
               font_size,            
               None,
               orders=orders,
               title=title,
               flag_linearity=flag_linearity,
               flag_fix=flag_fix,
               flag_replace=flag_replace,
               flag_optimal=flag_optimal)
        
       
        pl.savefig(output_fullpath)
        pl.close()
        
    if output_fullpath is None:
    
        # Display the image to the screen.

        doplot(spectra,
               norders,
               napertures,
               plot_apertures,
               ytype,
               colors,
               line_width,
               latex_xlabel,
               latex_ylabel,
               latex_yunits,
               yrange_buffer,
               (figure_size[0]*showscale,figure_size[1]*showscale),
               font_size*showscale,            
               plot_number,
               orders=orders,
               title=title,
               flag_linearity=flag_linearity,
               flag_fix=flag_fix,
               flag_replace=flag_replace,
               flag_optimal=flag_optimal)
               
        
        pl.show(block=showblock)
        if showblock is False: pl.pause(1)


def doplot(spectra:npt.ArrayLike,
           norders:int,
           napertures:int,
           plot_apertures:list,
           ytype:str,
           colors:str | list,
           line_width:int | float,
           latex_xlabel:str,
           latex_ylabel:str,
           latex_yunits:str,
           yrange_buffer:int | float,
           figure_size:tuple,
           font_size:int,
           plot_number:int,
           orders:list=None,
           title:str=None,
           flag_linearity:bool=False,
           flag_fix:bool=False,
           flag_replace:bool=False,
           flag_optimal:bool=False):

    """
    Actually does the plot independent of device.

    


    """

    #
    # Figure out the number of rows and columns
    #

    nplot = len(plot_apertures)

    nrows = 1
    ncols = 1

    if nplot >= 2:

        ncols = 2

    if nplot >=3 :

        nrows = 2

    #
    # Get the figure started
    #
        
    # Set the fonts

    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)

    if ytype == 'flux':

        ylabel = latex_ylabel

    elif ytype == 'uncertainty':

        ylabel = 'Uncertainty (' + latex_yunits + ')'

    elif ytype == 'snr':

        ylabel = 'Signal to Noise'


    elif ytype == 'flux and uncertainty':

        ylabel = latex_ylabel       
    
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
                       bottom=0.125,
                       top=0.9)

    for i in range(nplot):

        idx = np.arange(norders)*napertures+i
        
        wranges, yranges = get_ranges(spectra[idx,:,:], ytype,
                                      ybuffer_fraction=yrange_buffer)

        xrange = [np.nanmin(wranges), np.nanmax(wranges)]
        yrange = [np.nanmin(yranges), np.nanmax(yranges)]
        
        axe = pl.subplot(nrows, ncols, i+1)
        axe.plot([np.nan],[np.nan])

        axe.plot([np.nan], [np.nan])
        axe.set_xlim(xrange)
        axe.set_ylim(yrange)
        axe.set_xlabel(latex_xlabel)
        axe.set_ylabel(ylabel)

        axe.xaxis.set_minor_locator(AutoMinorLocator())    
        axe.tick_params(right=True, left=True, top=True, bottom=True,
                          which='both', direction='in', width=1.5)
        axe.tick_params(which='minor', length=3)
        axe.tick_params(which='major', length=5)
        axe.yaxis.set_minor_locator(AutoMinorLocator())
        
        
        # Get the coordinate transform ready for order number labels
    
        axis_to_data = axe.transAxes + axe.transData.inverted()
        data_to_axis = axis_to_data.inverted()

        for j in range(norders):
            
            # Get the plotting values
            
            xvalues = spectra[j*napertures+i, 0, :]
            
            if ytype == 'flux':
                
                yvalues = spectra[j*napertures+i,1,:]
                y2values = [np.nan]*len(yvalues)
                
            elif ytype == 'uncertainty':
                
                yvalues = spectra[j*napertures+i,2,:]
                y2values = [np.nan]*len(yvalues)
                
            elif ytype == 'snr':
                
                yvalues = spectra[j*napertures+i,1,:]/\
                    spectra[i*napertures+j,2,:]
                y2values = [np.nan]*len(yvalues)
                
            elif ytype == 'flux and uncertainty':
                
                yvalues = spectra[j*napertures+i,1,:]
                y2values = spectra[j*napertures+i,2,:]
                
            # Get the colors
                
            if isinstance(colors, list):
                    
                ncolors = len(colors)
                    
                if ncolors == 2:
                        
                    color = colors[0] if j % 2 == 0 else colors[1]
                        
                elif ncolors == 3:
                        
                    mod = i % 3
                    
                    if mod == 0: color = colors[0]
                    
                    if mod == 1: color = colors[1]
                    
                    if mod == 2:  color = colors[2]
                        
            else:
                
                color = colors
                            
            # Plot the spectrum
        
            axe.step(xvalues, yvalues, color=color, ls='-',
                     lw=line_width, where='mid')
            axe.step(xvalues, y2values, color='grey', lw=line_width,
                     where='mid')

            # Now label the order numbers

            if orders is not None:

                # Get the order number position

                ave_wave = np.mean(wranges[j, :])
                xnorm = data_to_axis.transform((ave_wave, 0))[0]
                
                # Print the text
                
                axe.text(xnorm, 1.01 + 0.01 * (j % 2), str(orders[j]),
                         color=color, transform=axe.transAxes)
                
                # Add the title
                
            if title is not None:
                
                axe.set_title(title+' - Aperture '+str(i+1), pad=20.0)
                
            else:

                axe.set_title('Aperture '+str(i+1), pad=20.0)

            np.set_printoptions(threshold=np.inf)                
            bitmask = spectra[j*napertures+i,3,:]
            bitmask[np.isnan(bitmask)] = 0

            if flag_linearity is True:
                
                mask = bit_set(bitmask.astype(int),0)
                z  = mask == 1
                axe.plot(xvalues[z], yvalues[z], marker='.', marksize=2,
                         color='red', linestyle='')
                
            if flag_replace is True:

                mask = bit_set(bitmask.astype(int),1)
                z  = mask == 1
                axe.plot(xvalues[z], yvalues[z], marker='.',
                         markersize=2, color='blue', linestyle='')
                

            
            
def get_ranges(spectra:npt.ArrayLike,
               ytype:str,
               ybuffer_fraction:int | float=0.05):

    """
    To determine the plot range of the spectra

    Parameters
    ----------
    spectra : ndarray
        An (norders*naps, 4, nwavelength) array of spectra.

    ytype : {'flux', 'uncertainty', 'snr', 'flux and uncertainty'}
        A string with the type of spectrum to plot:

    ybuffer_fraction : int or float, default 0.05
        The 

    Returns
    -------
    ndarray, ndarray
    
    
    """

    #
    # Check parameters
    #

    check_parameter('get_ranges', 'spectra', spectra, 'ndarray')

    check_parameter('get_ranges', 'ytype', ytype, 'str')

    check_parameter('get_ranges', 'ybuffer_fraction', ybuffer_fraction,
                    ['int','float'])        

    #
    # Get set up
    #
    
    norders = np.shape(spectra)[0]
    wranges = np.empty((norders, 2))
    yranges = np.empty((norders, 2))

    #
    # Start the loop over each order
    #
    
    for i in range(norders):

        # Do the wavelengths

        wave = spectra[i, 0, :]
        wranges[i, :] = [np.nanmin(wave), np.nanmax(wave)]

        ndat = len(wave)
        x_values = np.arange(ndat)

        # Now do the "flux"
        
        if ytype == 'flux':

#            array = spectra[i, 1, :]
#            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spectra_range(spectra[i, 1, :], robust=True,
                                              frac=ybuffer_fraction)

        if ytype == 'uncertainty':

#            array = spectra[i, 2, :]
#            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spectra_range(spectra[i, 2, :], robust=True,
                                              frac=ybuffer_fraction)

        if ytype == 'snr':

#            array = spectra[i, 1, :]/spectra[i, 2, :]
#            sg_array = robust_savgol(x_values, array, 11)['fit']
            
            yranges[i, :] = get_spectra_range(spectra[i, 1, :]/spectra[i, 2, :],
                                              robust=True,
                                              frac=ybuffer_fraction)

        if ytype == 'flux and uncertainty':

            array1 = spectra[i, 1, :]
            sg_array1 = robust_savgol(x_values, array1, 11)['fit']
            
            array2 = spectra[i, 2, :]
            sg_array2 = robust_savgol(x_values, array2, 11)['fit']            
                    
            yranges[i, :] = get_spectra_range(sg_array1, sg_array2,
                                           frac=ybuffer_fraction)
            
    return wranges, yranges
    
    
