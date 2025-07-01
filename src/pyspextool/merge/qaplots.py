import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)

from pyspextool.io.check import check_parameter

def plot_merges(plot_number:int,
                subplot_size:tuple,
                subplot_stackmax:int,
                font_size:int,
                scale:int | float,
                spectrum_linewidth:int | float,
                spine_linewidth:int | float,
                xlabel:str,
                spectra:npt.ArrayLike,
                orders:npt.ArrayLike,
                merged_spectrum:npt.ArrayLike):

    """
    To create a QA plot for order merging

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
      
    Returns
    -------
    None
    
    """

    #
    # Check parameters
    #

    check_parameter('plot_merges', 'plot_number', plot_number, 
                    ['int', 'NoneType'])

    check_parameter('plot_merges', 'subplot_size', subplot_size, 'tuple')

    check_parameter('plot_merges', 'subplot_stackmax', subplot_stackmax, 'int')

    check_parameter('plot_merges', 'font_size', font_size, 'int')

    check_parameter('plot_merges', 'scale', scale, ['int','float'])

    check_parameter('plot_merges', 'spectrum_linewidth', 
                    spectrum_linewidth, ['int','float'])        

    check_parameter('plot_merges', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        

    check_parameter('plot_merges', 'xlabel', xlabel, 'str')

    check_parameter('plot_merges', 'spectra', spectra, 'ndarray')

    check_parameter('plot_merges', 'orders', orders, 'ndarray')

    check_parameter('plot_merges', 'merged_spectrum', merged_spectrum, 
                    'ndarray')
    
    #
    # Get important information
    #

    # Get norders and napertures for the object spectra.
    
    norders = len(orders)

    #
    # Now start the plotting
    #
    
    # Set the fonts

    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size*scale}

    rc('font', **font)
    
    # Determine the plot size
    
    ncols = np.ceil(norders / subplot_stackmax).astype(int)

    nrows = np.min([norders,subplot_stackmax]).astype(int)

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
    for i in range(norders-1):

        print(i)



    

