import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)

from pyspextool.io.check import check_parameter

def plot_mergedorders(plot_number:int,
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

    check_parameter('plot_mergedorders', 'plot_number', plot_number,
                    ['int', 'NoneType'])

    check_parameter('plot_mergedorders', 'subplot_size', subplot_size, 'tuple')

    check_parameter('plot_mergedorders', 'subplot_stackmax', subplot_stackmax, 
                    'int')

    check_parameter('plot_mergedorders', 'font_size', font_size, 'int')

    check_parameter('plot_mergedorders', 'scale', scale, ['int','float'])

    check_parameter('plot_mergedorders', 'spectrum_linewidth', 
                    spectrum_linewidth, ['int','float'])        

    check_parameter('plot_mergedorders', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        

    check_parameter('plot_mergedorders', 'xlabel', xlabel, 'str')
    

