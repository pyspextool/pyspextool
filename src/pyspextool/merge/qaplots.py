import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
import matplotlib
from matplotlib import rc
from matplotlib.ticker import (AutoMinorLocator)

from pyspextool.io.check import check_parameter
from pyspextool.merge.core import get_spectra_position
from pyspextool.plot.limits import get_spectra_range

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
                aperture:int,
                merged_spectra:npt.ArrayLike,
                n_overlaps:float=1,
                reverse_order:bool=True):

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

    check_parameter('plot_merges', 'aperture', aperture, 'int')

    check_parameter('plot_merges', 'merged_spectra', merged_spectra, 
                    'list')    
    #
    # Get important information
    #

    # Get norders and napertures for the object spectra.
    
    norders = len(orders)

    shape = np.shape(spectra)

    napertures = int(shape[0]/norders)

    #
    # Now start the plotting
    #
    
    # Set the fonts
    # removed helvetica

    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size*scale}

    rc('font', **font)
    
    # Determine the plot size
    
    ncols = np.ceil((norders-1)*napertures / subplot_stackmax).astype(int)

    nrows = np.min([(norders-1)*napertures,subplot_stackmax]).astype(int)

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

    #
    # Start the loops over order and aperture.  The range is norders-1 because 
    # you are plotting pairs and don't need to loop over the last order.
    #

    m = 0
    for i in range(norders-1):

        # Find the correct order index

        if reverse_order is True:
            
            addorder_idx = norders*napertures-i*napertures-napertures
            anchororder_idx = addorder_idx-napertures

        else:
            
            addorder_idx = i*napertures
            anchororder_idx = addorder_idx+napertuers

        for j in range(napertures):

            addaperture_idx = addorder_idx + j
            anchoraperture_idx = anchororder_idx + j
            
            merged_wavelength = merged_spectra[j][0,:]
            merged_intensity = merged_spectra[j][1,:]
    
            # Get the two spectra

            add_wavelength = spectra[addaperture_idx,0,:]
            add_intensity = spectra[addaperture_idx,1,:]

            anchor_wavelength = spectra[anchoraperture_idx,0,:]
            anchor_intensity = spectra[anchoraperture_idx,1,:]

            pos, overlap = get_spectra_position(anchor_wavelength,
                                                add_wavelength)
            
            delta = overlap[1]-overlap[0]
            xrange = [overlap[0]-delta*0.2,overlap[1]+delta*0.2]

            z1 = np.where((add_wavelength > xrange[0]) & 
                          (add_wavelength < xrange[1]))[0]
            z2 = np.where((anchor_wavelength > xrange[0]) & 
                          (anchor_wavelength < xrange[1]))[0]


            yrange = get_spectra_range(add_intensity[z1],anchor_intensity[z2])

            axe = pl.subplot(nrows, ncols, plot_index[m])

            orders_addidx = norders-i-1 if reverse_order is True else i 
            orders_anchoridx = norders-i-2 if reverse_order is True else i-1 


            axe.step(anchor_wavelength[z2],
                     anchor_intensity[z2],
                     color='green',
                     where='mid',
                     lw=spectrum_linewidth,
                     label='Order '+str(orders[orders_anchoridx]))

            axe.step(add_wavelength[z1],
                     add_intensity[z1],
                     color='grey',
                     where='mid',
                     lw=spectrum_linewidth,
                     label='Order '+str(orders[orders_addidx]))

        # We first need to determine the index of the order number for 
        # labeling later.

            axe.set_title('Orders ' + str(orders[orders_addidx])+' & '+\
                          str(orders[orders_anchoridx])+', aperture '+\
                          str(j+1))
            axe.set_ylabel('Relative Intensity')
            axe.set_xlabel(xlabel)            




            axe.set_xlim(xrange)
            axe.set_ylim(yrange)

            axe.xaxis.set_minor_locator(AutoMinorLocator())    
            axe.tick_params(right=True, left=True, top=True, bottom=True,
                            which='both', direction='in', 
                            width=spine_linewidth)
            axe.tick_params(which='minor', length=3)
            axe.tick_params(which='major', length=5)
            axe.yaxis.set_minor_locator(AutoMinorLocator())

            # change all spines
            for axis in ['top','bottom','left','right']:
                axe.spines[axis].set_linewidth(spine_linewidth)
                
            axe.legend()

            m += 1

