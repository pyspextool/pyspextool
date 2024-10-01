import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import matplotlib
from matplotlib.ticker import (AutoMinorLocator)



from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.plot.limits import get_stack_range
from pyspextool.io.check import check_parameter


def plot_allorders(plot_number,
                   figure_size,
                   font_size,
                   plot_scalerange:tuple=False,
                   title:str=None):

    """
    To plot all the spectra

    Parameters
    ----------
    None

    Returns
    -------
    None
        

    """

    #
    # Check the parameters
    #



    check_parameter('plot_allorders', 'plot_number', plot_number, 'int')

    check_parameter('plot_allorders', 'figure_size', figure_size, 'tuple')

    check_parameter('plot_allorders', 'plot_scalerange', plot_scalerange,
                    'bool')

    check_parameter('plot_allorders', 'title', title, ['str','NoneType'])
    
    #
    # Rename variables for ease of use
    #

    wavelength = combine.state['wavelengths']
    intensity = combine.state['intensities']
    uncertainty = combine.state['uncertainties']
    snr = intensity/uncertainty
    bitmask = combine.state['bitmasks']

    #
    # Get the ranges
    #

    wavelength_range = [np.nanmin(wavelength),np.nanmax(wavelength)]

    wavelength_ranges = np.empty((combine.state['norders'], 2))    
    intensity_ranges = np.empty((combine.state['norders'], 2))
    uncertainty_ranges = np.empty((combine.state['norders'], 2))
    snr_ranges = np.empty((combine.state['norders'], 2))            

    # Do the loop

    for i in range(combine.state['norders']):

        wavelength_ranges[i,:] = np.array([np.nanmin(wavelength[0,i]),
                                           np.nanmax(wavelength[0,i])])

        intensity_ranges[i,:] = get_stack_range(intensity[0,i], savgol=True,
                                                frac=0.05)

        uncertainty_ranges[i,:] = get_stack_range(uncertainty[0,i], savgol=True,
                                                frac=0.05)        

        snr_ranges[i,:] = get_stack_range(snr[0,i], savgol=True, frac=0.05)
        
        
    # Find the minimum and maximum values

    intensity_range = [np.min(intensity_ranges), np.max(intensity_ranges)]
    snr_range = [np.min(snr_ranges), np.max(snr_ranges)]    

    #
    # Create the figure
    #

    #
    # Set the fonts
    #
    
    font = {'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    matplotlib.rc('font', **font)


    
    fig = pl.figure(num=plot_number,
                    figsize=figure_size)


    pl.clf()
    pl.subplots_adjust(left=0.1,
                       bottom=0.1, 
                       right=0.95, 
                       top=0.9, 
                       hspace=0.05)


    #
    # plot the flux
    #
        
    ax = fig.add_subplot(211)    
    ax.set_xlim(wavelength_range)
    ax.set_ylim(intensity_range)
#    ax.set_xlabel(combine.state['xlabel'])
    ax.set_ylabel(combine.state['ylabel'])

    ax.xaxis.set_minor_locator(AutoMinorLocator())    
    ax.tick_params(right=True, left=True, top=True, bottom=True,
                    which='both', direction='in',
                   width=setup.plots['spine_linewidth'])
    ax.tick_params(which='minor', length=3)
    ax.tick_params(which='major', length=5)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
        
    ax.get_xaxis().set_ticklabels([])


    
    axis_to_data = ax.transAxes + ax.transData.inverted()
    data_to_axis = axis_to_data.inverted()


    for i in range(combine.state['norders']):

        color = iter(cm.rainbow(np.linspace(0, 1, combine.state['nspectra'])))
        for j in range(combine.state['nspectra']):

            c = next(color)
            ax.step(wavelength[0,i,:], intensity[0,i,j,:],c=c,
                    lw=setup.plots['spectrum_linewidth'])                
#            pl.plot(wavelength[0,i,:], uncertainty[0,i,j,:],c='gray',
#                    lw=setup.plots['linewidth'])            

        # Get the order number position

        ave_wave = np.mean(wavelength_ranges[i, :])
        xnorm = data_to_axis.transform((ave_wave, 0))[0]

        # Print the text

        ax.text(xnorm, 1.01 + 0.01 * (i % 2), str(combine.state['orders'][i]),
                color='black', transform=ax.transAxes)

    # Add the title

    if title is not None:

        ax.set_title(title, pad=20.0)

    if plot_scalerange is True:

        ax.axvline(x=combine.state['scale_range'][0],color='0',ls='--')
        ax.axvline(x=combine.state['scale_range'][1],color='0',ls='--')        

    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(setup.plots['spine_linewidth'])


        
    #
    # Do the SNR plot
    #

        

    ax = fig.add_subplot(212)    
    ax.set_xlim(wavelength_range)
    ax.set_ylim(snr_range)
    ax.set_xlabel(combine.state['xlabel'])
    ax.set_ylabel('Signal-to-Noise Ratio')

    ax.xaxis.set_minor_locator(AutoMinorLocator())    
    ax.tick_params(right=True, left=True, top=True, bottom=True,
                    which='both', direction='in',
                   width=setup.plots['spine_linewidth'])
    ax.tick_params(which='minor', length=3)
    ax.tick_params(which='major', length=5)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
        



    
    axis_to_data = ax.transAxes + ax.transData.inverted()
    data_to_axis = axis_to_data.inverted()

    for i in range(combine.state['norders']):

        color = iter(cm.rainbow(np.linspace(0, 1, combine.state['nspectra'])))
        for j in range(combine.state['nspectra']):

            c = next(color)
            ax.step(wavelength[0,i,:], snr[0,i,j,:],c=c,
                    lw=setup.plots['spectrum_linewidth'])                
    
        # Get the order number position

        ave_wave = np.mean(wavelength_ranges[i, :])
        xnorm = data_to_axis.transform((ave_wave, 0))[0]


    if plot_scalerange is True:

        ax.axvline(x=combine.state['scale_range'][0],color='0',ls='--')
        ax.axvline(x=combine.state['scale_range'][1],color='0',ls='--')        

    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(setup.plots['spine_linewidth'])


        
