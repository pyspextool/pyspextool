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
                   spectrum_linewidth,
                   spine_linewidth,                   
                   spectra_labels,
                   scalerange:tuple=None,
                   title:str=None):

    """
    To plot all the spectra during the combine step in a device independent way

    Parameters
    ----------
    plot_number : int
        The plot number to pass to matplotlib

    figure_size : tuple
        A (2,) tuple giving the figure size to pass to matplotlib

    font_size : int
        An int giving the font size to pass to matplotlib

    spectrum_linewidth : int or float
        An int or float giving the spectrum line width to pass to matplotlib

    spine_linewidth : int or float
        An int or float giving the spine line width to pass to matplotlib
    
    spectrum_labels : list
        An (nspectra,) list giving the file names for the spectra.

    spectra_labels : list
        An (nspectra,) list giving the labels for each spectrum.
    
    scalerange : ndarray, default None
        If given, a (2,) ndarray giving the wavlength range over which the
        scale factors were determined.

    title : string, default None
        If given, a string giving the title of the plot.
    
    Returns
    -------
    None
        Plots 

    """

    #
    # Check the parameters
    #

    check_parameter('plot_allorders', 'plot_number', plot_number, 'int')

    check_parameter('plot_allorders', 'figure_size', figure_size, 'tuple')

    check_parameter('plot_allorders', 'font_size', font_size, ['int','float'])

    check_parameter('plot_allorders', 'spectrum_linewidth', spectrum_linewidth,
                    ['int','float'])        

    check_parameter('plot_allorders', 'spine_linewidth', spine_linewidth,
                    ['int','float'])        

    check_parameter('plot_allorders', 'spectra_labels', spectra_labels,
                    'list')        
    
    check_parameter('plot_allorders', 'scalerange', scalerange,
                    ['list','ndarray', 'NoneType'])

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
    
    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
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


    for i in range(combine.state['final_napertures']):

        #
        # plot the flux
        #
 
        idx = i+1
        ax = fig.add_subplot(2,combine.state['final_napertures']+1,i+1)    
        ax.set_xlim(wavelength_range)
        ax.set_ylim(intensity_range)

        if i == 0:
            ax.set_ylabel(combine.state['ylabel'])

        else:

            ax.get_yaxis().set_ticklabels([])        

            
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


        for j in range(combine.state['norders']):

            color = iter(cm.rainbow(np.linspace(0, 1,
                                                combine.state['nspectra'])))
            for k in range(combine.state['nspectra']):
                
                c = next(color)
                ax.step(wavelength[0,j,:], intensity[i,j,k,:],c=c,
                        lw=spectrum_linewidth)

            # Get the order number position
            
            ave_wave = np.mean(wavelength_ranges[j, :])
            xnorm = data_to_axis.transform((ave_wave, 0))[0]
            
            # Print the text

            ax.text(xnorm, 1.01 + 0.01 * (j % 2),
                    str(combine.state['orders'][j]),
                    color='black', transform=ax.transAxes)

        # change all spines
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(spine_linewidth)            
                
        ax.text(0.9, 0.9, 'Ap '+str(i+1), color='black', ha='right',
                transform=ax.transAxes)

        if scalerange is not None:
            
            ax.axvline(x=scalerange[0],color='0',ls='--')
            ax.axvline(x=scalerange[1],color='0',ls='--')
            
        

            

    for i in range(combine.state['final_napertures']):

        #
        # plot the uncertainty
        #

        idx = combine.state['final_napertures']+i+2
        ax = fig.add_subplot(2,combine.state['final_napertures']+1,idx)
        ax.set_xlim(wavelength_range)
        ax.set_ylim(intensity_range)

        if i == 0:
            ax.set_ylabel('Signal-to-Noise Ratio')

        else:

            ax.get_yaxis().set_ticklabels([])


        ax.set_xlim(wavelength_range)
        ax.set_ylim(snr_range)
        ax.set_xlabel(combine.state['xlabel'])

        
        ax.xaxis.set_minor_locator(AutoMinorLocator())    
        ax.tick_params(right=True, left=True, top=True, bottom=True,
                       which='both', direction='in',
                       width=setup.plots['spine_linewidth'])
        ax.tick_params(which='minor', length=3)
        ax.tick_params(which='major', length=5)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        
        
        axis_to_data = ax.transAxes + ax.transData.inverted()
        data_to_axis = axis_to_data.inverted()
        
        for j in range(combine.state['norders']):
            
            color = iter(cm.rainbow(np.linspace(0, 1,
                                                combine.state['nspectra'])))
            for k in range(combine.state['nspectra']):
                
                c = next(color)
                ax.step(wavelength[0,j,:], snr[i,j,k,:],c=c,
                        lw=spectrum_linewidth)                
                
                # Get the order number position
                
            ave_wave = np.mean(wavelength_ranges[j, :])
            xnorm = data_to_axis.transform((ave_wave, 0))[0]

        # change all spines
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(spine_linewidth)

        if scalerange is not None:
            
            ax.axvline(x=scalerange[0],color='0',ls='--')
            ax.axvline(x=scalerange[1],color='0',ls='--')


                      
            
#    fig.tight_layout()

    pl.suptitle(title)

    
    ax = fig.add_subplot(1,3,3)
    
    ax.get_yaxis().set_ticklabels([])
    ax.get_xaxis().set_ticklabels([])    

    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0)

        ax.tick_params(which='minor', length=0)
        ax.tick_params(which='major', length=0)

    colors = iter(cm.rainbow(np.linspace(0, 1,
                                        combine.state['nspectra'])))

        
    for i in range(combine.state['nspectra']):

        c = next(colors)

        ax.step(np.nan, np.nan,c=c,
                lw=1,label=spectra_labels[i])                

         
    leg = ax.legend(frameon=False)

    for color,text in zip(colors,leg.get_texts()):
        text.set_color(color)
    

        
