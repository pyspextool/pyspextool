import numpy as np
import matplotlib.pyplot as pl

from pyspextool.plot.limits import get_spec_range

def plot_profiles(profiles,slith_arc, doorders, apertures=None):

    norders = len(profiles)

    pl.figure(figsize=(8.5,8))
    pl.subplots_adjust(hspace=2)

    print(doorders)
    for i in range(norders):

        if doorders[i] == 0:
            plot_color='grey'
            profile_color = 'grey'
            aperture_color = 'grey'

        else:
            plot_color = 'black'
            profile_color = 'black'
            aperture_color = 'blue'            
        
        # Get the plot range

        yrange = get_spec_range(profiles[i]['p'], frac=0.2)

        # Plot the profile
        
        axe = pl.subplot(norders,1,norders-i)
        axe.plot(profiles[i]['y'], profiles[i]['p'],color=profile_color)        
        axe.spines['left'].set_color(plot_color)
        axe.spines['right'].set_color(plot_color)
        axe.spines['top'].set_color(plot_color)
        axe.spines['bottom'].set_color(plot_color)
        axe.xaxis.label.set_color(plot_color)
        axe.yaxis.label.set_color(plot_color)
        axe.tick_params(colors=plot_color, which='both')
        axe.set_xlim([0,slith_arc])
        axe.set_ylim(yrange)        
        axe.set_title('Order '+str(profiles[i]['order']),color=plot_color)

        # Plot the apertures
        
        if apertures is not None:

            axe.vlines(apertures[i,0],yrange[0], yrange[1],
                       color=aperture_color)
            axe.vlines(apertures[i,1],yrange[0], yrange[1],
                       color=aperture_color)            
        
    pl.show()

        
    
