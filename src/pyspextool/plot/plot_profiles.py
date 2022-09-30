import numpy as np
import matplotlib.pyplot as pl

from pyspextool.plot.limits import get_spec_range

def plot_profiles(profiles,slith_arc, apertures=None):

    norders = len(profiles)

    pl.figure(figsize=(8.5,8))
    pl.subplots_adjust(hspace=2)
    
    for i in range(norders):

        # Get the plot range

        yrange = get_spec_range(profiles[i]['p'], frac=0.2)

        # Plot the profile
        
        axe = pl.subplot(norders,1,norders-i)
        axe.plot(profiles[i]['y'], profiles[i]['p'],color='black')        
        axe.set_xlim([0,slith_arc])
        axe.set_ylim(yrange)        
        axe.set_title('Order '+str(profiles[i]['order']))

        # Plot the apertures
        
        if apertures is not None:

            axe.vlines(apertures[i,0],yrange[0], yrange[1], color='blue')
            axe.vlines(apertures[i,1],yrange[0], yrange[1], color='blue')            
        
    pl.show()

        
    
