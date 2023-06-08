import numpy as np
import matplotlib.pyplot as pl
import os

from pyspextool.plot.limits import get_spec_range
from pyspextool.extract.make_aperture_mask import make_aperture_mask


def plot_profiles(profiles, slith_arc, doorders, plot_number=None,
                  apertures=None, aperture_radii=None, psf_radius=None,
                  ps_bginfo=None, xs_bginfo=None, file_info=None,
                  plot_size=(6, 10)):
    """
    To plot the profiles.

    Parameters
    ----------

    file_info : dict, optional
        `"figsize"` : tuple
            (2,) tuple of the figure size (inches).

        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

    """

    #
    #  Check parameters
    #

    
    
    # Get basic information



    
    if file_info is None:

        # This is to the screen

        pl.ion()    
        plot_number = doplot(profiles, slith_arc, doorders, plot_size,
                             plot_number, apertures=apertures,
                             aperture_radii=aperture_radii,
                             psf_radius=psf_radius, ps_bginfo=ps_bginfo,
                             xs_bginfo=xs_bginfo)
        pl.show()
        pl.pause(1)
        return plot_number

    else:

        # This is to a file

        pl.ioff()
        doplot(profiles, slith_arc, doorders, plot_size, plot_number,
               apertures=apertures, aperture_radii=aperture_radii,
               psf_radius=psf_radius, ps_bginfo=ps_bginfo, xs_bginfo=xs_bginfo)

        pl.savefig(os.path.join(file_info['filepath'],
                                file_info['filename'] +
                                file_info['extension']))
        pl.close()
        return None
        

    
def doplot(profiles, slith_arc, doorders, plot_size, plot_number,
           apertures=None, aperture_radii=None, psf_radius=None,
           ps_bginfo=None, xs_bginfo=None):

    """
    Executes the plot indpendent of the device

    Parameters
    ----------


    Returns
    -------


    """

    # Get norders and napertures
    
    norders = len(profiles)

    if apertures is not None:

        if len(doorders) == 1:

            naps = len(apertures[0])

        else:

            naps = np.shape(apertures)[1]

    #
    # Make the figure
    #
    
    pl.figure(num=plot_number, figsize=plot_size)
    pl.subplots_adjust(hspace=2)

    for i in range(norders):

        if doorders[i] == 0:
            plot_color = 'grey'
            profile_color = 'grey'


        else:

            plot_color = 'black'
            profile_color = 'black'
            aperture_color = 'cyan'
            apradii_color = 'green'
            psfradius_color = 'blue'
            bg_color = 'red'

        # Get the plot range

        yrange = get_spec_range(profiles[i]['profile'], frac=0.2)

        # Plot the profile

        axe = pl.subplot(norders, 1, norders - i)
        axe.plot(profiles[i]['angle'], profiles[i]['profile'],
                 color=profile_color)
        axe.spines['left'].set_color(plot_color)
        axe.spines['right'].set_color(plot_color)
        axe.spines['top'].set_color(plot_color)
        axe.spines['bottom'].set_color(plot_color)
        axe.xaxis.label.set_color(plot_color)
        axe.yaxis.label.set_color(plot_color)
        axe.tick_params(colors=plot_color, which='both')
        axe.set_xlim([0, slith_arc])
        axe.set_ylim(yrange)
        axe.set_title('Order ' + str(profiles[i]['order']), color=plot_color)

        if doorders[i] == 0:
            continue

        if apertures is not None:

            # Now start the aperture loop

            for j in range(naps):

                axe.vlines(apertures[i, j], yrange[0], yrange[1],
                           color=aperture_color)

                if aperture_radii is not None:
                    mask = make_aperture_mask(profiles[i]['angle'],
                                              np.squeeze(apertures[i, :]),
                                              aperture_radii,
                                              psbginfo=ps_bginfo,
                                              xsbginfo=xs_bginfo)

                    # We want to plot the apertures, so set all other pixels
                    # to NaN

                    z = mask <= 0.0
                    tmp = profiles[i]['profile'] * 1
                    tmp[z] = np.nan

                    axe.plot(profiles[i]['angle'], tmp, color=apradii_color)

                    # We want to plot the background, so set all other pixels
                    # to NaN

                    z = mask != -1
                    tmp = profiles[i]['profile'] * 1
                    tmp[z] = np.nan
                    axe.plot(profiles[i]['angle'], tmp, color=bg_color)

                if psf_radius is not None:
                    axe.vlines(apertures[i, j] - psf_radius, yrange[0],
                               yrange[1], color=psfradius_color,
                               linestyle='dotted')
                    axe.vlines(apertures[i, j] + psf_radius, yrange[0],
                               yrange[1], color=psfradius_color,
                               linestyle='dotted')


    #
    # Get plot number
    #

    plot_number = pl.gcf().number
    return plot_number
