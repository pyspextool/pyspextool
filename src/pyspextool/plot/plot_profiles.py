import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
import numpy.typing as npt

from pyspextool.io.check import check_parameter
from pyspextool.plot.limits import get_spectra_range
from pyspextool.extract.extraction import make_aperture_mask


def plot_profiles(profiles:list,
                  slith_arc:int | float,
                  doorders:npt.ArrayLike,
                  aperture_positions:npt.ArrayLike=None,
                  aperture_signs:npt.ArrayLike=None,
                  aperture_radii:npt.ArrayLike=None,
                  psf_radius=None,
                  bg_regions:str=None,
                  bg_annulus:list=None,
                  profile_size:tuple=(6,4),
                  profilestack_max:int=4,                  
                  font_size:int=12,               
                  output_fullpath:str=None,
                  showblock:bool=False,
                  showscale:float | int=None,
                  plot_number:int=None):

    """
    To plot the profiles.

    Parameters
    ----------
    profiles : list
        A (norders, ) list where each element is a dictionary with the
        following keys:

        `"order"` : int
            The order number.
    
        `"angles"` : ndarray
            An (nangles, ) array of angles on the sky...

        `"profile"` : ndarray
            An (nangles, ) array giving the mean spatial profile.

    slith_arc : int or float
        The slit height in arcseconds.

    doorders : ndarray of int
        An (norders, ) array of order numbers.

    aperture_positions : ndarray
        An (norders, naps) array of aperture positions (in arcseconds).

    aperture_radii : ndarray
        An (norders, naps) array of aperture radii (in arcseconds).

    psf_radius : int or float
        The PSF radius used for optimal extraction (in arcseconds).

    """

    #
    #  Check parameters
    #

    check_parameter('plot_profiles', 'profiles', profiles, 'list')

    check_parameter('plot_profiles', 'slith_arc', slith_arc,['int', 'float'])

    check_parameter('plot_profiles', 'doorders', doorders, 'ndarray')
    
    if output_fullpath is not None:

        # Write the plot to disk.

        doplot(None,
               profile_size,
               profilestack_max,
               font_size,
               profiles,
               slith_arc,
               doorders,
               aperture_positions=aperture_positions,
               aperture_signs=aperture_signs,               
               aperture_radii=aperture_radii,
               psf_radius=psf_radius,
               bg_regions=bg_regions,
               bg_annulus=bg_annulus)

        pl.savefig(output_fullpath)
        pl.close()
    
    if output_fullpath is None:
    
        # Display the image to the screen.


        doplot(plot_number,
               (profile_size[0]*showscale,profile_size[1]*showscale),
               profilestack_max,
               font_size*showscale,
               profiles,
               slith_arc,
               doorders,
               aperture_positions=aperture_positions,
               aperture_signs=aperture_signs,                              
               aperture_radii=aperture_radii,
               psf_radius=psf_radius,
               bg_regions=bg_regions,
               bg_annulus=bg_annulus)               
        
        pl.show(block=showblock)
        if showblock is False: pl.pause(1)


       
def doplot(plot_number,
           profile_size,
           profilestack_max,
           font_size,
           profiles,
           slith_arc,
           doorders,
           aperture_positions=None,
           aperture_signs=None,                          
           aperture_radii=None,
           psf_radius=None,
           bg_regions=None,
           bg_annulus=None):
#
#    """
#    Executes the plot indpendent of the device
#
#    Parameters
#    ----------
#
#
#    Returns
#    -------
#
#
#    """
#

    # Set the fonts

    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)


    
# Determine the plot size
    
    norders = len(profiles)

    ncols = np.ceil(norders / profilestack_max).astype(int)

    nrows = np.min([norders,profilestack_max]).astype(int)

    plot_index = np.arange(1,nrows*ncols+1)
    
    plot_index = np.reshape(np.reshape(plot_index,(nrows,ncols)),
                            ncols*nrows,order='F')
    
    figure_size = (profile_size[0]*ncols, profile_size[1]*nrows)

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
      
    if aperture_positions is not None:

        naps = np.shape(aperture_positions)[1]

    for i in range(norders):

        profile = profiles[norders-i-1]

        

        if doorders[norders-i-1] == 0:
            plot_color = 'grey'
            profile_color = 'grey'


        else:

            plot_color = 'black'
            profile_color = 'black'
            aperture_color = 'violet'
            apradii_color = 'green'
            psfradius_color = 'blue'
            bg_color = 'red'

        # Get the plot range

        yrange = get_spectra_range(profile['profile'], frac=0.2)

        # Plot the profile

        axe = pl.subplot(nrows, ncols, plot_index[i])

        axe.step(profile['angles'],
                 profile['profile'],
                 color=profile_color,
                 where='mid')

        axe.spines['left'].set_color(plot_color)
        axe.spines['right'].set_color(plot_color)
        axe.spines['top'].set_color(plot_color)
        axe.spines['bottom'].set_color(plot_color)
        axe.xaxis.label.set_color(plot_color)
        axe.yaxis.label.set_color(plot_color)
        axe.tick_params(colors=plot_color, which='both')
        axe.set_xlim([0, slith_arc])
        axe.set_ylim(yrange)
        axe.set_ylabel('Relative Intensity')
        axe.set_xlabel('Slit Position (arcseconds)')        
        axe.set_title('Order ' + str(profiles[norders-i-1]['order']),
                      color=plot_color)

        axe.xaxis.set_minor_locator(AutoMinorLocator())    
        axe.tick_params(right=True, left=True, top=True, bottom=True,
                          which='both', direction='in', width=1.5)
        axe.tick_params(which='minor', length=3)
        axe.tick_params(which='major', length=5)
        axe.yaxis.set_minor_locator(AutoMinorLocator())
        

        if doorders[norders-i-1] == 0:
            continue

        if aperture_positions is not None:

            # Now start the aperture loop

            for j in range(naps):

                positions = aperture_positions[norders-i-1, :]
                axe.vlines(positions[j], yrange[0], yrange[1],
                           color=aperture_color)

                if aperture_radii is not None:

                    radii = aperture_radii[norders-i-1, :]
                    mask = make_aperture_mask(profile['angles'],
                                              positions,
                                              radii,
                                              bg_annulus=bg_annulus,
                                              bg_regions=bg_regions)

                    z = (mask > float(j)) & (mask <= float(j + 1))
                                       
                    # We want to plot the apertures, so set all other pixels
                    # to NaN

                    tmp = np.copy(profile['profile'])
                    tmp[~z] = np.nan

                    # Plot it, then fill it
                    
                    axe.step(profile['angles'], tmp,
                             color=apradii_color,
                             where='mid')


                    if aperture_signs is not None:

                        sign = aperture_signs[j]
                        y2 = np.nanmin(tmp) if sign == 1 else np.nanmax(tmp)

                        
                        axe.fill_between(profile['angles'],
                                         tmp,
                                         step='mid',
                                         y2=np.full(len(tmp),y2),
                                         alpha=0.4,
                                         color=apradii_color)
                    
                    # We want to plot the background, so set all other pixels
                    # to NaN

                    z = mask != -1
                    tmp = profile['profile'] * 1
                    tmp[z] = np.nan
                    axe.plot(profile['angles'], tmp, color=bg_color)

                    # Now add the psf_radius if given
                    
                if psf_radius is not None:

                    axe.vlines(positions[j] - psf_radius, yrange[0],
                               yrange[1], color=psfradius_color,
                               linestyle='dotted')
                    axe.vlines(positions[j] + psf_radius, yrange[0],
                               yrange[1], color=psfradius_color,
                               linestyle='dotted')

    # fix layout overlap
    
    pl.tight_layout()


