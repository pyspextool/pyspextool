import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate

from pyspextool.fit.fit_peak1d import *
from pyspextool.fit.polyfit import *
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_image import plot_image
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.arrays import find_index
from pyspextool.utils.loop_progress import loop_progress

def trace_spectrum_1dxd(image, order_mask, orders, wavecal, spatcal,
                        xranges, apertures, fit_degree=2, step_size=5,
                        summation_width=5, centroid_threshold=2, fwhm=0.8,
                        clupdate=True, qafileinfo=None):

    

    #
    # Check parameters
    #

    #
    # Get set up
    #

    nrows, ncols = np.shape(image)
    
    norders, naps = np.shape(apertures)

    coeffs = np.empty((naps*norders, fit_degree+1), dtype=float)

    half = int(summation_width/2.)

    xx, yy = make_image_indices(nrows, ncols)

#    fig = pl.gcf()

    plot_x = []
    plot_y = []
    plot_goodbad = []  
    plot_fits = []
    
    #
    # Start the loop over the orders
    #
    
    for i in range(norders):

        # Define the start and end points
        
        starts = xranges[i,0]
        stops = xranges[i,1]

        # Create an array of column numbers at which you will fit
        
        numstep = int((stops-starts)/step_size)+1
        columns = np.arange(numstep)*step_size+starts

        # Define arrays you fill in
        
        peaks_pix = np.full((numstep,naps),np.nan)
        peaks_arc = np.full((numstep,naps),np.nan)
#        fit_pix = np.full((numstep,naps),np.nan)                
#        fit_arc = np.full((numstep,naps),np.nan)                
        waves = np.full(numstep, np.nan)

        # Now loop over the columns
        
        for j in range(numstep):

            # Increase S/N by combining a bunch of columns
            
            colz = np.mean(image[:, max(columns[j]-half,0):
                                 min(columns[j]+half,ncols-1)],axis=1)

            omaskz = order_mask[:, columns[j]]
            wavez = wavecal[:, columns[j]]            

            # Grab just the pixels in the slit

            z = omaskz == orders[i]

            slitz = colz[z]
            slits = spatcal[z, columns[j]]
            slity = yy[z, columns[j]]
            
            waves[j] = wavez[z][0]

            
            # Now loop over each aperture

            for k in range(naps):

                guesss = apertures[i,k]
                
                f = interpolate.interp1d(slits,slitz)
                guessz = f(guesss)

                fit = fit_peak1d(slits, slitz,
                                 p0=[guessz, guesss, fwhm/2.354, 0])

                if np.abs(fit['parms'][1]-guesss) <= centroid_threshold:

                    peaks_arc[j,k] = fit['parms'][1]
                    f = interpolate.interp1d(slits,slity)                    
                    peaks_pix[j,k] = f(fit['parms'][1])
                    
#            pl.plot(slits, slitz)
#            pl.axvline(x=apertures[i,0],color='r', linestyle='dotted')
#            pl.axvline(x=peaks_arc[j,0])
#            pl.axvline(x=apertures[i,1],color='r', linestyle='dotted')
#            pl.axvline(x=peaks_arc[j,1])            
#            pl.pause(0.01)
#            pl.clf()

        for j in range(naps):

            # Generate index number to fill in results

            l = naps*i+j

            # Fit the trace in w/s space

            fit = poly_fit_1d(waves, peaks_arc[:,j], fit_degree,
                              robust={'thresh':4, 'eps':0.1})

            coeffs[l,:] = fit['coeffs']
#            fit_arc[:,j] = fit['yfit']
                        
            # Store the results for plotting
            
            plot_x = plot_x + list(columns)
            plot_y = plot_y + list(peaks_pix[:,j])
            plot_goodbad = plot_goodbad + list(fit['goodbad'])

        # Now do the coversion between arcseconds and pixels for plotting
            
#        for j in range(numstep):
#
#            omaskz = order_mask[:, columns[j]]
#            z = omaskz == orders[i]
#
#            slits = spatcal[z, columns[j]]
#            slity = yy[z, columns[j]]            
#
#            f = interpolate.interp1d(slits,slity)
#            fit_pix[j,:] =f(fit_arc[j,:])
#
#        for j in range(naps):
#
#            plot_fits.append(np.stack([columns, fit_pix[:,j]]))


            
        if clupdate is True:
            loop_progress(i, 0, norders, message='Tracing apertures...')            
    dictionary = {'coeffs':coeffs, 'x':np.array(plot_x), 'y':np.array(plot_y),
                  'goodbad':np.array(plot_goodbad, dtype=int)}

    return dictionary

