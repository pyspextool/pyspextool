import numpy as np
import matplotlib.pyplot as pl
import os

from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_fit_2d

def wavecal_solution_1d(orders, line_info, dispersion_degree,
                        xd=None, clupdate=True, qafileinfo=None):
    """
    To perform a 2D fit to the x position and order number of wavelength
    calibration lines.

    Parameters
    ----------
    orders:

    Returns
    -------
    numpy.ndarray
        An array of coefficients generated from poly_fit_2d.

    Notes
    -----
    The program determines a single wavelength solution for all orders
    in a spectrogram.   The technique is based on Hall et al. (1994, PASP,
    106, 315).  The ratio of two different wavelengths, lambda_1 and
    lambda_2, diffracted at the same angle but located in two different
    orders m_1 and m_2, is given by,

    lambda_1/lambda_2 = m_2/m_1



    Examples
    --------
    later

    """

    if xd is None:

        fit = poly_fit_1d(line_info['x'],
                          line_info['wavelength'].astype(np.float16),
                          dispersion_degree,
                          goodbad=line_info['goodbad'],
                          robust={'thresh': 4, 'eps': 0.1})

        if qafileinfo is not None:

            residuals = (line_info['wavelength'].astype(np.float16) -
                         fit['yfit']) * 1e4
                
            fig, axes1 = pl.subplots(1, figsize=qafileinfo['figsize'],
                                         constrained_layout=False)

            # Get the plot range 
            
            yrange = [-10 * fit['rms'] * 1e4, 10 * fit['rms'] * 1e4]
            
            # Get the number of good and bad points for the label

            zgood = fit['goodbad'] == 1
            ntot = np.sum(zgood)
            
            zbad = fit['goodbad'] == 0
            nbad = np.sum(zbad)

            # residuals as a function of column number 
                
            axes1.axhline(y=0, linestyle='dotted', color='black')
                
            axes1.scatter(line_info['x'][zgood], residuals[zgood],
                          color='red', edgecolors='black',
                          s=8 ** 2, alpha=0.8)

            if nbad !=0:
                axes1.plot(line_info['x'][zbad], residuals[zbad], 's',
                           markersize=13, markerfacecolor='none',
                           color='black')
                
            axes1.set(xlabel='Column (pixels)',
                      ylabel='Residual ($\mathrm{\AA}$)')
                
            axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])

            
            # Label the fit
                
            label = '$\lambda$ degree=' + str(dispersion_degree) + \
                    '\n RMS=' + "{:#.3g}".format(fit['rms'] * 1e4) + \
                    ' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
                    ', n$_\mathrm{bad}$=' + str(nbad)
                  
            axes1.text(0.02, 0.95, label, color='black', ha='left',
                       va='top', multialignment='left',
                       transform=axes1.transAxes)
                
            label = 'Not all bad data points may be shown.'
            axes1.text(0.98, 0.95, label, color='black', ha='right',
                       va='top', transform=axes1.transAxes)

            
            pl.savefig(os.path.join(qafileinfo['filepath'],
                                    qafileinfo['filename'] +
                                    '_residuals' +
                                    qafileinfo['extension']))            
        

    else:

        # Do the 2D fit.  First generate the scale factor for each line

        scales = line_info['order'] / xd['homeorder']
        scaled_wavelengths = line_info['wavelength'].astype(np.float16) * scales

        # Now do the fit
    
        fit = poly_fit_2d(line_info['x'], line_info['order'],
                          scaled_wavelengths, dispersion_degree, xd['orderdeg'],
                          goodbad=line_info['goodbad'],
                          robust={'thresh': 4, 'eps': 0.1})

        # Now do the qa plot

        if qafileinfo is not None:

            residuals = (scaled_wavelengths - fit['zfit']) * 1e4
                
            fig, (axes1, axes2) = pl.subplots(2, figsize=(8.5, 11),
                                            constrained_layout=False)
            # Get the colormap set up
            
            norders = len(orders)
            blue = np.linspace(0, 1, num=norders)
            red = 1 - blue
            
            # Get the plot range 
            
            yrange = [-10 * fit['rms'] * 1e4, 10 * fit['rms'] * 1e4]
            
            # Get the number of good and bad points for the label

            zgood = fit['goodbad'] == 1
            ntot = np.sum(zgood)
            
            zbad = fit['goodbad'] == 0
            nbad = np.sum(zbad)
            
            # Loop over each order so you can color them properly
            
            for i in range(norders):
                z = line_info['order'] == orders[i]
                
                # Residual as a function of order number
                
                axes1.axhline(y=0, linestyle='dotted', color='black')
                
                axes1.scatter(line_info['order'][zgood], residuals[zgood],
                            color=(red[i], 0, blue[i]), edgecolors='black',
                            s=8 ** 2, alpha=0.8)
                
#                bad = line_info['goodbad'] == 0
                
                axes1.plot(line_info['order'][zbad], residuals[zbad], 's',
                        markersize=13, markerfacecolor='none', color='black')
                
                axes1.set(xlabel='Order Number',
                        ylabel='Residual ($\mathrm{\AA}$)')
                
                axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
                # residuals as a function of column number 
                
                axes2.axhline(y=0, linestyle='dotted', color='black')
                
                axes2.scatter(line_info['x'][zgood], residuals[zgood],
                            color=(red[i], 0, blue[i]), edgecolors='black',
                            s=8 ** 2, alpha=0.8)
                
                axes2.plot(line_info['x'][zbad], residuals[zbad], 's',
                        markersize=13, markerfacecolor='none', color='black')
                
                axes2.set(xlabel='Column (pixels)',
                        ylabel='Residual ($\mathrm{\AA}$)')
                
                axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
                # Label the fit
                
                label = '$\lambda$ degree=' + str(dispersion_degree) + \
                  '\n order degree=' + str(xd['orderdeg']) + \
                  '\n RMS=' + "{:#.3g}".format(fit['rms'] * 1e4) + \
                  ' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
                  ', n$_\mathrm{bad}$=' + str(nbad)
                  
                axes1.text(0.02, 0.95, label, color='black', ha='left',
                           va='top', multialignment='left',
                           transform=axes1.transAxes)
                
                label = 'Not all bad data points may be shown.'
                axes1.text(0.98, 0.95, label, color='black', ha='right',
                           va='top', transform=axes1.transAxes)
                
                pl.savefig(os.path.join(qafileinfo['filepath'],
                                    qafileinfo['filename'] +
                                    '_residuals' +
                                    qafileinfo['extension']))
                
    solution = {'coeffs':fit['coeffs'], 'covar':fit['covar'], 'rms':fit['rms'],
                'nlines':ntot, 'ngood':ntot-nbad, 'nbad':nbad}
                
    return solution
