import numpy as np
import matplotlib.pyplot as pl
import os

from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_fit_2d

def wavecal_solution_1d(orders, line_info, dispersion_degree, xdinfo=None,
                        verbose=True, qa_plotsize=(5,8), qa_plot=None,
                        qa_fileinfo=None):
    """
    To calculate a wavelength solution in the pyspextool 1D or 1DXD case.

    Parameters
    ----------
    orders : ndarray
        A ndarray giving the order numbers.

    line_info : dict
        `"order"` : ndarray
             A (nlines,) ndarray giving the order number of each line.

        `"x"` : ndarray
             A (nlines,) ndarray giving the colunn of each line.

        `"wavelength"` : ndarray
             A (nlines,) ndarray giving the wavelengths of each line.

        `"goodbad"` : ndarray
             A (nlines,) ndarray.  0 = bad, 1=good.
        
    dispersion_degree : int
        The polynomial degree to fit in the dispersion direction.

    xdinfo : dict, optional
        `"homeorder"' : int
            The home order used to scale lines.  See Notes.

        '"orderdeg"' : int
            The polynomial order to fit in the order dimension.

    verbose : {None, True, False}, optional
        Set to True/False to override config.setup['verbose']

    qa_plot : {None, True, False}, optional
        Set to True/False to override config.setup['qa_plot'].  If set to True,
        quality assurance plots will be interactively generated.

    qa_plotsize : tuple, default=(5,8)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_fileinfo : dict, optional
        `"figsize"` : tuple
            A (2,) tuple giving the figure size.

        `"fullpath"` : str
            A string giving the path to write the file.

        '"filename"` : str
            The root of the filename.

        '"extension"` : str
            The file extension.

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

    """

    #
    # We just need to test whether this is a 1D or 1DXD fit
    #
    
    if xdinfo is None: # This is a 1D fit

        # Do the fit
        
        fit = poly_fit_1d(line_info['x'],
                          line_info['wavelength'].astype(np.float16),
                          dispersion_degree,
                          goodbad=line_info['goodbad'],
                          robust={'thresh': 4, 'eps': 0.1})

        # Do the QA plotting
        
        residuals = (line_info['wavelength'].astype(np.float16) -
                     fit['yfit']) * 1e4
            
        if qa_plot is True:

            pl.ion()
            do_1dplot(qa_plotsize, residuals, line_info['x'], fit['goodbad'],
                      fit['rms'], dispersion_degree)
            pl.show()
            pl.pause(1)
        
        if qa_fileinfo is not None:

            pl.ioff()
            do_1dplot(qa_fileinfo['figsize'], residuals, line_info['x'],
                      fit['goodbad'], fit['rms'], dispersion_degree)
            pl.savefig(os.path.join(qa_fileinfo['filepath'],
                                    qa_fileinfo['filename'] +
                                    '_residuals'+qa_fileinfo['extension']))
            pl.close()
        

    else: # This is a 1DXD fit

        # Scale the wavelengths to the home order

        scales = line_info['order'] / xdinfo['homeorder']
        scaled_wavelengths = line_info['wavelength'].astype(np.float16) * scales

        # Now do the fit
    
        fit = poly_fit_2d(line_info['x'], line_info['order'],
                          scaled_wavelengths, dispersion_degree,
                          xdinfo['orderdeg'], goodbad=line_info['goodbad'],
                          robust={'thresh': 4, 'eps': 0.1})
        
        # Now do the qa plot

        residuals = (scaled_wavelengths - fit['zfit']) * 1e4

        if qa_plot is True:

            pl.ion()
            do_1dxdplot(qa_plotsize, residuals, line_info['order'],
                        line_info['x'], orders, fit['goodbad'], fit['rms'],
                        dispersion_degree, xdinfo['orderdeg'])
            pl.show()
            pl.pause(1)

        if qa_fileinfo is not None:

            pl.ioff()
            do_1dxdplot(qa_fileinfo['figsize'], residuals, line_info['order'],
                        line_info['x'], orders, fit['goodbad'], fit['rms'],
                        dispersion_degree, xdinfo['orderdeg'])
          
            pl.savefig(os.path.join(qa_fileinfo['filepath'],
                                    qa_fileinfo['filename'] +
                                    '_residuals'+qa_fileinfo['extension']))
            pl.close()
                            
    #
    # Now get the results together and return the results
    #
    
    # Get the number of good and bad points

    zgood = fit['goodbad'] == 1
    ntot = np.sum(zgood)
            
    zbad = fit['goodbad'] == 0
    nbad = np.sum(zbad)

    # Create the dictionary
    
    solution = {'coeffs':fit['coeffs'], 'covar':fit['coeffs_covar'],
                'rms':fit['rms'], 'nlines':ntot, 'ngood':ntot-nbad, 'nbad':nbad}
                
    return solution



def do_1dxdplot(figsize, residuals, residual_orders, residual_columns, orders,
                goodbad, rms, dispersion_degree, order_degree):


    fig = pl.figure(figsize=figsize)

    axes1 = fig.add_subplot(211)
    axes2 = fig.add_subplot(212)        
    
#    fig, (axes1, axes2) = pl.subplots(2, figsize=figsize,
#                                      constrained_layout=False)

    # Get the colormap set up
            
    norders = len(orders)
    blue = np.linspace(0, 1, num=norders)
    red = 1 - blue
            
    # Get the plot range 
            
    yrange = [-10 * rms * 1e4, 10 * rms * 1e4]
            
    # Get the number of good and bad points for the label

    zgood = goodbad == 1
    ntot = np.sum(zgood)
            
    zbad = goodbad == 0
    nbad = np.sum(zbad)

    # Loop over each order so you can color them properly
            
    for i in range(norders):
        z = residual_orders == orders[i]
                
        # Residual as a function of order number
                
        axes1.axhline(y=0, linestyle='dotted', color='black')

        axes1.plot(residual_orders[z], residuals[z], 'o', 
                   markerfacecolor=(red[i], 0, blue[i]),
                   markeredgecolor='black',
                   markersize=8, alpha=0.8)                                                   
        bad = goodbad == 0
                
        axes1.plot(residual_orders[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes1.set(xlabel='Order Number', ylabel='Data-Fit ($\mathrm{\AA}$)')
                
        axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
        # residuals as a function of column number 
                
        axes2.axhline(y=0, linestyle='dotted', color='black')

        axes2.plot(residual_columns[z], residuals[z], 'o', 
                   markerfacecolor=(red[i], 0, blue[i]),
                   markeredgecolor='black',
                   markersize=8, alpha=0.8)                                                   
        axes2.plot(residual_columns[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes2.set(xlabel='Column (pixels)',
                  ylabel='Data-Fit ($\mathrm{\AA}$)')
                
        axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
        # Label the fit
                
        label = '$\lambda$ degree=' + str(dispersion_degree) + \
                '\n order degree=' + str(order_degree) + \
                '\n RMS=' + "{:#.3g}".format(rms * 1e4) + \
                ' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
                ', n$_\mathrm{bad}$=' + str(nbad)
                  
        axes1.text(0.02, 0.95, label, color='black', ha='left',
                   va='top', multialignment='left', transform=axes1.transAxes)
                
        label = 'Not all bad data points may be shown.'
        axes1.text(0.98, 0.95, label, color='black', ha='right',
                   va='top', transform=axes1.transAxes)    

    

def do_1dplot(figsize, residuals, residual_columns, goodbad, rms,
              dispersion_degree):

                
    fig, axes1 = pl.subplots(1, figsize=figsize, constrained_layout=False)

    # Get the plot range 
            
    yrange = [-10 * rms * 1e4, 10 * rms * 1e4]
            
    # Get the number of good and bad points for the label

    zgood = goodbad == 1
    ntot = np.sum(zgood)
    
    zbad = goodbad == 0
    nbad = np.sum(zbad)

    # residuals as a function of column number 
                
    axes1.axhline(y=0, linestyle='dotted', color='black')

    axes1.scatter(residual_columns[zgood], residuals[zgood],
                  color='red', edgecolors='black',
                  s=8 ** 2, alpha=0.8)
    
    if nbad !=0:
        axes1.plot(residual_columns[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes1.set(xlabel='Column (pixels)', ylabel='Residual ($\mathrm{\AA}$)')
                
        axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])

            
    # Label the fit
                
    label = '$\lambda$ degree=' + str(dispersion_degree) + \
            '\n RMS=' + "{:#.3g}".format(rms * 1e4) + \
            ' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
            ', n$_\mathrm{bad}$=' + str(nbad)
                  
    axes1.text(0.02, 0.95, label, color='black', ha='left',
               va='top', multialignment='left', transform=axes1.transAxes)
                
    label = 'Not all bad data points may be shown.'
    axes1.text(0.98, 0.95, label, color='black', ha='right',
               va='top', transform=axes1.transAxes)

            
    
