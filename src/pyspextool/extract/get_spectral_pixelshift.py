import numpy as np
import scipy
from scipy import signal
from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.plot.limits import get_spec_range
import matplotlib.pyplot as pl
import os


def get_spectral_pixelshift(xanchor, yanchor, xsource, ysource,
                            savitzky_golay=True, qafileinfo=None):
    """
    To determine the pixel shift between two spectra.

    Parameters
    ----------
    xanchor : array-like of int
        (nanchor,) array of pixel values.

    yanchor : array-like
        (nanchor,) array of intensity values.

    xsource : array-like of int
        (nsource,) array of pixel values.

    ysource : array-like
        (nsource,) array of intensity values.

    savitzky_golay : {False, True}, optional
        Set to smooth the `yanchor` and `ysource` with a savitzky-golay
        filter before cross correlation to protect against bad pixels.

    qafileinfo : dict, optional
        `"figsize"` : tuple
            (2,) tuple of the figure size (inches).

        `"filepath"` : str
            The directory to write the QA figure.

        `"filename"` : str
            The name of the file, sans suffix/extension.

        `"extension"` : str
            The file extension.  Must be compatible with the savefig
            function of matplotlib.

    Returns
    -------
    int
        The shift of `ysource` relative to `yanchor` in pixels.
        if `ysource` is shifted to the right of `yanchor` the returned
        value is positive.

    Notes
    -----
    Uses scipy.signal.correlate for the cross correlation.

    The savitzky-golay parameters are window_length=5 and polydeg=2.

    Examples
    --------
    later

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program xmc_corspec.pro.

    """

    # Convert to numpy arrays and do basic things

    xanchor = np.array(xanchor, dtype=int)
    yanchor = np.array(yanchor)
    xsource = np.array(xsource, dtype=int)
    ysource = np.array(ysource)

    ndat = len(xanchor)

    # Find the intersection

    junk, zanchor, zsource = np.intersect1d(xanchor, xsource,
                                            return_indices=True)

    # clip the results

    xanchor = xanchor[zanchor]
    yanchor = yanchor[zanchor]

    xsource = xsource[zsource]
    ysource = ysource[zsource]

    # Savitzky-Golay the results to protect against bad pixels

    if savitzky_golay is not False:

        sgyanchor = robust_savgol(xanchor, yanchor, 5)['fit']
        sgysource = robust_savgol(xsource, ysource, 5)['fit']

    else:

        sgyanchor = yanchor
        sgysource = ysource

    # Run the cross correlation

    xcor = scipy.signal.correlate(sgysource, sgyanchor,
                                  mode='same', method='fft')
    xcor = xcor / np.nanmax(xcor)

    lag = scipy.signal.correlation_lags(ndat, ndat, mode='same')

    # Now lets find the fit window

    # Determine the pixel when the derivative becomes positive again

    maxidx = np.argmax(xcor)

    dif = -1
    npix = 0
    while dif < 0 and maxidx + npix + 1 < ndat - 1:
        dif = xcor[maxidx + npix + 1] - xcor[maxidx + npix]
        npix += 1

    halfwin = np.around(npix * 1.5).astype('int')

    # Clip out the fit zone

    fitlag = lag[maxidx - halfwin:maxidx + halfwin]
    fitxcor = xcor[maxidx - halfwin:maxidx + halfwin]

    # Do the fit

    r = fit_peak1d(fitlag, fitxcor, nparms=4, positive=True)
    offset = r['parms'][1]

    if qafileinfo is not None:

        # Normalize the spectra for plotting purposes

        np.divide(ysource, np.median(ysource), out=ysource)
        np.divide(yanchor, np.median(yanchor), out=yanchor)

        if 'figsize' in qafileinfo.keys():

            figsize = qafileinfo['figsize']

        else:
            figsize = (8.5, 11)

        # A 3 panel figure

        #pl.rcParams['font.size'] = '8'
        #pl.rcParams['lines.linewidth'] = 0.75

        fig, (axes1, axes2, axes3) = pl.subplots(3, figsize=figsize,
                                                 constrained_layout=False)
        pl.subplots_adjust(hspace=0.5)

        # Plot the two spectra

        yrange = get_spec_range([yanchor, ysource], frac=0.1)

        axes1.margins(x=0)
        axes1.step(xanchor, yanchor, '#1f77b4')
        axes1.step(xsource, ysource, 'r')
        axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
        axes1.set(xlabel='Column Number', ylabel='Relative Intensity')

        axes1.text(0.95, 0.8, 'anchor', color='#1f77b4', ha='right',
                   transform=axes1.transAxes)

        axes1.text(0.95, 0.7, 'source', color='r', ha='right',
                   transform=axes1.transAxes)

        # Plot the entire cross correlation results

        yrange = get_spec_range(xcor, frac=0.1)

        axes2.margins(x=0)
        axes2.tick_params(axis='x')
        axes2.tick_params(axis='y')
        axes2.set_title('Cross Correlation')
        axes2.step(lag, xcor)
        axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
        axes2.set(xlabel='Lag (pixels)', ylabel='Relative Intensity')

        # Plot a zoom in of the cross correlation and the fit        

        yrange = get_spec_range([fitxcor, r['fit']], frac=0.1)

        axes3.margins(x=0)

        axes3.step(fitlag, fitxcor)
        axes3.step(fitlag, r['fit'], 'r')
        axes3.set_title('Fit of Cross Correlation')
        axes3.set_ylim(ymin=yrange[0], ymax=yrange[1])

        axes3.set(xlabel='Offset (pixels)', ylabel='Relative Intensity')
        axes3.axvline(x=offset, linestyle='dotted', color='r')

        axes3.text(0.95, 0.8, 'offset=' + "{:.1f}".format(offset) + ' pixels',
                   ha='right', c='r', transform=axes3.transAxes)

        # Save the figure

        pl.savefig(os.path.join(qafileinfo['filepath'],
                                qafileinfo['filename']) + \
                   '_shift'+qafileinfo['extension'])

    return offset
