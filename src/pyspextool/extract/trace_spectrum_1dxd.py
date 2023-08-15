from scipy import interpolate

from pyspextool.fit.fit_peak1d import *
from pyspextool.fit.polyfit import *
from pyspextool.io.check import check_parameter
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.loop_progress import loop_progress


def trace_spectrum_1dxd(image, order_mask, orders, wavecal, spatcal,
                        xranges, apertures, fit_degree=2, step_size=5,
                        summation_width=5, centroid_threshold=2, fwhm=0.8,
                        verbose=True):
    """
    To trace a spectrum.

    Parameters
    ----------
    image : numpy.ndarray
        A float (nrows, ncols) image with (cross-dispersed) spectral
        orders.  It is assumed that the dispersion direction is roughly
        aligned with the rows of `image` and the spatial axis is roughly
        aligned with the columns of `image`.  That is, orders go
        left-right and not up-down.

    order_mask : numpy.ndarray
        An int (nrows, cols) order mask.  A pixel is set to its
        order number.  Interorder pixels are set to zero.

    orders : numpy.ndarray
        An int (norders,) array of order numbers to be traced.

    wavecal : numpy.ndarray
        A float (nrows, cols) wavelength image.  A pixel is set to its
        wavelength.  Interorder pixels are set to NaN.

    spatcal : numpy.ndarray
        A float (nrows, cols) spatial image.  A pixel is set to its
        spatial position.  Interorder pixels are set to NaN.

    xranges : numpy.ndarray
        A float (norders, 2) array giving the column numbers over which
        to operate.  xranges[0,0] gives the starting column number for
        the order nearest the bottom of the image and xranges[0,1] gives
        the end column number for said order.

    apertures : numpy.ndarray
        A float (norders, naps) array of aperture positions.  By Spextool
        convention, apertures[0,0] is the aperture closest to the bottom of
        `image`.

    fit_degree : int, optional, default 2
        The polynomial degree for the fit.

    step_size : int, optional, default 5
        The step size as the function moves across the array identifying
        the peaks in the spatial dimension.

    summation_width : int, optional, default 5
        The number of columns to combine in order to increase the S/N
        of the data fit to identify the peaks in the spatial dimension.

    centroid_threshold : int, optional, default 2
        If (fit-guess) > centroid_threshold to peak identification is found
        to fail.

    fwhm : float, optional, default 0.8
        The approximate FWHM of the peaks in arcseconds.

    verbose : {True, False}, optional
        Set to True for command line updates during execution.

    Returns
    -------
    dictionary
        `"coeffs"` : numpy.ndarray
            A float (norders*naps,fit_degree+1) array giving the polynomial
            coefficients for each aperture.  By Spextool convention,
            dictionary['coeffs'][0,:] is the aperture closest to the bottom
            of `image`.

        `"x"` : numpy.ndarray
            A float (ndata,) array of all the x positions for the fits.

        `"y"` : numpy.ndarray
            A float (ndata,) array of all the y positions for the fits.

        `"goodbad"` : numpy.ndarray
            A int (ndata,) goodbad array.  1=good data, 0=data identified
            as an outlier by the fit.

    Notes
    -----
    Loops over each order and aperture and fits a gaussian to the slit
    at the aperture position in `apertures` with a FWHM of `fwhm`.
    Then fits the spatial positions (typically in arcseconds) as a function
    of wavelength to determine the trace coeffients stored in
    dictionary['tracecoeffs'].  The function also converts the spatial
    positions into (x,y) positions on the array for later plotting.

    Examples
    --------
    later

    """

    #
    # Check parameters
    #

    check_parameter('trace_spectrum_1dxd', 'image', image, 'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'order_mask', order_mask,
                    'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'orders', orders, 'ndarray', 1)

    check_parameter('trace_spectrum_1dxd', 'wavecal', wavecal, 'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'spatcal', spatcal, 'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'xranges', xranges, 'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'apertures', apertures, 'ndarray', 2)

    check_parameter('trace_spectrum_1dxd', 'fit_degree', fit_degree, 'int')

    check_parameter('trace_spectrum_1dxd', 'step_size', step_size, 'int')

    check_parameter('trace_spectrum_1dxd', 'summation_width', summation_width,
                    'int')

    check_parameter('trace_spectrum_1dxd', 'centroid_threshold',
                    centroid_threshold, 'int')

    check_parameter('trace_spectrum_1dxd', 'fwhm', fwhm, 'float')

    check_parameter('trace_spectrum_1dxd', 'verbose', verbose, 'bool')

    #
    # Get set up
    #

    nrows, ncols = np.shape(image)

    norders, naps = np.shape(apertures)
    coeffs = np.empty((naps * norders, fit_degree + 1), dtype=float)

    half = int(summation_width / 2.)

    xx, yy = make_image_indices(nrows, ncols)

    # intialize plotting stuff
    plot_x = []
    plot_y = []
    plot_goodbad = []

    #
    # Start the loop over the orders
    #

    #    pl.clf()

    for i in range(norders):

        # Define the start and end points

        starts = xranges[i, 0]
        stops = xranges[i, 1]

        # Create an array of column numbers at which you will fit

        numstep = int((stops - starts) / step_size) + 1
        columns = np.arange(numstep) * step_size + starts

        # Define arrays you fill in

        peaks_pix = np.full((numstep, naps), np.nan)
        peaks_arc = np.full((numstep, naps), np.nan)
        waves = np.full(numstep, np.nan)

        # Now loop over the columns

        for j in range(numstep):

            # Increase S/N by combining a bunch of columns

            colz = np.mean(image[:, max(columns[j] - half, 0):
                                    min(columns[j] + half, ncols - 1)], axis=1)

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

                guesss = apertures[i, k]

                f = interpolate.interp1d(slits, slitz)
                guessz = f(guesss)

                fit = fit_peak1d(slits, slitz,
                                 p0=[guessz, guesss, fwhm / 2.354, 0],
                                 ignore_optimizewarning=True)
               
                if not fit['goodfit']:
                    continue

                if np.abs(fit['parms'][1] - guesss) <= centroid_threshold:
                    peaks_arc[j, k] = fit['parms'][1]
                    f = interpolate.interp1d(slits, slity)
                    peaks_pix[j, k] = f(fit['parms'][1])

        #            title = 'Order ' + str(orders[i]) + ', Column ' + str(columns[j])
        #                    pl.plot(slits, slitz)
        #                    pl.axvline(x=apertures[i,0],color='r', linestyle='dotted')
        #                    pl.axvline(x=peaks_arc[j,0])
        #                    pl.axvline(x=apertures[i,1],color='r', linestyle='dotted')
        #                    pl.axvline(x=peaks_arc[j,1])
        #                    pl.title(title)
        #                    pl.pause(0.001)
        #                    pl.clf()

        for j in range(naps):
            # Generate index number to fill in results

            l = naps * i + j

            # Fit the trace in w/s space

            fit = poly_fit_1d(waves, peaks_arc[:, j], fit_degree,
                              robust={'thresh': 4, 'eps': 0.1})

            coeffs[l, :] = fit['coeffs']

            # Store the results for plotting

            plot_x = plot_x + list(columns)
            plot_y = plot_y + list(peaks_pix[:, j])
            plot_goodbad = plot_goodbad + list(fit['goodbad'])

        if verbose is True:
            loop_progress(i, 0, norders, message='Tracing apertures...')

    dictionary = {'coeffs': coeffs, 'x': np.array(plot_x),
                  'y': np.array(plot_y),
                  'goodbad': np.array(plot_goodbad, dtype=int)}

    return dictionary
