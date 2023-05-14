import numpy as np

from pyspextool.extract.make_aperture_mask import make_aperture_mask
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.arrays import trim_nan
from pyspextool.utils.loop_progress import loop_progress


def extract_extendedsource_1dxd(img, var, ordermask, orders, wavecal, spatcal,
                                appos, apradii, linmax_bitmask=None,
                                badpixel_mask=None, bginfo=None, verbose=True):
    """
    To extract and extended source.

    Parameters
    ----------
    img : numpy.ndarray
        An (nrows, ncols) image with spectral orders.  It is assumed
        that the dispersion direction is roughly aligned with the
        rows of `img` and the spatial axis is roughly aligned with
        the columns of `img`.  That is, orders go left-right and not
        up-down.

    var : numpy.ndarray
        (nrows, ncols) variance image for `img`.

    ordermask : numpy.ndarray
        (nrows, cols) image where each pixel is set to its
        order number.  Inter-order pixels are set to zero.

    orders : list of int
        (norders,) int array of the order numbers.  By Spextool convention,
        orders[0] is the order closest to the bottom of the array.

    wavecal : numpy.ndarray
        (nrows,ncols) image where each pixel is set to its wavelength.

    spatcal : numpy.ndarray
        (nrows, ncols) array where each pixel is set to its angular position
        on the sky (arcseconds).

    appos : numpy.ndarray
        (naps,) array of aperture positions (arcseconds).

    apradii : numpy.ndarray
        (naps,) array of aperture radii (arcseconds).

    linmax_bitmask : numpy.ndarray, optional
        (nrows, ncols) array where a pixel is set to unity if its value is
        beyond the linearity maximum.

    badpixel_mask : numpy.ndarray, optional
        (nrows, ncols) array where good pixels are set to unity and bad pixels
        are set to zero.

    bginfo : dict, optional


    verbose : {True, False}, optional
        Set to True for command line updates during execution.

    Returns
    -------
    dict

        spectrum : list
            A (norders,) list. Each entry is a (4, nwave) numpy.ndarray where:
            wave = (0,:)
            intensity = (1,:)
            uncertainty = (2,:)
            flags = (3,:)

        background : list
            A (norders,) list. Each entry is a (4, nwave) numpy.ndarray where:
            wave = (0,:)
            intensity = (1,:)
            uncertainty = (2,:)
            flags = (3,:)


        If bginfo is None, then background is set to None. 

    """

    # Do basic things

    nrows, ncols = img.shape

    norders = len(orders)

    # Convert to numpy arrays to deal with float/int versus list/array problem

    appos = np.array(appos, dtype='float', ndmin=1)
    apradii = np.array(apradii, dtype='float', ndmin=1)

    naps = len(apradii)

    # Deal with bginfo

    if bginfo is not None:

        bgregions = bginfo['regions']
        bgdeg = bginfo['degree']

    else:

        bgregions = None
        bgdeg = None

    # Create pixel coordinate arrays

    xx, yy = make_image_indices(nrows, ncols)

    #  Create linearity maximum mask

    if linmax_bitmask is None:
        linmax_bitmask = np.zeros((nrows, ncols)).astype(int)

    #  Create linearity maximum mask

    if badpixel_mask is None:
        badpixel_mask = np.ones((nrows, ncols)).astype(int)

        # Start the order loop

    spectrum_list = []
    background_list = []
    for i in range(norders):

        if verbose is True and i == 0:

            message = 'Extracting ' + str(naps) + ' apertures in ' + str(norders) + \
                      ' orders'

            if bgregions is not None:
                print(message + ' (with background subtraction)...')
            else:
                print(message + ' (without background subtraction)...')

        zordr = np.where(ordermask == orders[i])
        xmin = np.min(xx[zordr])
        xmax = np.max(xx[zordr])
        nwaves = xmax - xmin + 1

        #
        # Create the output arrays
        #

        spectrum_wave = np.full(nwaves, np.nan)
        spectrum_flux = np.full((naps, nwaves), np.nan)
        spectrum_unc = np.full((naps, nwaves), np.nan)
        spectrum_mask = np.zeros((naps, nwaves), dtype=int)

        if bgregions is not None:
            background_wave = np.full(nwaves, np.nan)
            background_flux = np.full((naps, nwaves), np.nan)
            background_unc = np.full((naps, nwaves), np.nan)
            background_mask = np.zeros((naps, nwaves), dtype=int)

            # Start the wavelength loop

        for j in range(nwaves):

            # get the slit mask

            colmask = ordermask[:, xmin + j]
            zslit = colmask == orders[i]

            # Store the wavelength

            spectrum_wave[j] = wavecal[zslit, xmin + j][0]

            # Carve out the slit values            

            slit_pix = yy[zslit, xmin + j]
            slit_arc = spatcal[zslit, xmin + j]
            slit_img = img[zslit, xmin + j]
            slit_var = var[zslit, xmin + j]
            slit_bpm = badpixel_mask[zslit, xmin + j]
            slit_lmm = linmax_bitmask[zslit, xmin + j]

            # Generate the aperture mask

            slitmask = make_aperture_mask(slit_arc, appos[i, :], apradii,
                                          xsbginfo=bgregions)

            #
            # Do the background subtraction
            #

            if bgregions is not None:
                # Fit the background

                z = (slitmask == -1)
                result = poly_fit_1d(slit_arc[z], slit_img[z], bgdeg,
                                     robust={'thresh': 4, 'eps': 0.1},
                                     silent=True)

                # Generate a background slit 

                slit_bg, slit_bg_var = poly_1d(slit_arc, result['coeffs'],
                                               covar=result['coeff_covar'])

                # Subtract the background and propagate

                slit_img = np.subtract(slit_img, slit_bg)
                slit_var = slit_var + slit_bg_var

            # Do the sum extraction

            for k in range(naps):

                # Find the apertures

                z = (slitmask > float(k)) & (slitmask <= float(k + 1))

                # Create partial array
                partial = slitmask[z] - float(k)

                #
                # Now do the sums
                #

                # Spectrum first

                spectrum_flux[k, j] = np.sum(slit_img[z] * partial)
                varval = np.sum(slit_var[z] * partial ** 2)
                spectrum_unc[k, j] = np.sqrt(np.abs(varval))

                # Background second
                if bgregions is not None:
                    background_flux[k, j] = np.sum(slit_bg[z] * partial)
                    varval = np.sum(slit_bg_var[z] * partial ** 2)
                    background_unc[k, j] = np.sqrt(np.abs(varval))
                    # Check the bitmask for linearity

            z = slit_lmm == 1
            if sum(z) is True:
                spectrum_mask[k, j] = 1

        # Store the results

        # Trim the NaNs if need be

        nonan = trim_nan(spectrum_wave, flag=2)

        # Generate the key

        for k in range(naps):
            spectrum = np.stack((spectrum_wave[nonan],
                                 spectrum_flux[k, nonan],
                                 spectrum_unc[k, nonan],
                                 spectrum_mask[k, nonan]))

            spectrum_list.append(spectrum)

            if bgregions is not None:
                background = np.stack((spectrum_wave[nonan],
                                       background_flux[k, nonan],
                                       background_unc[k, nonan],
                                       background_mask[k, nonan]))
                background_list.append(background)

        if verbose is True:
            loop_progress(i, 0, norders)

    if bgregions is not None:

        return {'spectra': spectrum_list, 'background': background_list}

    else:

        return {'spectra': spectrum_list, 'background': None}
