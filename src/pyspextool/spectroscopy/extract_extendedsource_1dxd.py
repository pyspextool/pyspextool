import numpy as np

from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.spectroscopy.make_aperture_mask import make_aperture_mask
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.arrays import trim_nan
from pyspextool.utils.loop_progress import loop_progress



def extract_extendedsource_1dxd(img, var, ordermask, orders, wavecal, spatcal,
                                appos, apradii, linmax_bitmask=None,
                                badpixel_mask=None, bginfo=None, clupdate=True):
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

    wavecal : numpy.ndarray
        (nrows,ncols) image where each pixel is set to its wavelength.

    spatcal : numpy.ndarray
        (nrows, ncols) array where each pixel is set to its angular position
        on the sky (arcseconds).

    appos : numpy.ndarray
        (naps,) array of aperture positions (arcseconds).

    apradius : numpy.ndarray
        (naps,) array of aperture radii (arcseconds).

    linmax_bitmask : numpy.ndarray, optional
        (nrows, ncols) array where a pixel is set to unity if its value is
        beyond the linearity maximum.

    badpixel_mask : numpy.ndarray, optional
        (nrows, ncols) array where good pixels are set to unity and bad pixels
        are set to zero.

    bginfo : later when you implement it.

    clupdate : {True, False}, optional
        Set to True for command line updates during execution.

    Returns
    -------
    list
        The value of each entry is a (4, nwave) numpy.ndarray where:
        wave = (0,:)
        intensity = (1,:)
        uncertainty = (2,:)
        flags = (3,:)

    Notes
    -----
    None

    Examples
    --------
    later?

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
        bgdeg = bginfo['deg']

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

    output_list = []
    for i in range(norders):

        if clupdate is not None and i == 0:
            print('Extracting apertures...')

        zordr = np.where(ordermask == orders[i])
        xmin = np.min(xx[zordr])
        xmax = np.max(xx[zordr])
        nwaves = xmax - xmin + 1

        # Create the output arrays

        owave = np.full(nwaves, np.nan)
        oflux = np.full((naps, nwaves), np.nan)
        ounc = np.full((naps, nwaves), np.nan)
        omask = np.zeros((naps, nwaves), dtype=int)

        # Start the wavelength loop

        for j in range(nwaves):

            # get the slit mask

            colmask = ordermask[:, xmin + j]
            zslit = colmask == orders[i]

            # Store the wavelength

            owave[j] = wavecal[zslit, xmin + j][0]

            # Carve out the slit values            

            slit_pix = yy[zslit, xmin + j]
            slit_arc = spatcal[zslit, xmin + j]
            slit_img = img[zslit, xmin + j]
            slit_var = var[zslit, xmin + j]
            slit_bpm = badpixel_mask[zslit, xmin + j]
            slit_lmm = linmax_bitmask[zslit, xmin + j]

            # Generate the aperture mask

            slitmask = make_aperture_mask(slit_arc, appos[i,:], apradii,
                                          xsbginfo=bgregions)
            
            # Do the background subtraction

            if bgregions is not None:

                z = (slitmask == -1)
                result = poly_fit_1d(slit_arc[z], slit_img[z],bgdeg,
                                     robust={'thresh':4, 'eps':0.1},
                                     silent=True)

                
                

            # Do the sum extraction

            for k in range(naps):

                # Find the apertures
                
                z = (slitmask > float(k)) & (slitmask <= float(k + 1))

                # Create partial array
                partial = slitmask[z] - float(k)

                # Now do the sums
                
                oflux[k, j] = np.sum(slit_img[z] * partial)
                varval = np.sum(slit_var[z] * partial**2)
                ounc[k, j] = np.sqrt(np.abs(varval))

            # Check the bitmask for linearity

            z = slit_lmm == 1
            if sum(z) is True:
                omask[k, j] = 1

        # Store the results

        # Trim the NaNs if need be

        nonan = trim_nan(owave, flag=2)

        # Generate the key

        for k in range(naps):
            arr = np.stack((owave[nonan], oflux[k, nonan], ounc[k, nonan],
                            omask[k, nonan]))

            output_list.append(arr)

        if clupdate is not None:
            loop_progress(i, 0, norders)

    return output_list
