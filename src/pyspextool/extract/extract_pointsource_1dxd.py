import numpy as np
from scipy.interpolate import interp1d

from pyspextool.extract.make_aperture_mask import make_aperture_mask
from pyspextool.extract.profiles import make_2d_profile
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.fit.polyfit import poly_1d
from pyspextool.io.check import check_parameter
from pyspextool.utils.math import moments
from pyspextool.utils.arrays import make_image_indices
from pyspextool.utils.arrays import trim_nan
from pyspextool.utils.loop_progress import loop_progress

def extract_pointsource_1dxd(image, variance, order_mask, orders, wavecal,
                             spatcal, trace_coefficients, aperture_radius,
                             aperture_sign, linmax_bitmask=None,
                             background_info=None, optimal_info=None,
                             badpixel_info=None, verbose=True):
    
    """
    To extract a point source.

    Parameters
    ----------
    image : numpy.ndarray
        An (nrows, ncols) image with spectral orders.  It is assumed
        that the dispersion direction is roughly aligned with the
        rows of `img` and the spatial axis is roughly aligned with
        the columns of `img`.  That is, orders go left-right and not
        up-down.

    variance : ndarray
        (nrows, ncols) variance image for `image`.

    order_mask : ndarray
        (nrows, cols) image where each pixel is set to its order number.  
        Inter-order pixels are set to zero.

    orders : ndarray of int
        (norders,) array of the order numbers.  By Spextool convention,
        orders[0] is the order closest to the bottom of the array.

    wavecal : ndarray
        (nrows,ncols) image where each pixel is set to its wavelength.

    spatcal : ndarray
        (nrows, ncols) array where each pixel is set to its angular position
        on the sky (arcseconds).

    trace_coefficients : ndarray
        An (norder*nap,) array of trace coefficients.

    aperture_radius : int or float
        The aperture radius.

    aperture_sign : ndarray
        An (naps,) array of aperture signs.  1=positive, -1=negative

    linmax_bitmask : ndarray or None, optional
        (nrows, ncols) array where a pixel is set to unity if its value is
        beyond the linearity maximum.

    background_info : dict or None, optional

        '`radius`: float or int
            The background radius value.

        '`width`: float or int
            The background width value.

        '`degree`: int
            The polynomial degree for the background fit.

    optimal_info : dict or None, optional
        '`images`' : ndarray
            An (norders,) array of rectified order dictionaries.  Each entry has
                '`image`': ndarray
                    An (nangles, nwavelengths) array of an rectified order.

                '`angle`': ndarray
                    An (nangles,) array of spatial angles.

                '`wavelength`': ndarray
                    An (nangles,) array of wavelengths.
                
        '`psfradius`' : int or float
            The radius for the PSF.

        '`usemeanprofile`' : {False, True}    
             Set to True to use the mean profile for all wavelengths.

        '`mask`' : ndarray
             An (nrows, ncols) array giving the location of bad pixels. 
             1 = good
             0 = bad

        '`thresh`' : float or int
             The sigma threshold to identify bad pixels.
    
    badpixel_info : dict or None, optional
        '`images`' : ndarray
            An (norders,) array of rectified order dictionaries.  Each entry has
                '`image`': ndarray
                    An (nangles, nwavelengths) array of an rectified order.

                '`angle`': ndarray
                    An (nangles,) array of spatial angles.

                '`wavelength`': ndarray
                    An (nangles,) array of wavelengths.

        '`usemeanprofile`' : {False, True}    
             Set to True to use the mean profile for all wavelengths.

        '`mask`' : ndarray
             An (nrows, ncols) array giving the location of bad pixels. 
             1 = good
             0 = bad

        '`thresh`' : float or int
             The sigma threshold to identify bad pixels.    

    verbose : {True, False}, optional
        Set to True for command line updates during execution.

    Returns
    -------
    list
        A (norders,) list. Each entry is a (4, nwave) numpy.ndarray where:
          wave = (0,:)
          intensity = (1,:)
          uncertainty = (2,:)
          flags = (3,:)

    """

    #
    # Check parameters
    #

    check_parameter('extract_pointsource_1d', 'image', image, 'ndarray')

    check_parameter('extract_pointsource_1d', 'variance', variance, 'ndarray')

    check_parameter('extract_pointsource_1d', 'order_mask', order_mask,
                    'ndarray')

    check_parameter('extract_pointsource_1d', 'orders', orders, 'ndarray')

    check_parameter('extract_pointsource_1d', 'wavecal', wavecal, 'ndarray')

    check_parameter('extract_pointsource_1d', 'spatcal', spatcal, 'ndarray')

    check_parameter('extract_pointsource_1d', 'trace_coefficients',
                    trace_coefficients, 'ndarray')

    check_parameter('extract_pointsource_1d', 'trace_coefficients',
                    trace_coefficients, 'ndarray')

    check_parameter('extract_pointsource_1d', 'aperture_radius',
                    aperture_radius, ['float','int'])

    check_parameter('extract_pointsource_1d', 'aperture_sign',
                    aperture_sign, 'ndarray')

    check_parameter('extract_pointsource_1d', 'linmax_bitmask',
                    linmax_bitmask, ['ndarray','NoneType'])

    check_parameter('extract_pointsource_1d', 'background_info',
                    background_info, ['dict','NoneType'])

    check_parameter('extract_pointsource_1d', 'optimal_info',
                    optimal_info, ['dict','NoneType'])

    check_parameter('extract_pointsource_1d', 'badpixel_info',
                    badpixel_info, ['dict','NoneType'])

    check_parameter('extract_pointsource_1d', 'verbose', verbose, 'bool')                                    

    #
    # Do basic things
    #
    
    nrows, ncols = image.shape

    norders = len(orders)

    # Convert to numpy arrays to deal with float/int versus list/array problem

    naps = len(aperture_sign)

    #
    # Deal with bginfo
    #
    
    if background_info is not None:

        psbginfo = (background_info['radius'], background_info['width'])
        bgdeg = background_info['degree']

    else:

        psbginfo = None
        bgdeg = None

    # Create pixel coordinate arrays

    xx, yy = make_image_indices(nrows, ncols)

    #  Set up linearity maximum mask

    if linmax_bitmask is None:
        linmax_bitmask = np.zeros((nrows, ncols)).astype(int)

    #  Set up bad pixel mask.  Use the oportunity to get the first part of the
    #  verbose message and set the use_profile variable.  

    if badpixel_info is None and optimal_info is None:

        badpixel_mask = np.ones((nrows, ncols)).astype(int)
        use_profile = False
        text = 'Sum extracting '

    else:

        use_profile = True

        if badpixel_info is not None:

            info = badpixel_info
            text = 'Sum extracting '
            
        if optimal_info is not None:

            info = optimal_info            
            text = 'Optimally extracting '
            
        badpixel_mask = info['mask']
        rectified_orders = info['images']
        atmosphere = info['atmosphere']
        thresh = info['thresh']
        use_mean_profile = info['usemeanprofile']
        
    # Start the order loop

    if verbose is True:
                
        message = text+str(naps)+' apertures in '+str(norders)+' orders'

        if psbginfo is not None:

            print(message + ' (with background subtraction)...')

        else:

            print(message + ' (without background subtraction)...')

        
    spectrum_list = []
    background_list = []

    for i in range(norders):

        if use_profile is True:
            
            r = make_2d_profile(rectified_orders[i],
                                trace_coefficients[i*naps:i*naps+naps],
                                np.full(naps, aperture_radius),
                                atmospheric_transmission=atmosphere,
                                use_mean_profile=False)

            profile_angle = r[0]
            profile_map = r[1]

        zordr = np.where(order_mask == orders[i])
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

        # Start the wavelength loop

        for j in range(nwaves):

            # get the slit mask

            colmask = order_mask[:, xmin + j]
            zslit = colmask == orders[i]

            # Store the wavelength

            spectrum_wave[j] = wavecal[zslit, xmin + j][0]

            # Carve out the slit values            

            slit_pix = yy[zslit, xmin + j]
            slit_arc = spatcal[zslit, xmin + j]
            slit_img = image[zslit, xmin + j]
            slit_var = variance[zslit, xmin + j]
            slit_bpm = badpixel_mask[zslit, xmin + j]
            slit_lmm = linmax_bitmask[zslit, xmin + j]

            if use_profile is True:

                function = interp1d(profile_angle, profile_map[:,i],
                                    fill_value=0.0, bounds_error=False)

                slit_prof = function(slit_arc)
                
            # Gernerate the slit positions using the tracecoeffs
        
            appos = np.empty(naps)

            for k in range(naps):
                l = i * naps + k
                wave = np.array(spectrum_wave[j], dtype='float', ndmin=1)
                appos[k] = poly_1d(wave, trace_coefficients[l])

            # Generate the aperture mask

            aperture_mask = make_aperture_mask(slit_arc, appos,
                                               np.full(naps, aperture_radius),
                                               psbginfo=psbginfo)

            # Generate the psf mask

            if optimal_info is not None:
            
                psf_mask = make_aperture_mask(slit_arc, appos,
                                    np.full(naps, optimal_info['psfradius']))
            
            #
            # Do the background subtraction
            #

            if psbginfo is not None:
                
                # Fit the background

                z_background = (aperture_mask == -1)
                result = poly_fit_1d(slit_arc[z_background],
                                     slit_img[z_background], bgdeg,
                                     robust={'thresh': 4, 'eps': 0.1},
                                     silent=True)

                # Generate a background slit 

                slit_bg, slit_bg_var = poly_1d(slit_arc, result['coeffs'],
                                               covar=result['coeffs_covar'])

                # Subtract the background and propagate the uncertainties

                slit_img = np.subtract(slit_img, slit_bg)
                slit_var = slit_var + slit_bg_var

            #
            # Scale the profile to the data
            #

            if use_profile is True:

                fit_img = poly_fit_1d(slit_prof, slit_img, 1, goodbad=slit_bpm,
                                      robust={'thresh':thresh, 'eps':0.1})

            #
            # Fix bad pixels if requested
            #

            if badpixel_info is not None:

                # We fit the variances.  Is this good?
                
                fit_var = poly_fit_1d(np.abs(slit_prof), slit_var, 1,
                                      goodbad=slit_bpm,
                                      robust={'thresh':thresh, 'eps':0.1})

                z_badpixels = fit_img['goodbad']*slit_bpm == 0

                if sum(z_badpixels) >= 1:
                    
                    slit_img[z_badpixels] = fit_img['yfit'][z_badpixels]
                    slit_var[z_badpixels] = fit_var['yfit'][z_badpixels]

            #
            # Do the extractions
            #
                    
            for k in range(naps):

                if optimal_info is not None:
                    
                    #
                    # Do optimal extraction
                    #

                    # Find the indices of the PSF pixels
                    
                    z_psf = (psf_mask > float(k)) & (psf_mask <= float(k + 1))

                    # Enforce positivity/negativity 

                    if aperture_sign[k] == 1:
                
                        slit_psf = slit_prof.clip(min=0.0)

                    else:

                        slit_psf = slit_prof.clip(max=0.0)

                    # Normalize the profile
                        
                    aperture_psf = aperture_sign[k]*\
                                   np.abs(slit_psf/np.nansum(slit_psf[z_psf]))

                    # Determine the pixels to actually use
                    
                    z_aperture = (aperture_mask > float(k)) & \
                                 (aperture_mask <= float(k + 1)) & \
                                 (fit_img['goodbad'] == 1) & \
                                 (aperture_psf != 0.0)

                    
                    if np.sum(z_aperture) > 0:

                        # Scale the data and variances
                        
                        vals = slit_img[z_aperture]/aperture_psf[z_aperture]
                        vars = slit_var[z_aperture]/aperture_psf[z_aperture]**2

                        # Compute the optimal estimate
                        
                        weights = 1/vars
                        wmean = np.sum(weights*vals)/np.sum(weights)
                        wvar = 1/np.sum(weights)

                        spectrum_flux[k, j] = wmean
                        spectrum_unc[k, j] = np.sqrt(wvar)

                else:

                    #
                    # Do the sum extraction
                    #
                        
                    # Find the indices for the aperture pixels

                    z_aperture = (aperture_mask > float(k)) & \
                                 (aperture_mask <= float(k + 1))

                    # Create partial pixel array

                    partial = aperture_mask[z_aperture] - float(k)

                    # Now do the sums

                    spectrum_flux[k, j] = np.sum(slit_img[z_aperture]*partial)*\
                                          aperture_sign[k]

                    spectrum_var = np.sum(slit_var[z_aperture] * partial ** 2)
                    spectrum_unc[k, j] = np.sqrt(np.abs(spectrum_var))

                # Check the bitmask for linearity

                z_linearity = slit_lmm[z_aperture] == 1
                if sum(z_linearity) is True:
                    spectrum_mask[k, j] = 1

                # Check bad pixels
                
                if badpixel_info is not None:

                    if sum(z_badpixels) == 1:
                        spectrum_mask[k,j] += 2

        #
        # Store the results
        #
    
        # Trim the NaNs if need be

        nonan = trim_nan(spectrum_wave, flag=2)

        # Generate the key

        for k in range(naps):

            spectrum = np.stack((spectrum_wave[nonan],
                                 spectrum_flux[k, nonan],
                                 spectrum_unc[k, nonan],
                                 spectrum_mask[k, nonan]))

            spectrum_list.append(spectrum)

        if verbose is True:
            loop_progress(i, 0, norders)

    return spectrum_list
