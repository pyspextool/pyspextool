import numpy as np
import numpy.typing as npt
from scipy.interpolate import interp1d
from astropy.io import fits
from os.path import basename
import matplotlib.pyplot as pl
from pyspextool.utils.interpolate import sinc_interpolation_fft, sinc_interpolation


from pyspextool.extract.profiles import make_2d_profile, make_aperture_mask
from pyspextool.fit.polyfit import polyfit_1d, poly_1d
from pyspextool.io.check import check_parameter
from pyspextool.utils.math import moments
from pyspextool.utils.split_text import split_text
from pyspextool.utils.arrays import make_image_indices, trim_nan, find_index
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.for_print import for_print

def extract_1dxd(image:npt.ArrayLike,
                 variance:npt.ArrayLike,
                 ordermask:npt.ArrayLike,
                 wavecal:npt.ArrayLike,
                 spatcal:npt.ArrayLike,
                 extract_orders:npt.ArrayLike,
                 trace_coefficients:npt.ArrayLike,
                 aperture_radii:npt.ArrayLike,
                 aperture_signs:npt.ArrayLike,
                 linmax_bitmask:npt.ArrayLike=None,
                 bg_annulus:list | npt.ArrayLike=None,
                 bg_regions:str=None,
                 bg_fitdegree:int=1,
                 optimal_info=None,
                 badpixel_info=None,
                 progressbar=True):
    
    """
    To extract spectra from a 1dxd formatted image.  

    Parameters
    ----------
    image : ndarray
        An (nrows, ncols) image with spectral orders.  It is assumed
        that the dispersion direction is roughly aligned with the
        rows of `img` and the spatial axis is roughly aligned with
        the columns of `img`.  That is, orders go left-right and not
        up-down.

    variance : ndarray
        (nrows, ncols) variance image for `image`.

    ordermask : ndarray
        (nrows, cols) image where each pixel is set to its order number.  
        Inter-order pixels are set to zero.

    wavecal : ndarray
        (nrows,ncols) image where each pixel is set to its wavelength.

    spatcal : ndarray
        (nrows, ncols) array where each pixel is set to its angular position
        on the sky (arcseconds).

    trace_coefficients : ndarray
        An (norder*nap,) array of trace coefficients.

    aperture_radii : ndarray
        An (norders, naps) array of aperture radii.  

    aperture_signs : ndarray
        An (naps,) array of aperture signs.  1=positive, -1=negative

    linmax_bitmask : ndarray, default None
        An (nrows, ncols) array where a pixel is set to unity if its value is
        beyond the linearity maximum.

    optimal_info : dict, default None
        '`images`' : ndarray
            An (norders,) array of rectified order dictionaries.  Each entry has

                '`image`': ndarray
                    An (nangles, nwavelengths) profile array of an rectified 
                    order.

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
                    An (nangles, nwavelengths) profile array of an rectified 
                    order.

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

    check_parameter('extract_1dxd', 'image', image, 'ndarray')

    check_parameter('extract_1dxd', 'variance', variance, 'ndarray')

    check_parameter('extract_1dxd', 'ordermask', ordermask, 'ndarray')

    check_parameter('extract_1dxd', 'wavecal', wavecal, 'ndarray')

    check_parameter('extract_1dxd', 'spatcal', spatcal, 'ndarray')

    check_parameter('extract_1dxd', 'extract_orders', extract_orders, 'ndarray')
    
    check_parameter('extract_1dxd', 'trace_coefficients', trace_coefficients,
                    'ndarray')

    check_parameter('extract_1dxd', 'aperture_radii', aperture_radii, 'ndarray')
    
    check_parameter('extract_1dxd', 'aperture_signs', aperture_signs, 'ndarray')

    check_parameter('extract_1dxd', 'linmax_bitmask', linmax_bitmask,
                    ['ndarray','NoneType'])

    check_parameter('extract_1dxd', 'bg_annulus', bg_annulus,
                    ['ndarray', 'list','NoneType'])    

    check_parameter('extract_1dxd', 'bg_regions', bg_regions,
                    ['str','NoneType'])    
    
    check_parameter('extract_1dxd', 'optimal_info', optimal_info,
                    ['dict','NoneType'])

    check_parameter('extract_1dxd', 'badpixel_info', badpixel_info,
                    ['dict','NoneType'])

    check_parameter('extract_1dxd', 'verbose',progressbar, 'bool')

    #
    # Get basic information
    #

    nrows, ncols = image.shape

    norders = len(extract_orders)

    naps = len(aperture_signs)

    background_subtraction = bg_annulus is not None or bg_regions is not None

    #
    # Get setup
    #

    # Create pixel coordinate arrays

    xx, yy = make_image_indices(nrows, ncols)
    
    #  Set up linearity maximum mask

    if linmax_bitmask is None:

        linmax_bitmask = np.zeros((nrows, ncols)).astype(int)

    #  Set up bad pixel mask.  Use the oportunity to get the first part of the
    #  verbose message, set `use_profile`, and store things

    if badpixel_info is None and optimal_info is None:

        badpixel_mask = np.ones((nrows, ncols)).astype(int)
        use_profile = False
        text = 'Sum extracting '

    else:

        # At least one of them is passed, so grab things you need

        use_profile = True

        if badpixel_info is not None:

            info = badpixel_info
            text = 'Sum extracting '
                                 
        if optimal_info is not None:

            if background_subtraction is False:

                message = ' Background subtraction is required for '+ \
                    'optimal extraction.  Please pass either `bg_annulus` '+ \
                    'or `bg_regions`.'
                raise pySpextoolError(message)

            
            info = optimal_info            
            text = 'Optimally extracting '
            
        badpixel_mask = info['mask']
        rectified_orders = info['images']
        atmosphere = info['atmosphere']
        thresh = info['thresh']
        use_meanprofile = info['usemeanprofile']

        
    #
    # Start the loop!
    #
                  
    spectrum_list = []
    background_list = []

    for i in range(norders):

        if use_profile is True:
            
            r = make_2d_profile(rectified_orders[i],
                                trace_coefficients[i*naps:i*naps+naps],
                                aperture_radii[i,:],
                                atmospheric_transmission=atmosphere,
                                use_mean_profile=False)
            
            profile_angles = r[0]
            profile_map = r[1]

        zordr = np.where(ordermask == extract_orders[i])
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

            # get the pixels in the slit

            zslit = ordermask[:, xmin + j] == extract_orders[i]

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

#                function = interp1d(profile_angles,
#                                    profile_map[:,j],
#                                    fill_value=0.0,
#                                    bounds_error=False,kind='quadratic')
#            
#                slit_prof = function(slit_arc)
                
                slit_prof = sinc_interpolation(profile_map[:,j],profile_angles, slit_arc)

            # Gernerate the aperture positions using the tracecoeffs
        
            aperture_positions = np.empty(naps)

            for k in range(naps):

                l = i * naps + k
                wave = np.array(spectrum_wave[j], dtype='float', ndmin=1)
                aperture_positions[k] = poly_1d(wave, trace_coefficients[l])[0]

            # Generate the aperture mask

            aperture_mask = make_aperture_mask(slit_arc,
                                               aperture_positions,
                                               aperture_radii[i,:],
                                               bg_annulus=bg_annulus,
                                               bg_regions=bg_regions)

            # Generate the psf mask

            if optimal_info is not None:

                psf_radii = np.full(naps, optimal_info['psfradius'])
                psf_mask = make_aperture_mask(slit_arc,
                                              aperture_positions,
                                              psf_radii)
            
            #
            # Do the background subtraction
            #

            if background_subtraction:
                
                # Fit the background

                z_background = (aperture_mask == -1)

                if bg_fitdegree == 0:

                    result = moments(slit_img[z_background],
                                     robust=4,
                                     goodbad=slit_bpm[z_background],
                                     silent=True)
                    slit_bg = np.full_like(slit_img, result['mean'])

                    slit_bgvar = np.full_like(slit_bg, result['stderr']**2)

                else:

                    result = polyfit_1d(slit_arc[z_background],
                                        slit_img[z_background],
                                        bg_fitdegree,
                                        goodbad=slit_bpm[z_background],
                                        yunc=np.sqrt(slit_var[z_background]),
                                        robust={'thresh': 4, 'eps': 0.1},
                                        silent=True)

                    # Generate a background slit 

                    slit_bg, slit_bgvar = poly_1d(slit_arc, result['coeffs'],
                                                   covar=result['coeffs_covar'],
                                                   talk=False)

                # Subtract the background and propagate the uncertainties

                slit_img -= slit_bg
                slit_var += slit_bgvar


#                if j > 150:
#                    print(j, slit_bpm)
#                    z = np.where(slit_bpm == 0)[0]
#                    fig = pl.figure(figsize=(15,10))    
#                    axes1 = pl.subplot(2,1,1)
#                    
#                    axes1.step(slit_arc, slit_img,where='mid',label='data')
#                    axes1.plot(slit_arc[z], slit_img[z],'or')
#                    axes1.step(slit_arc, slit_bg,where='mid',color='red',
#                               label='BG fits')
#                    axes1.set_title(str(j)+' '+str(spectrum_wave[j]))
#                    axes1.legend()
#                    pl.pause(0.1)
#                    pl.clf()
                    
#                pl.show()
                
                
            #
            # Scale the profile to the data
            #

            if use_profile is True:

#                if j > 192:
#
#                    silent = False
#                
#                else:
#
#                    silent = True

                fit_img = polyfit_1d(slit_prof,
                                     slit_img,
                                     1,
                                     goodbad=slit_bpm,
                                     robust={'thresh':thresh, 'eps':0.1},
                                     silent=True)

                z = np.where(fit_img['goodbad'] == 0)[0]

#                if j > 192:
#
#                    fig = pl.figure(figsize=(15,10))    
##                    fig = pl.figure(figsize=(15,10))    
##                    print(slit_bpm)
##                    print(z)
##                    print(xmin+j)
#                    axes1 = pl.subplot(2,1,1)
#
#                    scaled = poly_1d(slit_prof, fit_img['coeffs'])
#                    min = np.min([scaled,slit_img])
#                    max = np.max([scaled,slit_img])
#                    axes1.step(slit_arc, slit_img,where='mid',label='data')
#                    axes1.step(slit_arc, scaled,where='mid',color='green',
#                               label='profile')
#                    axes1.step(slit_arc[z], slit_img[z],'ro')
#                    axes1.set_ylim(min,max)
#                    axes1.set_title(str(j)+' '+str(spectrum_wave[j]))
#                    axes1.legend()
#                    pl.show()
##                    pl.draw()
##                    pl.pause(1)
##                    pl.clf()



            #
            # Fix bad pixels if requested
            #

            if badpixel_info is not None:

                # We fit the variances.  Is this good?
                
                fit_var = polyfit_1d(np.abs(slit_prof),
                                     slit_var, 1,
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

                    if aperture_signs[k] == 1:
                
                        slit_psf = slit_prof.clip(min=0.0)

                    else:

                        slit_psf = slit_prof.clip(max=0.0)

                    # Normalize the profile
                        
                    aperture_psf = aperture_signs[k]*\
                                   np.abs(slit_psf/np.nansum(slit_psf[z_psf]))
                    
                    # Determine the pixels to actually use
                    
                    z_aperture = (aperture_mask > float(k)) & \
                                 (aperture_mask <= float(k + 1)) & \
                                 (fit_img['goodbad'] == 1) & \
                                 (aperture_psf != 0.0)
                    
                    if np.sum(z_aperture) > 0:

                        # Scale the data and variances
                        
                        vals = slit_img[z_aperture]/aperture_psf[z_aperture]

                        vals2 = slit_img[z_aperture]

 #                       if j > 480 and k == 0:
 #
 #                           axes2 = pl.subplot(2,1,2,sharex=axes1)
 #                           axes2.plot(slit_arc[z_aperture],vals,'or')        
 #
 #                           #pl.draw()
 #                           #pl.pause(0.5)
 #                           #pl.clf()
                            

                        vars = slit_var[z_aperture]/aperture_psf[z_aperture]**2

                        # Compute the optimal estimate
                        
                        weights = 1/vars
                        wmean = np.sum(weights*vals)/np.sum(weights)
                        wvar = 1/np.sum(weights)

#                        if j > 480 and k == 0:
#
#                            print(j,wmean, wavecal[zslit, xmin + j][0]) 
#                         
#                            print('data', slit_img[z_aperture])
#                            print()
#                            print('var', slit_var[z_aperture])
#                            print()
#                            print('psf', aperture_psf[z_aperture])
#                            print()
#                            print('estimates', vals)
#                            print()
#                            print('estimates', vars)
#                            print()
#                            print()
#                            print()
#                            pl.show()

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
                                          aperture_signs[k]
                    
                    spectrum_var = np.sum(slit_var[z_aperture] * partial ** 2)

                    spectrum_unc[k, j] = np.sqrt(np.abs(spectrum_var))

                #
                # Deal with the flags
                #
                    
                # Check the bitmask for linearity

                z_linearity = slit_lmm[z_aperture] == 1

                if sum(z_linearity) > 0:
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

        if progressbar is True:
            loop_progress(i, 0, norders)

    return spectrum_list


def write_apertures_fits(spectra,
                         xranges,
                         aimage,
                         sky,
                         flat,
                         flat_fielded,
                         naps,
                         orders,
                         header_info,
                         aperture_positions,
                         aperture_radii,
                         plate_scale,
                         slith_pix,
                         slith_arc,
                         slitw_pix,
                         slitw_arc,
                         resolving_power,
                         xunits,
                         yunits,
                         latex_xunits,
                         latex_yunits,
                         latex_xlabel,
                         latex_ylabel,
                         latex_ulabel,                         
                         version,
                         output_fullpath,
                         wavecalinfo=None,
                         bg_annulus=None,
                         bg_regions=None,
                         bg_fitdegree=None,
                         optimal_info=None,
                         badpixel_info=None,
                         overwrite=True,
                         verbose=True):

    """
    To write a spextool spectral FITS file to disk

    Parameters
    ----------

    Returns
    -------
    None.  Writes FITS files to disk.


    """

    #
    # Check parameters
    #
    check_parameter('write_apertures_fits', 'spectra', spectra, 'list')

    check_parameter('write_apertures_fits', 'xranges', xranges, 'ndarray')

    check_parameter('write_apertures_fits', 'aimage', aimage, 'str')

    check_parameter('write_apertures_fits', 'flat', flat, 'str')

    check_parameter('write_apertures_fits', 'naps', naps, 'int')

    check_parameter('write_apertures_fits', 'orders', orders,
                    ['int', 'list', 'ndarray'])

    check_parameter('write_apertures_fits', 'header_info', header_info,
                    'dict')

    check_parameter('write_apertures_fits', 'aperture_positions',
                    aperture_positions, 'ndarray')

    check_parameter('write_apertures_fits', 'aperture_radii', aperture_radii,
                    ['int', 'float', 'ndarray'])

    check_parameter('write_apertures_fits', 'plate_scale', plate_scale, 'float')

    check_parameter('write_apertures_fits', 'slith_pix', slith_pix, 'float')

    check_parameter('write_apertures_fits', 'slith_arc', slith_arc, 'float')

    check_parameter('write_apertures_fits', 'slitw_pix', slitw_pix, 'float')

    check_parameter('write_apertures_fits', 'slitw_arc', slitw_arc, 'float')

    check_parameter('write_apertures_fits', 'resolving_power',
                    resolving_power, 'int')

    check_parameter('write_apertures_fits', 'xunits', xunits, 'str')

    check_parameter('write_apertures_fits', 'yunits', yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xunits', latex_xunits, 'str')

    check_parameter('write_apertures_fits', 'latex_yunits', latex_yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xlabel', latex_xlabel, 'str')

    check_parameter('write_apertures_fits', 'latex_ylabel', latex_ylabel, 'str')

    check_parameter('write_apertures_fits', 'version', version, 'str')

    check_parameter('write_apertures_fits', 'output_fullpath',
                    output_fullpath, 'str')

    check_parameter('write_apertures_fits', 'optimal_info', optimal_info,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'badpixel_info', badpixel_info,
                    ['NoneType', 'dict'])    

    check_parameter('write_apertures_fits', 'wavecalinfo', wavecalinfo,
                    ['dict', 'NoneType'])

    check_parameter('write_apertures_fits', 'overwrite', overwrite,
                    'bool')

    check_parameter('write_apertures_fits', 'verbose', verbose,
                    'bool')

    #    
    # Get set up
    #

    orders = np.asarray(orders)
    norders = len(orders)

    #
    # Get the list of spectra into an array format
    #

    # Determine the maximum size of each aperture array

    npixels = []
    for slice in spectra:
        npixels.append(np.shape(slice)[1])

    max_npixels = np.max(np.asarray(npixels))

    # Now create arrays into which the slices will be placed

    array = np.full((norders * naps, 4, max_npixels), np.nan)
    array[:,3,:] = 0.0

    # Now fill in the arrays

    l = 0
    for slice in spectra:
        array[l, :, 0:npixels[l]] = slice
        l += 1

    #
    # Now write the file(s) to disk
    #

    # Create the headers

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the original file header keywords.  Cut out any history first
    # if need be.

    is_history_present = 'HISTORY' in header_info

    if is_history_present is True:

        old_history = header_info['HISTORY']

        # remove it from the header_info dictionary

        header_info.pop('HISTORY')

    keys = list(header_info.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (header_info[keys[i]][0], header_info[keys[i]][1])

    # Add spextool keywords

    hdr['MODULE'] = ('extract', ' Creation module')
    hdr['VERSION'] = (version, ' pySpextool version')
    hdr['AIMAGE'] = (aimage, ' A image')
    hdr['SKYORDRK'] = (sky, ' Sky or dark image')
    hdr['FLAT'] = (flat, ' Flat field image')
    hdr['FLFIELD'] = (flat_fielded, ' Flat fielded?')    

    if wavecalinfo is not None:

        hdr['WAVECAL'] = (wavecalinfo['file'], ' Wavecal file')
        hdr['WCTYPE'] = (wavecalinfo['wavecaltype'],
                         ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (wavecalinfo['wavetype'], ' Wavelength type')
    
    hdr['NORDERS'] = (norders, ' Number of orders')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), ' Orders')
    hdr['NAPS'] = (naps, ' Number of apertures')

    hdr['PLTSCALE'] = (plate_scale, ' Plate scale (arcseconds pixel-1)')
    hdr['SLTH_ARC'] = (slith_arc, ' Nominal slit height (arcseconds)')
    hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
    hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit height (pixels)')
    hdr['SLTW_PIX'] = (slitw_pix, ' Slit width (pixels)')
    hdr['RP'] = (resolving_power, ' Average resolving power')

    # Add the aperture positions and radii

    for i in range(norders):

        name = 'APOSO' + str(orders[i]).zfill(3)
        comment = ' Aperture positions (arcseconds) for order ' + \
                  str(orders[i]).zfill(3)
        val = ','.join([str(round(elem, 4)) for elem in aperture_positions[i]])
        hdr[name] = (val, comment)

        name = 'ARADO' + str(orders[i]).zfill(3)
        comment = ' Aperture radii (arcseconds) for order ' + \
                  str(orders[i]).zfill(3)
        val = ','.join([str(round(elem, 3)) for elem in aperture_radii[i,:]])
        hdr[name] = (val, comment)

    # Deal with the background info

    if bg_annulus is not None:

        hdr['BGSUB'] = (True, ' Background subtraction?')        
        hdr['BGREGS'] = (None, ' Background regions (arcseconds)')
        hdr['BGSTART'] = (bg_annulus[0],
                          ' Background annulus start (arcseconds)')
        hdr['BGWIDTH'] = (bg_annulus[1],
                          ' Background annulus width (arcseconds)')        
        
        hdr['BGFITDEG'] = (bg_fitdegree, ' Background fit degree')        

    elif bg_regions is not None:

        hdr['BGSUB'] = (True, ' Background subtraction?')                
        hdr['BGREGS'] = (bg_regions,
                          ' Background regions (arcseconds)')        

        hdr['BGSTART'] = (None,
                          ' Background annulus start (arcseconds)')
        hdr['BGWIDTH'] = (None,
                          ' Background annulus width (arcseconds)')        
        
        hdr['BGFITDEG'] = (bg_fitdegree, ' Background fit degree')        
        
    else:

        hdr['BGSUB'] = (False, ' Background subtraction?')                
        hdr['BGREGS'] = (None,
                          ' Background regions (arcseconds)')        

        hdr['BGSTART'] = (None,
                          ' Background annulus start (arcseconds)')
        hdr['BGWIDTH'] = (None,
                          ' Background annulus width (arcseconds)')        
        
        hdr['BGFITDEG'] = (None, ' Background fit degree')        

    if optimal_info is not None:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')
        hdr['THRESH'] = (optimal_info['thresh'], ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (True, ' Optimal extraction?')
        hdr['SUMEXT'] = (False, 'Sum extraction?')
        hdr['PSFRAD'] = (optimal_info['psfradius'], ' PSF radius (arcseconds)')

                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        
        hdr['THRESH'] = (None, ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (False, ' Optimal extraction?')
        hdr['SUMEXT'] = (True, ' Sum extraction?') 
        hdr['PSFRAD'] = (None, ' PSF radius (arcseconds)')


    if badpixel_info is not None:

        hdr['BDPXFIX'] = (True, ' Bad pixels fixed?')        
        hdr['THRESH'] = (badpixel_info['thresh'], ' Bad pixel sigma threshold')
                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        

    # Add the S/N ratio 

    for i in range(norders):

        name = 'SNRO' + str(orders[i]).zfill(3)
        comment = ' Median S/N values for order ' + \
                  str(orders[i]).zfill(3)

        values = []
        for j in range(naps):

            idx = i*naps+j
        
            signal = array[idx,1,:]
            noise = array[idx,2,:]
            
            values.append(str(int(np.round(np.nanmedian(signal/noise)))))

        hdr[name] = (", ".join(values), comment)
                
    # Deal with units
        
    hdr['XUNITS'] = (xunits, ' Units of the x axis')
    hdr['YUNITS'] = (yunits, ' Units of the y axis')

    hdr['LXUNITS'] = (latex_xunits, ' LateX units of the x axis')
    hdr['LYUNITS'] = (latex_yunits, ' LateX units of the y axis')

    hdr['LXLABEL'] = (latex_xlabel, ' LateX x-axis label')
    hdr['LYLABEL'] = (latex_ylabel, ' LateX y-axis flux label')
    hdr['LULABEL'] = (latex_ulabel, ' LateX y-axis uncertainty label')    

    # Now add the HISTORY

    history = 'Spextool FITS files contain an array of size ' \
              '[nwaves,4,norders*naps]. The ith image (array[*,*,i]) ' \
              'contains the data for a single extraction aperture within ' \
              'an order such that, lambda=array[*,0,i], flux=array[*,1,i], ' \
              'uncertainty=array[*,2,i],flag=array[*,3,i].  The zeroth ' \
              'image (array[*,*,0]) contains the data for the aperture in ' \
              'the order closest to the bottom of the detector that is ' \
              'closest to the bottom of the slit (i.e. also closest to the ' \
              'bottom of the detector).  Moving up the detector, the FITS ' \
              'array is filled in with subsequent extraction apertures.  ' \
              'If no orders have been deselected in the extraction process, ' \
              'the contents of the ith aperture in order j can be found as ' \
              'follows: lambda=array[*,0,{j-min(orders)}*naps + (i-1)], ' \
              'flux=array[*,1,{j-min(orders)}*naps + (i-1)], ' \
              'uncertainty=array[*,2,{j-min(orders)}*naps + (i-1)], ' \
              'flag=array[*,3,{j-min(orders)}*naps + (i-1)].'

    history = split_text(history, length=65)

    if is_history_present is True:

        # Join the two histories

        history = old_history + [' '] + history
    
    for hist in history:
        hdr['HISTORY'] = hist

    # Set the file name for the spectra file and write

    hdr['FILENAME'] = (basename(output_fullpath), ' File name')

    fits.writeto(output_fullpath, array, hdr, overwrite=overwrite)


