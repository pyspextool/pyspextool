import numpy as np
from scipy import interpolate
import numpy.typing as npt

from pyspextool.io.check import check_parameter
from pyspextool.fit.polyfit import polyfit_1d, poly_1d
from pyspextool.utils.math import mean_data_stack, round
from pyspextool.fit.fit_peak1d import *
from pyspextool.utils.arrays import find_index, trim_nan
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.for_print import for_print

def find_peaks(profiles:list,
               method_info:dict,
               fwhm:int | float=0.8,
               ybuffer:int=1):
    
    """
    To determine the locations of peaks in spatial profiles


    Parameters
    ----------
    profiles : list
        An (norders, ) list where each element is a dictionary.  
        Each dictionary has two keys:

        `angles` : ndarray
            An (nangles, ) array of sky angles.

        `profile` : ndarray
            An (nangles, ) array of the profile.

    method_info : dict

        `method` : {'auto', 'guess', 'fixed'}

        `peaks` : int, float, list, ndarray

            If `method` is 'auto' then `peaks` is the number of peaks to find.

            If `method` is 'guess' then `peaks` are guess positions, either
            an int or float for a single peak, or a list/ndarray of int or
            float for multiple peaks.

            If `method` is 'fixed' then `peaks` are positions, either
            an int or float for a single peak, or a list/ndarray of int or
            float for multiple peaks.
    
    fwhm : float, default 0.8 (arcseconds)
        The approximate FWHM of the peak to be identified.
        Only used if `method` is 'auto' or 'guess' to estimate the standard
        deviation of the gaussian function fit to the profile.

    ybuffer : int, default 1
        The number of pixels on the edge of the orders to ignore.  Useful
        as sometimes there is a downturn that can mess with the finding
        routine.
    

    Returns
    -------
    ndarray, ndarray
    
    apertures : An (norders, naps) ndarray of aperture positions.  
                By pySpextool convention, apertures[0,:] are the apertures for
                the order closest to the bottom of the image.

    apsigns : An (norders, naps) ndarray of aperture signs 
              (1 = positive, -1=negative).  By pySpextool convention, 
              apsigns[0,:] are the aperture signs the order closest to the 
              bottom of the image.
    
    """

    #
    # Check parameters
    #

    check_parameter('find_peaks', 'profiles', profiles, 'list')

    check_parameter('find_peaks', 'method_info', method_info, 'dict')

    check_parameter('find_peaks', 'fwhm', fwhm, ['int', 'float'])

    check_parameter('find_peaks', 'ybuffer', ybuffer, 'int')            
    

    # Get the number of orders

    norders = len(profiles)

    #
    # Run things by the type of method requested
    #
    
    if method_info['method'] == 'auto':

        npeaks = method_info['peaks']

        # Create arrays

        peak_positions = np.empty((norders, npeaks), dtype=float)
        peak_signs = np.empty((norders, npeaks), dtype=int)

        # loop over the orders

        for i in range(norders):

            angles = profiles[i]['angles']
            profile = profiles[i]['profile']

            # Trim NaNs

            good = trim_nan(profile)
            angles = angles[good]
            profile = profile[good]

            # Remove `ybuffer` pixels

            angles = angles[ybuffer+1:-(ybuffer+1)]
            profile = profile[ybuffer+1:-(ybuffer+1)]            

            # Subtract the background

            profile -= np.median(profile)

            # Take the absolute values of the profile

            abs_profile = np.abs(profile)

            for j in range(npeaks):
                
                # Find the index number of the largest peak in the
                # absolute value of the profile

                idx = int(np.argmax(abs_profile))

                # Do the fit.  
                
                p0 = [abs_profile[idx],angles[idx], fwhm / 2.354, 0]
                               
                fit = fit_peak1d(angles, abs_profile, p0=p0)

                peak_positions[i, j] = fit['parms'][1]

                # Get the aperture sign relative to the original profile.
                # Find the index of the fitted position first.

                idx = int(find_index(angles, fit['parms'][1]))

                peak_signs[i, j] = 1 if profile[idx] > 0 else -1

                # Subtract the fit off of the absolute profile and repeat.

                abs_profile -= fit['fit']                

    elif method_info['method'] == 'guess':

        # Get apertures

        peaks_guess = np.array(method_info['peaks'])
        npeaks = np.size(peaks_guess) // norders

        # Create arrays

        peak_positions = np.empty((norders, npeaks), dtype=float)
        peak_signs = np.empty((norders, npeaks), dtype=int)

        # loop over the orders

        for i in range(norders):

            angles = profiles[i]['angles']
            profile = profiles[i]['profile']

            # Trim NaNs

            good = trim_nan(profile)
            angles = angles[good]
            profile = profile[good]

            # Subtract the median background off the profile

            profile -= np.median(profile)

            # Take the absolute value

            abs_profile = np.abs(profile)

            # Now loop over each aperture and do the fit

            for j in range(npeaks):

                # Get the index number for the guess position

                idx = int(find_index(angles, float(peaks_guess[i, j])))

                # Do the fit to the absolute value of the profile

                p0 = [abs_profile[idx], peaks_guess[i, j], fwhm / 2.354, 0]

                fit = fit_peak1d(angles, abs_profile, p0=p0)

                peak_positions[i, j] = fit['parms'][1]

                # Get the aperture sign relative to the original profile.
                # Find the index of the fitted position first.

                idx = int(find_index(angles, fit['parms'][1]))
                
                peak_signs[i, j] = 1 if profile[idx] > 0 else -1
                
    elif method_info['method'] == 'fixed':

        # Get apertures

        peak_guesses = np.array(method_info['peaks'])
        npeaks = np.size(peak_guesses) // norders

        # Create arrays

        peak_positions = np.empty((norders, npeaks), dtype=float)
        peak_signs = np.empty((norders, npeaks), dtype=int)

        # loop over the orders

        for i in range(norders):

            angles = profiles[i]['angles']
            profile = profiles[i]['profile']

            # Trim NaNs

            good = trim_nan(profile)
            angles = angles[good]
            profile = profile[good]

            # Subtract the median background level from the profile

            profile -= np.median(profile)

            #

            for j in range(npeaks):

                # Store results

                peak_positions[i, j] = peak_guesses[i, j]

                # Get the aperture sign

                idx = find_index(angles, float(peak_guesses[i, j]))
                peak_signs[i, j] = 1 if profile[int(idx)] > 0 else -1

    #
    # Now sort everything properly and return the results
    #
    
    for i in range(norders):
        idx = np.argsort(peak_positions[i, :])
        peak_positions[i, :] = peak_positions[i, idx]
        peak_signs[i, :] = peak_signs[i, idx]

    return peak_positions, peak_signs



def combine_aperturesigns(aperture_signs:npt.ArrayLike):

    """
    To compute the average aperture signs given by-order aperture signs.

    Parameters
    ----------
    aperture_signs : ndarray
        An (norders, naps) array of integers +1, -1, giving the aperture 
        signs for each order.
    
    Returns
    -------
    average_aperture_signs : ndarray
        An (naps,) array giving the average aperture signs

    label : str
        A string giving the apeture signs as a comma-separated string of 
        + and - signs. 

    """

    #
    # Check parameters
    #

    check_parameter('combine_aperturesigns', 'aperture_signs', aperture_signs,
                    'ndarray')

    #
    # Compute the average aperture sign.  If the answer is zero, it gets set
    # to +.
    #

    average = np.sum(aperture_signs, axis=0) / \
              np.sum(np.abs(aperture_signs), axis=0)

    average_aperturesigns = np.full_like(average,1,dtype=int)

    for i in range(len(average)):
        if average[i] < 0:
            average_aperturesigns[i] = -1

    #
    # Create string version of the aperture signs
    #

    label = ', '.join(list(average_aperturesigns.astype(str)))
    label = label.replace('-1', '-')
    label = label.replace('1', '+')

    return average_aperturesigns, label

    
def make_1d_profile(rectified_order:dict,
                    atmospheric_transmission:dict=None,
                    robust_threshold:int=5):

    """
    To create a mean spatial profile for a rectified order

    Parameters
    ----------
    rectified_order : dict

        'image': ndarray
            A (nangles, nwavelength) rectified order image.

        'wavelenths': ndarray
            A (nwavelength, ) array of wavelengths.

        '`angles`': ndarray
            A (nangles,) array of angles on the sky.

    atmospheric_transmission : dict, default None
        If not None,
        
        'wavelength': ndarray
            The wavelength array

        'transmission': ndarray
            The transmission array
         
    robust_threshold : int or float, default 5
        The threshold for identifying outliers when the mean profile is 
        created.  

    Returns
    -------
    ndarray, ndarray

        angles : ndarray
            The rectified_order['angles'] array

        mean_profile : ndarray
            An (nangles,) array giving the mean profile.

    """

    #
    # Check parameters
    #

    check_parameter('make_1d_profile', 'rectified_order', rectified_order,
                    'dict')

    check_parameter('make_1d_profile', 'atmospheric_transmission',
                    atmospheric_transmission, ['NoneType', 'dict'])

    check_parameter('make_1d_profile', 'robust_threshold', robust_threshold,
                    ['int', 'float'])        
    
    #
    # Build the profile
    #
    
    # Unpack the data
    
    img = rectified_order['image']
    angles = rectified_order['angles']
    wavelengths = rectified_order['wavelengths']
    
    nrows = len(angles)    

    # Subtract the background

    medbg = np.median(img, axis=0)

    bgimg = np.tile(medbg, (nrows, 1))

    np.subtract(img, bgimg, out=img)
    
    # Build the atmospheric transmission weight mask if requested

    if atmospheric_transmission is not None:

        # Do the interpolation of the atmosphere

        f = interpolate.interp1d(atmospheric_transmission['wavelength'],
                                 atmospheric_transmission['transmission'],
                                 fill_value=1)
        rt = f(wavelengths)

        # Clip the low points, (np.where works oddly in this mode!)
        
        rt = np.where(rt > 0.1, rt, 0.1)

        # Creat a weight array

        weights = np.tile((1 / rt) ** 2, (nrows, 1))
        weights = np.rot90(weights,3)
        
    else:
        
        weights = None
        
    # Collapse the profile using a mean, weighted by `weights`.
    
    mean, mvar, mask = mean_data_stack(np.fliplr(np.rot90(img, 3)),
                                       weights=weights,
                                       robust=robust_threshold)
    
    # Normalize by the total absolute flux and return the results

    mean /= np.sum(np.abs(mean))    

    return angles, mean



def make_2d_profile(rectified_order,
                    trace_coefficients,
                    aperture_radii,
                    atmospheric_transmission=None,
                    use_mean_profile=False,
                    robust_threshold=5):


    """
    To create a 2D spatial profile of a rectified order

    Parameters
    ----------
    rectified_order : dict

        '`image`': ndarray
            The 2D rectified order image.

        '`wavelength`': ndarray
            The wavelength array

        '`angle`': ndarray
            The spatial angle array

    trace_coefficients : ndarray
        A (naps,) ndarray of trace coefficients.  Each item in the list is a 
        (ncoeffs,) ndarray which when evaluated using 
        rectified_order['wavelength'] gives the aperture position in units of 
        rectified_order['angle'].

    aperture_radii : ndarray
        A (naps,) array of aperture radii values in units of 
        rectified_order['angle'].

    atmospheric_transmission : dict or None

        '`wavelength`': ndarray
            An (nwavelength,) array of wavelenths.

        '`transmission`': ndarray
            An (nwavelengths,) array giving the tranmission of the atmosphere.

    use_mean_profile : {False, True} 
        Set to True to use the mean profile for all wavelengths.  Useful when
        the S/N of the data is low.

    robust_threshold : int or float, default=5
        The threshold for identifying outliers when the mean profile is 
        created.  
    
    Returns
    -------

    Returns
    -------
    ndarray, ndarray

        angles : ndarray
            The rectified_order['angles'] array

        model_profile : ndarray
            An (nangles,nwavelengths) array giving the model profile.

    """
    
    #
    # Check parameters
    #
    
    check_parameter('make_2d_profile', 'rectified_order', rectified_order,
                    'dict')

    check_parameter('make_2d_profile', 'trace_coefficients', trace_coefficients,
                    'ndarray')

    check_parameter('make_2d_profile', 'aperture_radii', aperture_radii,
                    'ndarray')        

    check_parameter('make_2d_profile', 'atmospheric_transmission',
                    atmospheric_transmission, ['NoneType', 'dict'])

    check_parameter('make_2d_profile', 'use_mean_profile', use_mean_profile,
                    'bool')    

    check_parameter('make_2d_profile', 'robust_threshold', robust_threshold,
                    ['int', 'float'])


    #
    # Unpack the data and get useful information
    #
    
    img = rectified_order['image']
    angles = rectified_order['angles']
    wavelengths = rectified_order['wavelengths']

    nangles, nwavelengths = np.shape(img)

    napertures = len(aperture_radii)
    
    #
    # Setup coefficient array
    #

    coeffs = np.zeros((nangles,3))

    #
    # Create the mean profile and map
    #

    a_profile, mean_profile = make_1d_profile(rectified_order,
                              atmospheric_transmission=atmospheric_transmission,
                              robust_threshold=robust_threshold)

    mean_model = np.column_stack((mean_profile,)*nwavelengths)

    if use_mean_profile is True:

        return angles, mean_model

    else:

        # Set any zero values to NaN

        z = np.where(mean_model == 0)
        mean_model[z] = np.nan

        # Divide the mean map into the data

        ratio = np.divide(img,mean_model)

        # Determine the pixels to use to determine the scale factors based
        # on the trace and aperture radius.

        trace_arc = np.empty((nwavelengths, napertures))
        for i in range(napertures):

            trace_arc[:,i] = poly_1d(wavelengths,trace_coefficients[i])

        trace_arc = np.median(trace_arc,axis=0)

        aperture_mask = make_aperture_mask(angles, trace_arc, aperture_radii)

        z = aperture_mask == 0

        aperture_mask[z] = np.nan
        aperture_mask[~np.isnan(aperture_mask)] = 1

        aperture_mask_map = np.column_stack((aperture_mask,)*nwavelengths)

        normalization = np.nanmedian(np.multiply(ratio,aperture_mask_map),
                                     axis=0)
        z = normalization == 0.0
        normalization[z] = np.nan

    
        normalization_map = np.divide(img,\
                                      np.row_stack((normalization,)*nangles))

        # Build the atmospheric transmission weight mask if requested
    
        if atmospheric_transmission is not None:

            # Do the interpolation of the atmosphere

            f = interpolate.interp1d(atmospheric_transmission['wavelength'],
                                    atmospheric_transmission['transmission'],
                                    fill_value=1)
            rt = f(wavelengths)

            # Clip the low points, (np.where works oddly in this mode!)
            
            rt = np.where(rt > 0.1, rt, 0.1)
            
            # Create weight array

            yunc = 1/rt
                
        else:
            
            yunc = None

        for i in range(nangles):

            result = polyfit_1d(wavelengths,normalization_map[i,:],2,
                                robust={'thresh':3.5,'eps':0.01},
                                yunc=yunc)

            coeffs[i,:] = result['coeffs']
            
        #
        # Now create the spatial map using the coeffs array
        #

        model = np.empty((nangles,nwavelengths))
        for i in range(nangles):

            model[i,:] = poly_1d(wavelengths,coeffs[i,:])

        #
        # Return the results
        #

        return angles, model


def make_aperture_mask(angles:npt.ArrayLike,
                       aperture_positions:npt.ArrayLike,
                       aperture_radii:npt.ArrayLike,
                       bg_annulus:list=None,
                       bg_regions:str=None,
                       report=False):

    """
    To create an aperture mask for use in extraction.

    Parameters
    ----------
    angles : ndarray
        An array of angular positions along the slit (typically in arcseconds).

    aperture_positions : ndarray 
        An (norders, naps) array of aperture positions (same units as
        `angles`).
    
    aperture_radii : ndarray 
        An(norders, naps) array of aperture radii (same units as `slit_arc`).

    bg_annulus : list, default None
        If given, background regions will be defined as a 1D annulus.  That
        is, `bg_annulus`[0] is the start radius and `bg_annulus`[1] is the
        width of the 1D annulus.

    bg_regions : str, default None
        If given, background regions will be defined based on values passed
        by the user.  For example, `bg_regions`='1-2,14-15' will have a
        background region from 1-2 (units of `angles`) and 14-15 (units of
        `angles`).
    
    Returns
    -------
    ndarray
         An array of the aperture mask.  The pixel values denote their status:

         val = 0 - nothing
         val = -1 - a background pixel
         0 < val <= 1 - a pixel belonging to aperture 1
         1 < val <= 2 - a pixel belonging to aperture 2
         etc.

    Notes
    -----
    The function generates an array where each pixel value denotes its status.
    To construct the mask, the background regions are identified without
    regard to potential overlap.  If `psbginfo` is passed, then overlap is
    dealt with by clearing out from `aperture_positions`-`bgstart` to
    `aperture_positions`+`bgstart`.  Finally, the apertures are labeled by
    their aperture number, e.g. 1, 2, etc. Pixels at the edge of an aperture
    are adjust to reflect their fractional pixel value, e.g. 1 -> 0.5 or
    2 -> 1.6.


    """

    #
    # Check parameters
    #

    check_parameter('make_aperture_mask', 'angles',  angles, 'ndarray')

    check_parameter('make_aperture_mask', 'aperture_positions',
                    aperture_positions, 'ndarray')    

    check_parameter('make_aperture_mask', 'aperture_radii',  aperture_radii,
                    'ndarray')    

    check_parameter('make_aperture_mask', 'bg_annulus', bg_annulus,
                    ['list','NoneType'])
    
    check_parameter('make_aperture_mask', 'bg_regions', bg_regions,
                    ['str','NoneType'])

    #
    # Let's get started
    #
    
    # Now do basic things

    npix = len(angles)
    naps = len(aperture_positions)

    # Now check to see if the apertures overlap

    mask = np.zeros([npix])

    for i in range(naps):

        ap = [aperture_positions[i] - aperture_radii[i],
              aperture_positions[i] + aperture_radii[i]]

        idx = find_index(angles, ap)
        mask[int(idx[0]):int(idx[1])+1] += 1

    if report is True:

        print(idx)
        oldmask = mask
        
    test = mask > 1
    if sum(test) >= 1:

        exception = 'The apertures overlap.'
        raise pySpextoolError(exception)

    # Create the mask

    mask = np.zeros([npix])

    #
    # Let's deal with the background regions
    #

    if bg_regions is not None:

        regions = bg_regions.split(',')

        for bgrange in regions:

            left, right = bgrange.split('-')

            idx = find_index(angles, [float(left), float(right)])
            mask[int(idx[0]):int(idx[1])+1] = -1

    if bg_annulus is not None:

        bgstart = bg_annulus[0]
        bgwidth = bg_annulus[1]

        # Set all pixels that fall in the background region around
        # each aperture to -1 regardless of whether there is overlap.

        for i in range(naps):

            if bgstart <= aperture_radii[i]:

                error = '`bg_annulus`[0] (bgstart) must be greater than ' + \
                        'the aperture radius.'
                raise pySpextoolError(error)

            leftbg = [aperture_positions[i] - bgstart - bgwidth,
                      aperture_positions[i] - bgstart]

            rghtbg = [aperture_positions[i] + bgstart,
                      aperture_positions[i] + bgstart + bgwidth]
            
            leftidx = find_index(angles, leftbg)
            rghtidx = find_index(angles, rghtbg)
            
            mask[int(leftidx[0]):int(leftidx[1])] = -1
            mask[int(rghtidx[0]):int(rghtidx[1])] = -1

        # Now clear out from -bgstart to +bgstart around each aperture
        # position to deal with any overlap

        for i in range(naps):
            clr = [aperture_positions[i] - bgstart,
                   aperture_positions[i] + bgstart]
            clridx = find_index(angles, clr)
            mask[int(clridx[0]):int(clridx[1])+1] = 0

    #
    # Now fill in the apertures
    #
    
    for i in range(naps):

        ap = [aperture_positions[i] - aperture_radii[i],
              aperture_positions[i] + aperture_radii[i]]

        idx = find_index(angles, ap)
        mask[int(idx[0]):int(idx[1])+1] += i + 1

        # Fix end points to reflect fractional pixels


        if idx[0] - np.floor(idx[0]) >= 0.5:

            mask[int(idx[0])] = 0
            mask[int(idx[0]) + 1] = 0.5 + round(idx[0]) - idx[0] + i

        else:

            mask[int(idx[0])] = 0.5 - (idx[0] - np.floor(idx[0])) + i

        if idx[1] - np.floor(idx[1]) >= 0.5:

            mask[int(idx[1]) + 1] = 0.5 - (round(idx[1]) - idx[1]) + i

        else:

            mask[int(idx[1])] = 0.5 + (idx[1] - round(idx[1])) + i

    if report is True:

        for_print(np.arange(len(angles)),angles,oldmask, mask)

            
    return mask

        
   
    

