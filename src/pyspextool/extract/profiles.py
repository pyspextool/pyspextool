import numpy as np
from scipy import interpolate
#from pyspextool.plot.plot_image import plot_image
#import matplotlib.pyplot as pl
#from astropy.io import fits


from pyspextool.io.check import check_parameter
from pyspextool.extract.make_aperture_mask import make_aperture_mask
from pyspextool.fit.polyfit import poly_1d
from pyspextool.fit.polyfit import poly_fit_1d
from pyspextool.utils.math import mean_data_stack



def make_1d_profile(rectified_order, atmospheric_transmission=None,
                    robust_threshold=5):

    """
    To create a mean spatial profile for a rectified order

    Parameters
    ----------
    rectified_order : dict

        '`image`': ndarray
            A (nangles, nwavelength) rectified order image

        '`wavelenth`': ndarray
            A (nwavelength,) array of wavelengths

        '`angle`': ndarray
            A (nangles,) array of spatial angles

    atmospheric_transmission : dict or None

        '`wavelength`': ndarray
            The wavelength array

        '`transmission`': ndarray
            The transmission array

    robust_threshold : int or float, default=5
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
    angles = rectified_order['angle']
    wavelengths = rectified_order['wavelength']

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

        # Clip low points and create weight array

        rt = np.where(rt <= 0.1, rt, 0.1)

        weights = np.tile((1 / rt) ** 2, (nrows, 1))
        weights = np.rot90(weights, 3)

    else:

        weights = None

    # Collapse the profile using a mean, weighted by `weights`.

    mean, mvar, mask = mean_data_stack(np.fliplr(np.rot90(img, 3)),
                                       weights=weights, robust=robust_threshold)
    
    # Normalize by the total absolute flux and return the results

    mean /= np.sum(np.abs(mean))    

    return angles, mean



def make_2d_profile(rectified_order, trace_coefficients, aperture_radii,
                    atmospheric_transmission=None, use_mean_profile=False,
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
    angles = rectified_order['angle']
    wavelengths = rectified_order['wavelength']

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

            # Clip low points and create weight array

            rt = np.where(rt <= 0.1, rt, 0.1)

            yunc = 1/rt
                
        else:
            
            yunc = None

        for i in range(nangles):

            result = poly_fit_1d(wavelengths,normalization_map[i,:],2,
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
        

