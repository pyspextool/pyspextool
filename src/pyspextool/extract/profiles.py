import numpy as np
from scipy import interpolate


from pyspextool.io.check import check_parameter
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
        The threshold for identifying outliers.  See find_outliers.

    Returns
    -------
    ndarray, ndarray



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
    y = rectified_order['angle']
    w = rectified_order['wavelength']

    nrows = len(y)    

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
        rt = f(w)

        # Clip low points and create weight array

        rt = np.where(rt <= 0.1, rt, 0.1)

        weights = np.tile((1 / rt) ** 2, (nrows, 1))

    else:

        weights = None

    # Collapse the profile using a mean, weighted by `weights`.

    mean, mvar, mask = mean_data_stack(np.rot90(img, 3),
                                       weights=np.rot90(weights, 3),
                                       robust=robust_threshold)
    
    # Normalize by the total absolute flux and return the results

    mean /= np.sum(np.abs(mean))    

    return y, mean



def make_2d_profile(rectified_order, trace_coeffients, aperture_radii,
                    atmospheric_transmission=None, robust_threshold=5,
                    poly_degree=2):


    """
    To create a 2D spatial profile

    Parameters
    ----------
    rectified_order : dict

        '`image`': ndarray
            The 2D rectified order image.

        '`wavelength`': ndarray
            The wavelength array

        '`angle`': ndarray
            The spatial angle array

    trace_coefficients :

    aperture_radii :

    atmospheric_transmission : dict or None

        '`w`': ndarray
            The wavelength array

        '`t`': ndarray
            The transmission array

    robust_threshold : int or float, default=5
        The threshold for identifying outliers.  See find_outliers.
    
    poly_degree : int, default=2


    """
    
    #
    # Check parameters
    #

    check_parameter('make_2d_profile', 'rectified_order', rectified_order,
                    'dict')

    check_parameter('make_2d_profile', 'atmospheric_transmission',
                    atmospheric_transmission, ['NoneType', 'dict'])

    check_parameter('make_2d_profile', 'robust_threshold', robust_threshold,
                    ['int', 'float'])

    check_parameter('make_2d_profile', 'poly_degree', poly_degree, 'int')    
    
    #
    # Build the profile
    #    
    
    # Unpack the data
    
    img = rectified_order['image']
    a = rectified_order['angle']
    w = rectified_order['wavelength']

    na, nw = np.size(img)

#    # Creat the mean profile

    anglesmake_1d_profile(rectified_order, atmospheric_transmission=None,
                    robust_threshold=5)

    
