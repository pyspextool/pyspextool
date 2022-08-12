import numpy as np


def correct_uspex_amps(img):

    """
    To correct for bias voltage drift in an uSpeX FITS image

    Parameters
    ----------
    img : numpy array
        An uSpeX image

    Returns
    --------
    numpy.ndarray
        The uSpeX image with the bias variations "corrected".

    Procedure
    ---------
    There are 32 amplifiers that talk to 64 columns each.  The median
    intensity of the 64 reference pixels at the bottom of image are
    subtracted from all rows in the 64 columns.

    Example
    --------
    result = uspexampcor(img) (how do we do this properly?)


    Modification History
    --------------------
    2022-05-25 - Written by M. Cushing, University of Toledo.
                Based on the Spextool mc_uspexampcor.pro IDL program.
    """

    for i in range(0, 32):

        xl = 0+64*i
        xr = xl+63

        med = np.median(img[2044:2047+1, xl:xr+1])
        img[:, xl:(xr+1)] -= med

    return img
