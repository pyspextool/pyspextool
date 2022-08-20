import numpy as np
from scipy.interpolate import interp1d
from pyspextool.utils.arrays import trim_nan


def get_line_guess_position(spectra, orders, xranges, line_info):

    """
    To determine the (guess) column position of calibration lines

    Parameters
    ----------
    spectra :

    orders:

    xranges:

    line_info:





    Returns
    -------
    list
         A list of integers giving the individual file numbers


    Notes
    -----


    Examples
    --------

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program XS

    """
    # Get basic information

    norders = len(orders)
    nlines = len(line_info['strwave'])
    
    # Get output arrays set up

    x_cen = np.full(nlines, np.nan)
    x_left = np.full(nlines, np.nan)
    x_rght = np.full(nlines, np.nan)

    # Loop over each order

    for i in range(norders):

        # Find the lines for each order, if none, continue
        
        z  = line_info['order'] == orders[i]
        if np.sum(z) == 0:
            continue

        # Get set up for the interolation.  Get the stored wavecal
        # wavelengths, and their associated x values.

        w_stored = trim_nan(spectra[i,0,:],2,trim=True)
        x_stored = np.arange(xranges[i,0],xranges[i,1]+1)

        # Now let's get the wavelengths we want the positions of

        cen_w_want = line_info['strwave'][z].astype(float)        
        left_w_want = cen_w_want-line_info['leftwin'][z]
        rght_w_want = cen_w_want+line_info['rghtwin'][z]

        # Trim the left and right side to ensure they don't go below the
        # wavelength range of the order

        left_w_want = np.where(left_w_want < np.min(w_stored),
                             np.min(w_stored), left_w_want)

        rght_w_want = np.where(rght_w_want > np.max(w_stored),
                             np.max(w_stored), rght_w_want)

#        for_print(left_w_want,cen_w_want,rght_w_want)

        # Do the interplation 

        f = interp1d(w_stored,x_stored)

        left_x_want = f(left_w_want)
        cen_x_want = f(cen_w_want)
        rght_x_want = f(rght_w_want)                

        # Store the results

        x_cen[z] = cen_x_want
        x_left[z] = left_x_want
        x_rght[z] = rght_x_want

    # Add to line_info dictionary
 
    line_info['xleft'] = x_left
    line_info['xguess'] = x_cen
    line_info['xright'] = x_rght

    return line_info
        

        
