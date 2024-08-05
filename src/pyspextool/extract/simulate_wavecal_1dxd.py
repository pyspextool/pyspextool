import numpy as np
import numpy.typing as npt

from pyspextool.io.check import check_parameter


def simulate_wavecal_1dxd(ncols:int, nrows:int, edgecoeffs:npt.ArrayLike,
                          xranges:npt.ArrayLike, slith_arc:float):

    """
    To simulate Spextool wavecal and spatcal arrays.

    Will generate wavecal and spatcal files in the 1DXD case with the
    wavelengths replaced with the column numbers.


    Input Parameters
    ----------------
    ncols : int
        The number of columns of the image.

    nrows : int
        The number of rows of the image.

    edgecoeffs : ndarray
        (norders,`edgedeg`+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  

    xranges : ndarray
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives 
        the end column number for said order.

    slith_arc : float
        The nominal slit height (arcseconds).

    Returns
    -------
    wavecal : numpy.ndarray
        Wavecal (nrows,ncols) array where each pixel is set to its
        wavelength which in this case is the column number.

    spatcal : numpy.ndarray
        Spatcal (nrows,ncols) array where each pixel is set to its 
        angular position on the sky (in arcseconds).

    indices : list
        An (norders,) list where each element is a dictionary with the
        following keys:

        `"x"` : numpy.ndarray
            An (ncols,) array of x values (in pixels).

        `"y"` : numpy.ndarray
            An(nrows,) array of y values (in arcseconds).

        `"xidx"` : numpy.ndarray
            An (nrows, ncols) array of x indices.

        `"yidx"` : numpy.ndarray
            An (nrows, ncols) array of y indices.
        
    """

    #
    # Check parameters
    #

    check_parameter('simulate_wavecal_1dxd', 'ncols', ncols, 'int')

    check_parameter('simulate_wavecal_1dxd', 'nrows', nrows, 'int')

    check_parameter('simulate_wavecal_1dxd', 'edgecoeffs', edgecoeffs,
                    'ndarray')    

    check_parameter('simulate_wavecal_1dxd', 'xranges', xranges, 'ndarray')

    check_parameter('simulate_wavecal_1dxd', 'slith_arc', slith_arc,
                    ['int', 'float'])    

    #
    # Get basic info and do basic things
    #
    
    ndimen = edgecoeffs.ndim

    if ndimen == 2:
        norders = 1

    # Add a dimension for consistency with multi-order data

        edgecoeffs = np.expand_dims(edgecoeffs,axis=0)
        xranges = np.expand_dims(xranges,axis=0)        

        
    if ndimen == 3:
        norders = edgecoeffs.shape[0]

    # Create empty NaN arrays for the wavecal and spatcal arrays and an empty
    # list of the rectification indices
    
    wavecal = np.full([nrows, ncols], np.nan)
    spatcal = np.full_like(wavecal, np.nan)
    indices = []
    
    #
    # start the loop over order
    #
    
    y = np.arange(nrows)

    for i in range(norders):

        start = xranges[i, 0]
        stop = xranges[i, 1]

        x_pix = np.arange(stop - start + 1) + start
        nx = len(x_pix)

        # Get the top and bottom positions of the slit

        botedge = np.polynomial.polynomial.polyval(x_pix, edgecoeffs[i, 0, :])
        topedge = np.polynomial.polynomial.polyval(x_pix, edgecoeffs[i, 1, :])

        difference = topedge-botedge

        #
        # Create the rectification indices
        #

        # Do the x indices first
        
        ny = np.floor(np.min(difference)).astype(int)
        xidx = np.tile(x_pix, (ny, 1))
        
        # Now do the y indices

        y_pix = np.arange(ny)        
        ny = len(y_pix)
        yidx = np.tile(np.reshape(y_pix,(ny, 1)), (1, nx))

        # Get the linear transformation
        
        slope = difference/(ny-1)
        scale = np.tile(slope, (ny, 1))
        zpt = np.tile(botedge, (ny, 1))    

        yidx = yidx*scale+zpt

        y_arc = y_pix/y_pix[-1]*slith_arc

        # Store the results

        indices.append({'x': x_pix, 'y': y_arc, 'xidx': xidx, 'yidx': yidx})
        
        #
        # Now create the wavecal and spatcal arrays
        #
        
        # Creat the pixel to arcsecond transformation

        pixtoarc = np.empty([2, stop - start + 1])
        pixtoarc[1, :] = slith_arc / (difference)
        pixtoarc[0, :] = -1 * pixtoarc[1, :] * botedge

        # Fill things in

        for j in range(stop - start + 1):

            wavecal[np.floor(botedge[j]).astype('int'):
                    np.ceil(topedge[j]).astype('int'), x_pix[j]] = x_pix[j]

            # Create ysub to make things readable...

            ysub = y[np.floor(botedge[j]).astype('int'):
                     np.ceil(topedge[j]).astype('int')]

            spatcal[np.floor(botedge[j]).astype('int'):
                    np.ceil(topedge[j]).astype('int'), x_pix[j]] = \
                np.polynomial.polynomial.polyval(ysub, pixtoarc[:, j])

    return wavecal, spatcal, indices
