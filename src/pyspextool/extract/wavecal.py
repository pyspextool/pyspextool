import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as pl
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from astropy.io import fits


from pyspextool.fit.polyfit import polyfit_1d, polyfit_2d,poly_1d, poly_2d
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils.arrays import trim_nan
from pyspextool.plot.limits import get_spectra_range
from pyspextool.utils.loop_progress import loop_progress
from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.utils.arrays import idl_unrotate, idl_rotate
from pyspextool.pyspextoolerror import pySpextoolError


def find_lines_1dxd(spectra:dict,
                    orders:npt.ArrayLike,
                    line_info:dict,
                    pix_thresh:int,
                    qa_figuresize:tuple=(7,9),
                    qa_fontsize:tuple=12,                    
                    qa_fullpath:str=None,
                    verbose:bool=True):

    """
    To find the location of emission lines in cross-dispersed spectra.

    The "1dxd" label means the function can be used to identify lines 
    in spectra where the arrays columns are aligned with wavelength 
    (the 1d part) and can be crossed dispersed spectra (the xd part).


    Parameters
    ----------
    spectra : dict
        An (norder,) dictionary where each key has the form ORNNN_AP01
        where NNN is the order number.  The value of each entry is a
        (4, nwave) ndarray where:
        wavelength = (0,:)
        intensity = (1,:)
        uncertainty = (2,:)
        flags = (3,:)

    orders : ndarray
        An (norders,) array of int giving the order numbers.  By convention,
        orders[0] is the order closest to the bottom of the array.

    line_info : dict
        'order' : ndarray
            An (nlines, ) array of the order number of each line.

        'wavelength' : ndarray
            An (nlines, ) array of the wavelength of each line.

        'id' : ndarray of str
            An (nlines, ) array of the id of each line.

        'delta_wavelength_left' : ndarray 
            An (nlines, ) array of the delta lambda (in Angstroms) for the left
            edge of the fit window.

        'delta_wavelength_right' : ndarray 
            An (nlines, ) array of the delta lambda (in Angstroms) for the right
            edge of the fit window.
            
        'fit_type' : {'G', 'L', 'C'}
            Gaussian or Lorentzian or Centroid

        'num_parms' : int
            The number of parameters for the fit.  See fit_peak1d.

        'range_min_xguess' : ndarray
            An (nlines, ) array of the x value associated with 
            `'wavelength'` - `'delta_wavelength_left'`/1e4

        'xguess' : ndarray
            An (nlines, ) array of the x value associated with `'wavelength'`

        'range_max_xguess' : ndarray
            An (nlines, ) array of the x value associated with 
            `'wavelength'` + `'delta_wavelength_left'`/1e4

    pix_thresh: int
        The threshold (in pixels) beyond which an identification is deemed 
        bad.  That is, the fit is deemed bad if (abs(fit[1]-guess) > pix_thresh.


    Returns
    -------
        dict

        Adds four additional keys to the `line_info` dictionary.

        'x' : ndarray of float
            An (nlines,) array of the x position of each line.

        'fwhm' : ndarray of float
            An (nlines,) array of the fwhm of each line.

        'intensity' : ndarray float
            An (nlines,) array of the maximum value of each line.

        'goodbad' : ndarray of int 
            An (nlines,) goodbad array for the lines.

    """

    #
    # Check parameters
    #
    
    check_parameter('find_lines_1dxd', 'spectra', spectra, 'list')

    check_parameter('find_lines_1dxd', 'orders', orders, 'ndarray')

    check_parameter('find_lines_1dxd', 'line_info', line_info, 'dict')

    check_parameter('find_lines_1dxd', 'pix_thresh', pix_thresh,
                    ['float','int'])

    # Get basic information
    
    nlines = len(line_info['wavelength'])

    # Setup the output arrays

    line_xpos = np.full(nlines, np.nan)
    line_fwhm = np.full(nlines, np.nan)
    line_inten = np.full(nlines, np.nan)
    line_goodbad = np.zeros(nlines, dtype=np.uint8)

    # Setup the qaplot if asked

    if isinstance(qa_fullpath,str):

        pl.ioff()
        pl.rcParams['font.size'] = qa_fontsize
#        pl.rcParams['lines.linewidth'] = 0.75

        pdf = PdfPages(qa_fullpath)
    #                            
    # Loop over each line
    #
    
    for i in range(nlines):

        # Find the order associated with the line
        
        z = np.where(orders == line_info['order'][i])
        if np.size(z) == 0:
            continue

        # get spectra for this order and clip
    
        z = np.sum(z)
        x = trim_nan(spectra[z][0,:],2,trim=True)
        y = trim_nan(spectra[z][1,:],2,trim=True)        
        
        # Cut the line out

        zline = (x >=line_info['range_min_xguess'][i]) & \
                (x <=line_info['range_max_xguess'][i])
        if np.sum(zline) < line_info['num_parms'][i]:
            line_goodbad[i] = 0
            continue

        # Get fit type

        if line_info['fit_type'][i] == 'L':
            type = 'lorentzian'

        elif line_info['fit_type'][i] == 'G':
            type = 'gaussian'

        elif line_info['fit_type'][i] == 'C':
            type = 'centroid'

        else:
            print('Unknown `fittype`.')

        if type != 'centroid':

            offset = np.min(y[zline]) if line_info['num_parms'][i] == 3 else 0

            try:

                fit = fit_peak1d(x[zline],y[zline],
                                 nparms=line_info['num_parms'][i].item(),
                                 type=type,positive=True)

            except RuntimeError:

                fit = {'parms':[np.nan, np.nan, np.nan],
                       'fit':np.full(len(x[zline]), np.nan)}


            # Store the results
    
            line_xpos[i] = fit['parms'][1]
            line_fwhm[i] = fit['parms'][2]*2.354
            line_inten[i] = fit['parms'][0]
                    
        else:

            med = np.median((y[zline][0],y[zline][-1]))

            m1 = np.sum(x[zline]*(y[zline]-med)) / np.sum(y[zline]-med)
            
            m2 = np.sum((x[zline]-m1)**2*(y[zline]-med)) / np.sum(y[zline]-med)

            f = interp1d(x[zline],y[zline])
            inten = f(m1)

            line_xpos[i] = m1            
            line_fwhm[i] = np.sqrt(m2)*2.354
            line_inten[i] = inten

        # Now let's check to see whether it is a good find or not.


        if (line_xpos[i] <= line_info['xguess'][i]+pix_thresh) and \
           (line_xpos[i] >= line_info['xguess'][i]-pix_thresh) and \
           (line_fwhm[i] > 0) and (line_inten[i] > 0):
           line_goodbad[i] = 1

        if isinstance(qa_fullpath,str):
                                    
            fig, (axes1, axes2) = pl.subplots(2, figsize=qa_figuresize,
                                              constrained_layout=False)

            # Plot the entire spectrum with vertical line for the line
            
            yrange = get_spectra_range(y,frac=0.1)

            title = 'Order '+str(line_info['order'][i])+r': $\lambda$='+\
                line_info['wavelength'][i]+r' $\mu$m'
            
            axes1.step(x, y, 'black', where='mid')
            axes1.set_title(title)            
            axes1.set_ylim(ymin=0, ymax=yrange[1])
            axes1.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes1.axvline(x=line_info['xguess'][i], linestyle='solid',color='r')

            # Plot the line itself with fit.
            
            yrange = get_spectra_range(y[zline],frac=0.1)

            goodbad = 'Good Fit' if line_goodbad[i] == 1 else 'Bad Fit'
            
            axes2.step(x[zline], y[zline], 'black', where='mid')
            axes2.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes2.axvline(x=line_info['xguess'][i], linestyle='solid',
                          color='r')
            axes2.axvline(x=line_xpos[i], linestyle='solid', color='g')
            
            if type != 'centroid':
                axes2.step(x[zline], fit['fit'], 'g' , where='mid')

            axes2.set_title(goodbad)
            
            axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])

            pdf.savefig(fig)
            pl.close(fig)

            if verbose is not None:
                loop_progress(i, 0, nlines)

    if isinstance(qa_fullpath,str):
        pdf.close()
           
    # Add the results to the line_info dictionary

    line_info['x'] = line_xpos
    line_info['fwhm_pix'] = line_fwhm
    line_info['intensity'] = line_inten
    line_info['goodbad'] = line_goodbad

    return line_info



def get_line_guess_position(spectra:npt.ArrayLike,
                            orders:npt.ArrayLike | int,
                            xranges:npt.ArrayLike,
                            line_info:dict):

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

    """
    # Get basic information

    norders = len(orders)
    nlines = len(line_info['wavelength'])
    
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

        cen_w_want = line_info['wavelength'][z].astype(float)        
        left_w_want = cen_w_want-line_info['delta_wavelength_left'][z]
        rght_w_want = cen_w_want+line_info['delta_wavelength_right'][z]

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
 
    line_info['range_min_xguess'] = x_left
    line_info['xguess'] = x_cen
    line_info['range_max_xguess'] = x_rght

    return line_info





def get_spectral_pixelshift(xanchor:npt.ArrayLike,
                            yanchor:npt.ArrayLike,
                            xsource:npt.ArrayLike,
                            ysource:npt.ArrayLike,
                            smooth_source:bool=True,
                            qa_plotnumber:int=None,
                            qa_figuresize:tuple=(7,9),
                            qa_fontsize:int=12,                            
                            qa_show:bool=False,
                            qa_showscale:float | int=None,
                            qa_showblock:bool=False,
                            qa_fullpath:bool=None,
                            anchor_label:str='Anchor', 
                            source_label:str='Source'):

    """
    To determine the shift (in pixels) between two spectra.

    Parameters
    ----------
    xanchor : ndarray of int
        (nanchor,) array of pixel values.

    yanchor : ndarray
        (nanchor,) array of intensity values.

    xsource : ndarray of int
        (nsource,) array of pixel values.

    ysource : ndarray
        (nsource,) array of intensity values.

    smooth_source : {False, True}
        Set to True smooth the `yanchor` and `ysource` with a savitzky-golay
        filter before cross correlation to protect against bad pixels.
        Set to False to use the `yanchor` and `ysource` as given.

    output_fullpath : str, optional

    qa_show : {None, True, False}, optional
        Set to True/False to override config.setup['qa_show'].  If set to True,
        quality assurance plots will be interactively generated.

    Returns
    -------
    int
        The shift of `ysource` relative to `yanchor` in pixels.
        if `ysource` is shifted to the right of `yanchor` the returned
        value is positive.

    Notes
    -----
    Uses scipy.signal.correlate for the cross correlation.
    The savitzky-golay parameters are window_length=5 and polydeg=2.

    """

    #
    # Check parameters
    #

    check_parameter('get_spectral_pixelshift', 'xanchor', xanchor, 'ndarray')

    check_parameter('get_spectral_pixelshift', 'yanchor', yanchor, 'ndarray')

    check_parameter('get_spectral_pixelshift', 'xsource', xsource, 'ndarray')

    check_parameter('get_spectral_pixelshift', 'ysource', ysource, 'ndarray')

    check_parameter('get_spectral_pixelshift', 'smooth_source', smooth_source,
                    'bool')

    qa = check_qakeywords(show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock)

           
    ndat = len(xanchor)

    #
    # Find the intersection of the two x arrays.  That is, find the values
    # that are in common between the two arrays.

    junk, zanchor, zsource = np.intersect1d(xanchor, xsource,
                                            return_indices=True)

    # clip the results

    xanchor = xanchor[zanchor]
    yanchor = yanchor[zanchor]

    xsource = xsource[zsource]
    ysource = ysource[zsource]


    #
    # Savitzky-Golay the results to protect against bad pixels
    #

    if smooth_source is not False:

        sgyanchor = robust_savgol(xanchor, yanchor, 5)['fit']
        sgysource = robust_savgol(xsource, ysource, 5)['fit']

    else:

        sgyanchor = yanchor
        sgysource = ysource

    #
    # Check for NaNs
    #

    z = np.isnan(sgyanchor)
    if np.sum(z) != 0:

        sgyanchor[z] = 0

    z = np.isnan(sgysource)
    if np.sum(z) != 0:

        sgysource[z] = 0
               
    #     
    # Run the cross correlation
    #
    
    xcor = correlate(sgysource, sgyanchor, mode='same', method='direct')
    
    xcor = xcor / np.nanmax(xcor)

    lag = correlation_lags(ndat, ndat, mode='same')

    # Assume the peak occurs within 100 pixels of zero and clip the results.

    z = np.logical_and(lag > -100, lag < 100)
    xcor = xcor[z]
    lag = lag[z]
    
    #
    # Now lets find the fit window.  We don't want to fit the entire thing, just
    # around the peak of the cross correlation.  So work your way awake from
    # the peak until the derivative changes sign.
    #

    maxidx = np.argmax(xcor)    

    dif = -1
    npix = 0
    while dif < 0 and maxidx + npix + 1 < ndat - 1:
        dif = xcor[maxidx + npix + 1] - xcor[maxidx + npix]
        npix += 1
        
    halfwin = np.around(npix).astype('int')

    # Clip out the fit zone

    fitlag = lag[maxidx - halfwin:maxidx + halfwin]
    fitxcor = xcor[maxidx - halfwin:maxidx + halfwin]

    # Do the fit
    
    r = fit_peak1d(fitlag, fitxcor, nparms=4, positive=True)
    offset = r['parms'][1]
    
    #
    # Do the QA plotting
    #
    
    if qa['show'] is True:
        
        plot_spectral_pixelshift(qa_plotnumber,
                                 (qa_figuresize[0]*qa['showscale'],
                                  qa_figuresize[1]*qa['showscale']),
                                 qa_fontsize*qa['showscale'],
                                 xanchor,
                                 yanchor,
                                 xsource,
                                 ysource,
                                 lag,
                                 xcor,
                                 fitlag,
                                 fitxcor,
                                 r['fit'],
                                 offset,
                                 anchor_label=anchor_label,
                                 source_label=source_label)
        
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False: 
            
            pl.pause(1)
    
    if isinstance(qa_fullpath,str):

        plot_spectral_pixelshift(None,
                                 qa_figuresize,
                                 qa_fontsize,
                                 xanchor,
                                 yanchor,
                                 xsource,
                                 ysource,
                                 lag,
                                 xcor,
                                 fitlag,
                                 fitxcor,
                                 r['fit'],
                                 offset,
                                 anchor_label=anchor_label,
                                 source_label=source_label)


        pl.savefig(qa_fullpath)
        pl.close()

    return offset


def make_interp_indices_1d(edgecoeffs:npt.ArrayLike,
                           xranges:npt.ArrayLike,
                           slith_arc:int | float,
                           array_output:bool=False):

    """
    To generate 1D indices for resampling of a spectral order.

    The program generates indices that can be used with 
    scipy.interpolate.RegularGridInterpolator to "straighten" a spectral order.
    The "1D" refers to the fact that it is assumed that each column of the 
    detector represents a single wavelength.  


    Parameters
    ----------
    edgecoeffs : ndarray
        A (`edgedeg`+1, 2) float array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,:]
        gives the coefficients for the bottom of the order and 
        edgecoeffs[1,:] gives the coefficients for the top of the order.  

    xranges : ndarray
        An (2, ) array giving the column numbers over which to 
        operate.  xranges[0] gives the starting column number for the 
        order and xranges[1] gives the end column number for the order.

    slith_arc : float
        The nominal slit height (arcseconds).

    array_output : {False, True}, optional
        Set to pack the results into a single array (see Returns).

    Returns
    -------
    array_output=False
        wavemap : ndarray
        The (ny, nx) array of columns associated with the resampled order
        spatmap : ndarray
        The (ny, nx) spatial positions associated with the resampled order 
        in arcseconds.
        xidx : ndarray
        The (ny, nx) array of x coordinates.
        yidx : ndarray
        The (ny, nx) array of y coordinates.

    array_output=True
    
    ndarray 
        A (2, nx+1, ny+1) float array.          
        [0,0,1:] = x_pix
        [0,1:,0] = y_arc
        [0,:,:] = the array of x coordinates.
        [1,0,1:] = x_pix
        [1,1:,0] = y_arc
        [1,:,:] = the array of y coordinates.    

    Examples
    --------
    later

    """

    # Check the parameters

    check_parameter('make_interp_indices_1d', 'edgecoeffs', edgecoeffs,
                    'ndarray', 2)

    check_parameter('make_interp_indices_1d', 'xranges', xranges, 'ndarray', 1)

    check_parameter('make_interp_indices_1d', 'slith_arc', slith_arc,
                    ['int', 'float'])
        
    # Generate the the positions at the top and bottom of the slit
    
    x_pix = np.arange(xranges[0], xranges[1]+1)
    top = poly_1d(x_pix, edgecoeffs[1,:])
    bot = poly_1d(x_pix, edgecoeffs[0,:])    

    # Find the minimum number of pixels spanned by the slit
    
    dif = top-bot
    ny = np.floor(np.min(dif)).astype(int)

    #
    # Generate the rectification grid
    # 

    # Build the xidx array
    
    xidx = np.tile(x_pix, (ny, 1))

    # Now do the yidx array.

    y_pix = np.arange(ny)

    slope = dif/(ny-1)
    nx = len(x_pix)

    yimg = np.tile(np.reshape(y_pix,(ny, 1)), (1, nx))
    
    scale = np.tile(slope, (ny, 1))
    zpt = np.tile(bot, (ny, 1))    

    yidx = yimg*scale+zpt

    y_arc = y_pix/y_pix[-1]*slith_arc

    #
    # Now return the results
    #
    
    if array_output is False:

        wavemap = np.tile(x_pix,(ny,1))
        spatmap = np.rot90(np.tile(y_arc, (nx, 1)), k=3)
        
        return xidx, yidx, wavemap, spatmap

    else:
        
        indices = np.full((2, ny+1, nx+1),np.nan)
        
        indices[0,1:, 1:] = xidx
        indices[0,0,1:] = x_pix
        indices[0,1:,0] = y_arc
        
        indices[1,1:, 1:] = yidx
        indices[1,0,1:] = x_pix
        indices[1,1:,0] = y_arc
        
        return indices



def mix_orders(image1:npt.ArrayLike,
               image1_mask:npt.ArrayLike,
               image2:npt.ArrayLike,
               image2_mask:npt.ArrayLike,
               order_mask:npt.ArrayLike,
               orders:npt.ArrayLike | list,
               orders_from_image1:npt.ArrayLike | list,
               orders_from_image2:npt.ArrayLike | list,
               interorder_value=0.0):

    """
    To create an image with orders from two different images.

    Parameters
    ----------
    image1 : ndarray
        An (nrows, ncols) array.

    image1_mask : ndarray
        An (nrows, ncols) bit mask.

    image2 : ndarray
        An (nrows, ncols) array.

    image2_mask : ndarray
        An (nrows, ncols) bit mask.

    order_mask : ndarray
        An (nrows, ncols) array where each pixel is set to its order number.

    orders : ndarray
        An (norders,) array of the orders numbers in order_mask.

    orders_from_image1: ndarray
        A array of orders numbers to grab from image 1.

    orders_from_image2: ndarray
        A array of orders numbers to grab from image 2.

    interorder_value : float, optional
        The value to use for inter-order pixels.  The default is np.nan.
             
    Returns
    -------
    ndarray

       
    """

    #
    # Check parameters
    #

    check_parameter('mix_orders', 'image1', image1, 'ndarray')

    check_parameter('mix_orders', 'image1_mask', image1_mask, 'ndarray')

    check_parameter('mix_orders', 'image2', image2, 'ndarray')

    check_parameter('mix_orders', 'image2_mask', image2_mask, 'ndarray')

    check_parameter('mix_orders', 'order_mask', order_mask, 'ndarray')

    check_parameter('mix_orders', 'orders', orders, ['ndarray','list'])

    check_parameter('mix_orders', 'orders_from_image1', orders_from_image1,
                    ['ndarray','list'])                

    check_parameter('mix_orders', 'orders_from_image2', orders_from_image2,
                    ['ndarray','list'])

    check_parameter('mix_orders', 'interorder_value', interorder_value,
                    ['int','float'])

    #
    # Do checking to make sure the requested orders exist.
    #

    # Check the orders from image 1

    if np.sum(np.in1d(orders_from_image1, orders)) != len(orders_from_image1):

        message = 'Order number in parameter `orders_from_image1` is not in ' \
            '`orders`.'
        raise pySpextoolError(message)
    
    # Check the orders from image 2

    if np.sum(np.in1d(orders_from_image2, orders)) != len(orders_from_image2):

        message = 'Order number in parameter `orders_from_image2` is not in ' \
            '`orders`.'
        raise pySpextoolError(message)

    # Check whether threre is overlap in orders_from_image1 and
    # orders_from_image2

    if np.sum(np.in1d(orders_from_image2, orders_from_image1)) != 0:

        message = 'At least one order number is in both parameters ' \
            '`orders_from_image1` and `orders_from_image2`.'
        raise pySpextoolError(message)
       
    #
    # Loop over each request
    #

    output_image = np.full_like(image1, interorder_value, dtype=float)
    output_mask = np.full_like(image1_mask, 0)

    for order in orders_from_image1:

        z = np.where(order_mask == order)
        output_image[z] = image1[z]
        output_mask[z] = image1_mask[z]

    for order in orders_from_image2:
       
        z = np.where(order_mask == order)
        output_image[z] = image2[z]
        output_mask[z] = image2_mask[z]

    #
    # Return the new array
    #

    return output_image, output_mask



def plot_1d_residuals(figsize,
                      residuals:npt.ArrayLike,
                      residual_columns:npt.ArrayLike,
                      goodbad:npt.ArrayLike,
                      rms:float | int,
                      dispersion_degree:int):

    """
    To plot residuals for a 1D fit in a device independent way

    residuals : ndarray
        An (npoints, ) array of residuals

    """

    
                
    fig, axes1 = pl.subplots(1, figsize=figsize, constrained_layout=False)

    # Get the plot range 
            
    yrange = [-10 * rms * 1e4, 10 * rms * 1e4]
            
    # Get the number of good and bad points for the label

    zgood = goodbad == 1
    ntot = np.sum(zgood)
    
    zbad = goodbad == 0
    nbad = np.sum(zbad)

    # residuals as a function of column number 
                
    axes1.axhline(y=0, linestyle='dotted', color='black')

    axes1.scatter(residual_columns[zgood], residuals[zgood],
                  color='red', edgecolors='black',
                  s=8 ** 2, alpha=0.8)
    
    if nbad !=0:
        axes1.plot(residual_columns[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes1.set(xlabel='Column (pixels)', ylabel=r'Residual ($\mathrm{\AA}$)')
                
        axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])

            
    # Label the fit
                
    label = r'$\lambda$ degree=' + str(dispersion_degree) + \
        '\n RMS=' + "{:#.3g}".format(rms * 1e4) + \
        r' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
        r', n$_\mathrm{bad}$=' + str(nbad)
                  
    axes1.text(0.02, 0.95, label, color='black', ha='left',
               va='top', multialignment='left', transform=axes1.transAxes)
                
    label = 'Not all bad data points may be shown.'
    axes1.text(0.98, 0.95, label, color='black', ha='right',
               va='top', transform=axes1.transAxes)

    
def plot_1dxd_residuals(plot_number:int,
                        figure_size:tuple,
                        font_size:float | int,
                        residuals:npt.ArrayLike,
                        residual_orders:npt.ArrayLike,
                        residual_columns:npt.ArrayLike,
                        orders:npt.ArrayLike,
                        goodbad:npt.ArrayLike,
                        rms:float,
                        dispersion_degree:float | int,
                        order_degree:float | int):

    #
    # Set the fonts
    #
    
    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)

    
    #
    # Create the figure
    #

    fig = pl.figure(num=plot_number, figsize=figure_size)

    pl.subplots_adjust(left=0.1,
                       bottom=0.1, 
                       right=0.95, 
                       top=0.9, 
                       hspace=0.2)

    axes1 = fig.add_subplot(211)
    axes2 = fig.add_subplot(212)        

    # Get the colormap set up
            
    norders = len(orders)
    blue = np.linspace(0, 1, num=norders)
    red = 1 - blue
            
    # Get the plot range 
            
    yrange = [-10 * rms * 1e4, 10 * rms * 1e4]
            
    # Get the number of good and bad points for the label

    zgood = goodbad == 1
    ntot = np.sum(zgood)
            
    zbad = goodbad == 0
    nbad = np.sum(zbad)

    # Loop over each order so you can color them properly
            
    for i in range(norders):
        z = residual_orders == orders[i]
                
        # Residual as a function of order number
                
        axes1.axhline(y=0, linestyle='dotted', color='black')

        axes1.plot(residual_orders[z], residuals[z], 'o', 
                   markerfacecolor=(red[i], 0, blue[i]),
                   markeredgecolor='black',
                   markersize=8, alpha=0.8)
        
        axes1.plot(residual_orders[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes1.set(xlabel='Order Number', ylabel=r'Data-Fit ($\mathrm{\AA}$)')
                
        axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
        # residuals as a function of column number 
                
        axes2.axhline(y=0, linestyle='dotted', color='black')

        axes2.plot(residual_columns[z], residuals[z], 'o', 
                   markerfacecolor=(red[i], 0, blue[i]),
                   markeredgecolor='black',
                   markersize=8, alpha=0.8)
        
        axes2.plot(residual_columns[zbad], residuals[zbad], 's',
                   markersize=13, markerfacecolor='none', color='black')
                
        axes2.set(xlabel='Column (pixels)',
                  ylabel=r'Data-Fit ($\mathrm{\AA}$)')
                
        axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
                
        # Label the fit
                
        label = r'$\lambda$ degree=' + str(dispersion_degree) + \
                '\n order degree=' + str(order_degree) + \
                '\n RMS=' + "{:#.3g}".format(rms * 1e4) + \
                r' $\mathrm{\AA}$\n n$_\mathrm{tot}$=' + str(ntot) + \
                r', n$_\mathrm{bad}$=' + str(nbad)
                  
        axes1.text(0.02, 0.95, label, color='black', ha='left',
                   va='top', multialignment='left', transform=axes1.transAxes)
                
        label = 'Not all bad data points may be shown.'
        axes1.text(0.98, 0.95, label, color='black', ha='right',
                   va='top', transform=axes1.transAxes)    

    # Deal with the tickmarks
        

    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())

    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())

    
def plot_spectral_pixelshift(plot_number:int,
                             figure_size:tuple,
                             font_size:tuple,
                             xanchor:npt.ArrayLike,
                             yanchor:npt.ArrayLike,
                             xsource:npt.ArrayLike,
                             ysource:npt.ArrayLike,
                             lag:npt.ArrayLike,
                             xcorrelation:npt.ArrayLike,
                             fit_lag:npt.ArrayLike,
                             fit_xcorrelation:npt.ArrayLike,
                             fit:npt.ArrayLike,
                             offset:int | float,
                             anchor_label:str='Anchor', 
                             source_label:str='Source'):

    """
    Plot the cross correlation results in a device independent way

    Parameters
    ----------
    plot_number : int
        The plot number to pass to matplotlib.  Use None if unknown.

    figure_size : tuple
        A (2, ) tuple giving the plot size in inches.

    font_size : int or float
        The font size.

    xanchor : ndarray
        A (npixel1, ) array of pixel values for the anchor spectrum

    yanchor : ndarray
        A (npixel1, ) array of intensity values for the anchor spectrum

    xsource : ndarray
        A (npixel2, ) array of pixel values for the source spectrum

    ysource : ndarray
        A (npixel2, ) array of intensity values for the source spectrum

    lag : ndarray
        A (npixel3, ) array of lag values for the cross correlation

    xcorrelation : ndarray
        A (npixel3, ) array of cross correlation values.

    fit_lag : ndarray
        A (npixel4, ) array of lag values that were fit

    fit_xcorrelation : ndarray
        A (npixel4, ) array of cross correlation values that were fit

    fit : ndarray
        A (npixel4, ) of fitted cross correlation values

    offset : float or int
        The mean value of the fit


    Returns
    -------
    None
    
    """

    #
    # Normalize each spectrum 
    #
    
    np.divide(ysource, np.median(ysource), out=ysource)
    np.divide(yanchor, np.median(yanchor), out=yanchor)

    #
    # Make the figure
    #

    # Set the fonts

    # removed helvetica - problem for windows OS
    font = {
    #'family' : 'helvetica',
            'weight' : 'normal',
            'size'   : font_size}

    rc('font', **font)

    # Start the figure
    
    fig = pl.figure(num=plot_number, figsize=figure_size)
    pl.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    hspace=0.5)
    
    # Create the spectral plot
    
    axes1 = fig.add_subplot(311)
    
    yrange = get_spectra_range([yanchor, ysource], frac=0.1)

    axes1.margins(x=0)
    axes1.step(xanchor, yanchor, 'black', where='mid')
    axes1.step(xsource, ysource, 'r', where='mid')
    axes1.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes1.set(xlabel='Column Number', ylabel='Relative Intensity')

    axes1.text(0.95, 0.8, anchor_label, color='black', ha='right',
               transform=axes1.transAxes)

    axes1.text(0.95, 0.7, source_label, color='r', ha='right',
               transform=axes1.transAxes)

    # Deal with the tick marks
    
    axes1.xaxis.set_minor_locator(AutoMinorLocator())    
    axes1.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes1.tick_params(which='minor', length=3)
    axes1.tick_params(which='major', length=5)
    axes1.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Plot the entire cross correlation results

    axes2 = fig.add_subplot(312)    

    yrange = get_spectra_range(xcorrelation, frac=0.1)

    axes2.margins(x=0)
    axes2.tick_params(axis='x')
    axes2.tick_params(axis='y')
    axes2.set_title('Cross Correlation')
    axes2.step(lag, xcorrelation,color='black', where='mid')
    axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])
    axes2.set(xlabel='Lag (pixels)', ylabel='Relative Intensity')

    # Deal with the tick marks
    
    axes2.xaxis.set_minor_locator(AutoMinorLocator())    
    axes2.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes2.tick_params(which='minor', length=3)
    axes2.tick_params(which='major', length=5)
    axes2.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Plot a zoom in of the cross correlation and the fit

    axes3 = fig.add_subplot(313)    

    yrange = get_spectra_range([fit_xcorrelation, fit], frac=0.1)

    axes3.margins(x=0)
    axes3.step(fit_lag, fit_xcorrelation, color='black', where='mid')
    axes3.step(fit_lag, fit, 'r', where='mid')
    axes3.set_title('Fit of Cross Correlation')
    axes3.set_ylim(ymin=yrange[0], ymax=yrange[1])

    axes3.set(xlabel='Offset (pixels)', ylabel='Relative Intensity')
    axes3.axvline(x=offset, linestyle='dotted', color='r')

    axes3.text(0.95, 0.8, 'offset=' + "{:.1f}".format(offset) + ' pixels',
               ha='right', c='r', transform=axes3.transAxes)

    # Deal with the tick marks
    
    axes3.xaxis.set_minor_locator(AutoMinorLocator())    
    axes3.tick_params(right=True, left=True, top=True, bottom=True,
                      which='both', direction='in', width=1.5)
    axes3.tick_params(which='minor', length=3)
    axes3.tick_params(which='major', length=5)
    axes3.yaxis.set_minor_locator(AutoMinorLocator())






def read_line_list(filename:str,
                   delta_to_microns:bool=False):
    
    """
    To read a Spextool line list into memory.

    Parameters
    ----------
    filename : str
        The fullpath to a Spextool line list, e.g. "ShortXD_lines.dat".

    Returns
    -------
    dict
        'order_number' : ndarray
            (nlines, ) int array of the order number

        'wavelength' : ndarray
            (nlines, ) str array of the wavelengths (microns).

        'id' : ndarray
            (nlines, ) str array with the line identification.

        'delta_wavelength_left' : ndarray
            (nlines,) float array with the delta_lambda value
            for the lower wavelength of the fit range.

        'delta_wavelength_right' : ndrray
            (nlines, ) float array with the delta_lambda value
            for the upper wavelength of the fit range.

        'fit_type' : ndarray
            (nlines, ) str array giving the type of fit to use.
            "G"=gaussian, "L"=Lorentzian, "C"=centroid

        'num_params' : ndarray
            (nlines, ) int array giving he number of terms in a
            gaussian or lorentzian fit.  see fitpeak1d for details.

    Notes
    -----
    The wavelengths are given in microns and the delta_wavelenths are
    given in Angstroms unless delta_to_microns is set to True.

    """

    #
    # Check parameters
    #

    check_parameter('read_line_list', 'filename', filename, 'str')

    check_parameter('read_line_list', 'delta_to_microns', delta_to_microns,
                    'bool')

    # set up empty lists

    order = []
    swave = []
    lineid = []
    lwin = []
    rwin = []
    fittype = []
    nterms = []

    # Load them in

    f = open(filename, 'r')
    for line in f:

        if line[0] == '#':
            continue

        vals = line.strip().split('|')
        order.append(int(vals[0]))
        swave.append(str(vals[1]).strip())
        lineid.append(str(vals[2]).strip())
        lwin.append(float(vals[3]))
        rwin.append(float(vals[4]))
        fittype.append(str(vals[5]).strip())
        nterms.append(int(vals[6]))

    # Convert to numpy array and dictionary-ify

    lineinfo = {'order': np.array(order), 'wavelength': np.array(swave),
                'id': np.array(lineid),
                'delta_wavelength_left': np.array(lwin),
                'delta_wavelength_right': np.array(rwin),
                'fit_type': np.array(fittype), 'num_parms': np.array(nterms)}

    if delta_to_microns:
        lineinfo['delta_wavelength_left'] = \
            lineinfo['delta_wavelength_left'] / 1e4
        lineinfo['delta_wavelength_right'] = \
            lineinfo['delta_wavelength_right'] / 1e4

    return lineinfo


def read_wavecal_file(file:str):

    """
    To read a pySpextool wavecal calibration file.

    Parameters
    ----------
    file : str
        The fullpath to a calibration file.

        Typically, the file will end in "_wavecalinfo.fits".

    Returns
    -------
    dict
        `"spectra"` : ndarray
            (4,nwave) flat array where:

            spectra[0,:] = wave
            spectra[1,:] = ''flux''
            spectra[2,:] = uncertainty
            spectra[3,:] = bit-set mask

        `"norders"` : int
            The number of orders

        `"orders"` : ndarray of int
            (norders,) int array of the order numbers.  By Spextool convention,
            orders[0] is the order closest to the bottom of the array.

        `"wcaltype"` : str
            '1D' - A single spectrum where wavelength aligns with the columns
                   of the array.
            '1DXD' - Cross-dispersed spectra where wavelength aligns with the
                   columns of the array.
            '2D' - A single spectrum where wavelength does not align with
                   the columns of the array.
            '2DXD' - Cross-dispersed spectra where wavelength does not align
                   with the columns of the array.

        `"arctype"`: str
            'on' - Images taken with lights "on"
            'on-off' - Images taken with lights "on" and "off"    

        `"arcorders"` : ndarray of int
            An array of order numbers to grab from the "arc" image.

        `"skyorders"` : ndarray of int
            An array of order numbers to grab from the "sky" image.
    
        '"linelist"` : str
            The name of the file with the line list.

        `"xranges"` : numpy.ndarray
            An (norders,2) array giving the column numbers over which to
            operate.  xranges[0,0] gives the starting column number for the
            order nearest the bottom of the image and xranges[0,1] gives
            the end column number for said order.

        `"apradius"` : float
            The extraction aperture radius (in arcseconds).

        `"xcororder"` : int
            The order number to be used for the cross-correlation.

        `"dispdeg"` : int
            The polynomial order for the pixel-to-wavelength coefficients.

        if `wcaltype`== '1DXD'

        `"ordrdeg"` : int, optional
            If `"wcaltype" == '1DXD', the polynomial order for the
            2D pixel-to-wavelength coefficients.

        `"p2wcoeffs"` : numpy.ndarray
            Float array of the polynomial coefficients to convert a pixel
            number to a wavelength.

        `"homeorder"` : int
            The order the other orders are scaled to for the fit.

        `"wavefmt"` : str
            The format string for wavelengths.

        `"spatfmt"` : str
            The format string for spatial angles.

    """

    #
    # Check the parameters
    #

    check_parameter('read_wavecal_file', 'file', file, 'str')

    # Open the file
    try:
        hdul = fits.open(file)
    except OSError as e:
        msg = (
            f"Could not open the wavecal file:  {file} \n"
            "Please check that the file downloaded from Git LFS properly. "
            "It should be 3-30 MB. Try `git lfs pull` or"
            "download directly from <URL TBD>"
        )
        logger.error(msg)
        raise e

    # Store the spectrum
    spectra = hdul[0].data
    if spectra.ndim == 2:
        spectra = np.expand_dims(spectra, 0)

    result = {'spectra': spectra}

    # Clean the header and grab important keywords

    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems

    norders = hdul[0].header['NORDERS']
    result.update({'norders': norders})

    val = hdul[0].header['ORDERS'].split(',')
    orders = np.array([int(x) for x in val])
    result.update({'orders': orders})

    wcaltype = hdul[0].header['WCALTYPE']
    result.update({'wcaltype': wcaltype.strip()})

    arctype = hdul[0].header['ARCTYPE']
    result.update({'arctype': arctype.strip()})    

    val = hdul[0].header['ARCORDRS'].split(',')
    arcorders = np.array([int(x) for x in val])
    result.update({'arcorders': arcorders})

    val = hdul[0].header['SKYORDRS'].strip()
    if val == '':

        result.update({'skyorders': []})        

    else:

        val = hdul[0].header['SKYORDRS'].split(',')
        skyorders = np.array([int(x) for x in val])
        result.update({'skyorders': skyorders})        

        
    slits = hdul[0].header['SLITS'].split(',')
    slits = np.array([float(x) for x in slits])
    result.update({'slits': slits})

    usestored = hdul[0].header['USESTORE'].split(',')
    usestored = np.array([eval(x) for x in usestored])
    result.update({'usestored': usestored})

    linelist = hdul[0].header['LINELIST']
    result.update({'linelist': linelist.strip()})

    xranges = np.empty((norders, 2), dtype=int)

    for i in range(norders):
        val = hdul[0].header['OR' + str(orders[i]).zfill(3) + '_XR'].split(',')
        val = [int(x) for x in val]
        xranges[i, :] = val

    result.update({'xranges': xranges})

    val = hdul[0].header['EXTAP']
    result.update({'apradius': val})

    # Get the cross correlation spectrum

    xcororder = hdul[0].header['XCORORDR']
    result.update({'xcororder': xcororder})

    z = orders == xcororder
    xcorspec = np.squeeze(spectra[z, :, :])
    xcorspec[0, :] = np.arange(np.squeeze(xranges[z, 0]),
                               np.squeeze(xranges[z, 1]) + 1)

    result.update({'xcorspec': xcorspec})

    val = hdul[0].header['NLINES']
    result.update({'nlines': val})

    val = hdul[0].header['NGOOD']
    result.update({'ngood': val})

    val = hdul[0].header['NBAD']
    result.update({'nbad': val})

    val = hdul[0].header['RMS']
    result.update({'rms': val})

    dispdeg = hdul[0].header['DISPDEG']
    result.update({'dispdeg': dispdeg})

    if wcaltype == '1dxd':

        val = hdul[0].header['HOMEORDR']
        result.update({'homeorder': val})

        ordrdeg = hdul[0].header['ORDRDEG']
        result.update({'ordrdeg': ordrdeg})

        ncoeffs = (dispdeg + 1) * (ordrdeg + 1)

    elif wcaltype == '1d':

        ncoeffs = (dispdeg + 1)

    # Now get the pixel to wavelength coefficients

    p2wcoeffs = np.empty(ncoeffs)

    for i in range(ncoeffs):
        key = 'P2W_C' + str(i).zfill(2)
        p2wcoeffs[i] = hdul[0].header[key]

    result.update({'coeffs': p2wcoeffs})

    # Get the covariance matrix

    p2wcovar = np.empty((ncoeffs, ncoeffs))
    for i in range(ncoeffs):

        for j in range(ncoeffs):
            key = 'COV_' + str(i).zfill(2) + str(j).zfill(2)
            p2wcovar[j, i] = hdul[0].header[key]

    result.update({'covar': p2wcovar})

    val = hdul[0].header['WAVEFMT']
    result.update({'wavefmt': val.strip()})

    val = hdul[0].header['SPATFMT']
    result.update({'spatfmt': val.strip()})

    hdul.close()

    return result



def read_wavecal_fits(fullpath:str,
                      rotate:bool=True):

    """
    Reads a pyspextool "wavecal" calibration FITS file.

    Parameters
    ----------
    fullpath : str
        The fullpath to a pySpextool wavecal file.

    rotate : {True, False}
        Set to True to rotate the images give the ROTATION keyword.
    
    Returns
    -------
    dict
       
    """

    #
    # Check the parameters
    #

    check_parameter('read_wavecal_fits', 'fullpath', fullpath, 'str')

    check_parameter('read_wavecal_fits', 'rotate', rotate, 'bool')    

    #
    # Read the data 
    #

    hdul = fits.open(fullpath)
    hdul[0].verify('silentfix')

    hdr = hdul[0].header

    if rotate is True:

        rotation = hdr['ROTATION']

    else:

        rotation = 0

    wavecal = idl_rotate(hdul[1].data, rotation)
    spatcal = idl_rotate(hdul[2].data, rotation)


    orders = hdr['ORDERS'].split(',')
    orders = np.array([int(o) for o in orders])
    
    # Grab the indices

    indices = []

    for i in range(hdr['NORDERS']):
        tmp = hdul[i + 3].data
        x = tmp[0, 0, 1:]
        y = tmp[0, 1:, 0]
        xidx = tmp[0, 1:, 1:]
        yidx = tmp[1, 1:, 1:]

        indices.append({'order':orders[i], 'w': x, 'a': y, 'xidx': xidx,
                        'yidx': yidx})

    hdul.close()

    #
    # Now create the wavecalinfo dictionary
    #

    wavecalinfo = {'wavecal': wavecal}
    wavecalinfo.update({'spatcal': spatcal})
    wavecalinfo.update({'rectindices': indices})

    # Do the header info

    wavecalinfo.update({'flatfile': hdr['FLATFILE']})
    wavecalinfo.update({'rotation': hdr['ROTATION']})
    wavecalinfo.update({'norders': hdr['NORDERS']})
    wavecalinfo.update({'orders': hdr['ORDERS']})
    wavecalinfo.update({'exttype': hdr['EXTTYPE']})
    wavecalinfo.update({'wctype': hdr['WCTYPE']})
    wavecalinfo.update({'dispdeg': hdr['DISPDEG']})

    if wavecalinfo['wctype'] in ['1DXD', '2DXD']:

        wavecalinfo.update({'orderdeg': hdr['ORDRDEG']})
        wavecalinfo.update({'homeorder': hdr['HOMEORDR']})

        ncoeffs = (hdr['ORDRDEG'] + 1) * (hdr['DISPDEG'] + 1)

    else:

        ncoeffs = hdr['DISPDEG'] + 1

    wavecalinfo.update({'rms': hdr['RMS']})
    wavecalinfo.update({'nlines': hdr['NLINES']})
    wavecalinfo.update({'ngood': hdr['NGOOD']})
    wavecalinfo.update({'nbad': hdr['NGOOD']})
    wavecalinfo.update({'stored': hdr['STORED']})
    wavecalinfo.update({'version': hdr['VERSION']})

    # Get the xranges, coeffs, and covariance matrix

    xranges = np.empty([hdr['NORDERS'], 2], dtype=int)
    coeffs = np.empty(ncoeffs, dtype=float)
    covar = np.empty([ncoeffs, ncoeffs], dtype=float)

    for i in range(hdr['NORDERS']):

        name = 'OR' + str(orders[i]).zfill(3) + '_XR'
        xranges[i, :] = [int(x) for x in hdr[name].split(',')]

    wavecalinfo.update({'xranges':xranges})
                                                      
    for i in range(ncoeffs):
        name = 'COEFF_' + str(i).zfill(2)
        coeffs[i] = hdr[name]

    for i in range(ncoeffs):

        for j in range(ncoeffs):
            name = 'COV_' + str(i).zfill(2) + str(j).zfill(2)
            covar[j, i] = hdr[name]

    wavecalinfo.update({'coefficients':coeffs})
    wavecalinfo.update({'covariance':covar})                           

                           
    return wavecalinfo


def simulate_wavecal_1dxd(ncols:int,
                          nrows:int,
                          edgecoeffs:npt.ArrayLike,
                          xranges:npt.ArrayLike,
                          slith_arc:float):

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
    wavecal : ndarray
        Wavecal (nrows,ncols) array where each pixel is set to its
        wavelength which in this case is the column number.

    spatcal : ndarray
        Spatcal (nrows,ncols) array where each pixel is set to its 
        angular position on the sky (in arcseconds).

    indices : list
        An (norders,) list where each element is a dictionary with the
        following keys:

        'x' : ndarray
            An (ncols,) array of x values (in pixels).

        'y' : ndarray
            An(nrows,) array of y values (in arcseconds).

        'xidx' : ndarray
            An (nrows, ncols) array of x indices.

        'yidx' : ndarray
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




def wavecal_solution_1d(orders:npt.ArrayLike,
                        line_info:dict,
                        dispersion_degree:float | int,
                        xd_info:dict=None,
                        verbose:bool=False,
                        qa_plotnumber:int=None,
                        qa_figuresize:tuple=(7,9),
                        qa_fontsize:tuple=12,                        
                        qa_show:bool=False,
                        qa_showscale:float | int=1,
                        qa_showblock:bool=False,
                        qa_fullpath:str=None):

    """
    To calculate a wavelength solution in the pySpextool 1D or 1DXD case.

    Parameters
    ----------
    orders : ndarray
        An (norders,) array with the order numbers.

    line_info : dict
        `"order"` : ndarray
             A (nlines,) array giving the order number of each line.

        `"x"` : ndarray
             A (nlines,) array giving the x position of each line.

        `"wavelength"` : ndarray
             A (nlines,) array giving the wavelength of each line.

        `"goodbad"` : ndarray
             A (nlines,) ndarray.  0 = a bad fit, 1= a good fit.
        
    dispersion_degree : int
        The polynomial degree to fit in the dispersion direction.

    xdinfo : dict, optional
        `"homeorder"' : int
            The home order used to scale lines.  See Notes.

        '"orderdeg"' : int
            The polynomial degree to fit in the order dimension.

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].

    qa_size : tuple
        The size of the QA plot window in inches.  The default is (7,9).  
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_scale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does NOT affect plots written
        to disk.  Set to None to default to setup.state['qa_scale'].
        
    qa_block : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_writeinfo : dict, optional
        `"figsize"` : tuple
            A (2,) tuple giving the figure size.

        `"fullpath"` : str
            A string giving the path to write the file.

        '"filename"` : str
            The root of the filename.

        '"extension"` : str
            The file extension.

    Returns
    -------
    numpy.ndarray
        An array of coefficients generated from poly_fit_2d.

    Notes
    -----
    The program determines a single wavelength solution for all orders
    in a spectrogram.   The technique is based on Hall et al. (1994, PASP,
    106, 315).  The ratio of two different wavelengths, lambda_1 and
    lambda_2, diffracted at the same angle but located in two different
    orders m_1 and m_2, is given by,

    lambda_1/lambda_2 = m_2/m_1

    """

    #
    # Check parameters and keywords
    #

    check_parameter('wavecal_solution_1d', 'orders', orders, 'ndarray')

    check_parameter('wavecal_solution_1d', 'line_info', line_info, 'dict')

    check_parameter('wavecal_solution_1d', 'dispersion_degree',
                    dispersion_degree, ['float','int'])        
    
    check_parameter('wavecal_solution_1d', 'xd_info', xd_info, 
                    ['dict', 'NoneType'])

    check_parameter('wavecal_solution_1d', 'verbose', verbose, 'bool')

#    check_parameter('wavecal_solution_1d', 'qa_size', qa_size, 'tuple')
#    
#    check_parameter('wavecal_solution_1d', 'qa_show', qa_show, 'bool')
#
#    check_parameter('wavecal_solution_1d', 'qa_scale', qa_scale,
#                    ['float','int'])            
#
#    check_parameter('wavecal_solution_1d', 'qa_block', qa_block, 'bool')
#
#    check_parameter('wavecal_solution_1d', 'qa_writeinfo', qa_writeinfo,
#                    ['dict','NoneType'])    
    
    
    check_qakeywords(verbose=verbose)
    
    #
    # We just need to test whether this is a 1D or 1DXD fit
    #
    
    if xd_info is None: # This is a 1D fit

        # Do the fit
        
        fit = polyfit_1d(line_info['x'],
                         line_info['wavelength'].astype(np.float16),
                         dispersion_degree,
                         goodbad=line_info['goodbad'],
                         robust={'thresh': 4, 'eps': 0.1})
        
        # Do the QA plotting
        
        residuals = (line_info['wavelength'].astype(np.float16) -
                     fit['yfit']) * 1e4
            
        if qa_show is True:

            pl.ion()
            plot_1d_residuals(qa_showsize,
                              residuals,
                              line_info['x'],
                              fit['goodbad'],
                              fit['rms'],
                              dispersion_degree)
            pl.show()
            pl.pause(1)
        
#        if qa_writeinfo is not None:
#
#            pl.ioff()
#            plot_1d(qa_writeinfo['figsize'], residuals, line_info['x'],
#                      fit['goodbad'], fit['rms'], dispersion_degree)
#            pl.savefig(join(qa_writeinfo['filepath'],
#                            qa_writeinfo['filename'] +
#                            '_residuals'+qa_writeinfo['extension']))
#            pl.close()
        

    else: # This is a 1DXD fit

        # Scale the wavelengths to the home order


        scales = line_info['order'] / xd_info['homeorder']

        scaled_wavelengths = line_info['wavelength'].astype(np.float64) * scales
        
        # Now do the fit

#        f = open('python.dat', 'w')
#
#        for k in range(len(line_info['x'])):
#                       
#            f.write('%s %s %s %s %s\n' % (line_info['order'][k], \
#                                       line_info['x'][k], \
#                                       line_info['wavelength'][k],
#                                       scaled_wavelengths[k], \
#                                       line_info['goodbad'][k]))
#
#        f.close()

        
        fit = polyfit_2d(line_info['x'],
                         line_info['order'],
                         scaled_wavelengths,
                         dispersion_degree,
                         xd_info['orderdeg'],                          
                         goodbad=line_info['goodbad'],
                         robust={'thresh': 4, 'eps': 0.1})
        
        # Now do the qa plot

        residuals = (scaled_wavelengths - fit['zfit']) * 1e4

        if qa_show is True:
                        
            plot_1dxd_residuals(qa_plotnumber,
                                (qa_figuresize[0]*qa_showscale,
                                 qa_figuresize[1]*qa_showscale),
                                qa_fontsize*qa_showscale,
                                residuals,
                                line_info['order'],
                                line_info['x'],
                                orders,
                                fit['goodbad'],
                                fit['rms'],
                                dispersion_degree,
                                xd_info['orderdeg'])

            pl.show(block=qa_showblock)
            if qa_showblock is False:  
            
                pl.pause(1)

        if isinstance(qa_fullpath,str):

            plot_1dxd_residuals(None,
                                qa_figuresize,
                                qa_fontsize,
                                residuals,
                                line_info['order'],
                                line_info['x'],
                                orders,
                                fit['goodbad'],
                                fit['rms'],
                                dispersion_degree,
                                xd_info['orderdeg'])

          
            pl.savefig(qa_fullpath)
            pl.close()
                            
    #
    # Now get the results together and return the results
    #
    
    # Get the number of good and bad points

    zgood = fit['goodbad'] == 1
    ntot = np.sum(zgood)
            
    zbad = fit['goodbad'] == 0
    nbad = np.sum(zbad)

    # Create the dictionary
    
    solution = {'coeffs':fit['coeffs'],
                'covar':fit['coeffs_covar'],
                'rms':fit['rms'],
                'nlines':ntot,
                'ngood':ntot-nbad,
                'nbad':nbad}
                
    return solution


def write_wavecal1d_fits(order_mask:npt.ArrayLike,
                         xranges:npt.ArrayLike,
                         coeffs:npt.ArrayLike,
                         covar:npt.ArrayLike,
                         dispersion_degree:int,
                         rms:float | int,
                         nlines:int,
                         ngood:int,
                         nbad:int,
                         wavecal_pixels:npt.ArrayLike,
                         stored_solution_offset: int | float,
                         spatcal:npt.ArrayLike,
                         indices,
                         rotate:int,
                         flatname:str,
                         skysnames:str,
                         oname:str,
                         version:str,
                         xdinfo:dict=None,
                         stored_solution:bool=False):
    
    """
    To write a Spextool 1D (or 1DXD) wavecal file to disk.

    Parameters
    ----------    
    orders : ndarray of int
        An (norders,) int array of the order numbers.  By Spextool convention,
        orders[0] is the order closest to the bottom of the array.


    
    flatname : str
        A str giving the name of the pySpextool flat field file.
    
    stored_solution : {False, True}
        Set to True to indicate that a stored solution was used.
        Set to False to indicate that the solution was generated.
    
    """


    
    
    #
    # Check parameters
    #


    
    # Get basic things

    nrows, ncols = np.shape(wavecal_pixels)    
    orders = np.unique(order_mask)[1:]    
    norders = len(orders)

    # Create the primary HDU

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Generate the wavelengths

    wavecal = np.full_like(wavecal_pixels,np.nan)
    if xdinfo is None:

        # This is the pure 1D case.

        for i in range(norders):

            #
            # Create the 2D wavecal image
            #
            
            # Find the pixels in the order
            
            z = np.where(order_mask == orders[i])
            
            # Evaluate the polynomial at these values

            wgrid = poly_1d(wavecal_pixels[z]+stored_solution_offset, coeffs)
            
            # Now store the results
                        
            wavecal[z] = wgrid

            #
            # compute the wavelengths for the indices
            #

            wgrid = poly_1d(indices[i][0,0,1:]+stored_solution_offset, coeffs)

            indices[i][0,0,1:] = wgrid
            
            
        ncoeffs = dispersion_degree + 1
        wctype = '1D'

    else:

        # This is the 1DXD case.

        for i in range(norders):

            #
            # Create the 2D wavecal image
            #
            
            # Find the pixels in the order and get the scale factor
            
            z = np.where(order_mask == orders[i])
            scale = xdinfo['homeorder'] / orders[i]

            # Evaluate the polynomial at these values            

            wgrid = poly_2d(wavecal_pixels[z]+stored_solution_offset,
                            order_mask[z],
                            dispersion_degree,
                            xdinfo['orderdeg'],
                            coeffs)*scale
        
            # Now store the results
                        
            wavecal[z] = wgrid            

            #
            # compute the wavelengths for the indices
            #
        
            order_array = np.full_like(indices[i][0,0,1:],orders[i])
            
            wgrid = poly_2d(indices[i][0,0,1:]+stored_solution_offset,
                            order_array,
                            dispersion_degree,
                            xdinfo['orderdeg'],
                            coeffs) * scale

            indices[i][0,0,1:] = wgrid
            
        ncoeffs = (dispersion_degree + 1) * (xdinfo['orderdeg'] + 1)

        wctype = '1DXD'

    # Fill in header things
        
    hdr['FLATFILE'] = (flatname, ' Associated flat-field image')
    hdr['SKYFILES'] = (skysnames, ' Sky images')    
    hdr['ROTATION'] = (rotate, ' IDL rotate value')
    hdr['NORDERS'] = (norders, ' Number of orders identified')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), 'Orders identified')
    hdr['EXTTYPE'] = ('1D', ' Extraction type')
    hdr['WCTYPE'] = (wctype, ' Wavelength calibration type')
    hdr['DISPDEG'] = (dispersion_degree, ' Dispersion fit degree')

    if xdinfo is not None:
        hdr['ORDRDEG'] = (xdinfo['orderdeg'], ' Order fit degree')
        hdr['HOMEORDR'] = (xdinfo['homeorder'], ' Home Order')

    hdr['RMS'] = (rms, 'RMS of fit in Angstroms')
    hdr['NLINES'] = (nlines, ' Number of lines in the fit')
    hdr['NGOOD'] = (ngood, ' Number of good points')
    hdr['NBAD'] = (nbad, 'Number of bad points')
    hdr['STORED'] = (stored_solution, ' Use stored solution?')
    hdr['VERSION'] = (version, 'pySpextool version')

    # Write the extract ranges

    for i in range(norders):
        name = 'OR' + str(orders[i]).zfill(3) + '_XR'
        comment = ' Extraction range for order ' + str(orders[i]).zfill(3)
        hdr[name] = (','.join(str(x) for x in xranges[i, :]), comment)

    # Write the coefficients

    for i in range(ncoeffs):
        name = 'COEFF_' + str(i).zfill(2)
        comment = ' c' + str(i) + ' coefficient for solution'
        hdr[name] = (coeffs[i], comment)

    # Write the covariance 

    for i in range(ncoeffs):

        for j in range(ncoeffs):
            name = 'COV_' + str(i).zfill(2) + str(j).zfill(2)
            comment = str(i) + ',' + str(j) + \
                      ' (col,row) element of the covariance'
            hdr[name] = (covar[j, i], comment)

    # Write the results

    wavecal = np.float32(wavecal)
    spatcal = np.float32(spatcal)
            
    wavecal_hdu = fits.ImageHDU(idl_unrotate(wavecal, rotate))
    spatcal_hdu = fits.ImageHDU(idl_unrotate(spatcal, rotate))

    list_hdu = [phdu, wavecal_hdu, spatcal_hdu]

    # Add the indices

    for i in range(norders):

        idx_hdu = fits.ImageHDU(np.float32(indices[i]))
        list_hdu.append(idx_hdu)

    hdu = fits.HDUList(list_hdu)
    hdu.writeto(oname, overwrite=True)




    
        
    

