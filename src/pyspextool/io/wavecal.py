import numpy as np
from astropy.io import fits

from pyspextool.fit.polyfit import poly_2d
from pyspextool.spectroscopy.make_order_mask import make_order_mask
from pyspextool.utils.arrays import idl_unrotate


def read_line_list(file, delta_to_microns=False):
    """
    To read a Spextool line list into memory.

    Parameters
    ----------
    file : str
        The fullpath to a Spextool line list, e.g. "ShortXD_lines.dat".

    Returns
    -------
    dict
        `"order_number"` : numpy.ndarray
            (nlines,) int array of the order number

        `"wavelength"` : numpy.ndarray
            (nlines,) str array of the wavelengths (microns).

        `"id"` : numpy.ndarray
            (nlines,) str array with the line identification.

        `"delta_wavelength_left"` : numpy.ndarray
            (nlines,) float array with the delta_lambda value
            for the lower wavelength of the fit range.

        `"delta_wavelength_right"` : numpy.ndrray
            (nlines,) float array with the delta_lambda value
            for the upper wavelength of the fit range.

        `"fit_type"` : numpy.ndarray
            (nlines,) str array giving the type of fit to use.
            "G"=gaussian, "L"=Lorentzian, "C"=centroid

        `"num_params"` : numpy.ndarray
            (nlines,) int array giving he number of terms in a
            gaussian or lorentzian fit.  see fitpeak1d for details.

    Notes
    -----
    The wavelengths are given in microns and the delta_wavelenths are
    given in Angstroms.

    Examples
    --------
    later

    """

    # set up empty lists

    order = []
    swave = []
    lineid = []
    lwin = []
    rwin = []
    fittype = []
    nterms = []

    # Load them in

    f = open(file, 'r')
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


def read_wavecal_file(file):
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
        `"spectra"` : numpy.ndarray
            (4,nwave) flat array where:

            spectra[0,:] = wave
            spectra[1,:] = ''flux''
            spectra[2,:] = uncertainty
            spectra[3,:] = bit-set mask

        `"norders"` : int
            The number of orders

        `"orders"` : numpy.ndarray of int
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

    Examples
    --------
    later

    """

    # Open the file
    
    hdul = fits.open(file)

    # Store the spectrum 

    spectra = hdul[0].data
    if spectra.ndim == 2:
        spectra = np.expand_dims(spectra,0)

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
    xcorspec[0, :] = np.arange(xranges[z, 0], xranges[z, 1] + 1)

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

    p2wcovar = np.empty((ncoeffs,ncoeffs))
    for i in range(ncoeffs):

        for j in range(ncoeffs):
            key = 'COV_' + str(i).zfill(2) + str(j).zfill(2)
            p2wcovar[j,i] = hdul[0].header[key]
           
    result.update({'covar': p2wcovar})
   
    val = hdul[0].header['WAVEFMT']
    result.update({'wavefmt': val.strip()})

    val = hdul[0].header['SPATFMT']
    result.update({'spatfmt': val.strip()})

    hdul.close()

    return result


def write_wavecal_1dxd(ncols, nrows, orders, edgecoeffs, xranges, coeffs,
                       covar, dispersion_degree, order_degree, home_order, rms,
                       nlines, ngood, nbad, wavecal, spatcal, rotate, flatname,
                       oname, version, stored_solution=False, overwrite=True):
    # Get basic things

    norders = len(orders)

    # Create the primary HDU

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Create the order mask

    omask = make_order_mask(ncols, nrows, edgecoeffs, xranges, orders)

    # Generate the wavelengths

    for i in range(norders):
        z = np.where(omask == orders[i])
        scale = home_order / orders[i]
        wavecal[z] = poly_2d(wavecal[z], omask[z], dispersion_degree,
                             order_degree, coeffs) * scale

    hdr['FLATFILE'] = (flatname, ' Assocaited flat-field image')
    hdr['NORDERS'] = (norders, ' Number of orders identified')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), 'Orders identified')
    hdr['EXTTYPE'] = ('1D', ' Extraction type')
    hdr['WCTYPE'] = ('1DXD', ' Wavelength calibration type')
    hdr['STORED'] = (stored_solution, ' Use stored solution?')
    hdr['RMS'] = (rms, 'RMS of 1DXD wavecal fit in Angstroms')
    hdr['NLINES'] = (nlines, ' Number of lines in the fit')
    hdr['NGOOD'] = (ngood, ' Number of good points')
    hdr['NBAD'] = (nbad, 'Number of bad points')
    hdr['VERSION'] = (version, 'pySpextool version')

    # Write the extract ranges
    
    for i in range(norders):
        name = 'OR' + str(orders[i]).zfill(3) + '_XR'
        comment = ' Extraction range for order ' + str(orders[i]).zfill(3)
        hdr[name] = (','.join(str(x) for x in xranges[i, :]), comment)

    # Write the coefficients

    ncoeffs = (dispersion_degree+1) * (order_degree+1)

    for i in range(ncoeffs):

        name = 'COEFF_' + str(i).zfill(2)
        comment = ' c'+str(i)+' coefficient for solution'
        hdr[name] = (coeffs[i], comment)

    # Write the covariance 

    for i in range(ncoeffs):

        for j in range(ncoeffs):        

            name = 'COV_'+str(i).zfill(2)+str(j).zfill(2)
            comment = str(i)+','+str(j)+\
            ' (col,row) element of the covariance'
            hdr[name] = (covar[j,i], comment)


        
    # Write the results

    wavecal_hdu = fits.ImageHDU(idl_unrotate(wavecal, rotate))
    spatcal_hdu = fits.ImageHDU(idl_unrotate(spatcal, rotate))

    hdu = fits.HDUList([phdu, wavecal_hdu, spatcal_hdu])
    hdu.writeto(oname, overwrite=overwrite)

    # fits.writeto('wavecal.fits', wavecal, overwrite=True)
