import numpy as np
from astropy.io import fits

def read_line_list(file):

    """
    To read a Spextool line list into memory.


    Parameters
    ----------
    file : str
        The fullpath to a Spextool line list.


    Returns
    -------
    dict
        ``"order"``
            (numpy.ndarray of int) ->  (nlines,) the order number

        ``"swave"``
            (numpy.ndarray of str) ->  (nlines,) the wavelengths (microns)

        ``"lineid"``
            (numpy.ndarray of str) ->  (nlines,) the line identification

        ``"leftwin"``
            (numpy.ndarray of float) ->  (nlines,) the width of the window 
            to the left of line center to fit (arcseconds)

        ``"rghtwin"``
            (numpy.ndarray of float) ->  (nlines,) the width of the window 
            to the right of line center to fit (arcseconds)

        ``"fittype"``
            (numpy.ndarray of str) -> (nlines,) G=gaussian, L=Lorentzian, 
            C=centroid

        ``"nterms"``
            (numpy.ndarray of int) -> (nlines,) The number of terms in a
            gaussian or lorentzian fit.  see fitpeak1d.

    Notes
    -----
    None

    Examples
    --------
    later

    Modification History
    --------------------
    2022-07-05 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program readlinelist.pro

    """
    # Test for `file` being a string

    #if isinstance(file,str) is False:
    #    raise TypeError(file must be a string.)

    # Check to make sure file exists.

    # set up empty strings
    
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
        swave.append(vals[1])
        lineid.append(vals[2])
        lwin.append(float(vals[3]))
        rwin.append(float(vals[4]))
        fittype.append(vals[5])
        nterms.append(int(vals[6]))

# Convert to numpy array and dictionary-ify 
        
    lineinfo = {'order': np.array(order), 'strwave': np.array(swave),
                'lindid': np.array(lineid), 'leftwin': np.array(lwin),
                'rghtwin': np.array(rwin), 'fittype': np.array(fittype),
                'nterms': np.array(nterms)}

    return lineinfo


def read_wavecal_file(file):
    """
    To read a pySpextool wavecal calibration file.

    Input Parameters
    ----------------
    file : str
        The fullpath to a calibration file.

        Typicallt the file will end in "_wavecalinfo.fits".

    Returns
    -------
    dict
        A dictionary with the following keywords:

        spectra : numpy.ndarray of float
            (4,n_wave) array
            spectra[0,:] = wave
            spectra[1,:] = ''flux''
            spectra[2,:] = uncertainty
            spectra[3,:] = bit-set mask

        norders : int
            The number of orders

        orders : numpy.ndarray of int
            The order numbers.  orders[0] is the order closest to the bottom
            of the array.

        wcaltype : str
            '1D' - A single spectrum where wavelength aligns with the columns
                   of the array.
            '1DXD' - Cross-dispersed spectra where wavelength aligns with the
                   columns of the array.
            '2D' - A single spectrum where wavelength does not align with
                   the columns of the array.
            '2DXD' - Cross-dispersed spectra where wavelength does not align
                   with the columns of the array.

        linelist : str
            The name of the file with the line list.

        xranges : array_like of float
            An (norders,2) array giving the column numbers over which to
            operate.  xranges[0,0] gives the starting column number for the
            order nearest the bottom of the image and xranges[0,1] gives
            the end column number for said order.

        apradius : float
            The extraction aperture radius (in arcseconds).

        xcororder : int
            The order number to be used for the cross-correlation.

        dispdeg : int
            The polynomial order for the pixel-to-wavelength coefficients.

        if `wcaltype`== '1DXD'

        ordrdeg : int
            The polynomial order for the 2D pixel-to-wavelength coefficients.
        p2wcoeffs : numpy.ndarray of float
            The polynomial coefficients to convert a pixel number to a
            wavelength.

        homeorder : int
            The order the other orders are scaled to for the fit.

        wavefmt : str
            The format string for wavelengths.

        spatfmt : str
            The format string for spatial angles.

    Notes
    -----


    Examples
    --------
    later

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_readwavecalinfo.pro

    """

    hdul = fits.open(file)
    spectra = hdul[0].data

    # Clean the header and grab important keywords

    hdul[0].verify('silentfix')  # this was needed for to correct hdr problems

    result = {'spectra': spectra}

    #    naps = hdul[0].header['NAPS']
    #    result.update({'naps':naps})

    norders = hdul[0].header['NORDERS']
    result.update({'norders': norders})

    val = hdul[0].header['ORDERS'].split(',')
    orders = np.array([int(x) for x in val])
    result.update({'orders': orders})

    wcaltype = hdul[0].header['WCALTYPE']
    result.update({'wcaltype': wcaltype})

    linelist = hdul[0].header['LINELIST']
    result.update({'linelist': linelist})

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

    if wcaltype == '1DXD':

        dispdeg = hdul[0].header['DISPDEG']
        result.update({'dispdeg': dispdeg})

        ordrdeg = hdul[0].header['ORDRDEG']
        result.update({'ordrdeg': ordrdeg})

        # Now get the pixel to wavelength coefficients

        ncoeffs = (dispdeg + 1) * (ordrdeg + 1)
        p2wcoeffs = np.empty(ncoeffs)

        for i in range(ncoeffs):
            key = 'P2W_C' + str(i).zfill(2)
            p2wcoeffs[i] = hdul[0].header[key]

        result.update({'p2wcoeffs': pswcoeffs})

        val = hdul[0].header['HOMEORDR']
        result.update({'homeorder': val})

    val = hdul[0].header['WAVEFMT']
    result.update({'wavefmt': val})

    val = hdul[0].header['SPATFMT']
    result.update({'spatfmt': val})

    hdul.close()

    return result
