import numpy as np
from astropy.io import fits

from pyspextool.io.check import check_parameter
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.fitsheader import get_headerinfo


def read_spectra_fits(file:str):

    """
    To read a pyspextool FITS file and return ndarray and dictionary.

    Parameters
    ----------
        file : str
            The fullpath to a pyspextool spectral FITS file.

    Returns
    -------
        tuple (ndarray, dict)

            tuple(0) : ndarray
                An (norders*napertures, 4, nwavelength) array of spectra.

            tuple(1) : dict

                `'header'` : astropy.io.fits.header.Header

                `'instr'` : str
                    A str giving the name of the instrument

                `'obsmode'` : str
                    A str giving the obseving mode.

                `'norders'` : int
                    An int giving the number of orders.
                
                `'orders'` : ndarray
                    An (norders,) array giving the order numbers.

                `'xunits'` : str
                    A str giving the xunits, e.g. 'um'.
                
                `'yunits'` : str
                    A str giving the yunits, e.g. 'W m-2 um-1'
                
                `'slith_pix'` : int
                    An int giving the nominal slit height in pixels.
                
                `'slith_arc'` : float
                    A float giving the slit height in arcseconds.
                
                `'slitw_pix'` : int
                    An int giving the slit width in pixels.
                
                `'slitw_arc'`: float
                    A float giving the slight width in arcseconds.
                
                `'module'` : str 
                    A str giving the pyspextool module that created the file.
                
                `'history'` : astropy.io.fits.header._HeaderCommentaryCards
                    The astropy header.

                `'wranges'` : list
                    A (norders*naps) list of (2,) arrays giving the minimum 
                    and maximum wavelength value for each order/aperture.

    """

    #
    # Check parameters
    #

    check_parameter('read_spectra_fits', 'file', file, 'str')

    #
    # Read the file
    #

    hdul = fits.open(file)
    hdul[0].verify('silentfix')  # this was needed to correct hdr problems

    spectra = hdul[0].data
    header = hdul[0].header
    hdul.close()
    
    #
    # Check to see if it is a (py)Spextool FITS file.  Assume yes if it has 
    # these three keywords that it is.
    #

    try:

        header['NAPS']
        header['ORDERS']
        header['NORDERS']

    except:

        message = file + ' is not a pySpextool FITS file.'
        raise pySpextoolError(message)

    #
    # Now add minimum information to output dictionary
    #

    val = header['ORDERS'].split(',')
    orders = np.array([int(x) for x in val])

    dictionary = {'astropyheader': header,
                  'napertures': header['NAPS'],
                  'norders': header['NORDERS'],
                  'orders':orders}
    
    #
    # Now add optional things
    #

    # First set up a mapping between FITS keywords and dictionary keys.  

    keys = [('instr', 'INSTR'),
            ('obsmode', 'MODE'),
            ('norders', 'NORDERS'),
            ('xunits', 'XUNITS'),
            ('yunits', 'YUNITS'),
            ('lxunits', 'LXUNITS'),
            ('lyunits', 'LYUNITS'),
            ('lxlabel', 'LXLABEL'),
            ('lylabel', 'LYLABEL'),
            ('slith_pix', 'SLTH_PIX'),
            ('slith_arc', 'SLTH_ARC'),
            ('slitw_pix', 'SLTW_PIX'),
            ('slitw_arc', 'SLTW_ARC'),
            ('module', 'MODULE'),
            ('resolving_power', 'RP'),            
            ('history','HISTORY'),
            ]
            
    # Now go through and search for these FITS keywords.

    for key in keys:

        try:

            dictionary[key[0]] = header[key[1]]

        except:

            dictionary[key[0]] = None

    #
    # Obtain the wavelength ranges of each order/aperture
    #

    wavelength_ranges = []

    for i in range(header['NORDERS']):

        for j in range(header['NAPS']):
            
            idx = i * header['NAPS'] + j
            min = np.nanmin(spectra[idx, 0, :])
            max = np.nanmax(spectra[idx, 0, :])
            
        wavelength_ranges.append(np.array([min, max]))

    dictionary['wranges'] = wavelength_ranges
    
    #
    # return results
    #

    return spectra, dictionary

