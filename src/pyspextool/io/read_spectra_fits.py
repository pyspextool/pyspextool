import numpy as np
from astropy.io import fits

from pyspextool.io.check import check_parameter


def read_spectra_fits(file):
    """
    To read a pyspextool FITS file and keywords.

    Parameters
    ----------
        file : str
            The fullpath to a pyspextool spectral FITS file.

    Returns
    -------
        tuple (ndarray, dict)

            tuple(0) : ndarray
                A (norders*napertures, 4, nwavelength) array.

            tuple(1) : dic

                `'header'` : astropy.io.fits.header.Header

                `'instr'` : str

                `'obsmode'` : str

                `'norders'` : int
                
                `'orders'` : ndarray

                `'dispersions'` : ndarray
                
                `'xunits'` : str
                
                `'yunits'` : str
                
                `'slith_pix'` : int
                
                `'slith_arc'` : float
                
                `'slitw_pix'` : int
                
                `'slitw_arc'`: float
                
                `'creationmodule'` : str 
                
                `'history'` : astropy.io.fits.header._HeaderCommentaryCards

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
    # Adjust if it is an old Spextool file
    #

    if np.ndim(spectra) == 2:

        spectra = np.expand_dims(spectra, 0)
    
    #
    # Check to see if it is a pySpextool  FITS file.
    #

    try:

        header['NAPS']

    except:

        message = file + ' is not a pySpextool FITS file.'
        raise ValueError(message)

    #
    # Start pulling the keywords
    #

    dictionary = {'header': header, 'instr': header['INSTR'],
                  'obsmode': header['MODE'], 'norders': header['NORDERS']}

    val = header['ORDERS'].split(',')
    orders = np.array([int(x) for x in val])

#    dispersions = np.empty(dictionary['norders'])

#    for i in range(dictionary['norders']):

#        key = 'OR'+str(orders[i]).zfill(3)+'_DP'
#        dispersions[i] = header[key]

    

    add = {'orders': orders,
           'napertures': header['NAPS'],
#           'dispersions':dispersions,
           'xunits': header['XUNITS'],
           'yunits': header['YUNITS'],
           'lxunits': header['LXUNITS'],
           'lyunits': header['LYUNITS'],
           'lxlabel': header['LXLABEL'],
           'lylabel': header['LYLABEL'],
           'lulabel': header['LULABEL'],           
           'slith_pix': header['SLTH_PIX'],
           'slith_arc': header['SLTH_ARC'],
           'slitw_pix': header['SLTW_PIX'],
           'slitw_arc': header['SLTW_ARC'],
           'resolvingpower':header['RP'],
           'module': header['MODULE'],
           'history': header['HISTORY']}

    dictionary.update(add)

    #
    # Return the results
    #

    return spectra, dictionary
