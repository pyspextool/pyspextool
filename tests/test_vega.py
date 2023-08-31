import numpy as np
from astropy.io import fits
import glob
import os

import pyspextool as ps
from pyspextool import config as setup
from pyspextool.telluric.vega import vega_xcorrelate
from pyspextool.io.read_spectra_fits import read_spectra_fits
import pyspextool.fit.polyfit as polyfit


def test_vega_xcorrelate():

    """
    Not nice, just a place holder

    """

    
    ps.pyspextool_setup(qa_path='tests/test_data/misc_data/',
                        qa_extension='.png')

    objspectra, hdrdict = read_spectra_fits('tests/test_data/misc_data/spectra11-18.fits')

    object_wavelength = objspectra[3,0,:]
    object_flux = objspectra[3,1,:]
    
    
    hdul = fits.open(setup.state['package_path']+'/data/Vega5000.fits') 
    data  = hdul[1].data
    
    vega_wavelength = data['wavelength']
    vega_flux = data['flux density']
    vega_continuum = data['continuum flux density']
    

    hdul.close()

    #
    # Get the wavelengths for the datda
    #

    z1 = (object_wavelength > 0.975)
    z2 = (object_wavelength < 0.99)
    

    zleft = np.logical_and(z1,z2)


    z1 = (object_wavelength > 1.02)
    z2 = (object_wavelength < 1.04)
    

    zright = np.logical_and(z1,z2)

    z = np.logical_or(zleft,zright)

    result = polyfit.poly_fit_1d(object_wavelength[z],object_flux[z],2)
    object_continuum = polyfit.poly_1d(object_wavelength,result['coeffs'])
    
    object_normalized = object_flux/object_continuum

    z1 = (object_wavelength > 0.97)
    z2 = (object_wavelength < 1.05)


    z = np.logical_and(z1,z2)

    qafileinfo = {'figsize': (8.5, 11), 'filepath': setup.state['qa_path'],
                  'filename': 'test', 'extension': setup.state['qa_extension']}
 
    vega_xcorrelate(object_wavelength[z], object_normalized[z],
                    vega_wavelength, vega_flux, vega_continuum,
                        resolving_power=2000, qa_show=False,
                        qa_fileinfo=qafileinfo)

    png_files = glob.glob("tests/test_data/misc_data/test_xcorrelate.png")

    assert len(png_files) == 1
    
    # CLEANUP
    # remove generated files
    for files in png_files:
        os.remove(files)


