import numpy as np
from astropy.io import fits
import glob
import os
import pytest

import pyspextool as ps
from pyspextool import config as setup
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.utils.spectra import normalize_continuum

def test_normalize_continuum():

    
    ps.pyspextool_setup(qa_path='tests/test_data/misc_data/',
                        qa_extension='.png')

    objspectra, hdrdict = read_spectra_fits('tests/test_data/misc_data/spectra11-18.fits')

    wavelength = objspectra[3,0,:]
    flux = objspectra[3,1,:]
    

    with pytest.raises(ValueError):
        normalize_continuum(wavelength, flux, np.array([1,10]), 4)

    with pytest.raises(ValueError):
        normalize_continuum(wavelength, flux, np.array([0,1.1]), 4)

    with pytest.raises(ValueError):
        normalize_continuum(wavelength, flux, np.array([1.1,1]), 4)        

    qafileinfo = {'filepath':'tests/test_data/misc_data/',
                  'filename':'RV',
                  'extension':'.png'}
    normalize_continuum(wavelength, flux, np.array([0.98,0.995, 1.02,1.04]), 4,
                        qa_fileinfo=qafileinfo)        
        
    png_files = glob.glob("tests/test_data/misc_data/RV.png")

    # CLEANUP
    # remove generated files
    for files in png_files:
        os.remove(files)



    
    
    #
    # Get the wavelengths for the datda
    #

#    z1 = (object_wavelength > 0.975)
#    z2 = (object_wavelength < 0.99)
#    
#
#    zleft = np.logical_and(z1,z2)
#
#
#    z1 = (object_wavelength > 1.02)
#    z2 = (object_wavelength < 1.04)
#    
#
#    zright = np.logical_and(z1,z2)
#
#    z = np.logical_or(zleft,zright)
#
#    result = polyfit.poly_fit_1d(object_wavelength[z],object_flux[z],2)
#    object_continuum = polyfit.poly_1d(object_wavelength,result['coeffs'])
#    
#    object_normalized = object_flux/object_continuum
#
#    z1 = (object_wavelength > 0.97)
#    z2 = (object_wavelength < 1.05)
#
#
#    z = np.logical_and(z1,z2)
#
#    qafileinfo = {'figsize': (8.5, 11), 'filepath': setup.state['qa_path'],
#                  'filename': 'test', 'extension': setup.state['qa_extension']}
# 
#    vega_xcorrelate(object_wavelength[z], object_normalized[z],
#                    vega_wavelength, vega_flux, vega_continuum,
#                        resolving_power=2000, qa_show=False,
#                        qa_fileinfo=qafileinfo)
#
#    png_files = glob.glob("tests/test_data/misc_data/test_xcorrelate.png")

#    assert len(png_files) == 1
    
    # CLEANUP
    # remove generated files
#    for files in png_files:
#        os.remove(files)


