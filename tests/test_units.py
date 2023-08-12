import numpy as np
from pyspextool.utils.units import *

def test_convert_wavelength():

    wavelength = np.array([1.1,2])

    result = convert_wavelength(wavelength, 'um', 'A')
    np.testing.assert_array_equal(result,np.array([1.1e4,2e4]))

    result = convert_wavelength(wavelength, 'um', 'um')
    np.testing.assert_array_equal(result,np.array([1.1,2]))

    result = convert_wavelength(wavelength, 'um', 'nm')
    np.testing.assert_array_equal(result,np.array([1.1e3,2e3]))    

def test_convert_fluxdensity():

    wavelength = np.array([1.0,2.0])
    flux_density = np.array([1.0,1.0])

    # These first tests will ensure that the internal conversion between
    # erg s-1 cm-1 A-1 to the user desired units worked.
    
    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'W m-2 um-1')
    np.testing.assert_array_equal(result,np.array([1.0,1.0]))

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'W m-2 Hz-1')
    np.testing.assert_allclose(result,np.array([3.335640951981e-15,
                                                1.334256380792e-14]))

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'erg s-1 cm-2 A-1')
    np.testing.assert_allclose(result,np.array([0.1,0.1]))

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'erg s-1 cm-2 Hz-1')
    np.testing.assert_allclose(result,np.array([3.335640951981e-12,
                                                1.334256380792e-11]))

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'Jy')
    np.testing.assert_allclose(result,np.array([333564095198.15210,
                                                1334256380792.6084]))

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'mJy')
    np.testing.assert_allclose(result,np.array([333564095198.15210,
                                                1334256380792.6084])*1e3)

    result = convert_fluxdensity(wavelength, flux_density, 'um', 'W m-2 um-1',
                                 'uJy')
    np.testing.assert_allclose(result,np.array([333564095198.15210,
                                                1334256380792.6084])*1e6)                                
    # These second tests ensure that the intial conversion to ergs s-1 cm-1 A-1
    # works properly

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'W m-2 um-1')
    np.testing.assert_allclose(result,np.array([10.0,10.0]))

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'W m-2 Hz-1')
    np.testing.assert_allclose(result,np.array([3.335640951981e-14,
                                                1.334256380792e-13]))

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'erg s-1 cm-2 A-1')
    np.testing.assert_array_equal(result,np.array([1.0,1.0]))        

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'erg s-1 cm-2 Hz-1')
    np.testing.assert_allclose(result,np.array([3.335640951981e-11,
                                                1.334256380792e-10]))

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'Jy')
    np.testing.assert_allclose(result,np.array([3335640951981.5205,
                                                13342563807926.082]))

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'mJy')
    np.testing.assert_allclose(result,np.array([3335640951981.5205,
                                                13342563807926.082])*1e3)

    result = convert_fluxdensity(wavelength, flux_density, 'um',
                                 'erg s-1 cm-2 A-1', 'uJy')
    np.testing.assert_allclose(result,np.array([3335640951981.5205,
                                                13342563807926.082])*1e6)                

    

    
