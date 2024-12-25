import numpy as np
from pyspextool.utils.irplanck import irplanck
from astropy.modeling.models import BlackBody
import astropy.units as u

def test_irplank():
    # 1.2 microns, 6000 K
    value = irplanck(1.2, 6000)

    # 1.2 microns, 6000 K, RETURS 7506318  W m-2 um-1 sr-1
    np.testing.assert_allclose(value, 7506318.497767141, 0.001)

    wavelength = 1.2 * u.micron
    temperature = 6000 * u.K
    scale = 1 * u.W / (u.m**2 * u.um * u.sr)

    bb = BlackBody(temperature=temperature, scale=scale)
    bb_result = bb(wavelength)
    np.testing.assert_allclose(bb_result.value, 7506318., 0.001)