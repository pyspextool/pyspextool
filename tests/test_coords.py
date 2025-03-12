import numpy as np
from pyspextool.utils.coords import ten, sixty
import math 
import pytest
from astropy.coordinates import Angle
import astropy.units as u

@pytest.mark.parametrize('input_ten, input_angle,expected', 
    [
        ('-00:00:40.04424', '-00:00:40.04424',-0.0111234), # ten('-00:00:40.04424')
        ([-0.0,0.0,41.04424], '-0 0 41.04424',-0.011401178), # ten([-0.0,0.0,40.04424])
        ([0.0,0.0,-42.04424], '-42.04424s', -0.011678955556), # ten([0.0,0.0,-40.04424])
        (np.array([-0.0,0.0,43.04424]), '-0 0 43.04424', -0.0119567334), # ten(np.array([-0.0,0.0,40.04424]))
        (np.array([0.0,0.0,-44.04424]), "-44.04424″", -0.012234511)  #ten(np.array([0.0,0.0,-40.04424]))
    ])
def test_ten(input_ten, input_angle,expected):
    result = ten(input_ten)
    assert np.isclose(result,expected)
    print('ten', result)

    angle = Angle(input_angle, unit=u.deg)
    print('angle', angle)    
    assert np.isclose(angle.value, expected), f"angle provided: {input_angle}, angle returned: {angle}"
    assert angle.unit.is_equivalent(u.degree)

@pytest.mark.parametrize('input_sixty, output_sixty,plus, trailsign,expected',
    [
        (-0.0111234, [0,0,-40.04424], 0, False, '-00:00:40.04'),
        (-0.0111234, [-0,0,40.04424], 1, True, '-00:00:40.04'),
        (0.0111234, [0,0,40.04424], 0, False,'00:00:40.04'),
        (0.0111234, [0,0,40.04424], 1, True,'+00:00:40.04')

    ])
def test_sixty(input_sixty, output_sixty, plus, trailsign, expected):

    result = sixty(input_sixty, trailsign=trailsign)
    np.testing.assert_array_equal(result,output_sixty)

    result = sixty(input_sixty, colons={'dec':2,'plus':plus}, trailsign=trailsign)
    assert result == expected

    angle = Angle(input_sixty, unit=u.deg)
    angle_string = angle.to_string(unit=u.degree, sep=':', alwayssign=trailsign,pad=True, precision=2)      
    assert angle_string == expected
