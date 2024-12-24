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

def test_sixty():

    result = sixty(-0.0111234)
    np.testing.assert_array_equal(result,[0,0,-40.04424])

    result = sixty(-0.0111234,trailsign=True)
    np.testing.assert_array_equal(result,[-0,0,40.04424])

    result = sixty(-0.0111234, colons={'dec':2,'plus':1}, trailsign=True)
    assert result == '-00:00:40.04'

    result = sixty(0.0111234)
    np.testing.assert_array_equal(result,[0,0,40.04424])

    result = sixty(0.0111234,trailsign=True)
    np.testing.assert_array_equal(result,[0,0,40.04424])

    result = sixty(0.0111234, colons={'dec':3,'plus':0}, trailsign=True)
    assert result == '00:00:40.044'

    result = sixty(0.0111234, colons={'dec':2,'plus':1}, trailsign=True)
    assert result == '+00:00:40.04'        
    
