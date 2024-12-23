import numpy as np
from pyspextool.utils.coords import ten, sixty
import math 
import pytest
from astropy.coordinates import SkyCoord
import astropy.units as u

@pytest.mark.parametrize('value,expected', 
    [
        ('-00:00:40.04424',-0.0111234), # ten('-00:00:40.04424')
        ([-0.0,0.0,40.04424],-0.0111234), # ten([-0.0,0.0,40.04424])
        ([0.0,0.0,-40.04424],-0.0111234), # ten([0.0,0.0,-40.04424])
        (np.array([-0.0,0.0,40.04424]), -0.0111234), # ten(np.array([-0.0,0.0,40.04424]))
        (np.array([0.0,0.0,-40.04424]), -0.0111234)  #ten(np.array([0.0,0.0,-40.04424]))
    ])
def test_ten(value, expected):
    result = ten(value)
    assert math.isclose(result,expected)

    c = SkyCoord(ra=1, dec=value, unit=('hourangle', 'deg'))
    assert math.isclose(c.dec.value,expected), f"dec provided: {value}, dec returned: {c.dec}"
    assert c.dec.unit.is_equivalent(u.degree)

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
    
