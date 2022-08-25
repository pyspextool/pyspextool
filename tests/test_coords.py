import numpy as np
from pyspextool.utils.coords import *
import math 

def test_ten():

    result = ten('-00:00:40.04424')
    assert math.isclose(result,-0.0111234)

    result = ten([-0.0,0.0,40.04424])
    assert math.isclose(result,-0.0111234)

    result = ten([0.0,0.0,-40.04424])
    assert math.isclose(result,-0.0111234)

    result = ten(np.array([-0.0,0.0,40.04424]))
    assert math.isclose(result,-0.0111234)

    result = ten(np.array([0.0,0.0,-40.04424]))
    assert math.isclose(result,-0.0111234)

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
    
