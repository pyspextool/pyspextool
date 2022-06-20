##
# contains tests for findorders
##

from pyspextool.segdeg import ten, sixty
import numpy as np

def test_ten():
	x= '-00:00:40.04424'
	assert ten(x)== -0.0111234

	x= [-0.0,0.0,40.04424]
	assert ten(x)==-0.0111234

	x= np.array([0.0,0.0,-40.04424])
	assert ten(x)==-0.0111234


def test_sixty():

	x = -0.0111234
	assert (sixty(x)==np.array([0.0, 0.0, -40.04424])).all()

	assert (sixty(x,trailsign=True)==np.array([-0.0, 0.0, 40.04424])).all()

	assert (sixty(x,colons={'dec':2,'plus':1},trailsign=True)=='-00:00:40.04')

	x = 0.0111234
	assert (sixty(x,colons={'dec':2,'plus':1},trailsign=True)=='+00:00:40.04')