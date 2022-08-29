#################
## contains tests for image.py
###############
import numpy as np
from pyspextool.utils.image import *
import pytest

def test_scalestack():
	spec1 = np.array([1,1,1,1])
	spec2 = np.array([2,2,2,2])
	spec3 = np.array([4,4,4,4])
	specstack = np.stack((spec1,spec2,spec3))
	scaledstack, var, scales = scalestack(specstack,None)
	assert (scaledstack == np.array([[2., 2. ,2., 2.], [2. , 2. , 2.,  2.], [2., 2. ,2. , 2.]])).all()
	assert (scales == [2.,  1.,  0.5]).all()

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

