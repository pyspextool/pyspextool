##
# contains tests for combflagstack
##
import numpy as np
from pyspextool.combflagstack import combflagstack

def test_combflagstack():
	
	spec1 = np.array([0,4,4,2])
	spec2 = np.array([0,0,3,1])
	stack = np.stack((spec1,spec2))

	assert (combflagstack(stack) == np.array([0, 4, 7 ,3])).all()

	img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
	img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
	stack = np.stack((img1,img2))
	assert (combflagstack(stack) == np.array([[1, 2, 0], [3, 0, 4], [0, 0, 0]])).all()
