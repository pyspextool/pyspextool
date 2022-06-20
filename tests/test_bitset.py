import numpy as np
from pyspextool.bitset import bitset

def test_biset_mask():
	#tests if bitset mask is working correctly
	mask = bitset(np.array([3,4,1]),[0])
	assert (mask==np.array([1, 0, 1])).all()