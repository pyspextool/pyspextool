##
# contains tests for mkimgidxs.py
##

from pyspextool.mkimgidxs import mkimgidxs
import numpy as np

def test_mkimgidxs():
	ximg, yimg = mkimgidxs(3,5)

	assert (ximg[0]==[0,1, 2 ,3 , 4]).all()
	assert (ximg== ximg[0]).all()

	assert (yimg[0]==0).all()
	assert (yimg[1]==1).all()
	assert (yimg[2]==2).all()