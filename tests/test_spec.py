#################
## contains tests for spec.py
###############
import numpy as np
from pyspextool.utils.spec import *
import pytest


def test_combflagstack():
	spec1 = np.array([0,4,4,2])
	spec2 = np.array([0,0,3,1])
	stack = np.stack((spec1,spec2))

	assert (combflagstack(stack) == np.array([0, 4, 7 ,3])).all()

	img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
	img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
	stack = np.stack((img1,img2))
	assert (combflagstack(stack) == np.array([[1, 2, 0], [3, 0, 4], [0, 0, 0]])).all()

@pytest.mark.skip(reason="example not implemented in docstring")
def test_extr_xs1dxd():
	assert  0 

@pytest.mark.skip(reason="example not implemented in docstring")
def test_findorders():
	assert  0 

@pytest.mark.skip(reason="example not implemented in docstring")
def test_find_top_bot():
	assert  0 

@pytest.mark.skip(reason="example not implemented in docstring")
def test_getspecpixshift():
	assert  0 

def test_mkapmask():
	slit_arc = np.arange(100)
	res=mkapmask(slit_arc,[50,60],[4,4.1],psbginfo=[11,20])
	assert (res == np.array(['0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.\
        0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.\
       -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.\
        0.   0.   0.   0.   0.5  1.   1.   1.   1.   1.   1.   1.   0.5  0.\
        1.6  2.   2.   2.   2.   2.   2.   2.   1.6  0.   0.   0.   0.   0.\
        0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.\
       -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.   0.   0.   0.   0.\
        0.   0.'.split()]).astype(float)).all()

@pytest.mark.skip(reason="example not implemented in docstring")
def test_simwavecal1dxd():
	assert  0 

