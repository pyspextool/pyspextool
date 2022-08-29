#################
## contains tests for fitsutils.py
###############
import numpy as np
from pyspextool.utils.fitsutils import *
import pytest

@pytest.mark.skip(reason="example not implemented in docstring")
def test_avehdrs():
	assert  0 

@pytest.mark.skip(reason="example not implemented in docstring")
def test_gethdrinfo():
	assert  0 

@pytest.mark.skip(reason="example not implemented in docstring")
def test_getimgrange():
	assert  0 

def test_idlrotate():
	img = np.array([[1,2,3],[4,5,6],[7,8,9]])
	rotated= idlrotate(img,6)
	assert (rotated== np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])).all() 

def test_medcomb():
	ss= np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
	msk= np.ones((3,4),dtype=int)
	med,unc=medcomb(ss,mask=msk)

	assert (med==np.array([1., 5., 4., 6.])).all()
	assert len(unc)==len(med)