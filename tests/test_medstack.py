##
# contains tests for medstack.py
##

from pyspextool.medstack import medcomb
import numpy as np

def test_medstack():
	ss= np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
	msk= np.ones((3,4),dtype=int)
	med,unc=medcomb(ss,mask=msk)

	assert (med==np.array([1., 5., 4., 6.])).all()
	assert len(unc)==len(med)