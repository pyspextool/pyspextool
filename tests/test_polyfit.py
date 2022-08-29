#################
## contains tests for polyfit.py
###############
import numpy as np
from pyspextool.fit.polyfit import *
import pytest

@pytest.mark.skip(reason="example not implemented in docstring")
def test_poly2d():
	pass

def test_polyfit1d():
	x=  np.array([10. ,20. ,30. ,40. ,50., 60. ,70. ,80. ,90.,100.])
	y=  np.array([0.37,0.58,0.83,1.15,1.36,1.62,1.90,2.18,2.45,np.nan])
	yerr = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
	res=  polyfit1d(x,y,1,yunc=yerr)
	#test coefficients
	assert np.isclose(res['coeffs'][1],  0.026,  rtol=1e-2)
	assert np.isclose(res['coeffs'][0],  0.071,  rtol=1e-2)
	#test other stats?

@pytest.mark.skip(reason="example not implemented in docstring")
def test_polyfit2d():
	pass

@pytest.mark.skip(reason="example not implemented in docstring")
def test_robustpolyfit1d():
	pass
