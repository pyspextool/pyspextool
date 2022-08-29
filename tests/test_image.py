#################
## contains tests for image.py
###############
import numpy as np
from pyspextool.fit.image import *
import pytest

def test_imgpoly():
	img= np.array([[1,1,1],[2,2,2],[3,3,3]])
	coeffs0=np.array([[2,2,2],[2,2,2],[2,2,2]])
	coeffs1= np.array([[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]])
	coeffs= np.stack((coeffs0,coeffs1))

	res=imgpoly(img,coeffs)

	assert (res[0]== np.array([2.5, 2.5, 2.5])).all()
	assert (res[1]== np.array([3.,  3. , 3. ])).all()
	assert (res[-1]== np.array([3.5, 3.5, 3.5])).all()