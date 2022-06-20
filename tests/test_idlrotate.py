###
# contains test for idl rotate function
### 
import numpy as np
from pyspextool.idlrotate import  idlrotate

def test_idlrotate_zero():
	##test idl rotation using a 2d-matrix 
	img = np.array([[1,2,3],[4,5,6],[7,8,9]])
	rotated= idlrotate(img,6)
	assert (rotated== np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])).all() 