##
# contains tests for mkimgidxs.py
##

from pyspextool.polyfit1d import polyfit1d
import numpy as np

def test_polyfit1d():
	x=  np.array([10. ,20. ,30. ,40. ,50., 60. ,70. ,80. ,90.,100.])
	y=  np.array([0.37,0.58,0.83,1.15,1.36,1.62,1.90,2.18,2.45,np.nan])
	yerr = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
	res=  polyfit1d(x,y,1,yunc=yerr)
	pass
	#will implement test later