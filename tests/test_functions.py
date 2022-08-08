#################
## contains tests for utils in pyspextool
###############

import numpy as np
from pyspextool.calibration import *
from pyspextool.fit import *
from pyspextool.io import *
from pyspextool.plotting import *
from pyspextool.uspex import *
from pyspextool.utils import *

import pytest
xfail = pytest.mark.xfail #taking care of tests that are expected to fial

@xfail(reason='examples not implemented')
def test_avehdrs():
	pass

@xfail(reason='examples not implemented')
def test_bicucof():
	pass

@xfail(reason='examples not implemented')
def test_bicuval():
	pass

@xfail(reason='examples not implemented')
def test_findorders():
	pass

def test_biset_mask():
	#tests if bitset mask is working correctly
	mask = bitset(np.array([3,4,1]),[0])
	assert (mask==np.array([1, 0, 1])).all()

def test_combflagstack():
	spec1 = np.array([0,4,4,2])
	spec2 = np.array([0,0,3,1])
	stack = np.stack((spec1,spec2))

	assert (combflagstack(stack) == np.array([0, 4, 7 ,3])).all()

	img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
	img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
	stack = np.stack((img1,img2))
	assert (combflagstack(stack) == np.array([[1, 2, 0], [3, 0, 4], [0, 0, 0]])).all()

def test_dictaddentry():
    dicta = {'HA':1,'PA':2,'MJD': 3}
    res=dictaddentry(dicta,'MJD','before','new',4)
    assert res['HA']== 1
    assert res['PA']== 2
    assert res['new']== 4
    assert res['MJD']== 3

    dicta = {'HA':1,'PA':2,'MJD':3}
    res=dictaddentry(dicta,'HA','after','new',(3,4))
    assert res['HA']== 1
    assert res['new']== (3, 4)
    assert res['PA']== 2
    assert res['MJD'] ==3

def test_fsextract():
	assert (fsextract('1-3,5,7,10-12','index')==np.array([1, 2, 3, 5, 7, 10, 11, 12])).all()
	res=  fsextract('spc00001.a.fits,spc00002.a.fits','filename')
	assert res[0]=='spc00001.a.fits'
	assert res[1]=='spc00002.a.fits'
	assert len(res)==2

def test_idlrotate():
	##test idl rotation using a 2d-matrix 
	img = np.array([[1,2,3],[4,5,6],[7,8,9]])
	rotated= idlrotate(img,6)
	assert (rotated== np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])).all() 

def test_imgpoly():
	img= np.array([[1,1,1],[2,2,2],[3,3,3]])
	coeffs0=np.array([[2,2,2],[2,2,2],[2,2,2]])
	coeffs1= np.array([[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]])
	coeffs= np.stack((coeffs0,coeffs1))

	res=imgpoly(img,coeffs)

	assert (res[0]== np.array([2.5, 2.5, 2.5])).all()
	assert (res[1]== np.array([3.,  3. , 3. ])).all()
	assert (res[-1]== np.array([3.5, 3.5, 3.5])).all()

def test_medstack():
	ss= np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
	msk= np.ones((3,4),dtype=int)
	med,unc=medcomb(ss,mask=msk)

	assert (med==np.array([1., 5., 4., 6.])).all()
	assert len(unc)==len(med)

def test_mkimgidxs():
	ximg, yimg = mkimgidxs(3,5)

	assert (ximg[0]==[0,1, 2 ,3 , 4]).all()
	assert (ximg== ximg[0]).all()

	assert (yimg[0]==0).all()
	assert (yimg[1]==1).all()
	assert (yimg[2]==2).all()

def test_ten():
	x= '-00:00:40.04424'
	assert ten(x)== -0.0111234

	x= [-0.0,0.0,40.04424]
	assert ten(x)==-0.0111234

	x= np.array([0.0,0.0,-40.04424])
	assert ten(x)==-0.0111234


def test_sixty():
	x = -0.0111234
	assert (sixty(x)==np.array([0.0, 0.0, -40.04424])).all()

	assert (sixty(x,trailsign=True)==np.array([-0.0, 0.0, 40.04424])).all()

	assert (sixty(x,colons={'dec':2,'plus':1},trailsign=True)=='-00:00:40.04')

	x = 0.0111234
	assert (sixty(x,colons={'dec':2,'plus':1},trailsign=True)=='+00:00:40.04')


def test_polyfit1d():
	x=  np.array([10. ,20. ,30. ,40. ,50., 60. ,70. ,80. ,90.,100.])
	y=  np.array([0.37,0.58,0.83,1.15,1.36,1.62,1.90,2.18,2.45,np.nan])
	yerr = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
	res=  polyfit1d(x,y,1,yunc=yerr)
	pass
	#will implement test later