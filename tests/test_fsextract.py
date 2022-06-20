##
# contains tests for fsextract.py
##

from pyspextool.fsextract import fsextract
import numpy as np


def test_fsextract():
	assert (fsextract('1-3,5,7,10-12','index')==np.array([1, 2, 3, 5, 7, 10, 11, 12])).all()

	res=  fsextract('spc00001.a.fits,spc00002.a.fits','filename')
	assert res[0]=='spc00001.a.fits'
	assert res[1]=='spc00002.a.fits'
	assert len(res)==2
