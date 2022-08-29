#################
## contains tests for ioutils.py
###############
import numpy as np
from pyspextool.utils.ioutils import *
import pytest

def test_biset():
	#tests if bitset mask is working correctly
	mask = bitset(np.array([3,4,1]),[0])
	assert (mask==np.array([1, 0, 1])).all()

@pytest.mark.skip(reason="will implement later once we figure out path navigation")
def test_check_dir():
	pass

@pytest.mark.skip(reason="will implement later once we figure out path navigation")
def test_check_file():
	pass

def test_fsextract():
	assert fsextract('1-3,5,7,10-12','index') ==[1, 2, 3, 5, 7, 10, 11, 12]
	assert fsextract('spc00001.a.fits,spc00002.a.fits','filename') == ['spc00001.a.fits', 'spc00002.a.fits']

@pytest.mark.xfail(reason="will implement later once we figure out path navigation")
def test_mkfullpath():
	files = '1-5'
	dirct = '../../uSpeXdata/raw/'
	res=mkfullpath(dirct,files,indexinfo={'nint':5,'prefix':'spc-','suffix':'.[ab].fits'},exist=True)
	assert res[0]== '../../uSpeXdata/raw/spc-00001.a.fits'
	assert res[1]=='../../uSpeXdata/raw/spc-00002.b.fits'
	assert res[-1]=='../../uSpeXdata/raw/spc-00005.a.fits'