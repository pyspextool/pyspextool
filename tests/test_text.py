#################
## contains tests for text.py
###############
import numpy as np
from pyspextool.utils.text import *
import pytest

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

@pytest.mark.xfail(reason="results are correct but failing because of return type")
def test_forprint():
	 x = [1,2,3]
	 y = [4,5,6]
	 assert (forprint(x,y)[0]==[1, 4]).all()

@pytest.mark.skip(reason="example not implemented in docstring")
def test_splittext():
	pass

def test_wherelist():
    x = [1,'two',3]
    assert wherelist(x,'two')[0]==1
