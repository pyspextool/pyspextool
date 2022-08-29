#################
## contains tests for plot.py
###############
import numpy as np
from pyspextool.plotting.plot import *
import pytest


def test_bufrange():
	res=bufrange((1,11),frac=0.1)
	assert res==(0.0, 12.0)

@pytest.mark.skip(reason="example not implemented in docstring")
def test_getyrange():
	assert  0 