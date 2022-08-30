import pytest
from pyspextool.utils.math import *


def test_round():
    result = round(np.array([-3.5,-2.5,2.5,3.5]))
    assert result[0] == -4.0
    assert result[1] == -3.0
    assert result[2] == 3.0
    assert result[3] == 4.0
    