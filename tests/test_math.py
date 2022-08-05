import numpy as np
from pyspextool.utils.math import *


def test_findidx():
    array = [1, 2.5, 3, 5.5]
    test_points = [1.5, 3.1, 6]
    result = findidx(array, test_points)
    assert result[0] == 1./3.
    assert result[1] == 2.04
    assert np.isnan(result[2])
