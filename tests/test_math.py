import pytest
import numpy as np
from pyspextool.utils.math import *


def test_find_index():
    array = [1, 2.5, 3, 5.5]
    points = [1.5, 3.1, 6]
    result = find_index(array, points)
    assert result[0] == 1./3.
    assert result[1] == 2.04
    assert np.isnan(result[2])

    single_point = 1.5
    result_single = find_index(array, single_point)
    assert result_single == 1./3.


def test_find_index_errors():
    # input not monotonic
    with pytest.raises(ValueError):
        array_bad = [1, 2.5, 2, 5.5]
        test_points = [1.5, 3.1, 6]
        result_bad = find_index(array_bad, test_points)

    # inputs should be array_like
    with pytest.raises(ValueError):
        array_bad2 = ('1234', 4.0, [1, 2])
        result_bad2 = find_index(array_bad2, test_points)

    with pytest.raises(ValueError):
        array = [1., 5, 10]
        test_points_bad = (1, 2.0, [1, 2])
        result_bad3 = find_index(array, test_points_bad)


def test_round():
    result = round(np.array([-3.5,-2.5,2.5,3.5]))
    assert result[0] == -4.0
    assert result[1] == -3.0
    assert result[2] == 3.0
    assert result[3] == 4.0
    