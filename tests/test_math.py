import pytest
from pyspextool.utils.math import *


def test_findidx():
    array = [1, 2.5, 3, 5.5]
    points = [1.5, 3.1, 6]
    result = findidx(array, points)
    assert result[0] == 1./3.
    assert result[1] == 2.04
    assert np.isnan(result[2])

    single_point = 1.5
    result_single = findidx(array, single_point)
    assert result_single == 1./3.

def test_findidx_errors():
    # input not monotonic, should raise ValueError
    with pytest.raises(ValueError):
        array_bad = [1, 2.5, 2, 5.5]
        test_points = [1.5, 3.1, 6]
        result_bad = findidx(array_bad, test_points)

    # inputs should be array_like
    with pytest.raises(ValueError):
        array_bad2 = ('1234', 4.0, [1, 2])
        result_bad2 = findidx(array_bad2, test_points)