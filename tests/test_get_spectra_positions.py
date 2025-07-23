import numpy as np
from pyspextool.merge.core import get_spectra_position


def test_get_spectra_position():


    anchor = np.array([5,6,7,8,9,10])
    add = np.array([1,2,3])

    result = get_spectra_position(anchor,add)
    assert result[0] == 'left'
    assert result[1] == [3,5]

    result = get_spectra_position(add, anchor)
    assert result[0] == 'right'
    assert result[1] == [3,5]

    anchor = np.array([5,6,7,8,9,10])
    add = np.array([1,2,3,4,5.5])

    result = get_spectra_position(anchor,add)
    assert result[0] == 'left+overlap'
    assert result[1] == [5,5.5]

    result = get_spectra_position(add, anchor)
    assert result[0] == 'right+overlap'
    assert result[1] == [5,5.5]

    anchor = np.array([5,6,7,8,9,10])
    add = np.array([5.5,6,8])

    result = get_spectra_position(anchor,add)
    assert result[0] == 'inside'
    assert result[1] == [5.5,8]

    result = get_spectra_position(add, anchor)
    assert result[0] == 'encompass'
    assert result[1] == [5.5,8]

