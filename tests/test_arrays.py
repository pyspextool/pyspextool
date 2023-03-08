from pyspextool.utils.arrays import *
import numpy as np
import pytest


def test_find_index():
    array = [1, 2.5, 3, 5.5]
    points = [1.5, 3.1, 6]
    result = find_index(array, points, ends_to_nan=True)
    assert result[0] == 1./3.
    assert result[1] == 2.04
    assert np.isnan(result[2])

    single_point = 1.5
    result_single = find_index(array, single_point)
    assert result_single == 1./3.

    x = [1, 2.5, 3, 5.5]
    x_want = [-0.1, 1.5, 3.1, 6]
    result = find_index(x, x_want, ends_to_nan=True)
    np.testing.assert_allclose(result, [np.nan, 0.3333333333, 2.04, np.nan])


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


def test_make_image_indices():

    ximg, yimg = make_image_indices(3, 5)
    np.testing.assert_array_equal(ximg, np.array([[0, 1, 2, 3, 4],
                                                  [0, 1, 2, 3, 4],
                                                  [0, 1, 2, 3, 4]]))

    np.testing.assert_array_equal(yimg, np.array([[0, 0, 0, 0, 0],
                                                  [1, 1, 1, 1, 1],
                                                  [2, 2, 2, 2, 2]]))


def test_trim_nan():

    x = np.array([np.nan, 2, 3, 4, np.nan, 7, 89, 90, np.nan])

    result = trim_nan(x, flag=0)
    np.testing.assert_array_equal(result, [True, True, True, True,
                                           True, True, True, True, False])

    result = trim_nan(x, flag=0, trim=True)
    np.testing.assert_array_equal(result, [np.nan, 2, 3, 4, np.nan, 7, 89, 90])

    result = trim_nan(x, flag=1)
    np.testing.assert_array_equal(result, [False, True, True, True,
                                           True, True, True, True, True])

    result = trim_nan(x, flag=1, trim=True)
    np.testing.assert_array_equal(result, [2, 3, 4, np.nan, 7, 89, 90, np.nan])

    result = trim_nan(x, flag=2)
    np.testing.assert_array_equal(result, [False, True, True, True,
                                           True, True, True, True, False])

    result = trim_nan(x, flag=2, trim=True)
    np.testing.assert_array_equal(result, [2, 3, 4, np.nan, 7, 89, 90])

    result = trim_nan(x, flag=3)
    np.testing.assert_array_equal(result, [False, True, True, True,
                                           False, True, True, True, False])

    result = trim_nan(x, flag=3, trim=True)
    np.testing.assert_array_equal(result, [2, 3, 4, 7, 89, 90])


def test_idl_rotate():

    img = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    
    result = idl_rotate(img, 0)
    np.testing.assert_array_equal(result, img)

    result = idl_rotate(img, 1)
    np.testing.assert_array_equal(result,np.array([[7, 4, 1],
                                                   [8, 5, 2],
                                                   [9, 6, 3]]))

    result = idl_rotate(img, 2)
    np.testing.assert_array_equal(result, np.array([[9, 8, 7],
                                                    [6, 5, 4],
                                                    [3, 2, 1]]))

    result = idl_rotate(img, 3)
    np.testing.assert_array_equal(result, np.array([[3, 6, 9],
                                                    [2, 5, 8],
                                                    [1, 4, 7]]))

    result = idl_rotate(img, 4)
    np.testing.assert_array_equal(result, np.array([[1, 4, 7],
                                                    [2, 5, 8],
                                                    [3, 6, 9]]))

    result = idl_rotate(img, 5)
    np.testing.assert_array_equal(result, np.array([[3, 2, 1],
                                                    [6, 5, 4],
                                                    [9, 8, 7]]))

    result = idl_rotate(img, 6)
    np.testing.assert_array_equal(result, np.array([[9, 6, 3],
                                                    [8, 5, 2],
                                                    [7, 4, 1]]))

    result = idl_rotate(img, 7)
    np.testing.assert_array_equal(result, np.array([[7, 8, 9],
                                                    [4, 5, 6],
                                                    [1, 2, 3]]))


def test_idl_unrotate():

    home = np.array([[1, 2, 3],
                     [4, 5, 6],
                     [7, 8, 9]])
    
    result = idl_unrotate(home, 0)
    np.testing.assert_array_equal(result, home)
    
    img = np.array([[7, 4, 1],
                    [8, 5, 2],
                    [9, 6, 3]])
    result = idl_unrotate(img, 1)
    np.testing.assert_array_equal(result, home)

    img = np.array([[9, 8, 7],
                    [6, 5, 4],
                    [3, 2, 1]])
    result = idl_unrotate(img, 2)
    np.testing.assert_array_equal(result, home)

    img = np.array([[3, 6, 9],
                    [2, 5, 8],
                    [1, 4, 7]])
    result = idl_unrotate(img, 3)
    np.testing.assert_array_equal(result, home)

    img = np.array([[1, 4, 7],
                    [2, 5, 8],
                    [3, 6, 9]])
    result = idl_unrotate(img, 4)
    np.testing.assert_array_equal(result, home)

    img = np.array([[3, 2, 1],
                    [6, 5, 4],
                    [9, 8, 7]])
    result = idl_unrotate(img, 5)
    np.testing.assert_array_equal(result, home)

    img = np.array([[9, 6, 3],
                    [8, 5, 2],
                    [7, 4, 1]])
    result = idl_unrotate(img, 6)
    np.testing.assert_array_equal(result, home)

    img = np.array([[7, 8, 9],
                    [4, 5, 6],
                    [1, 2, 3]])
    result = idl_unrotate(img, 7)
    np.testing.assert_array_equal(result, home)

    
