import numpy as np
import pyspextool as ps
from pyspextool.utils import interpolate

def test_interpolate():

    #
    # Test for _noxnan_linearinterp1d
    # 

    input_x  = np.array([4.5, 5.5, 6.5, 7])
    input_y  = np.array([1. , 1  , 1  , 1])
    input_u  = np.array([2. , 2  , 2  , 2])

    output_x  = np.array([1,2,3,4,5,6])

    result = interpolate._noxnan_linearinterp1d(input_x,
                                                input_y,
                                                output_x,
                                                input_u = input_u)

    np.testing.assert_array_equal(result['y'],
                                  np.array([np.nan,np.nan,np.nan,np.nan,1,1]))

    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([np.nan,np.nan,np.nan,np.nan,
                                            2.44948974,2.44948974]))

    #
    # Test for linear_interp1d
    # 

    
#    result = interpolate.linear_interp1d(input_x,
#                                         input_y,
#                                         output_x,
#                                         input_u = input_u)
#    
#
#    np.testing.assert_array_equal(result[0],
#                                  np.array([np.nan,np.nan,np.nan,np.nan,1,1]))
#
#    np.testing.assert_allclose(result[1],
#                                  np.array([np.nan,np.nan,np.nan,np.nan,
#                                            2.44948974,2.44948974]))

