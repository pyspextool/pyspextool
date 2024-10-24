import numpy as np
from pyspextool.plot.limits import *

def test_buffer_range():

    range = (1.5,2.5)
    result = buffer_range(range,0.1)
    assert result == (1.4, 2.6)


def test_get_image_range():

    assert 1 == 1
    

    
def test_get_spectra_range():

    spec1 = np.array([1,2,3])
    spec2 = np.array([4,5,6])

    range = get_spectra_range(spec1,spec2)
    assert range[0] == 1.0
    assert range[1] == 6.0

    range = get_spectra_range(spec1,spec2,frac=0.1)
    assert range[0] == 0.5
    assert range[1] == 6.5


    
