import numpy as np
from unittest import TestCase
from pyspextool.utils.irplanck import irplanck


value = irplanck(1.2,6000)

np.testing.assert_allclose(value,7506318.497767141,0.001)
