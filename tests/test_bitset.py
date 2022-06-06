import numpy as np
from pyspextool.mc_bitset import bitset


mask = bitset(np.array([3,4,1]),[0])
assert (mask==np.array([1, 0, 1])).all()