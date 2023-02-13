"""
PySpexTool is a package intended to contain core functionality and some
common tools needed for performing reduction for SpeX instrument in
Python.
"""

from . import extract
from .setup_utils import pyspextool_setup as pyspextool_setup
from .setup_utils import set_paths as set_paths
from .setup_utils import set_instrument as set_instrument

pyspextool_setup()



