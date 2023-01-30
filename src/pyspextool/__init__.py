"""
PySpexTool is a package intended to contain core functionality and some
common tools needed for performing reduction for SpeX instrument in
Python.
"""

import os
import sys
from .cl.setup_utils import setup
from .cl.make_flat import make_flat
from .cl.make_wavecal import make_wavecal
from .cl.load_image import load_image
from .cl.set_extraction_type import set_extraction_type
from .cl.make_spatial_profiles import make_spatial_profiles
from .cl.locate_aperture_positions import locate_aperture_positions
from .cl.select_orders import select_orders
from .cl.trace_apertures import trace_apertures
from .cl.define_aperture_parameters import define_aperture_parameters
from .cl.extract_apertures import extract_apertures

print('run the __init__ in pyspextool/')
