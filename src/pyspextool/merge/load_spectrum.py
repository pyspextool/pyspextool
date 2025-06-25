import os
import numpy as np
import logging

from pyspextool import config as setup
from pyspextool.merge import config as config
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import make_full_path
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.load_atmosphere import load_atmosphere
from pyspextool.utils.interpolate import linear_interp1d

def load_spectrum(
    file: str,
    outputfile_root:str,
    merge_aperture:int=1,
    verbose:bool = None):

    """

    Parameters
    ----------
    file : str
        The name of the pySpextool FITS file containing the multi-order spectra.

    outputfile_root : str
        A string giving the root of the output file name, e.g. 'Wolf359'. 

    merge_aperture : int, default 1
        The aperture number to merge.

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].

    Returns
    -------
    None
    Loads data into memory


    """

    #
    # Check the parameters and keywords
    #

    check_parameter("load_spectrum", "file", file, "str")

    check_parameter("load_spectrum", "outputfile_root", outputfile_root, "str")

    check_parameter("load_spectrum", "merge_aperture", merge_aperture, "int")

    check_parameter("load_spectrum", "verbose", verbose, ["NoneType", "bool"])

    check_qakeywords(verbose=verbose)

    #
    # Clear the state variables
    #

    config.state.clear()

    #
    # Store user inputs
    #

    config.state["file"] = file
    config.state["outputfile_root"] = outputfile_root

    logging.info(" Order Merging\n--------------------\n")
    logging.info(" Loading the spectrum.")

    fullpath = make_full_path(setup.state["proc_path"], 
                              config.state["file"], exist=True)

    spectra, data = read_spectra_fits(fullpath)

    spectra = spectra.astype(np.float64)


    hdrinfo = get_headerinfo(data["header"], 
                             keywords=setup.state["telluric_keywords"])

    #
    # Store the results
    #
    config.state["spectra"] = spectra
    config.state["data"] = data
    config.state["hdrinfo"] = hdrinfo
    config.state["norders"] = data["norders"]
    config.state["napertures"] = data["napertures"]
    config.state["orders"] = data["orders"]
    config.state["xlabel"] = hdrinfo["LXLABEL"][0]
    config.state["merge_aperture"] = merge_aperture

    #
    # Set done variables
    #

    config.state["load_done"] = True
    config.state["merge_done"] = False

#

#    #
#    # Now get the object ranges
#    #
#
#    wavelength_ranges = []
#
#    for i in range(tc.state["object_norders"]):
#
#        idx = i * tc.state["object_napertures"]
#        min = np.nanmin(tc.state["object_spectra"][idx, 0, :])
#        max = np.nanmax(tc.state["object_spectra"][idx, 0, :])
#
#        wavelength_ranges.append(np.array([min, max]))
#
#    # Store the results
#
#    tc.state["object_wavelengthranges"] = wavelength_ranges

