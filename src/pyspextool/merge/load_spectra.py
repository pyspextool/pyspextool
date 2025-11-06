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

def load_spectra(
    file: str,
    outputfile_root:str,
    verbose:bool = None):

    """

    Parameters
    ----------
    file : str
        The name of the pySpextool FITS file containing the multi-order spectra.

    outputfile_root : str
        A string giving the root of the output file name, e.g. 'Wolf359'. 

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

    check_parameter("load_spectra", "file", file, "str")

    check_parameter("load_spectra", "outputfile_root", outputfile_root, "str")

    check_parameter("load_spectra", "verbose", verbose, ["NoneType", "bool"])

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
    logging.info(" Loading the spectra.")

    fullpath = make_full_path(setup.state["proc_path"], 
                              config.state["file"], exist=True)

    spectra, data = read_spectra_fits(fullpath)

    spectra = spectra.astype(np.float64)


#    print(data.keys())
    hdrinfo = get_headerinfo(data["astropyheader"])
    # commented this next line out as there are no "telluric_keywords" 
#                             keywords=setup.state["telluric_keywords"])

    #
    # Store the results
    #
    config.state["rawspectra"] = spectra
    config.state["data"] = data
    config.state["hdrinfo"] = hdrinfo
    config.state["norders"] = data["norders"]
    config.state["napertures"] = data["napertures"]
    config.state["orders"] = data["orders"]
    config.state["xlabel"] = hdrinfo["LXLABEL"][0]

    if data["norders"] == 1:

        message = ' The file has only 1 order.  No merging required.'
        raise pySpextoolError(message)

    logging.info(" There are "+str(data["napertures"])+\
                 " apertures in this file.")

    #
    # Set done variables
    #

    config.state["load_done"] = True
    config.state["merge_done"] = False

