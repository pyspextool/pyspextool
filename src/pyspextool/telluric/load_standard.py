import os
import numpy as np
import logging
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.telluric import config

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import make_full_path
from pyspextool.io.load_atmosphere import load_atmosphere
from pyspextool.io.query_simbad import query_simbad
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io import read
from pyspextool.io.sptype2teff import sptype2teff
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.setup_utils import mishu # pooch class for caching files


logger = logging.getLogger(__name__)

def load_standard(
    standard_fullfilename: str,
    standard_info: list | str | dict,
    output_filename: str,
    correction_type: list = "A0 V",
    verbose: bool = None):

    """
    Loads standard spectra and ancillary things.

    Parameters
    ----------
    standard_fullfilename : str
        The full file name (including the suffix '.fits') of the pySpextool 
        FITS file containing the standard spectra.

    standard_info : str, list, dict

        If a string is passed, it is assumed to be the name of the standard.
        SIMBAD is queried for the spectral type and B- and V-band magnitudes.

        If a list is passed, it is assumed to contain the coordinates of the
        standard.  standard_info[0] = RA, standard_info[1] = Dec.  Each element
        of the list is a string which gives the sexigesimal coordidates,
        separated by a either a ' ', :, or h/m/s, d/m/s.  SIMBAD is then
        queried for the id, spectral type, and B- and V-band magnitudes.

        If a dict is passed, it is assumed to contain a dictionary with the
        standard star information.

        `"id"` : str
            The name of the standard.

        `"sptype"` : str
            The standard star spectral type.

        `"bmag"` : int or float
            The standard star B-band magnitude

        `"vmag"` : int or float
            The standard star V-band magnitude.

    output_filename : str
        A string giving the root of the output file name, e.g. 'Wolf359'.

    correction_type : {'A0 V', 'reflectence', 'basic'}
        The type of telluric correction requested.
        type = 'A0 V' -> standard correction using an A0 V standard star
        type = 'reflectence' -> ratio the object and standard.
        type = 'basic' -> ratio the object and standard and multiply by a
               Planck function of the appropriate temperature.

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].

    Returns
    -------
    None
    Loads data into memory:

        config.state['object_file']
        config.state['standard_fullfilename']
        config.state['standard_info']
        config.state['output_filename']
        config.state['reflectance']
        config.state['ew_scale']
        config.state['standard_rv']
        config.state['standard_redshift']
        config.state['1+z']
        config.state['rms_deviation']
        config.state['max_deviation']
        config.state['load_done']
        config.state['prepare_done']
        config.state['rv_done']
        config.state['kernel_done']

    """

    #
    # Check the parameters and keywords
    #

    check_parameter("load_standard", "standard_fullfilename", 
                    standard_fullfilename, "str")

    check_parameter("load_standard", "standard_info", 
                    standard_info,  ["str", "list", "dict"])

    check_parameter("load_standard", "output_filename", 
                    output_filename, "str")

    check_parameter("load_standard", "correction_type", 
                    correction_type, "str",
                    possible_values=setup.state["telluric_correctiontypes"])

    check_parameter("load_standard", "verbose", 
                    verbose, ["NoneType", "bool"])

    check_qakeywords(verbose=verbose)

    #
    # Log the process
    #
    
    logging.info(" Telluric Correction\n--------------------------\n")

    #
    # Start the output dictionary
    #

    output = {}

    output["standard_fullfilename"] = standard_fullfilename
    output["standard_info"] = standard_info
    output["telluric_output_filename"] = output_filename
    output["correction_type"] = correction_type.replace(" ", "")

    #
    # Load and store the standard file.
    #

    logging.info(" Loading the standard spectrum.")

    fullpath = make_full_path(
        setup.state["proc_path"], 
        output["standard_fullfilename"], 
        exist=True)

    spectra, dict = read_spectra_fits(fullpath)

    # Ensure the standard has only one aperture.

    if dict["napertures"] != 1:

        message = "The standard must have one aperture."
        raise pySpextoolError(message)

    # Store results and useful information

    output['standard_spectra'] = spectra
    output['standard_dictionary'] = dict
    output['standard_astropyheader'] = dict['astropyheader']
    output["slitw_pix"] = dict["slitw_pix"]
    output["resolving_power"] = dict["resolving_power"]
    output["standard_orders"] = dict["orders"]
    output["latex_xlabel"] = dict["lxlabel"]
    output["instrument_mode"] = dict["obsmode"]
    output["standard_norders"] = dict["norders"]
    
    # Compute the minimum and maximum wavelengths for each standard order 
    #(This may go away eventually once you implement mikeline_convolve).
    
    wavelength_ranges = []

    for i in range(dict["norders"]):

        # Do min max first

        min = np.nanmin(spectra[i, 0, :])
        max = np.nanmax(spectra[i, 0, :])

        wavelength_ranges.append(np.array([min, max]))

    # Store the results

    output["standard_wavelengthranges"] = wavelength_ranges
    
    #
    # Load and store the standard's B mag, V mag, and effective temperature.
    #

    result = query_simbad(standard_info)

    output['standard_name'] = result['name']
    output['standard_sptype'] = result['sptype']
    output['standard_vmag'] = result['vmag']
    output['standard_bmag'] = result['bmag']
    output['standard_teff'] = sptype2teff(result['sptype'])

    # Set default things

    output['standard_rv'] = 0.0
    output['vega_pixelshift'] = 0.0

    #
    # Load and store the hydrogen lines
    #
    
    result = read.read_hlines_file()

    output["H_wavelengths"] = result[0]
    output["H_latex_ids"] = result[1]
    output["H_ascii_ids"] = result[1]
    
    #
    # Load the IP coefficients for the specific slit
    #

    file = os.path.join(setup.state["instrument_path"], "IP_coefficients.dat")

    # Read the file

    slitw_arc, c0, c1, c2 = np.loadtxt(file, comments="#", unpack=True)

    # Find the right slit width and store the coefficients

    z = np.where(slitw_arc == dict["slitw_arc"])[0]

    ip_coefficients = np.squeeze(np.array((c0[z], c1[z], c2[z])))

    output["ip_coefficients"] = ip_coefficients

    config.state["ip_coefficients"] = ip_coefficients

    #
    # Load and store the atmospheric transmission
    #

    atmospheric_transmission = np.copy(spectra)
    wavelengths, transmission = load_atmosphere(dict['resolving_power'])

    for i in range(dict['norders']):

        atmosphere = linear_interp1d(wavelengths, transmission, spectra[i, 0, :])
        atmospheric_transmission[i, 1, :] = atmosphere

    output["atmospheric_transmission"] = atmospheric_transmission

    #
    # Get and store the mode telluric information
    #

    fullpath = os.path.join(setup.state["instrument_path"], 
        output["instrument_mode"]+"_telluric.dat")

    result = read.read_telluric_file(fullpath)
    keys = list(result.keys())
    vals = list(result.values())


    for i in range(len(keys)):

        output[keys[i]] = vals[i]

    #
    # Get and store the Vega model 
    #

    root = "Vega" + output["telluric_vegamodel"] + ".fits"
        
    fullpath = mishu.fetch(root)

    hdul = fits.open(fullpath)
    data = hdul[1].data

    hdul.close()

    # Store the results

    output["vega_wavelength"] = data["wavelength"]
    output["vega_fluxdensity"] = data["flux density"]
    output["vega_continuum"] = data["continuum flux density"]
    output["vega_fitted_continuum"] = data["fitted continuum flux density"]

    normalized = data["flux density"] / data["continuum flux density"]

    output["vega_normalized_fluxdensity"] = normalized

    #
    # Clear the state variable and store user inputs
    #

    config.state.clear()

    keys = list(output.keys())
    vals = list(output.values())

    for i in range(len(keys)):

        config.state[keys[i]] = vals[i]

    #
    # Set the done variables
    #

    config.state["load_done"] = True
    config.state["normalize_done"] = False
    config.state["rv_done"] = False
    config.state["kernel_done"] = False

    return output

