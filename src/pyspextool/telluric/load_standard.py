import os
import numpy as np
import logging
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.telluric import config

from pyspextool.fit.polyfit import polyfit_1d
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import make_full_path
from pyspextool.io.load_atmosphere import load_atmosphere
from pyspextool.io.query_simbad import query_simbad
from pyspextool.io.read_spectra_fits import read_spectra_fits
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
    verbose: bool = None,
    new=False):

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

    check_parameter("load_standard", "standard_info", standard_info, 
        ["str", "list", "dict"])

    check_parameter("load_standard", "output_filename", output_filename, "str")

    check_parameter("load_standard", "correction_type", correction_type,
        "str", possible_values=setup.state["telluric_correctiontypes"])

    check_parameter("load_standard", "verbose", verbose, ["NoneType", "bool"])

    check_qakeywords(verbose=verbose)

    #
    # Log the process
    #
    
    logging.info(" Telluric Correction\n--------------------------\n")

    #
    # Clear the state variable and store user inputs
    #

    config.state.clear()

    config.state["standard_fullfilename"] = standard_fullfilename
    config.state["standard_info"] = standard_info
    config.state["telluric_output_filename"] = output_filename
    config.state["correction_type"] = correction_type.replace(" ", "")

    #
    # Load and store the standard file.
    #

    logging.info(" Loading the standard spectrum.")

    fullpath = make_full_path(setup.state["proc_path"], 
                              config.state["standard_fullfilename"], 
                              exist=True)

    spectra, dict = read_spectra_fits(fullpath)

    # Ensure the standard has only one aperture.

    if dict["napertures"] != 1:

        message = "The standard must have one aperture."
        raise pySpextoolError(message)

    # Store results and useful information

    config.state['standard_spectra'] = spectra
    config.state['standard_dictionary'] = dict
    config.state['standard_astropyheader'] = dict['astropyheader']

    config.state["slitw_pix"] = dict["slitw_pix"]
    config.state["resolving_power"] = dict["resolving_power"]
    config.state["standard_orders"] = dict["orders"]
    config.state["latex_xlabel"] = dict["lxlabel"]
    config.state["instrument_mode"] = dict["obsmode"]
    config.state["standard_norders"] = dict["norders"]
    
    # Compute the minimum and maximum wavelengths and dispersions for each
    # standard order (This may go away eventually once you implement
    # mikeline_convolve).
    
    nwavelengths = np.shape(spectra)[-1]
    pixels = np.arange(nwavelengths)
    dispersions = np.empty(dict["norders"])
    wavelength_ranges = []

    for i in range(dict["norders"]):

        # Do min max first

        min = np.nanmin(spectra[i, 0, :])
        max = np.nanmax(spectra[i, 0, :])

        wavelength_ranges.append(np.array([min, max]))

        # Now do the dispersion

        fit = polyfit_1d(pixels, spectra[i, 0, :], 1)
        dispersions[i] = fit["coeffs"][1]

    # Store the results

    config.state["standard_wavelengthranges"] = wavelength_ranges
    config.state["standard_dispersions"] = dispersions
    config.state["standard_fwhm"] = dispersions * dict["slitw_pix"]
    
    #
    # Load and store the standard star's B mag, V mag, and effective temperature.
    #

    result = query_simbad(standard_info)

    config.state['standard_name'] = result['name']
    config.state['standard_sptype'] = result['sptype']
    config.state['standard_vmag'] = result['vmag']
    config.state['standard_bmag'] = result['bmag']
    config.state['standard_teff'] = sptype2teff(result['sptype'])
    config.state['standard_rv'] = 0.0

    #
    # Load and store the hydrogen lines
    #
    
    fullpath = os.path.join(setup.state["package_path"], "data", "HI.dat")

    wavelength, lineid = np.loadtxt(fullpath, comments="#", unpack=True, dtype="str", 
                                    delimiter="|")

    config.state["H_wavelengths"] = np.array(wavelength).astype(float)
    config.state["H_ids"] = lineid
    
    #
    # Load the IP coefficients for the specific slit
    #

    file = os.path.join(setup.state["instrument_path"], "IP_coefficients.dat")

    # Read the file

    slitw_arc, c0, c1, c2 = np.loadtxt(file, comments="#", unpack=True)

    # Find the right slit width and store the coefficients

    z = np.where(slitw_arc == dict["slitw_arc"])[0]

    ip_coefficients = np.squeeze(np.array((c0[z], c1[z], c2[z])))

    config.state["ip_coefficients"] = ip_coefficients

    #
    # Load and store the atmospheric transmission
    #

    atmospheric_transmission = np.copy(spectra)
    wavelengths, transmission = load_atmosphere(dict['resolving_power'])

    for i in range(dict['norders']):

        atmosphere = linear_interp1d(wavelengths, transmission, spectra[i, 0, :])
        atmospheric_transmission[i, 1, :] = atmosphere

    config.state["atmospheric_transmission"] = atmospheric_transmission

    #
    # Get and store the mode information
    #

    fullpath = os.path.join(setup.state["instrument_path"], "telluric_modeinfo.dat")

    result = _load_modeinfo(fullpath, dict["obsmode"])

    config.state["kernel_method"] = result['kernel_method'].replace(" ", "")
    config.state["vega_model"] = result['model']
    config.state["normalization_order"] = result['normalization_order']
    config.state["normalization_line"] = result['normalization_line']
    config.state["normalization_window"] = result['normalization_window']
    config.state["normalization_degree"] = result['normalization_degree']
    config.state["normalization_fittype"] = result['normalization_fittype']
    config.state["radialvelocity_nfwhm"] = result['radialvelocity_nfwhm']
    config.state["deconvolution_nfwhm"] = result['deconvolution_nfwhm']

    #
    # Get and store the Vega model 
    #

#    if new is True:
#
#        root = "Vega" + result["model"] + "_new.fits"
#
#    else:

    root = "Vega" + result["model"] + ".fits"
        
#    fullpath = os.path.join(setup.state["package_path"], "data", root)


#    root = "Vega" + result["model"] + ".fits"
#    print(root)


    fullpath = mishu.fetch(root)


    hdul = fits.open(fullpath)
    data = hdul[1].data

    hdul.close()

    # Store the results

    config.state["vega_wavelength"] = data["wavelength"]
    config.state["vega_fluxdensity"] = data["flux density"]
    config.state["vega_continuum"] = data["continuum flux density"]
    config.state["vega_fitted_continuum"] = data["fitted continuum flux density"]

    normalized = data["flux density"] / data["continuum flux density"]
    config.state["vega_normalized_fluxdensity"] = normalized

    # Compute the dispersions over the standard order wavelengths

    vega_dispersions = np.empty(dict["norders"])
    wavelength_ranges = wavelength_ranges

    for i in range(dict["norders"]):


        zleft = data["wavelength"] > wavelength_ranges[i][0]
        zright = data["wavelength"] < wavelength_ranges[i][1]

        zselection = np.logical_and(zleft, zright)
        pixels = np.arange(np.sum(zselection))

        fit = polyfit_1d(pixels, data["wavelength"][zselection], 1)
        vega_dispersions[i] = fit["coeffs"][1]

    # Store the results

    config.state["vega_dispersions"] = vega_dispersions

    #
    # Set the done variables
    #

    config.state["load_done"] = True
    config.state["prepare_done"] = False
    config.state["rv_done"] = False
    config.state["kernel_done"] = False




def _load_modeinfo(
    fullpath:str,
    mode:str):    

    """
    Load the specifics for a specific mode.

    Parameters
    ----------
    fullpath : str
        The path to the modeinfo file.
        Expected to be  pyspextool/instruments/<spex | uspex>/telluric_modeinfo.dat

    mode : str
        The mode name.  
    
    Returns
    -------
    dict
        `"method"` : str
    

    Loads data into memory:

        tc.state['kernel_method']
        tc.state['model']
        tc.state['normalized_order']
        tc.state['normalization_line']
        tc.state['normalization_window']
        tc.state['normalization_degree']
        tc.state['radialvelocity_nfwhm']
        tc.state['deconvolution_nfwhm']


    """

    # 
    # Check parameters
    #

    check_parameter('_load_modeinfo', 'fullpath', fullpath, "str")

    check_parameter('_load_modeinfo', 'mode', mode, "str")

    #
    # Read the file
    #

    values = np.loadtxt(fullpath, comments="#", delimiter="|", dtype="str")

    # Deal with the fact that there might only be one mode

    if np.ndim(values) == 1:

        values = np.expand_dims(values, 0)


    # Figure out which mode we are dealing with

    modes = list(values[:, 0])
    modes = np.array([x.strip() for x in modes])

    # Grab the values for that mode

    z = np.where(modes == mode)[0]
    
    method = str(values[z, 1][0])
    model = str(values[z, 2][0]).strip()

    order = None if str(values[z, 3][0]).strip() == "" else int(values[z, 3][0])

    line = None if str(values[z, 4][0]).strip() == "" else float(values[z, 4][0])

    if str(values[z, 5][0]).strip() == "":

        window = None

    else:

        window = str(values[z, 5][0]).split()
        window = [float(x) for x in window]

    degree = None if str(values[z, 6][0]).strip() == "" else int(values[z, 6][0])

    fittype = (
        None if str(values[z, 7][0]).strip() == "" else str(values[z, 7][0]).strip()
    )

    rv_nfwhm = None if str(values[z, 8][0]).strip() == "" else float(values[z, 8][0])

    dc_nfwhm = None if str(values[z, 9][0]).strip() == "" else float(values[z, 9][0])

    #
    # Save and return the results
    #
    
    dict = {'kernel_method':method,
            'model':model,
            'normalization_order':order,
            'normalization_line':line,            
            'normalization_window':window,            
            'normalization_degree':degree,            
            'normalization_fittype':fittype,
            'radialvelocity_nfwhm':rv_nfwhm,
            'deconvolution_nfwhm':dc_nfwhm}

    return dict

            



    

    








