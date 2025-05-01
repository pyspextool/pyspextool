import os
import numpy as np
import logging
from astropy.coordinates import SkyCoord, angular_separation
from astropy.io import fits
from astropy.table.table import Table
import astropy.units as u
from astroquery.simbad import Simbad

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.utils.coords import ten
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import make_full_path
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.fit.polyfit import polyfit_1d
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.load_atmosphere import load_atmosphere
from pyspextool.io.sptype2teff import sptype2teff
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.setup_utils import mishu

logger = logging.getLogger(__name__)

def load_spectra(
    object_file: str,
    standard_file: str,
    standard_info: list | str | dict,
    output_filename: str,
    correction_type: list = "A0 V",
    verbose: bool = None,
    new=False,
):
    """

    Parameters
    ----------
    object_file : str
        The name of the pySpextool FITS file containing the object spectra.

    standard_file : str
        The name of the pySpextool FITS file containing the standard spectra.

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

        tc.state['object_file']
        tc.state['standard_file']
        tc.state['standard_info']
        tc.state['output_filename']
        tc.state['reflectance']
        tc.state['ew_scale']
        tc.state['standard_rv']
        tc.state['standard_redshift']
        tc.state['1+z']
        tc.state['rms_deviation']
        tc.state['max_deviation']
        tc.state['load_done']
        tc.state['prepare_done']
        tc.state['rv_done']
        tc.state['kernel_done']



    """

    #
    # Check the parameters and keywords
    #

    check_parameter("load_spectra", "object_file", object_file, "str")

    check_parameter("load_spectra", "standard_file", standard_file, "str")

    check_parameter(
        "load_spectra", "standard_info", standard_info, ["str", "list", "dict"]
    )

    check_parameter("load_spectra", "output_filename", output_filename, "str")

    check_parameter(
        "load_spectra",
        "correction_type",
        correction_type,
        "str",
        possible_values=setup.state["telluric_correctiontypes"],
    )

    check_parameter("load_spectra", "verbose", verbose, ["NoneType", "bool"])

    check_qakeywords(verbose=verbose)

    #
    # Clear the state variables
    #

    tc.state.clear()

    #
    # Store user inputs
    #

    tc.state["object_file"] = object_file
    tc.state["standard_file"] = standard_file
    tc.state["standard_info"] = standard_info
    tc.state["output_filename"] = output_filename
    tc.state["type"] = correction_type.replace(" ", "")

    logging.info(" Telluric Correction\n--------------------------\n")

    #
    # Load the files into memory
    #

    # load_data()
    _load_object_data()

    #
    # Load the standard star information
    #
    _load_H_lines()
    _load_standard_data()
    _load_standard_info()

    # Compare the object and standard spectra to make sure they
    # are compatible and to determine the angular
    # separation and airmass difference
    _compare_object_standard()

    #
    # Get the mode info and store pertinent information
    #

    _load_modeinfo()

    #
    # Load the instrument profile parameters
    #

    _load_ipcoefficients()

    #
    # Load Vega Model if required
    #

    if tc.state["type"] == "A0V":
        _load_vegamodel(new=new)

    #
    # Set variables that need to be set to a default value.
    #

    tc.state["ew_scale"] = 1.0
    tc.state["standard_rv"] = 0.0
    tc.state["standard_redshift"] = 0.0
    tc.state["1+z"] = 0.0
    tc.state["rms_deviation"] = np.nan
    tc.state["max_deviation"] = np.nan
    tc.state["shift"] = np.zeros(
        (tc.state["object_norders"], tc.state["object_napertures"])
    )

    #
    # Set the done variables
    #

    tc.state["load_done"] = True
    tc.state["prepare_done"] = False
    tc.state["rv_done"] = False
    tc.state["kernel_done"] = False


def _load_standard_data(standard_filename: str = None):
    """
    Loads the standard observations.

    Parameters
    ----------
    None

    Uses
    - standard_filename: str
        Filename of a standard spectrum file.

    Returns
    -------
    None

    Loads data into memory.

        tc.state['standard_spectra']
        tc.state['standard_data']
        tc.state['standard_hdrinfo']
        tc.state['mode']
        tc.state['resolving_power']
        tc.state['standard_norders']
        tc.state['standard_orders']
        tc.state['standard_wavelengthranges']
        tc.state['standard_dispersions']
        tc.state['standard_fwhm']
        tc.state['slitw_pix']
        tc.state['atmospheric_transmission']

    """

    logging.info(" Loading the standard spectrum.")

    # Read in the standard file

    if standard_filename is None:
        fullpath = make_full_path(
            setup.state["proc_path"], tc.state["standard_file"], exist=False
        )
    else:
        fullpath = make_full_path(
            setup.state["proc_path"], standard_filename, exist=False
        )
    
    print(fullpath)

    standard_spectra, standard_data = read_spectra_fits(fullpath)

    standard_hdrinfo = get_headerinfo(
        standard_data["header"], keywords=setup.state["telluric_keywords"]
    )

    #
    # Store the results
    #

    tc.state["standard_spectra"] = standard_spectra
    tc.state["standard_data"] = standard_data
    tc.state["standard_hdrinfo"] = standard_hdrinfo
    tc.state["mode"] = standard_hdrinfo["MODE"][0]
    tc.state["resolving_power"] = standard_hdrinfo["RP"][0]
    tc.state["standard_norders"] = standard_data["norders"]
    tc.state["standard_orders"] = standard_data["orders"]

    # Does the standard have only one aperture?

    if standard_data["napertures"] != 1:

        message = "The standard has more than one aperture and can only " + "have one."
        raise pySpextoolError(message)

    #
    # Compute the minimum and maximum wavelengths and dispersions for each
    # standard order (This may go away eventually once you implement
    # mikeline_convolve).
    #

    nwavelengths = np.shape(tc.state["standard_spectra"])[-1]
    pixels = np.arange(nwavelengths)
    dispersions = np.empty(tc.state["standard_norders"])
    wavelength_ranges = []

    for i in range(tc.state["standard_norders"]):

        # Do min max first

        min = np.nanmin(tc.state["standard_spectra"][i, 0, :])
        max = np.nanmax(tc.state["standard_spectra"][i, 0, :])

        wavelength_ranges.append(np.array([min, max]))

        # Now do the dispersion

        fit = polyfit_1d(pixels, tc.state["standard_spectra"][i, 0, :], 1)
        dispersions[i] = fit["coeffs"][1]

    # Store the results

    tc.state["standard_wavelengthranges"] = wavelength_ranges
    tc.state["standard_dispersions"] = dispersions
    tc.state["standard_fwhm"] = dispersions * standard_data["slitw_pix"]
    tc.state["slitw_pix"] = standard_data["slitw_pix"]

    #
    # Load the atmospheric transmission
    #

    atmospheric_transmission = np.copy(tc.state["standard_spectra"])
    wavelengths, transmission = load_atmosphere(2000)

    for i in range(tc.state["standard_norders"]):

        atmosphere = linear_interp1d(
            wavelengths, transmission, tc.state["standard_spectra"][i, 0, :]
        )
        atmospheric_transmission[i, 1, :] = atmosphere

    # Store the results

    tc.state["atmospheric_transmission"] = atmospheric_transmission


def _load_H_lines(fullpath: str = None):
    """
    #
    # Load the Hydrogen lines
    #
    Parameters
    ----------
    fullpath : Path
        The path to the H1.dat file.
        Expected to be  pyspextool/data/H1.dat

    Sets:
    - tc.state['H_wavelengths']
    - tc.state['H_ids']
    """
    if fullpath is None:
        fullpath = os.path.join(setup.state["package_path"], "data", "HI.dat")

    wavelength, lineid = np.loadtxt(
        fullpath, comments="#", unpack=True, dtype="str", delimiter="|"
    )

    # Store the results

    tc.state["H_wavelengths"] = np.array(wavelength).astype(float)
    tc.state["H_ids"] = lineid


def _load_object_data(fullpath: str = None):
    """
    Loads the object spectrum.

    Parameters
    ----------
    None

    Uses
    - fullpath
        Path to the object spectrum file.

    Returns
    -------
    None

    Loads data into memory.
        tc.state['object_spectra']
        tc.state['object_data']
        tc.state['object_hdrinfo']
        tc.state['object_norders']
        tc.state['object_napertures']
        tc.state['object_orders']
        tc.state['xlabel']
        tc.state['H_wavelengths']
        tc.state['H_ids']


    """
    logging.info(" Loading the object spectrum.")

    # Object first

    if fullpath is None:
        fullpath = make_full_path(
            setup.state["proc_path"], tc.state["object_file"], exist=True
        )

    object_spectra, object_data = read_spectra_fits(fullpath)

    object_hdrinfo = get_headerinfo(
        object_data["header"], keywords=setup.state["telluric_keywords"]
    )

    # Check whether it was created by the extract or combine module -
    # temporary fix

    if object_hdrinfo["AVE_AM"][0] is None:

        object_hdrinfo2 = get_headerinfo(
            object_data["header"], keywords=setup.state["combine_keywords"]
        )
        object_hdrinfo["AVE_AM"] = [object_hdrinfo2["AM"][0], "Average airmass"]

    #        object_hdrinfo2 = get_headerinfo(header,keywords=setup.state['combine_keywords'])
    #        object_hdrinfo['AVE_AM'] = [object_hdrinfo2['AM'][0],'Average airmass']

    #
    # Store the results
    #
    tc.state["object_spectra"] = object_spectra
    tc.state["object_data"] = object_data
    tc.state["object_hdrinfo"] = object_hdrinfo
    tc.state["object_norders"] = object_data["norders"]
    tc.state["object_napertures"] = object_data["napertures"]
    tc.state["object_orders"] = object_data["orders"]
    tc.state["xlabel"] = object_hdrinfo["LXLABEL"][0]

    #
    # Now get the object ranges
    #

    wavelength_ranges = []

    for i in range(tc.state["object_norders"]):

        idx = i * tc.state["object_napertures"]
        min = np.nanmin(tc.state["object_spectra"][idx, 0, :])
        max = np.nanmax(tc.state["object_spectra"][idx, 0, :])

        wavelength_ranges.append(np.array([min, max]))

    # Store the results

    tc.state["object_wavelengthranges"] = wavelength_ranges


def _compare_object_standard(
    standard_fullpath: str = None, object_fullpath: str = None
):
    """
    Compare the object and standard spectra to determine the angular
    separation and airmass difference.

    Parameters
    ----------
    standard_full_path : str
        The full path to the standard spectrum file.
    object_full_path : str
        The full path to the object spectrum file.

    Returns
    -------
    None

    Loads data into memory:
        tc.state['delta_angle']
        tc.state['delta_airmass']

    """
    if standard_fullpath is not None:
        standard_data = read_spectra_fits(standard_fullpath)(1)
        standard_hdrinfo = get_headerinfo(
            standard_data["header"], keywords=setup.state["telluric_keywords"]
        )
    else:
        standard_data = tc.state["standard_data"]
        standard_hdrinfo = tc.state["standard_hdrinfo"]

    if object_fullpath is not None:
        object_data = read_spectra_fits(object_fullpath)(1)
        object_hdrinfo = get_headerinfo(
            object_data["header"], keywords=setup.state["telluric_keywords"]
        )
    else:
        object_data = tc.state["object_data"]
        object_hdrinfo = tc.state["object_hdrinfo"]

    # Does the standard have at least the order numbers of the object?

    intersection = np.intersect1d(object_data["orders"], standard_data["orders"])

    if np.size(intersection) != object_data["norders"]:

        message = "The standard lacks an order the object has."
        raise pySpextoolError(message)

    # Determine the angular separation between the object and standard

    obj_long = 15 * ten(object_hdrinfo["RA"][0], toradians=True)
    obj_lat = ten(object_hdrinfo["DEC"][0], toradians=True)
    std_long = 15 * ten(standard_hdrinfo["RA"][0], toradians=True)
    std_lat = ten(standard_hdrinfo["DEC"][0], toradians=True)

    angle = angular_separation(obj_long, obj_lat, std_long, std_lat) * 180 / np.pi

    # Determine the airmass difference

    delta_airmass = object_hdrinfo["AVE_AM"][0] - standard_hdrinfo["AVE_AM"][0]

    #
    # Store the results
    #

    tc.state["delta_angle"] = angle
    tc.state["delta_airmass"] = delta_airmass


def _load_vegamodel(standard_filename: str = None, new=False):
    """
    Loads the proper Vega model given the observing setup. 
    The observing setup is determined by the standard spectrum.

    Parameters
    ----------
    standard_filename : str
        The filename of a standard spectrum.

    Returns
    -------
    dict
        A dictionary containing the Vega model data.
        "vega_file" : str
            The name of the Vega model file.

        `"vega_wavelength"` : ndarray
            The wavelengths of the Vega model in microns.

        `"vega_fluxdensity"` : ndarray
            The flux density of the Vega model in erg s- cm-2 A-1

        `"vega_continuum"` : ndarray
            The continuum of the Vega model in erg s- cm-2 A-1

        `"vega_fitted_continuum"` : ndarray
            The fitted continuum of the Vega model in erg s- cm-2 A-1

        `"vega_normalized_fluxdensity"` : ndarray
            The normalized flux density of the Vega model.
            = `"vega_fluxdensity"/ `"vega_fitted_continuum"`

        `"vega_dispersions"` : ndarray
            The average dispersion (microns pixel-1) of the vega model
            over the wavelength range of each order.

    """

    logging.info(" Loading the Vega model.")

    #
    # Determine which Vega model to use and load
    #
    try:
        mode = tc.state["mode"]
    except KeyError:
        if standard_filename is not None:
            _load_standard_data(standard_filename)
        else:
            _load_standard_data()
        mode = tc.state["mode"]

    logging.debug(f"Mode: {mode}")

    try:
        model = tc.state["model"]
    except KeyError:        
        _load_modeinfo()
        model = tc.state["model"]

    logging.debug(f"Model: {model}")

    if new is True:

        root = "Vega" + model + "_new.fits"

    else:

        root = "Vega" + model + ".fits"

    vega_file = os.path.join(setup.state["package_path"], "data", root)

    if os.path.isfile(vega_file) is False:
        vega_file = mishu.fetch(root)

    hdul = fits.open(vega_file)
    data = hdul[1].data

    vega_wavelength = data["wavelength"]
    vega_fluxdensity = data["flux density"]
    vega_continuum = data["continuum flux density"]
    vega_fitted_continuum = data["fitted continuum flux density"]

    hdul.close()

    #
    # Compute the dispersions over the order wavelengths
    #

    vega_dispersions = np.empty(tc.state["standard_norders"])
    wavelength_ranges = tc.state["standard_wavelengthranges"]
    norders = len(vega_dispersions)

    for i in range(norders):

        zleft = vega_wavelength > wavelength_ranges[i][0]
        zright = vega_wavelength < wavelength_ranges[i][1]

        zselection = np.logical_and(zleft, zright)

        pixels = np.arange(np.sum(zselection))

        fit = polyfit_1d(pixels, vega_wavelength[zselection], 1)
        vega_dispersions[i] = fit["coeffs"][1]

    #
    # Create the output and return
    #

    normalized = vega_fluxdensity / vega_fitted_continuum

    result = {
    "vega_file": root,
    "vega_wavelength": vega_wavelength,
    "vega_fluxdensity": vega_fluxdensity,
    "vega_continuum": vega_continuum,
    "vega_fitted_continuum": vega_fitted_continuum,
    "vega_normalized_fluxdensity": normalized,
    "vega_dispersions": vega_dispersions,
    }

    tc.state["vega_wavelength"] = result["vega_wavelength"]
    tc.state["vega_fluxdensity"] = result["vega_fluxdensity"]
    tc.state["vega_continuum"] = result["vega_continuum"]
    tc.state["vega_fitted_continuum"] = result["vega_fitted_continuum"]
    tc.state["vega_normalized_fluxdensity"] = result["vega_normalized_fluxdensity"]
    tc.state["vega_dispersions"] = result["vega_dispersions"]
    
    return result


def _load_modeinfo(file: str = None):
    """
    Load the mode info given

    Parameters
    ----------
    path : str
        The path to the modeinfo file.
        Expected to be  pyspextool/instruments/<spex | uspex>/telluric_modeinfo.dat

    Returns
    -------
    None

    Loads data into memory:

        tc.state['method']
        tc.state['model']
        tc.state['normalized_order']
        tc.state['normalization_line']
        tc.state['normalization_window']
        tc.state['normalization_degree']
        tc.state['radialvelocity_nfwhm']
        tc.state['deconvolution_nfwhm']


    """
    if file is None:
        file = os.path.join(setup.state["instrument_path"], "telluric_modeinfo.dat")

    values = np.loadtxt(file, comments="#", delimiter="|", dtype="str")

    # Deal with the fact that there might only be one mode

    if np.ndim(values) == 1:

        values = np.expand_dims(values, 0)

    # Figure out which mode we are dealing with

    modes = list(values[:, 0])
    modes = np.array([x.strip() for x in modes])

    # Grab the values for that mode

    z = np.where(modes == tc.state["mode"])[0]

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

    # Save the results

    tc.state["method"] = method
    tc.state["model"] = model
    tc.state["normalization_order"] = order
    tc.state["normalization_line"] = line
    tc.state["normalization_window"] = window
    tc.state["normalization_degree"] = degree
    tc.state["normalization_fittype"] = fittype
    tc.state["radialvelocity_nfwhm"] = rv_nfwhm
    tc.state["deconvolution_nfwhm"] = dc_nfwhm


def _load_ipcoefficients():
    """
    Load the mode info given

    Parameters
    ----------
    None

    Returns
    -------
    None
    Loads data into memory:

        tc.state['ip_coefficients']

    """

    # Get the file name

    file = os.path.join(setup.state["instrument_path"], "IP_coefficients.dat")

    # Read the file

    slitw_arc, c0, c1, c2 = np.loadtxt(file, comments="#", unpack=True)

    # Find the right slit width and store the coefficients

    z = np.where(slitw_arc == tc.state["standard_hdrinfo"]["SLTW_ARC"][0])[0]

    ip_coefficients = np.squeeze(np.array((c0[z], c1[z], c2[z])))

    tc.state["ip_coefficients"] = ip_coefficients


def _load_standard_info():
    """
    Load the standard star spectral type, B-, and V-band magnitudes.

    Parameters
    ----------
    None

    Returns
    -------
    None
    Loads data into memory:

        tc.state['standard_name']
        tc.state['standard_sptype']
        tc.state['standard_vmag']
        tc.state['standard_bmag']
        tc.state['standard_teff']

    """

    if isinstance(tc.state["standard_info"], str):

        #
        # user has passed the name of the standard
        #

        # Get SIMBAD information of the standard

        logging.info(
            " Querying SIMBAD for standard star "
            + tc.state["standard_info"]
            + " information."
        )

        Simbad.add_votable_fields("sptype", "flux(B)", "flux(V)")
        table = Simbad.query_object(tc.state["standard_info"])

        if isinstance(table, Table):
            standard_name = table["MAIN_ID"][0]
            standard_sptype = table["SP_TYPE"][0]
            standard_vmag = float(table["FLUX_V"][0])
            standard_bmag = float(table["FLUX_B"][0])

        else:

            message = (
                'Standard name "{}"" was not found in SIMBAD; provide '
                "the correct name, correct coordinates, or a dictionary "
                'containing keys "id", "sptype", "bmag", '
                'and "vmag".'.format(tc.state["standard_name"])
            )

            raise pySpextoolError(message)

    if isinstance(tc.state["standard_info"], list):

        #
        # user has passed the coordinates of the standard
        #

        logging.info(" Querying SIMBAD for standard star information...")

        Simbad.add_votable_fields("id", "sptype", "flux(B)", "flux(V)")

        c = SkyCoord(
            tc.state["standard_info"][0],
            tc.state["standard_info"][1],
            unit=(u.hourangle, u.deg),
        )

        table = Simbad.query_region(c, radius="0d1m0s")

        if isinstance(table, Table):
            standard_name = table["MAIN_ID"][0]
            standard_sptype = table["SP_TYPE"][0]
            standard_vmag = float(table["FLUX_V"][0])
            standard_bmag = float(table["FLUX_B"][0])

        else:

            message = (
                'Standard coordiantes "{}"" were not found in SIMBAD; '
                "provide the correct name, correct coordinates, or a "
                'dictionary containing keys "id", "sptype", "bmag", '
                'and "vmag".'.format(tc.state["standard_info"])
            )

            raise pySpextoolError(message)

    if isinstance(tc.state["standard_info"], dict):

        #
        # user has passed the standard information
        #

        standard_name = tc.state["standard_info"]["id"]
        standard_sptype = tc.state["standard_info"]["sptype"]
        standard_vmag = float(tc.state["standard_info"]["vmag"])
        standard_bmag = float(tc.state["standard_info"]["bmag"])

    #
    # Store the results
    #

    tc.state["standard_name"] = standard_name
    tc.state["standard_sptype"] = standard_sptype
    tc.state["standard_vmag"] = standard_vmag
    tc.state["standard_bmag"] = standard_bmag
    tc.state["standard_teff"] = sptype2teff(standard_sptype)
