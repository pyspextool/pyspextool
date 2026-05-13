from os.path import basename, join
import logging
import numpy as np
from datetime import date, datetime
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool import utils

__all__ = ["convert_to_fits", 
           "spectrum_isplottable",
           "add_missing_keywords",
           "add_wavelength_keywords",
           "add_observation_date",
           "check_header",]

logger = logging.getLogger("pyspextool")


def spectrum_isplottable(spectrum_path, raise_error=True, show_plot=False):
    """
    Check if spectrum is plottable
    """
    # load the spectrum and make sure it's a Spectrum1D object

    try:
        # spectrum: Spectrum1D = load_spectrum(spectrum_path) #astrodbkit2 method
        spectrum = Spectrum1D.read(spectrum_path)
    except Exception as e:
        msg = (
            str(e) + f"\n{spectrum_path}: \n" "unable to load file as Spectrum1D object"
        )
        if raise_error:
            logger.error(msg)
            raise IOError(msg)
        else:
            logger.warning(msg)
            return False

    # checking spectrum has good units and not only NaNs
    try:
        wave: np.ndarray = spectrum.spectral_axis.to(u.micron).value
        flux: np.ndarray = spectrum.flux.value
    except AttributeError as e:
        msg = str(e) + f"{spectrum_path}: unable to parse spectral axis"
        if raise_error:
            logger.error(msg)
            raise IOError(msg)
        else:
            logger.warning(msg)
            return False
    except u.UnitConversionError as e:
        msg = e + f"{spectrum_path}: unable to convert spectral axis to microns"
        if raise_error:
            logger.error(msg)
            raise IOError(msg)
        else:
            logger.warning(msg)
            return False
    except ValueError as e:
        msg = e + f"{spectrum_path}: Value error"
        if raise_error:
            logger.error(msg)
            raise IOError(msg)
        else:
            logger.warning(msg)
            return False

    # check for NaNs
    nan_check: np.ndarray = ~np.isnan(flux) & ~np.isnan(wave)
    wave = wave[nan_check]
    flux = flux[nan_check]
    if not len(wave):
        msg = f"{spectrum_path}: spectrum is all NaNs"
        if raise_error:
            logger.error(msg)
            raise IOError(msg)
        else:
            logger.warning(msg)
            return False

    if show_plot:
        plt.plot(wave, flux)
        plt.show()

    return True


def convert_to_fits(input_file, output_path="."):
    """
    Convert a pyspextool file to a specutils compatibile FITS file
    TODO: write tests
    TODO: think about header keywords and Spectrum DM 1.2

    Parameters
    ----------
    input_file : str
    output_path : str

    Returns
    -------
    None

    Writes converted FITS file to output_path

    """

    # expect files with individual orders
    # input file may have multiple orders and multiple aperatures (for extended sources)
    spectra, header_dict = read_spectra_fits(input_file)
    #print(header_dict)
    #print(header_dict["XUNITS"])
    header = header_dict["astropyheader"]

    spectrum_data_out = make_table_of_spectra(spectra, header)
    
    new_header = do_header_things(header, input_file)

    # replace spaces in object name with underscores
    object_name = header["OBJECT"].replace(" ", "_").replace("/", "_")
    
    # convert obsdate to format YYYYMonDD
    obs_date = datetime.strptime(header["AVE_DATE"], "%Y-%m-%d").strftime("%Y%b%d")

    # Make the HDUs
    hdu1 = fits.BinTableHDU(data=spectrum_data_out)
    hdu1.header["EXTNAME"] = "SPECTRUM"
    hdu1.header.set("OBJECT", object_name, "Object Name")
    hdu0 = fits.PrimaryHDU(header=new_header)

    # Write the MEF with the header and the data
    spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

    fits_filename = f"{object_name}_{obs_date}.fits"
    fits_path = join(output_path, fits_filename)
    spectrum_mef.writeto(fits_path, overwrite=True, output_verify="exception")
    # TODO: think about overwrite
    logger.info(f"Wrote {fits_path}")

    #print(header_dict)

    return


def make_table_of_spectra(spectra, header):
    # TODO: ADD UCD for FluxAxis.ucd = phot.flux.density;em.wl 

    x_units = header["XUNITS"]
    y_units = header["YUNITS"]
    if y_units == "DN s-1":
        y_units = "count s-1"
        # TODO: fluxAxis.ucd = arith.rate;phot.count 

    n_orders     = header["NORDERS"] # each order has a different wavelength range
    n_aperatures = header["NAPS"]

    for order in range(n_orders):
        for aperture in range(n_aperatures):
            idx = order * n_aperatures + aperture
            non_nans = utils.arrays.trim_nan(spectra[idx, 0, :],flag=2) # flag=2 removes leading and trailing NaNs
            wavelength = spectra[idx, 0, :][non_nans] # floating point
            flux = spectra[idx, 1, :][non_nans] # floating point
            flux_unc = spectra[idx, 2, :][non_nans] # floating point
            mask = spectra[idx, 3, :][non_nans].astype(np.uint8) # read in as floating point, but can be converted to 8-bit array

            if idx == 0:
                spectrum_data_out = Table(
                    {
                        "wavelength": wavelength * u.Unit(x_units),
                        "flux": flux * u.Unit(y_units),
                        "flux_uncertainty": flux_unc * u.Unit(y_units),
                        "mask": mask * u.Unit("1"),  # mask is a bitmask, so unit is 1
                        "order": np.zeros(len(mask)) + order,  # keep track of the orders
                    }
                )
            else:
                # Create a new table and append to spectrum data out
                spectrum_data_new = Table([wavelength * u.Unit(x_units), flux * u.Unit(y_units), 
                                           flux_unc * u.Unit(y_units), mask * u.Unit("1"), np.zeros(len(mask)) + order], 
                                           names=('wavelength', 'flux', 'flux_uncertainty', 'mask', 'order'))
                spectrum_data_out = vstack( [spectrum_data_out, spectrum_data_new] )
                

    # Sort the order of the wavelengths
    if n_orders == 1:

        spectrum_data_out.sort('wavelength')

    else: # conserve the list of orders by wavelength order

        # find unique orders
        orders = np.unique(spectrum_data_out['order'])

        # compute the minimum wavelength for each order
        min_wave_per_order = {o: np.min(spectrum_data_out[spectrum_data_out['order'] == o]['wavelength']) for o in orders}

        # sort orders by their minimum wavelength
        sorted_orders = sorted(orders, key=lambda o: min_wave_per_order[o])

        # stack tables in the desired order, keeping wavelength sorted within each order
        table_list = []
        for o in sorted_orders:
            spectrum_order = spectrum_data_out[spectrum_data_out['order'] == o]
            spectrum_order.sort('wavelength')
            table_list.append(spectrum_order)

        # replace the spectrum data with the new table    
        spectrum_data_out = vstack(table_list)


    return spectrum_data_out


def do_header_things(header, input_file):
    # Dealing with the primary header HDU0
    header["HISTORY"] = (
        f"Converted {basename(input_file)} using pyspextool.io.convert_to_fits"
    )

    # Update Date the file was written
    header["DATE"] = date.today().strftime("%Y-%m-%d")

    # Add in the missing keywords (if any)
    new_header = add_missing_keywords(header=header)

    return new_header










##############################################################################################################













def add_missing_keywords(header=None, format="ivoa-spectrum-dm-1.2", keywords=None):
    """Finds the keywords that are missing from a header and adds them with blank values

    Inputs
    -------
    header: fits.Header
        a fits header object or dictionary of header values

    format: string
        header data models to enforce. Options, 'simple-spectrum' and 'ivoa-spectrum-dm-1.2`

    keywords: list
        a list of keywords to check for if format is not specified

    Returns
    -------
    FITS header object

    Examples
    --------
    Returns a header with keywords but blank values
    >>> new_header = add_missing_keywords()

    Adds missing keywords (with blank values) to an existing header
    using the simple spectrum data model
    >>> new_header = add_missing_keywords(old_header, format='simple-spectrum')
    """

    # If no header was provided, start with a blank one
    if header is None:
        header = fits.Header()

    if keywords is None and format is None:
        format = "ivoa-spectrum-dm-1.2"

    keywords = get_keywords(format)

    if format == "ivoa-spectrum-dm-1.2":
        header["VOCLASS"] = "Spectrum-1.2"

    # Add in missing values from current header
    coords = SkyCoord('%s %s'%(header["RA"],header["DEC"]), unit=(u.hourangle, u.deg))
    header["RA_TARG"]  = coords.ra.degree
    header["DEC_TARG"] = coords.dec.degree
    header["INSTRUME"] = header["INSTR"]
    header["OBSERVAT"] = "Maunakea Observatories"
    header["DATE-OBS"] = header["DATE"]
    header["TMID"]     = header["AVE_MJD"]
    header["TELAPSE"]  = header["TOTITIME"]

    missing_keywords = []
    # Loop through the original header and add keywords with blank values to the new header
    for keyword, comment in keywords:
        value = header.get(keyword)
        if value is None:
            header.set(keyword, None, comment)
            missing_keywords.append((keyword, comment))

    # Loop over missing keywords and print for copy and paste purposes
    print(
        "COPY AND PASTE THE FOLLOWING COMMANDS INTO YOUR SCRIPT \n"
        "Replace <value> with the appropriate value for your dataset \n"
        "If you're not sure of the correct value, use None \n"
        "If you started with a header object not called `header`, \n"
        "replace 'header' with the name of your header object \n"
        "Use the `astrodb_utils.fits.add_wavelength_keywords` function \n"
        "to add the SPEC_VAL, SPEC_BW, TDMID1, TDMAX1, and SPECBAND keywords \n"
    )

    for keyword, comment in missing_keywords:
        print(f"header.set('{keyword}', \"<value>\")")  # {comment}")

    return header


def add_wavelength_keywords(header=None, wavelength_data=None):
    """Uses wavelength array to generate header keywords

    Inputs
    -------
    wavelength_data: astropy.units.Quantity
        an array of wavelengths. should include units.

    header_dict: header
        a Header object

    Returns
    -------
    None

    Examples
    --------
    >>> wavelength = np.arange(5100, 5300)*u.AA
    >>> add_wavelength_keywords(header=new_header, wavelength_data = wavelength)

    """

    # Make new, blank header
    if header is None:
        header = fits.Header()

    # Use wavelength data to calculate header values
    w_min = min(wavelength_data).astype(np.single)
    w_max = max(wavelength_data).astype(np.single)
    width = (w_max - w_min).astype(np.single)
    w_mid = ((w_max + w_min) / 2).astype(np.single)
    bandpass = assign_ucd(w_mid)

    header.set("SPECBAND", bandpass)
    header.set("SPEC_VAL", w_mid.value, f"[{w_mid.unit}] Characteristic spec coord")
    header.set("SPEC_BW", width.value, f"[{width.unit}] Width of spectrum")
    header.set("TDMIN1", w_min.value, f"[{w_min.unit}] Starting wavelength")
    header.set("TDMAX1", w_max.value, f"[{w_max.unit}] Ending wavelength")
    header["HISTORY"] = (
        "Wavelength keywords added by astrodb_utils.fits.add_wavelength_keywords"
    )

    # return header


def add_observation_date(header=None, date=None):
    """Adds the observation date to the header

    Inputs
    -------
    header: fits.Header
        a fits header object or dictionary of header values

    date: string
        the date of the observation

    Returns
    -------
    None

    Examples
    --------
    >>> add_observation_date(header, date='2021-06-01')
    """

    if header is None:
        header = fits.Header()

    if date is None:
        raise ValueError("Date of observation is required")

    try:
        obs_date = dateparser.parse(date)
        logger.debug(f"Date of observation: {date} parsed to {obs_date}")
        if obs_date is not None:
            obs_date_short = obs_date.strftime("%Y-%m-%d")
            obs_date_long = obs_date.strftime("%b %d, %Y")
            header.set("DATE-OBS", obs_date_short, "date of the observation")
            print(f"Date of observation: {obs_date_long}")
            print(f"DATE-OBS set to : {obs_date_short}.")
        else:
            raise ValueError(f"Date could not be parsed by dateparser.parse: {date}")
    except Exception as e:
        raise e


def check_header(header=None, format="ivoa-spectrum-dm-1.2", ignore_simbad=False):
    """
    Check the header of a FITS file for required keywords and other properties.

    Parameters
    ----------
    header : astropy.io.fits.Header
        The header object to be checked.
    format : str, optional
        The data model of the FITS file. Default is 'ivoa-spectrum-dm-1.2'.
    ignore_simbad : bool, optional
        Whether to ignore checking SIMBAD coordinates. Default is False.

    Returns
    -------
    bool
        True if the header passes all checks, False otherwise.

    Raises
    ------
    ValueError
        If the header is not provided.

    """
    # TODO: Check DOI

    result = True

    if header is None:
        raise ValueError("Header is required")

    # check for missing keywords
    keywords = get_keywords(format)
    missing_keywords = []
    for keyword, comment in keywords:
        value = header.get(keyword)
        if value is None:
            if len(missing_keywords) == 0:
                print("The following keywords are not set in the header:")
            missing_keywords.append((keyword, comment))
            print(f"{keyword} : {comment}")

    coord = make_skycoord(header)
    if coord is None:
        print("Coordinates are required")
        result = False
        logger.debug(f"Coords: {coord}, result: {result}")

    if ignore_simbad is False:
        name_check = check_simbad_name(header)

        if name_check is True:
            print("make sure SIMBAD coords match header coords")
            # check_ra_dec_simbad(simbad_name_results):
        else:
            result = False
        logger.debug(f"name check returns {result}")

    if check_date(header) is False:
        result = False
        logger.debug(f"Date check returns {result}")

    logger.debug(f"Header check returns {result}")
    return result


def get_keywords(format):

    formats = ["simple-spectrum", "ivoa-spectrum-dm-1.2"]
    if format not in formats:
        msg = f"(Format must be one of these: {formats})"
        raise ValueError(msg)

    if format == "simple-spectrum":
        keywords = [
            ("OBJECT", "Name of observed object"),
            ("RA_TARG", "[deg] target position"),
            ("DEC_TARG", "[deg] target position"),
            ("INSTRUME", "Instrument name"),
            ("TELESCOP", "Telescope name"),
            ("DATE-OBS", "Date of observation"),
            ("TELAPSE", "[s] Total elapsed time (s)"),
            ("AUTHOR", "Authors of original dataset"),
            ("VOREF", "URL, DOI, or bibcode of original publication"),
            ("VOPUB", "Publisher"),  # SIMPLE
            ("CONTRIB1", "Contributor who generated this header"),
            ("SPEC_VAL", "[angstrom] Characteristic spectral coordinate"),
            ("SPEC_BW", "[angstrom] width of spectrum"),
            ("TDMIN1", "Start in spectral coordinate"),
            ("TDMAX1", "Stop in spectral coordinate"),
            ("SPECBAND", "SED.bandpass"),
            ("APERTURE", "[arcsec] slit width"),
        ]
    elif format == "ivoa-spectrum-dm-1.2":
        keywords = [
            ("OBJECT", "Name of observed object"),
            ("RA_TARG", "[deg] target position"),
            ("DEC_TARG", "[deg] target position"),
            ("INSTRUME", "Instrument ID"),
            ("TELESCOP", "Telescope"),        # not in DM 1.2
            ("OBSERVAT", "Observatory name"), # not in DM 1.2
            ("VOCLASS", "Data model name and version"),
            ("VOPUB", "Publisher"),
            ("VOREF", "URL, DOI, or bibcode of original publication"),
            ("TITLE", "Dataset title "),
            ("AUTHOR", "Authors of the original dataset"),
            ("CONTRIB1", "Contributor who generated this file"),  # optional
            ("DATE-OBS", "Date of observation"),  # optional
            ("TMID", "[d] MJD of exposure mid-point"),
            ("TELAPSE", "[s] Total elapsed time (s)"),
            ("SPEC_VAL", "[angstrom] Characteristic spectral coordinate"),
            ("SPEC_BW", "[angstrom] width of spectrum"),
            ("TDMIN1", "Start in spectral coordinate"),
            ("TDMAX1", "Stop in spectral coordinate"),
            ("SPECBAND", "SED.bandpass"),
            ("APERTURE", "[arcsec] slit width"),
        ]

    return keywords


def make_skycoord(header):
    """
    Make a SkyCoord object from the RA and Dec in the header
    Look for coordinates in either RA_TARG and DEC_TARG or RA and DEC keywords
    """

    if header.get("RA_TARG") is None and header.get("RA") is None:
        print("RA_TARG or RA is required")
        ra = None
        return None

    #print(header.get("RA_TARG"))
    #print(header.get("RA"))
    #print(header.get("DEC_TARG"))
    #print(header.get("DEC"))

    if header.get("RA_TARG") is not None:
        ra = float(header.get("RA_TARG"))
    elif header.get("RA") is not None:
        ra = float(header.get("RA"))

    if ra > 360:
        print("RA_TARG does not appear to be in degrees")
        print(f"RA_TARG: {ra}")
        print("RA_TARG should be in degrees")
        return None

    logger.debug(f"RA: {ra}")

    if header.get("DEC_TARG") is None and header.get("DEC") is None:
        print("DEC_TARG or DEC is required")
        dec = None
        return None

    if header.get("DEC_TARG") is not None:
        dec = float(header.get("DEC_TARG"))
    elif header.get("DEC") is not None:
        dec = float(header.get("DEC"))

    if dec > 90 or dec < -90:
        print("DEC_TARG value is out of the expected range.")
        print(f"DEC_TARG: {dec}")
        print("DEC_TARG should be in degrees in the range -90 to 90")
        return None

    logger.debug(f"DEC: {dec}")

    # Check if ra and dec could be read into SkyCoord object and converted to sexagesimal
    try:
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        coord_str = coord.to_string("hmsdms")
        print(f"coordinates converted to sexagesimal: {coord_str}")
        return coord
    except Exception as e:
        print(
            f"coordinates ({ra},{dec}) could not be converted to Skycoord object: {e}"
        )
        return None


def check_simbad_name(header):
    # search SIMBAD for object name
    object_name = header.get("OBJECT")

    if object_name is None:
        print("OBJECT name is required")
        return False

    simbad_name_results = Simbad.query_object(object_name)
    logger.debug(f"SIMBAD results for object name {object_name}: {simbad_name_results}")
    print(f"SIMBAD results for object name {object_name}: {simbad_name_results}")
    coord = make_skycoord(header)

    if simbad_name_results is None and coord is None:
        print(f"Object name {object_name} not found in SIMBAD")
        return False

    # search for object using RA/Dec coordinates
    elif simbad_name_results is None and coord is not None:
        coord_str = coord.to_string("hmsdms")
        print(f'Object name {object_name} not found in SIMBAD, trying 60" coord search')
        simbad_coord_results = Simbad.query_region(coord, radius="0d0m60s")
        if simbad_coord_results is not None:
            print(f'SIMBAD objects within 60" search of {coord_str}:')
            print(simbad_coord_results)
            return False

    if len(simbad_name_results) == 1:
        print(f"Object name {object_name} found in SIMBAD")
        print(simbad_name_results)
        return True
    else:
        print(f"Multile object name {object_name} matches found in SIMBAD")
        print(simbad_name_results)
        return False


def check_ra_dec_simbad(simbad_name_results):
    # TODO: make this a skycoord object  comparison and figure out correct tolerance
    result = True
    # check RA and  Dec agree with SIMBAD
    # ra_simbad = simbad_name_results['RA']
    # dec_simbad = simbad_name_results['DEC']

    # ra_check = np.isclose(ra, ra_simbad, atol=0.1) # check if ra is close to simbad ra
    # dec_check = np.isclose(dec, dec_simbad, atol=0.1) # check if dec
    # if not ra_check or not dec_check:
    #    print("RA_TARG and DEC_TARG do not match SIMBAD coordinates")
    #    print(f"RA_TARG: {ra}, DEC_TARG: {dec}")
    #    print(f"SIMBAD RA: {ra_simbad}, DEC: {dec_simbad}")
    #     result = False

    return result


def check_date(header):
    # check date can be turned to dateTime object
    date = header.get("DATE-OBS")
    if date is None:
        print("DATE-OBS is not set in header")
        result = False
    else:
        try:
            obs_date = dateparser.parse(date)
            obs_date_long = obs_date.strftime("%b %d, %Y")
            print(f"DATE-OBS set to : {date}.")
            print(f"Date of observation: {obs_date_long}")
            result = True
        except Exception as e:
            print(
                f"Date ({date})could not be converted to Python DateTime object \n {e}"
            )
            result = False

    logger.debug(f"Date of observation: {date} returns {result}")
    return result

