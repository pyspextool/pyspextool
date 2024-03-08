from os.path import basename, join
import logging
import numpy as np
from datetime import date, datetime
from astropy.io import fits
import astropy.units as u
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from astropy.table import Table
from pyspextool.io.read_spectra_fits import read_spectra_fits

__all__ = ["convert_to_fits", "spectrum_isplottable"]

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

    spectra, header_dict = read_spectra_fits(input_file)

    wavelength = spectra[0, 0, :]
    flux = spectra[0, 1, :]
    flux_unc = spectra[0, 2, :]
    # mask?

    # header = compile_header(wavelength, **spectrum_info_all)
    header = header_dict["header"]
    header["HISTORY"] = (
        f"Converted {basename(input_file)} using pyspextool.io convert_to_fits"
    )

    # replace spaces in object name with underscores
    object_name = header["OBJECT"].replace(" ", "_").replace("/", "_")

    # convert obsdate to format YYYYMonDD
    obs_date = datetime.strptime(header["AVE_DATE"], "%Y-%m-%d").strftime("%Y%b%d")

    # Update Date the file was written
    header["DATE"] = date.today().strftime("%Y-%m-%d")

    x_units = header["XUNITS"]
    y_units = header["YUNITS"]
    if y_units == "DN s-1":
        y_units = "count s-1"

    spectrum_data_out = Table(
        {
            "wavelength": wavelength * u.Unit(x_units),
            "flux": flux * u.Unit(y_units),
            "flux_uncertainty": flux_unc * u.Unit(y_units),
        }
    )

    # Make the HDUs
    hdu1 = fits.BinTableHDU(data=spectrum_data_out)
    hdu1.header["EXTNAME"] = "SPECTRUM"
    hdu1.header.set("OBJECT", object_name, "Object Name")
    hdu0 = fits.PrimaryHDU(header=header)

    # Write the MEF with the header and the data
    spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

    fits_filename = f"{object_name}_{obs_date}.fits"
    fits_path = join(output_path, fits_filename)
    spectrum_mef.writeto(fits_path, overwrite=True, output_verify="exception")
    # TODO: think about overwrite
    logger.info(f"Wrote {fits_path}")

    return
