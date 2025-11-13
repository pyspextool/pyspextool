from os.path import basename, join
import logging
import numpy as np
from datetime import date, datetime
from astropy.io import fits
import astropy.units as u
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool import utils

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

    return header