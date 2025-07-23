from os.path import basename, join
import logging
import numpy as np
from datetime import date, datetime
from astropy.io import fits
import astropy.units as u
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io.votable import from_table, writeto
from astropy.io.votable.tree import VOTableFile, Resource, Param
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io import convert_fits

__all__ = ["convert_to_votable"]

logger = logging.getLogger("pyspextool")


def convert_to_votable(input_file, output_path="."):
    """
    Convert a pyspextool file to a VO Table
    Note: this is a work in progress, not all features are implemented yet.

    Parameters
    ----------
    input_file : str
    output_path : str

    Returns
    -------
    None

    Writes converted VOTable file to output_path

    """
    # Most of this is being developed in the convert_fits.py file
    # If we want to convert to VOTable, we should use the make_table_of_spectra function to start

    make_table_of_spectra = convert_fits.make_table_of_spectra(input_file)

    # Add in the spectrum data
    init_table['flux'].unit        = y_units
    init_table['flux'].description = 'Flux'
    init_table['flux'].meta['ucd'] = 'phot.flux.density;em.opt'

    init_table['flux_uncertainty'].unit        = y_units
    init_table['flux_uncertainty'].description = 'Flux uncertainty'
    init_table['flux_uncertainty'].meta['ucd'] = 'stat.error;phot.flux.density'

    init_table['wavelength'].unit        = x_units
    init_table['wavelength'].description = 'Wavelength'
    init_table['wavelength'].meta['ucd'] = 'em.wl'


    # Take all the header data and make a Table
    # Create the VOTable structure manually
    votable  = VOTableFile()
    resource = Resource()
    votable.resources.append(resource)

    # Convert Astropy table to VOTable Table
    votable_table = from_table(init_table).resources[0].tables[0]
    resource.tables.append(votable_table)

    # Add FITS header keywords as PARAMs
    for key in header:
        param = Param(votable,
                      name=key,
                      datatype='char',
                      arraysize='*',
                      value=str(header[key]))
        resource.params.append(param)

        # Optionally add comment as description
        if key in header.comments:
            param.description = header.comments[key]
        resource.params.append(param)


    # Write to file
    vo_filename = f"{object_name}_{obs_date}.xml"
    vo_path = join(output_path, vo_filename)
    writeto(votable, vo_path)


    # TODO: think about overwrite
    logger.info(f"Wrote {vo_path}")

    return
