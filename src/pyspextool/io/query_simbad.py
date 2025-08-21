from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table.table import Table
import astropy.units as u
from astroquery.simbad import Simbad
from astroquery.exceptions import NoResultsWarning
import warnings

from pyspextool.io.check import check_parameter
from pyspextool.pyspextoolerror import pySpextoolError

def query_simbad(
    info:str | list | dict):

    """
    Query SIMBAD for spectral type, B magnitude, and V magnitudes.  
    
    Parameters
    ----------
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

    Returns
    -------
    dict
        `"name"` : str
            A string giving the name of a star.
        
        `"sptype"` : str
            A string giving the spectral type of a star.

        `"bmag"` : float
            A float giving the B-band magnitude of a star.

        `"vmag"` : float
            A float giving the V-band magnitude of a star.

    """

    #
    # Check parameters
    #

    check_parameter("_query_simbad", "info", info, ["str", "list", "dict"])

    warnings.filterwarnings('error')

   #
   # Load the standard information based on the type of input
   # 

    if isinstance(info, str):

        #
        # user has passed the name of the standard
        #

        Simbad.add_votable_fields("sp_type", "B", "V")

        try:

            table = Simbad.query_object(info)

        except NoResultsWarning as w:

            message = (
                f'Standard name {info} was not found in SIMBAD; provide '
                "the correct name, correct coordinates, or a dictionary "
                'containing keys "id", "sptype", "bmag", '
                f'and "vmag": {w}')

            raise pySpextoolError(message)


        name = table["main_id"][0]
        sptype = table["sp_type"][0]
        vmag = float(table["V"][0])
        bmag = float(table["B"][0])


    if isinstance(info, list):

        #
        # user has passed the coordinates of the standard
        #

        Simbad.add_votable_fields("main_id", "sp_type", "B", "V")

        c = SkyCoord(info[0], info[1],unit=(u.hourangle, u.deg))

        try:

            table = Simbad.query_region(c, radius="0d1m0s")


        except NoResultsWarning as w:

            message = (
                f"Standard coordinates {info} were not found in SIMBAD; "
                "provide the correct name, correct coordinates, or a "
                'dictionary containing keys "id", "sptype", "bmag", '
                f'and "vmag": {w}'
            )

            raise pySpextoolError(message)

        name = table["main_id"][0]
        sptype = table["sp_type"][0]
        vmag = float(table["V"][0])
        bmag = float(table["B"][0])

    if isinstance(info, dict):

        #
        # user has passed the standard information as needed
        #

        name = info["id"]
        sptype = info["sptype"]
        vmag = float(info["vmag"])
        bmag = float(info["bmag"])

    #
    # Store the results and return
    #

    dictionary = {"name":name,
                  "sptype":sptype,
                  "vmag":float(f'{vmag:3f}'),
                  "bmag":float(f'{bmag:3f}')}

    warnings.resetwarnings()
    return dictionary
