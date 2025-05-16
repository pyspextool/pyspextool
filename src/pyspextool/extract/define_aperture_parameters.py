import numpy as np
import logging
from os.path import join

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.pyspextoolerror import pySpextoolError


def define_aperture_parameters(
    aperture_radii: int | float | list,
    bg_annulus: list = None,
    bg_regions: str = None,
    bg_fit_degree: int = 0,
    psf_radius: int | float = None,
    verbose: bool = None,
    qa_show: bool = None,
    qa_showscale: float | int = None,
    qa_showblock: bool = None,
    qa_write: bool = None,
):
    """
    To define the extraction parameters

    Parameters
    ----------
    aperture_radii : int, float, or list

        If len(`aperture_radii`) is 1, then this value will be used
        for all apertures.  if len(`aperture_radii`) is not 1, then
        it must be of the same length as the number of apertures previously
        defined.

    bg_annulus : list, default None
        A (2, ) list giving the start radius of the 1D annulus and the
        width of the 1D annulus.

    bg_regions : str, default None
        A string giving the background regions, e.g. '1-2, 14-15'.

    bg_fit_degree : int, default 1
        The polynomial degree of the fit to the background.

    psf_radius : int or float, default None
        The radius used to to normalize the PSF as part of optimal extraction.

    Returns
    -------
    None



    """

    #
    # Check if we can proceed
    #

    if extract.state["trace_done"] is False:

        message = "Previous steps not complete.  Please run trace_apertures.py"
        raise pySpextoolError(message)

    #
    # Check common parameters
    #

    check_parameter(
        "define_aperture_parameters",
        "aperture_radii",
        aperture_radii,
        ["int", "float", "list"],
    )

    check_parameter(
        "define_aperture_parameters", "bg_annulus", bg_annulus, ["list", "NoneType"]
    )

    check_parameter(
        "define_aperture_parameters", "bg_regions", bg_regions, ["str", "NoneType"]
    )

    check_parameter(
        "define_aperture_parameters", "bg_fit_degree", bg_fit_degree, ["int"], [1, 2]
    )

    check_parameter(
        "define_aperture_parameters",
        "psf_radius",
        psf_radius,
        ["int", "float", "NoneType"],
    )

    check_parameter(
        "define_aperture_parameters", "verbose", verbose, ["NoneType", "bool"]
    )

    check_parameter(
        "define_aperture_parameters", "qa_write", qa_write, ["NoneType", "bool"]
    )

    check_parameter(
        "define_aperture_parameters", "qa_show", qa_show, ["NoneType", "bool"]
    )

    check_parameter(
        "define_aperture_parameters",
        "qa_showscale",
        qa_showscale,
        ["int", "float", "NoneType"],
    )

    check_parameter(
        "define_aperture_parameters", "qa_showblock", qa_showblock, ["NoneType", "bool"]
    )

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write,
    )

    logging.info(" Defining aperture parameters.")

    #
    # Check user inputs against a bunch of constraints
    #

    if bg_annulus is not None and bg_regions is not None:

        message = "Cannot pass `bg_annulus` and `bg_regions` at the same time."
        raise pySpextoolError(message)

    #
    # Get aperture_radii into the right format and test to make sure the
    # user is passing the right number of radii.
    #

    if isinstance(aperture_radii, list):

        # The user is passing multiple radii.

        aperture_radii = np.tile(aperture_radii, (extract.state["norders"], 1))

        naps = np.shape(aperture_radii)[1]

        if naps != extract.state["naps"]:

            message = (
                "Number of radii in `aperture_radii` must equal the "
                + "number of apertures.  And "
                + str(naps)
                + " .ne. "
                + str(extract.state["naps"])
                + "."
            )
            raise pySpextoolError(message)

    else:

        # The user passed a single aperture.  So use for all apertures.

        aperture_radii = np.full(extract.state["naps"], aperture_radii)

    # Now tile the results in order to have a (norders, naps) array.

    aperture_radii = np.tile(aperture_radii, (extract.state["norders"], 1))

    #
    # Now check the psf_radius

    bgtest = bg_annulus is None and bg_regions is None
    if psf_radius is not None and bgtest:

        message = (
            "`bg_annulus` or `bg_region` is required if " + "`psf_radius` is passed."
        )
        raise pySpextoolError(message)

    if psf_radius is not None:

        # Make sure it is larger than aperture_radii

        if np.sum(psf_radius < aperture_radii) != 0:
            message = "`psf_radius` must be >= `aperture_radii`."
            raise pySpextoolError(message)

    #
    # Store the results
    #

    extract.state["aperture_radii"] = aperture_radii
    extract.state["bg_annulus"] = bg_annulus
    extract.state["bg_regions"] = bg_regions
    extract.state["bg_fitdegree"] = bg_fit_degree
    extract.state["psf_radius"] = psf_radius

    #
    # Do the QA plots
    #

    if qa["show"] is True:

        plot_profiles(
            extract.state["profiles"],
            extract.state["slith_arc"],
            extract.state["doorders"],
            aperture_positions=extract.state["aperture_positions"],
            aperture_signs=extract.state["average_aperturesigns"],
            aperture_radii=aperture_radii,
            psf_radius=psf_radius,
            bg_regions=bg_regions,
            bg_annulus=bg_annulus,
            plot_number=setup.plots["profiles"],
            profilestack_max=setup.plots["profilestack_max"],
            profile_size=setup.plots["profile_size"],
            font_size=setup.plots["font_size"],
            showscale=qa["showscale"],
            showblock=qa["showblock"],
        )

    if qa["write"] is True:

        filename = (
            extract.state["qafilename"]
            + "_extractedprofiles"
            + setup.state["qa_extension"]
        )
        fullpath = join(setup.state["qa_path"], filename)

        plot_profiles(
            extract.state["profiles"],
            extract.state["slith_arc"],
            extract.state["doorders"],
            aperture_positions=extract.state["aperture_positions"],
            aperture_signs=extract.state["average_aperturesigns"],
            aperture_radii=aperture_radii,
            psf_radius=psf_radius,
            bg_regions=bg_regions,
            bg_annulus=bg_annulus,
            profilestack_max=setup.plots["profilestack_max"],
            profile_size=setup.plots["profile_size"],
            font_size=setup.plots["font_size"],
            output_fullpath=fullpath,
        )

    #
    # Set the done variable
    #

    extract.state["parameters_done"] = True
