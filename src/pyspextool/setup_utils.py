import os
from astropy.io import fits
import numpy as np
import logging
from pyspextool import config as setup
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter, check_path, check_file
from importlib.resources import files  # Python 3.10+

# TODO:  test logging works as expected. run some commands in the REPL

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

def pyspextool_setup(
    instrument=setup.state["instruments"][0],
    raw_path: str = None,
    cal_path: str = None,
    proc_path: str = None,
    verbose: bool = False,
    qa_show: bool = None,
    qa_write: bool = None,
    qa_path: str = None,
    qa_extension: str = None,
):
    """
    Set the pyspextool instrument, paths, and quality assurance settings

    Parameters
    ----------
    instrument : str, optional
        Default: 'uspex'
        The name of the instrument.
        Must be one of options provided in config.setup['instruments'].
        Currently, only 'spex' and 'uspex' are supported.

    raw_path : str, optional.
        The path to the raw data.
        Default: current working directory

    cal_path : str, optional
        The path where calibration files will be written.
        Default: current working directory

    proc_path : str, optional
        The path to where the processed files will be written.
        Default: current working directory

    verbose : bool, default = True
        verbose = True sets the logging level to DEBUG
            Lots of information will be printed to the screen.
        verbose = False sets the logging level to INFO
            Only important information will be printed to the screen.

    qa_show : {True, False}, optional
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.
        Default: False

    qa_write : {True, False}, optional
        True: Save the quality assurance plots to disk.
        False: Do not save the quality assurance plots to disk.
        Set the path with `qa_path`.
        Set extension with `qa_extension`.
        Default: True

    qa_path : str, optional
        If qa_write is True, this is the path where the files will be written.
        Default: current working directory

    qa_extension : {'pdf', 'png'}, optional
        if qa_write is True, this is the file extension used for the files.
        Default: pdf
        Options: pdf, png

    Returns
    -------
    None

    """
    # Set up verbose scale and logging

    if verbose is True:
        logging.getLogger().setLevel(logging.INFO)
        setup.state["verbose"] = True
    elif verbose is False:
        logging.getLogger().setLevel(logging.ERROR)
        setup.state["verbose"] = False

    logging.info(f"Verbose set to {setup.state['verbose']}")
    # Set the instrument

    if instrument is not None:
        set_instrument(instrument)

    logging.info(f"Instrument set to {setup.state['instrument']}")
    # Set the paths

    paths = {
        "raw_path": raw_path,
        "cal_path": cal_path,
        "proc_path": proc_path,
    }

    set_paths(paths)

    logging.info("Paths set")
    print("paths: ", paths)

    # Set the quality assurance settings

    set_qa_state(
        qa_show=qa_show, qa_write=qa_write, qa_path=qa_path, qa_extension=qa_extension
    )

    logging.info("QA settings set")

    msg = f"""
    Pyspextool Setup
    ----------------
    Instrument: {setup.state['instrument']}

    Rawpath: {setup.state['raw_path']}
    Calpath: {setup.state['cal_path']}
    Procpath: {setup.state['proc_path']}

    QA Plot: {setup.state['qa_show']}

    QA Write: {setup.state['qa_write']}
    QA path: {setup.state['qa_path']}
    QA Extension: {setup.state['qa_extension']}

    Verbose: {setup.state['verbose']}
    """

    logging.debug(msg)

    return  # setup.state


def set_paths(paths: dict = None):
    """
    Set the pyspextool paths

    Parameters
    ----------
    raw_path : str, optional
        The path to the raw directory.

    cal_path : str, optional
        The path to the calibration directory.

    proc_path : str, optional
        The path to the processed directory.


    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter(
        "set_parameters", "raw_path", paths["raw_path"], ["str", "NoneType"]
    )

    check_parameter(
        "set_parameters", "cal_path", paths["cal_path"], ["str", "NoneType"]
    )

    check_parameter(
        "set_parameters", "proc_path", paths["proc_path"], ["str", "NoneType"]
    )

    #
    # Load the .pyspextool file if it exists.
    #

    # Change to use the current working directory rather than HOME
    # home_path = os.path.expanduser('~')
    # file_name = '.pyspextool_' + setup.state['instrument'] + '.dat'
    # full_path = os.path.join(home_path, file_name)

    # Does the file exist?

    # if os.path.isfile(full_path) is True:

    #     # Yes.  Load the previous paths.

    #     f = open(full_path, 'r')
    #     paths = []

    #     for line in f:
    #         paths.append(line.strip())

    #     setup.state['raw_path'] = paths[0]
    #     setup.state['cal_path'] = paths[1]
    #     setup.state['proc_path'] = paths[2]
    #     setup.state['qa_path'] = paths[3]

    # else:

    #     # No.  Use the current working directory

    # Should only use the current working directory
    # if the user has not defined the paths
    cwd = os.path.abspath(os.getcwd())

    # #
    # Now let's modify the paths based on the user requests.
    #

    if paths["raw_path"] is not None:
        paths["raw_path"] = check_path(paths["raw_path"], make_absolute=False)
        setup.state["raw_path"] = paths["raw_path"]
        logging.debug(f"Set raw_path to {paths['raw_path']}")
    else:
        setup.state["raw_path"] = cwd

    if paths["cal_path"] is not None:
        try:
            paths["cal_path"] = check_path(paths["cal_path"], make_absolute=False)
            logging.debug(f"Set cal_path to {paths['cal_path']}")
        except ValueError as e:
            # os.mkdir(paths["cal_path"])
            # paths["cal_path"] = check_path(paths["cal_path"], make_absolute=True)
            logging.error(f"Can't verify cal_path directory {paths['cal_path']}")
            raise (e)

        setup.state["cal_path"] = paths["cal_path"]
    else:
        setup.state["cal_path"] = cwd

    if paths["proc_path"] is not None:
        paths["proc_path"] = check_path(paths["proc_path"], make_absolute=False)
        setup.state["proc_path"] = paths["proc_path"]
    else:
        setup.state["proc_path"] = cwd

    #
    # Now write the paths to the user home directory
    #

    # dat_file_name = f".pyspextool_{setup.state['instrument']}.dat"
    # f = open(os.path.join(home_path, dat_file_name), 'w')
    # f.write('%s \n' % setup.state['raw_path'])
    # f.write('%s \n' % setup.state['cal_path'])
    # f.write('%s \n' % setup.state['proc_path'])
    # f.write('%s \n' % setup.state['qa_path'])
    # f.close()

    # logging.info(f'Created {dat_file_name} in {home_path}')

    return


def set_instrument(instrument_name: str):
    """
    Set the instrument.

    Parameters
    ----------
    instrument_name : str
        Default = 'uspex'
        The name of the instrument.  Must be one of
        config.setup['instruments'].
        Currently, only 'spex' and 'uspex' are supported.

    Returns
    -------
    None

    """

    if instrument_name is None:
        instrument_name = setup.state["instruments"][0]

    #
    # Check parameter and store results
    #

    check_parameter("set_instrument", "instrument_name", instrument_name, "str")

    setup.state["instrument"] = instrument_name

    #
    # Set the package path
    #

    setup.state["package_path"] = str(files("pyspextool"))

    #
    # Check to make sure the instrument path exists
    #

    instrument_data_path = os.path.join(
        setup.state["package_path"], "instrument_data", instrument_name + "_dir"
    )
    data_path = os.path.join(setup.state["package_path"], "data")

    check_path(instrument_data_path)
    check_path(data_path)

    setup.state["instrument_path"] = instrument_data_path

    #
    # Now get the instrument file and load
    #

    instrument_info_file = os.path.join(instrument_data_path, instrument_name + ".dat")

    check_file(instrument_info_file)

    instrument_info = read_instrument_file(instrument_info_file)

    if instrument_name in ["uspex", "spex"]:
        setup.state["irtf"] = True

    # Fill out the state variables

    setup.state["suffix"] = instrument_info["SUFFIX"]
    setup.state["nint"] = instrument_info["NINT"]
    setup.state["extract_keywords"] = instrument_info["XSPEXTOOL_KEYWORDS"]
    setup.state["combine_keywords"] = instrument_info["COMBINE_KEYWORDS"]
    setup.state["telluric_keywords"] = instrument_info["TELLURIC_KEYWORDS"]

    # Now store linearity numbers

    setup.state["lincormax"] = instrument_info["LINCORMAX"]
    setup.state["linearity_info"] = {"max": setup.state["lincormax"], "bit": 0}

    # Get the bad pixel mask

    bad_pixel_mask_file = os.path.join(
        instrument_data_path, setup.state["instrument"] + "_bdpxmk.fits"
    )

    check_file(bad_pixel_mask_file)

    setup.state["raw_bad_pixel_mask"] = fits.getdata(bad_pixel_mask_file)

    #
    # Grab the Spextool keywords
    #

    keywords_path = os.path.join(data_path, "pyspextool_keywords.dat")

    keywords = np.loadtxt(keywords_path, comments="#", dtype="str").tolist()

    setup.state["pyspextool_keywords"] = keywords

    msg = f"""
    Instrument Setup
    ----------------
    Instrument: {setup.state['instrument']}
    State: {setup.state['suffix']}
    NINT: {setup.state['nint']}
    Extract Keywords: {setup.state['extract_keywords']}
    Combine Keywords: {setup.state['combine_keywords']}
    Telluric Keywords: {setup.state['telluric_keywords']}
    Linearity Max: {setup.state['lincormax']}
    Bad Pixel Mask: {bad_pixel_mask_file}
    """

    logging.debug(msg)

    return


def set_qa_state(
    qa_show: bool = None,
    qa_write: bool = None,
    qa_path: str = None,
    qa_extension: str = None,
):
    """
    qa_show : {True, False}, optional
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.
        Default: False

    qa_write : {True, False}, optional
        True: Save the quality assurance plots to disk.
        False: Do not save the quality assurance plots to disk.
        Set the path with `qa_path`.
        Set extension with `qa_extension`.
        Default: True

    qa_path : str, optional
        If qa_write is True, this is the path where the files will be written.
        Default: current working directory

    qa_extension : {'pdf', 'png'}, optional
        if qa_write is True, this is the file extension used for the files.
        Default: pdf
        Options: pdf, png

    Returns
    -------
    None

    """
    check_parameter(
        "set_parameters",
        "qa_extension",
        qa_extension,
        ["NoneType", "str"],
        possible_values=setup.state["qa_extensions"],
    )

    check_parameter("set_parameters", "qa_show", qa_show, ["NoneType", "bool"])

    check_parameter("set_parameters", "qa_write", qa_write, ["NoneType", "bool"])

    cwd = os.path.abspath(os.getcwd())

    if qa_path is not None:
        try:
            qa_path = check_path(qa_path, make_absolute=True)
            logging.debug(f"Set qa_path to {qa_path}")
        except ValueError as e:
            # os.mkdir(qa_path)
            # check_path(qa_path, make_absolute=True)
            logging.error(f"Can't verify qa_path directory {qa_path}")
            raise (e)

        setup.state["qa_path"] = qa_path
    else:
        setup.state["qa_path"] = cwd

    # Set the qa extension filetype

    setup.state["qa_extension"] = qa_extension

    #
    # Check defaults
    #

    if qa_show is not None:
        setup.state["qa_show"] = qa_show

    if qa_write is not None:
        setup.state["qa_write"] = qa_write

    if qa_extension is not None:
        setup.state["qa_extension"] = qa_extension

    else:
        setup.state["qa_extension"] = setup.state["qa_extensions"][0]

    msg = f"""
    QA Setup
    ----------------
    QA Show: {setup.state['qa_show']}

    QA Write: {setup.state['qa_write']}
    QA Path: {setup.state['qa_path']}
    QA Extension: {setup.state['qa_extension']}
    """
    logging.debug(msg)

    return
