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
    raw_path:str=None,
    cal_path:str=None,
    proc_path:str=None,
    qa_path:str=None,
    verbose:bool=False,
    qa_show:bool=False,
    qa_scale:float=1.0,
    qa_block:bool=False,
    qa_write:bool=False,
    qa_extension:str=None):

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

    qa_path : str, optional
        If qa_write is True, this is the path where the files will be written.
        Default: current working directory
    
    verbose : bool, default = True
        verbose = True sets the logging level to DEBUG
            Lots of information will be printed to the screen.
        verbose = False sets the logging level to INFO
            Only important information will be printed to the screen.

    qa_show : {False, True}, optional
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.

    qa_scale : float, default=1.0
        A scale factor that resizes the plots shown on the screen by a factor
        of 'qashow_scale`. 
    
    qa_block : {False, True}
        Set to True to stop the workflow after each plot is shown.
    
    qa_write : {True, False}, optional
        True: Save the quality assurance plots to disk.
        False: Do not save the quality assurance plots to disk.
        Set the path with `qa_path`.
        Set extension with `qa_extension`.

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

    set_paths(raw_path, cal_path, proc_path, qa_path)

    logging.info("Paths set")

    # Set the quality assurance settings

    set_qa_state(qa_show, qa_scale, qa_block, qa_write, qa_extension)

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


def set_paths(raw_path, cal_path, proc_path, qa_path):

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

    qa_path : str, optional
        The path to the quality assurance plots directory.
  
    Returns
    -------
    None

    """

    #
    # Check parameters
    #

    check_parameter("set_parameters", "raw_path", raw_path,
                    ['NoneType', "str"])

    check_parameter("set_parameters", "cal_path", cal_path, 
                    ['NoneType', "str"])
                    
    check_parameter("set_parameters", "proc_path", proc_path, 
                    ['NoneType', "str"])
                    
    check_parameter("set_parameters", "qa_path", qa_path, 
                    ['NoneType', "str"])
                    
    #
    # Load the paths
    #

    # Get the current working directory in case it is needed.
    
    cwd = os.path.abspath(os.getcwd())
    
    # Modify the paths based on the user requests.
    
    if raw_path is not None:

        raw_path = check_path(raw_path, make_absolute=True)
        setup.state["raw_path"] = raw_path
        logging.debug(f"Set raw_path to {raw_path}")

    else:

        setup.state["raw_path"] = cwd

    if cal_path is not None:

        cal_path = check_path(cal_path, make_absolute=True)
        setup.state["cal_path"] = cal_path
        logging.debug(f"Set cal_path to {cal_path}")

    else:

        setup.state["cal_path"] = cwd
                                        
    if proc_path is not None:

        proc_path = check_path(proc_path, make_absolute=True)
        setup.state["proc_path"] = proc_path
        logging.debug(f"Set proc_path to {proc_path}")
                                        
    else:

        setup.state["proc_path"] = cwd

    if qa_path is not None:

        qa_path = check_path(qa_path, make_absolute=True)
        setup.state["qa_path"] = qa_path
        logging.debug(f"Set qa_path to {qa_path}")
                                        
    else:

        setup.state["qa_path"] = cwd

        
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

    check_parameter("set_instrument", "instrument_name", instrument_name,
                    "str", possible_values=["spex", "uspex"])

    setup.state["instrument"] = instrument_name

    #
    # Set the package path
    #

    setup.state["package_path"] = str(files("pyspextool"))

    #
    # Check to make sure the instrument path exists
    #

    instrument_data_path = os.path.join(
        setup.state["package_path"], "instrument_data",
        instrument_name + "_dir")

    data_path = os.path.join(setup.state["package_path"], "data")

    check_path(instrument_data_path)
    check_path(data_path)

    setup.state["instrument_path"] = instrument_data_path

    #
    # Now get the instrument file and load
    #

    instrument_info_file = os.path.join(instrument_data_path,
                                        instrument_name + ".dat")

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


def set_qa_state(qa_show:bool,
                 qa_scale:float,
                 qa_block:bool,
                 qa_write:bool,
                 qa_extension:str,):

    """
    To set the quality assurance plot settings.
    
    
    Parameters
    ----------
    
    qa_show : {False, True}
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.

    qashow_scale : float, default=1.0
        A scale factor that resizes the plots shown on the screen by a factor
        of 'qashow_scale`. 

    qashow_block : {False, True}
        Set to True to stop the workflow after each plot is shown.
    
    qa_write : {True, False}
        True: Save the quality assurance plots to disk.
        False: Do not save the quality assurance plots to disk.
        Set the path with `qa_path`.
        Set extension with `qa_extension`.
        Default: True

    qa_extension : {'pdf', 'png'}, optional
        if qa_write is True, this is the file extension used for the files.
        Default: pdf
        Options: pdf, png

    Returns
    -------
    None

    """

    #
    # Check parameters
    #
    
    check_parameter("set_qa_state", "qa_show", qa_show, ["NoneType", "bool"])

    check_parameter("set_qa_state", "qa_scale", qa_scale, 'float')
    
    check_parameter("set_qa_state", "qa_block", qa_block, ["NoneType", "bool"])
    
    check_parameter("set_qa_state", "qa_write", qa_write, ["NoneType", "bool"])

    check_parameter("set_qa_state", "qa_extensioan", qa_extension,
        ["NoneType", "str"], possible_values=setup.state["qa_extensions"])

    #
    # Set the values
    #
    
    setup.state["qa_show"] = qa_show

    setup.state["qa_block"] = qa_block

    setup.state["qa_scale"] = qa_scale

    setup.state["qa_write"] = qa_write

    if qa_extension is not None:

        setup.state["qa_extension"] = qa_extension

    else:

        setup.state["qa_extension"] = setup.state["qa_extensions"][0]

#
#
#
#    msg = f"""
#    QA Setup
#    ----------------
#    QA Show: {setup.state['qa_show']}
#
#    QA Write: {setup.state['qa_write']}
#    QA Path: {setup.state['qa_path']}
#    QA Extension: {setup.state['qa_extension']}
#    """
#    logging.debug(msg)

    return
