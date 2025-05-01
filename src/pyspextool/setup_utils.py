import os
from astropy.io import fits
import numpy as np
import logging
import pooch

from pyspextool import config as setup
from pyspextool.io.read_instrument_file import read_instrument_file
from pyspextool.io.check import check_parameter, check_path, check_file
from pyspextool.pyspextoolerror import pySpextoolError
from importlib.resources import files  # Python 3.10+
from importlib.metadata import version, PackageNotFoundError

# TODO:  test logging works as expected. run some commands in the REPL

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

try:
    __version__ = version("pyspextool")
except PackageNotFoundError:
    # package is not installed
    pass


mishu = pooch.create(
        path=pooch.os_cache("pyspextool"),
        base_url="https://pyspextool.s3.us-east-1.amazonaws.com/",
        registry={
            "uspex_lincorr.fits": "9ba8c54dc9de08aab81a67cd37ee7d5c6aaad2aec6a13537cc8c412b896aca58",
            "uspex_bias.fits": "d1dbbffc882123de5f3e877ca14dbc6740e2e751d1c31d0647583305c6163cc6",
            "spex_lincorr.fits": "47fcbd6b854f1b80fc65615978dffdcfa793af24d12b3ce3b199efae6d78040f",
            "Vega50000.fits": "517f38feaaabe35443fcbd9a2670085b61af0e7dfd05a28a6c3c4f79ed7d7737",
        },
    )

def set_version():
    setup.state["version"] = __version__
    logger.debug(f"Version set to {setup.state['version']}")


def pyspextool_setup(instrument=setup.state["instruments"][0],
                     raw_path:str=None,
                     cal_path:str=None,
                     proc_path:str=None,
                     qa_path:str=None,
                     search_extension:str=setup.state["search_extensions"][0],
                     verbose:bool=True,
                     qa_show:bool=False,
                     qa_showscale:float=1.0,
                     qa_showblock:bool=False,
                     qa_write:bool=False,
                     qa_extension:str=setup.state["qa_extensions"][0]):

    """
    Set the pyspextool instrument, paths, and quality assurance settings

    Parameters
    ----------
    instrument : setup.state['instruments']
        The name of the instrument.
        Must be one of options provided in config.setup['instruments'].

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

    search_extension : setup.state['search_extensions']
        The file extension used to search for files in `raw_path`.
    
    verbose : bool, default True
        verbose = True sets the logging level to INFO
            Lots of information will be printed to the screen.
        verbose = False sets the logging level to ERROR
            Only important information will be printed to the screen.

    qa_show : {False, True}, optional
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.

    qa_showscale : float, default 1.0
        A scale factor that resizes the plots shown on the screen by a factor
        of 'qashow_showscale`. 
    
    qa_showblock : {False, True}
        Set to True to stop the workflow after each plot is shown.
    
    qa_write : {True, False}, optional
        True: Save the quality assurance plots to disk.
        False: Do not save the quality assurance plots to disk.
        Set the path with `qa_path`.
        Set extension with `qa_extension`.

    qa_extension : setup.state['qa_extensions']
        if qa_write is True, this is the file extension used for the files.

    Returns
    -------
    None
    Loads data into:

        setup.state["raw_path"]
        setup.state["cal_path"]
        setup.state["proc_path"]
        setup.state["qa_path"]
        setup.state["package_path"]
        setup.state["instrument_path"]
        setup.state["irtf"]
        setup.state["suffix"]
        setup.state["nint"]
        setup.state["extract_keywords"]
        setup.state["combine_keywords"]
        setup.state["telluric_keywords"]
        setup.state["lincormax"]
        setup.state["linearity_info"]
        setup.state["raw_bad_pixel_mask"]
        setup.state["pyspextool_keywords"]
        setup.state["qa_show"]
        setup.state["qa_showblock"]
        setup.state["qa_showscale"]
        setup.state["qa_write"]
        setup.state["qa_extension"]
       
    """

    #
    # Check parameters
    #

    check_parameter('pyspextool_setup', 'instrument', instrument, 'str',
                    possible_values=setup.state['instruments'])

    check_parameter('pyspextool_setup', 'raw_path', raw_path,
                    ['NoneType', 'str'])

    check_parameter('pyspextool_setup', 'cal_path', cal_path,
                    ['NoneType', 'str'])

    check_parameter('pyspextool_setup', 'proc_path', proc_path, 
                    ['NoneType', 'str'])

    check_parameter('pyspextool_setup', 'qa_path', qa_path, 
                    ['NoneType', 'str'])

    check_parameter('pyspextool_setup', 'search_extension', search_extension, 
                    'str', possible_values=setup.state['search_extensions'])

    check_parameter('pyspextool_setup', 'qa_show', qa_show, 'bool')

    check_parameter('pyspextool_setup', 'qa_showscale', qa_showscale,
                    ['float','int'])

    check_parameter('pyspextool_setup', 'qa_showblock', qa_showblock, 'bool')

    check_parameter('pyspextool_setup', 'qa_write', qa_write, 'bool')

    check_parameter('pyspextool_setup', 'qa_extensioan', qa_extension,
        'str', possible_values=setup.state['qa_extensions'])

    #
    # Report what you are doing
    #

    message = ' pySpextool Setup'
    logging.info(message+'\n'+'-'*(len(message)+5)+'\n')

    #
    # Store the search extension
    #

    setup.state['search_extension'] = search_extension

    #
    # Set up verbose scale and logging
    #

    if verbose is True:

        logging.getLogger().setLevel(logging.INFO)
        setup.state["verbose"] = True

    elif verbose is False:

        logging.getLogger().setLevel(logging.ERROR)
        setup.state["verbose"] = False

    logging.info(
        f" Verbose set to {setup.state['verbose']}. \n"
        f" Logging level set to {logging.getLogger().getEffectiveLevel()}"
    )

    #
    # Set the instrument
    #

    set_instrument(instrument)

    logging.info(f" Instrument set to {setup.state['instrument']}")

    #
    # Set the paths
    #

    set_paths(raw_path, cal_path, proc_path, qa_path)

    logging.info(" Paths set")

    # Set the quality assurance settings

    set_qa_state(qa_show, qa_showscale, qa_showblock, qa_write, qa_extension)

    logging.info(" QA settings set")

    # Set the version number
    set_version()

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

    Version: {setup.state['version']}
    """

    logging.debug(msg)

    return  # setup.state


def set_paths(raw_path:str,
              cal_path:str,
              proc_path:str,
              qa_path:str):

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
    # Load the paths
    #

    # Get the current working directory in case it is needed.
    
#    cwd = os.path.abspath(os.getcwd())
    
    # Modify the paths based on the user requests.

    if raw_path is not None:

        raw_path = check_path(raw_path, make_absolute=True)

    setup.state["raw_path"] = raw_path
    logging.debug(f"Set raw_path to {raw_path}")

    if cal_path is not None:

        cal_path = check_path(cal_path, make_absolute=True)

    setup.state["cal_path"] = cal_path
    logging.debug(f"Set cal_path to {cal_path}")
                                    
    if proc_path is not None:

        proc_path = check_path(proc_path, make_absolute=True)
                                        
    setup.state["proc_path"] = proc_path
    logging.debug(f"Set proc_path to {proc_path}")


    if qa_path is not None:

        qa_path = check_path(qa_path, make_absolute=True)

    setup.state["qa_path"] = qa_path
    logging.debug(f"Set qa_path to {qa_path}")
    
    #
    # Now ensure that `raw_path` does not equal any of the other paths
    #

    if setup.state['raw_path'] is not None:

        if setup.state['raw_path'] == setup.state['cal_path']:
            
            message = 'The parameter `raw_path` cannot be the same as the '+\
                'parameter `cal_path`.'
            raise pySpextoolError(message)

        if setup.state['raw_path'] == setup.state['proc_path']:

            message = 'The parameter `raw_path` cannot be the same as the '+\
                'parameter `proc_path`.'
            raise pySpextoolError(message)

        if setup.state['raw_path'] == setup.state['qa_path']:
            
            message = 'The parameter `raw_path` cannot be the same as the '+\
                'parameter `qa_path`.'
            raise pySpextoolError(message)
        
    return


def set_instrument(instrument_name: str):

    """
    Set the instrument.

    Parameters
    ----------
    instrument_name : setup.state['instruments']
        The name of the instrument. 

    Returns
    -------
    None

    """

    setup.state["instrument"] = instrument_name

    #
    # Set the package path
    #

    setup.state["package_path"] = str(files("pyspextool"))

    #
    # Check to make sure the instrument path exists
    #

    instrument_data_path = os.path.join(
        setup.state["package_path"], "instruments",instrument_name)

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
    setup.state["extract_keywords"] = instrument_info["EXTRACT_KEYWORDS"]
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
                 qa_showscale:float,
                 qa_showblock:bool,
                 qa_write:bool,
                 qa_extension:str,):

    """
    To set the quality assurance plot settings.
    
    
    Parameters
    ----------
    
    qa_show : {False, True}
        True: Display the quality assurance plots to the screen.
        False: Do not display the quality assurance plots to the screen.

    qa_showscale : float, default=1.0
        A scale factor that resizes the plots shown on the screen by a factor
        of 'qashow_scale`. 

    qa_showblock : {False, True}
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
    # Set the values
    #
    
    setup.state["qa_show"] = qa_show

    setup.state["qa_showblock"] = qa_showblock

    setup.state["qa_showscale"] = qa_showscale

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