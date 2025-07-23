from pyspextool.merge.load_spectra import load_spectra
from pyspextool.merge.merge_orders import merge_orders
from pyspextool.merge.write_spectrum import write_spectrum
from pyspextool.io.check import check_parameter, check_qakeywords


def merge(
        file: str,
        outputfile_root:str,
        merge_apertures:int=None,
        verbose:bool=None,
        qa_show:bool=None,
        qa_showscale:float=None,
        qa_showblock:bool=None,
        qa_write:bool=None):

    """
    Parameters
    ----------
    file : str
        The name of the pySpextool FITS file containing the multi-order spectra.

    outputfile_root : str
        A string giving the root of the output file name, e.g. 'Wolf359'. 

    merge_apertures : int, list, default None
        The aperture number or a list of aperture numbers to merge.  
        The apetures are index starting with 1.  

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].

    Returns
    -------
    None
        Writes a pySpextool FITS file to disk.

    """

    #
    # Check parameters and qa keywords
    #

    check_parameter("merge", "file", file, "str")

    check_parameter("merge", "outputfile_root", outputfile_root, "str")

    check_parameter("merge", "merge_apertures", merge_apertures, 
                    ["NoneType", "int", "list"])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Load the spectrum
    #

    load_spectra(file,
                 outputfile_root,
                 verbose=qa['verbose'])

    # 
    # Merge the orders
    #
    
    merge_orders(merge_apertures=merge_apertures,
                 verbose=qa['verbose'],
                 qa_show=qa['show'],
                 qa_showscale=qa['showscale'],
                 qa_showblock=qa['showblock'],
                 qa_write=qa['write'])

    #
    # Write the file to disk
    #

    write_spectrum(verbose=qa['verbose'],
                   qa_show=qa['show'],
                   qa_showscale=qa['showscale'],
                   qa_showblock=qa['showblock'],
                   qa_write=qa['write'])
