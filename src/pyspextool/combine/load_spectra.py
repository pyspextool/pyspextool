from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.io.check import check_parameter, check_qakeywords, check_file

from pyspextool.io.files import files_to_fullpath
from pyspextool.io.read_spectra_fits import read_spectra_fits


def load_spectra(files:str | list,
                 output_name:str,
                 input_extension:str='.fits',
                 verbose:bool=None,
                 qa_show:bool=None,
                 qa_showscale:float | int=None,
                 qa_showblock:bool=None,
                 qa_write:bool=None):


    """
    To load spectra to later be combined.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.

        e.g. ['spc', '1-2']

    output_name : str
        The output file name sans the suffix.

    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does NOT affect plots written
        to disk.  Set to None to default to setup.state['qa_scale'].    
    


    Returns
    -------

    """

    #
    # Check the parameters and QA keywords
    #

    check_parameter('load_spectra', 'files', files, ['str', 'list'])

    check_parameter('load_spectra', 'output_name', output_name, 'str')    

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Create the file names
    #

    results = files_to_fullpath(setup.state['proc_path'],
                                files,
                                setup.state['nint'],
                                '',
                                input_extension)

    input_files = results[0]
    file_readmode = results[1]

    check_file(input_files)

    combine.state['input_files'] = input_files
    combine.state['nfiles'] = len(input_files)

    #
    # Read the first file and store useful things
    #

    first_spectra, info = read_spectra_fits(combine.state['input_files'][0])

    combine.state['npixels'] = np.size(first_spectra[0,0,:])
    combine.state['module'] = info['module']    
    combine.state['orders'] = info['orders']
    combine.state['norders'] = len(info['orders'])
    combine.state['napertures'] = info['napertures']
    combine.state['xlabel'] = info['lxlabel']
    combine.state['xunits'] = info['xunits']    
    
    
    
