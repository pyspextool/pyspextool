import importlib
from pyspextool.extract import config # sets initial state dictionary

def make_flat(files, output_name, prefix='flat-', suffix='.[ab]',
              extension='.fits*', input_method='index', normalize=True,
              verbose=True, overwrite=True, qafile=False):

    """
    To create a (normalized) pyspextool flat field file.

    A wrapper to create a flat field for a specific instrument.

    Parameters
    ----------
    files : str or list of str.
        See input_method parameter for more information.
        
    output_name : str
        The filename of the flat field image to written to disk.      

    prefix : str, default='flat-', optional
        The prefix of the FITS file name (i.e. the stuff before the 
        file numbers.)

    suffix : str, default='.[ab]', optional
        The suffix of the FITS file name (i.e. the stuff after the
        file numbers but before the extension.)

    extension : str, default='.fits*', optional
        The file extension

    input_method : {'index', 'filename'}, optional
        `index`: `files` is a str of file numbers, e.g., '1-3,5-10'.

        `filename`: a str or list of strings with full file names, e.g.
            ['flat-00001.a.fits',''flat-00002.a.fits']

    normalize : {True, False}, optional
        Set to True to normalize the orders.

    verbose : {True, False}, optional
        Set to True for command line updates during execution. 

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.

    qafile : {True, False}, optional
        Set to create the quality assurance plot.

    Returns
    -------
    None
        Writes a FITS file to disk.

    Examples
    --------
    later

    """

    #
    # We do not need to check parameters because the make_flat functions do
    # that.
    #

    #
    # Load the instrument module
    #
    

    module = 'pyspextool.instrument_data.'+config.state['instrument_name']+\
        '_dir.'+config.state['instrument_name']

    instr = importlib.import_module(module)

    #
    # Make the flat
    #
    
    instr.make_flat(files, output_name, prefix=prefix, suffix=suffix,
                    extension=extension, input_method=input_method,
                    normalize=normalize, verbose=verbose, overwrite=overwrite,
                    qafile=qafile)

