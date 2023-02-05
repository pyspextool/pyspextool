import importlib
from pyspextool.extract import config # sets initial state dictionary

def make_wavecal(files, flat_file, output_name, prefix='arc-', suffix='.[ab]',
                 extension='.fits*', input_method='index',
                 use_stored_solution=False, verbose=True, qafile_shift=True,
                 qafile_findlines=True, qafile_fitlines=True, overwrite=True):

    """
    TBD
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
    # Make the wavecal file
    #
    
    instr.make_wavecal(files, flat_file, output_name, prefix=prefix,
                       suffix=suffix, extension=extension,
                       input_method=input_method,
                       use_stored_solution=use_stored_solution,
                       qafile_shift=qafile_shift, verbose=verbose, 
                       qafile_findlines=qafile_findlines,
                       qafile_fitlines=qafile_fitlines, overwrite=overwrite)

