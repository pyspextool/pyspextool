import logging

from pyspextool import config as setup
from pyspextool.telluric import config as telluric
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.check import check_parameter, check_keywords
from pyspextool.telluric.load_spectra import load_spectra
from pyspextool.telluric.prepare_line import prepare_line
from pyspextool.telluric.get_radialvelocity import get_radialvelocity
from pyspextool.telluric.get_kernels import get_kernels
from pyspextool.telluric.make_telluric_spectra import make_telluric_spectra
from pyspextool.telluric.correct_spectra import correct_spectra
from pyspextool.telluric.write_spectra import write_spectra

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def telluric_correction(object_file:str,
                        standard_file:str,
                        standard_info:list | str | dict,
                        output_filename:str,
                        output_units:str='W m-2 um-1',
                        reflectance:bool=False,
                        verbose:bool=None,
                        overwrite:bool=True,
                        write_telluric_spectra:bool=True,
                        write_model_spectra:bool=False,                        
                        qa_show:bool=None,
                        qa_scale:float=None,
                        qa_write:bool=None,
                        qa_block:bool=None):
                        
    """
    To correct spectra for telluric absorption and flux calibrate

    Parameters
    ----------
    object_file : str 
        The name of the pySpextool FITS file containing the object spectra.

    standard_file : str 
        The name of the pySpextool FITS file containing the standard star
        spectra.

    standard_info : str, list, dict
        If a string is passed, it is assumed to be the name of the standard.
        SIMBAD is queried for the spectral type and B- and V-band magnitudes.

        If a list is passed, it is assumed to contain the coordinates of the
        standard.  standard_info[0] = RA, standard_info[1] = Dec.  Each element
        is a string which gives the sexigesimal coordidates, separated by a
        blank space, :, or hms, dms.  SIMBAD is then queried for the id,
        spectral type, and B- and V-band magnitudes.

        If a dict is passed, it is assumed to contain the standard information.

        `"id"` : str
            The name of the standard.
    
        `"sptype"` : str
            The standard star spectral type.

        `"bmag"` : int or float
            The standard star B-band magnitude

        `"vmag"` : int or float
            The standard star V-band magnitude.
      
    output_filename : str
        The output file name sans the suffix, e.g. 'Wolf359'.

    qa_show : {None, True, False}
        Set to True to show a QA plot to the screen.
        Set to False to not show a QA plot to the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_block : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
  
    output_units : {'W m-2 um-1', 'erg s-1 cm-2 A-1', 'W m-2 Hz-1', 
                    'ergs s-1 cm-2 Hz-1', 'Jy', 'mJy', 'uJy'}
        The flux density units of the output.

        If reflectance=True, this has no effect.

    refelectence : {False, True} 
        Set to simply divide by the standard and return a relative 
        reflectance spectrum.


    write_model_spectra : {False, True}
        Set to True to write the modified Vega model to disk.

    write_telluric_spectra : {False, True}
        Set to True to write the telluric correction spectrum to disk.
    

    verbose : {None, True, False}
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file.  

    overwrite: {True, False}, optional
        Set to False to to not overwrite a file already on disk.

    Returns
    -------
    None
        Writes pyspextool FITS file to disk.

    """
    
    #
    # Check the parameters and keywords
    #

    check_parameter('telluric_correction', 'object_file', object_file, 'str')

    check_parameter('telluric_correction', 'standard_file', standard_file,
                    'str')

    check_parameter('telluric_correction', 'standard_info', standard_info,
                    ['str','list','dict'])
        
    check_parameter('telluric_correction', 'output_filename', output_filename,
                    'str')

    check_parameter('telluric_correction', 'reflectance', reflectance, 'bool')

    check_parameter('telluric_correction', 'write_telluric_spectra',
                    write_telluric_spectra, 'bool')

    check_parameter('telluric_correction', 'write_model_spectra',
                    write_model_spectra, 'bool')
    
    check_parameter('telluric_correction', 'output_units', output_units, 'str')
        
    check_parameter('telluric_correction', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'qa_scale', qa_scale,
                    ['NoneType', 'int','float'])

    check_parameter('telluric_correction', 'qa_block', qa_block,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('telluric_correction', 'overwrite', overwrite, 'bool')    

    keywords = check_keywords(verbose=verbose, qa_show=qa_show,
                              qa_scale=qa_scale, qa_block=qa_block,
                              qa_write=qa_write)

    #
    # Load the spectra
    #

    load_spectra(object_file,
                 standard_file,
                 standard_info,
                 output_filename,  
                 reflectance=reflectance,
                 verbose=keywords['verbose'])


    #
    # Are we doing a solar system object?
    #

    if reflectance is False:

        # Nope.  So do the full correction.

    
        #
        # Prepare line
        #

        if telluric.state['normalization_order'] != None:
    
            prepare_line(telluric.state['normalization_order'],
                         telluric.state['normalization_window'],
                         telluric.state['normalization_fittype'],
                         telluric.state['normalization_degree'],
                         verbose=keywords['verbose'],
                         qa_show=keywords['qa_show'],
                         qa_scale=keywords['qa_scale'],
                         qa_block=keywords['qa_block'],
                         qa_write=keywords['qa_write'])


        #
        # Measure radial velocity
        #

        if telluric.state['radialvelocity_nfwhm'] != None:    

            get_radialvelocity(telluric.state['radialvelocity_nfwhm'],
                               verbose=keywords['verbose'],
                               qa_show=keywords['qa_show'],
                               qa_scale=keywords['qa_scale'],
                               qa_block=keywords['qa_block'],
                               qa_write=keywords['qa_write'])

    
        #
        # Get kernels
        #

        if telluric.state['deconvolution_nfwhm'] != None:
    
            get_kernels(telluric.state['deconvolution_nfwhm'],
                        verbose=keywords['verbose'],
                        qa_show=keywords['qa_show'],
                        qa_scale=keywords['qa_scale'],
                        qa_block=keywords['qa_block'],
                        qa_write=keywords['qa_write'])
            
        else:

            get_kernels(verbose=keywords['verbose'],
                        qa_show=keywords['qa_show'],
                        qa_scale=keywords['qa_scale'],
                        qa_block=keywords['qa_block'],
                        qa_write=keywords['qa_write'])
        

    #
    # Construct the telluric correction spectra
    #

    make_telluric_spectra(intensity_unit=output_units)

    #
    # Correct the spectra for telluric absorption and flux calibrate
    #

    correct_spectra()

    #
    # Write the corrected spectra to disk
    #

    write_spectra(write_model_spectra=write_model_spectra,
                  write_telluric_spectra=write_telluric_spectra,
                  verbose=keywords['verbose'],
                  qa_show=keywords['qa_show'],
                  qa_scale=keywords['qa_scale'],
                  qa_block=keywords['qa_block'],
                  qa_write=keywords['qa_write'],
                  overwrite=overwrite)
    
    
