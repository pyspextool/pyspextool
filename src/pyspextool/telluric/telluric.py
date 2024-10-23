import logging

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.telluric.load_spectra import load_spectra
from pyspextool.telluric.prepare_line import prepare_line
from pyspextool.telluric.get_radialvelocity import get_radialvelocity
from pyspextool.telluric.get_kernels import get_kernels
from pyspextool.telluric.make_telluric_spectra import make_telluric_spectra
from pyspextool.telluric.correct_spectra import correct_spectra
from pyspextool.telluric.write_spectra import write_spectra

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def telluric(object_file:str,
             standard_file:str,
             standard_info:list | str | dict,
             output_filename:str,
             output_units:str='W m-2 um-1',
             reflectance:bool=False,
             write_telluric_spectra:bool=True,
             write_model_spectra:bool=False,                        
             verbose:bool=None,
             qa_show:bool=None,
             qa_showscale:float=None,
             qa_showblock:bool=None,
             qa_write:bool=None):

                        
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
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    

    Returns
    -------
    None
        Writes pyspextool FITS file to disk.

    """
    
    #
    # Check the parameters and keywords
    #

    check_parameter('telluric', 'object_file', object_file, 'str')

    check_parameter('telluric', 'standard_file', standard_file, 'str')

    check_parameter('telluric', 'standard_info', standard_info,
                    ['str','list','dict'])
        
    check_parameter('telluric', 'output_filename', output_filename, 'str')

    check_parameter('telluric', 'reflectance', reflectance, 'bool')

    check_parameter('telluric', 'write_telluric_spectra',
                    write_telluric_spectra, 'bool')

    check_parameter('telluric', 'write_model_spectra',
                    write_model_spectra, 'bool')
    
    check_parameter('telluric', 'output_units', output_units, 'str')

    check_parameter('telluric', 'verbose', verbose, ['NoneType', 'bool'])
    
    check_parameter('telluric', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('telluric', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('telluric', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('telluric', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])


    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    #
    # Load the spectra
    #

    load_spectra(object_file,
                 standard_file,
                 standard_info,
                 output_filename,  
                 reflectance=reflectance,
                 verbose=qa['verbose'])

    #
    # Are we doing a solar system object?
    #

    if reflectance is False:

        # Nope.  So do the full correction.
    
        #
        # Prepare line
        #

        if tc.state['normalization_order'] != None:
    
            prepare_line(tc.state['normalization_order'],
                         tc.state['normalization_window'],
                         tc.state['normalization_fittype'],
                         tc.state['normalization_degree'],
                         verbose=qa['verbose'],
                         qa_show=qa['show'],
                         qa_showscale=qa['showscale'],
                         qa_showblock=qa['showblock'],
                         qa_write=qa['write'])
    
        #
        # Measure radial velocity
        #

        if tc.state['radialvelocity_nfwhm'] != None:    

            get_radialvelocity(tc.state['radialvelocity_nfwhm'],
                               verbose=qa['verbose'],
                               qa_show=qa['show'],
                               qa_showscale=qa['showscale'],
                               qa_showblock=qa['showblock'],
                               qa_write=qa['write'])
    
        #
        # Get kernels
        #

        if tc.state['deconvolution_nfwhm'] != None:
    
            get_kernels(tc.state['deconvolution_nfwhm'],
                        verbose=qa['verbose'],
                        qa_show=qa['show'],
                        qa_showscale=qa['showscale'],
                        qa_showblock=qa['showblock'],
                        qa_write=qa['write'])
            
        else:

            get_kernels(verbose=qa['verbose'],
                        qa_show=qa['show'],
                        qa_showscale=qa['showscale'],
                        qa_showblock=qa['showblock'],
                        qa_write=qa['write'])

    #
    # Construct the telluric correction spectra
    #

    make_telluric_spectra(intensity_unit=output_units,
                          verbose=qa['verbose'])

    #
    # Correct the spectra for telluric absorption and flux calibrate
    #

    correct_spectra()

    #
    # Write the corrected spectra to disk
    #

    write_spectra(write_model_spectra=write_model_spectra,
                  write_telluric_spectra=write_telluric_spectra,
                  verbose=qa['verbose'],
                  qa_show=qa['show'],
                  qa_showscale=qa['showscale'],
                  qa_showblock=qa['showblock'],
                  qa_write=qa['write'])
    
    
