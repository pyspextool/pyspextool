from pyspextool import config as setup
from pyspextool.telluric import config 

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.telluric.load_standard import load_standard
from pyspextool.telluric.prepare_line import prepare_line
from pyspextool.telluric.get_radialvelocity import get_radialvelocity
from pyspextool.telluric.get_kernels import get_kernels
from pyspextool.telluric.adjust_ews import adjust_ews
from pyspextool.telluric.make_telluric_spectra import make_telluric_spectra
from pyspextool.telluric.write_telluric_spectra import write_telluric_spectra

def make_correction_spectra(
    standard_filename:str,
    standard_info:list | str | dict,
    correction_filename:str,
    correction_type:str='A0 V',
    output_units:str='W m-2 um-1',
    write_model_spectra:bool=False,                        
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    To create a "telluric correction" spectra and write a FITS file to disk.

    Parameters
    ----------
    standard_filename : str 
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
      
    correction_filename : str
        The output file name sans the suffix, e.g. 'Wolf359'.

    correction_type : {'A0 V', 'reflectance', 'basic'}
        The type of telluric correction requested.
        type = 'A0 V' -> standard correction using an A0 V standard star
        type = 'reflectance' -> ratio the object and standard.
        type = 'basic' -> ratio the object and standard and multiply by a
               Planck function of the appropriate temperature.
    
    output_units : {'W m-2 um-1', 'erg s-1 cm-2 A-1', 'W m-2 Hz-1', 
                    'ergs s-1 cm-2 Hz-1', 'Jy', 'mJy', 'uJy'}
        The flux density units of the output.

        If `correction_type`=True, this has no effect.
       
    write_model_spectra : {False, True}
        Set to True to write the modified Vega model to disk.

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
        Writes pyspextool FITS file(s) to disk.

    """

    #
    # Check the parameters and keywords
    #

    check_parameter('make_correction_spectra', 'standard_filename', 
                    standard_filename, 'str')

    check_parameter('make_correction_spectra', 'standard_info', standard_info,
                    ['str','list','dict'])
        
    check_parameter('make_correction_spectra', 'correction_filename', 
                    correction_filename, 'str')

    check_parameter('make_correction_spectra', 'correction_type', correction_type, \
                    'str', possible_values=setup.state['telluric_correctiontypes'])


    check_parameter('make_correction_spectra', 'write_model_spectra',
                    write_model_spectra, 'bool')
    
    check_parameter('make_correction_spectra', 'output_units', output_units, 'str')

    check_parameter('make_correction_spectra', 'verbose', verbose, ['NoneType', 'bool'])
    
    check_parameter('make_correction_spectra', 'qa_write', qa_write, \
                    ['NoneType', 'bool'])

    check_parameter('make_correction_spectra', 'qa_show', qa_show, \
                    ['NoneType', 'bool'])

    check_parameter('make_correction_spectra', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('make_correction_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Load and store the standard spectra
    #

    load_standard(
        standard_filename,
        standard_info,
        correction_filename,
        correction_type=correction_type,
        verbose=qa['verbose'])

    #
    # Are we doing a solar system object?
    #

    if config.state['correction_type'] == 'A0V':

        # Nope.  So do the full correction. <-- but missing RV
    
        #
        # Prepare line
        #

        if config.state['normalization_order'] is not None:
    
            prepare_line(
                config.state['normalization_order'],
                config.state['normalization_line'],
                config.state['resolving_power'],                         
                config.state['normalization_window'],
                config.state['normalization_fittype'],
                config.state['normalization_degree'],
                verbose=qa['verbose'],
                qa_show=qa['show'],
                qa_showscale=qa['showscale'],
                qa_showblock=qa['showblock'],
                qa_write=qa['write'])

        #
        # Measure radial velocity
        #

        if config.state['radialvelocity_nfwhm'] is not None:    

            get_radialvelocity(
                config.state['radialvelocity_nfwhm'],
                verbose=qa['verbose'],
                qa_show=qa['show'],
                qa_showscale=qa['showscale'],
                qa_showblock=qa['showblock'],
                qa_write=qa['write'])
            
        #
        # Get kernels
        #

        if config.state['deconvolution_nfwhm'] is not None:
    
            get_kernels(
                config.state['deconvolution_nfwhm'],
                verbose=qa['verbose'],
                qa_show=qa['show'],
                qa_showscale=qa['showscale'],
                qa_showblock=qa['showblock'],
                qa_write=qa['write'])
            
        else:

            get_kernels(
                verbose=qa['verbose'],
                qa_show=qa['show'],
                qa_showscale=qa['showscale'],
                qa_showblock=qa['showblock'],
                qa_write=qa['write'])

        #
        # Get EW scales
        #

        adjust_ews(
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])

    #
    # Make the telluric correction spectra
    # 

    make_telluric_spectra(
        intensity_unit=output_units,
        verbose=qa['verbose'])

    #
    # Write the results to disk
    #

    write_telluric_spectra(
        write_model_spectra=write_model_spectra,
        verbose=qa['verbose'],
        qa_show=qa['show'],
        qa_showscale=qa['showscale'],
        qa_showblock=qa['showblock'],
        qa_write=qa['write'])

    


