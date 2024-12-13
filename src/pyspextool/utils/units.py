import numpy as np
import numpy.typing as npt
from astropy import units as u
from pyspextool.io.check import check_parameter

def convert_fluxdensity(wavelength:npt.ArrayLike,
                        flux_density:npt.ArrayLike,
                        input_wunit:str,
                        input_funit:str,
                        output_funit:str):

    """
    To convert between flux density units

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths

    flux_density : ndarray
        A (nwave,) array of flux densities

    input_wunit : {'um', 'A', 'nm'}
        The units of `wavelength`
        'um'=microns
        'A'=Angstroms
        'nm'=nanometers

    input_funit : {'W m-2 um-1', 'erg s-1 cm-1 A-1', 'Jy', 'mJy', 'uJy',
                   'W m-2 Hz-1', 'erg s-1 cm-2 Hz-1'}
        The units of 'fluxdensity'

    output_funit : {'W m-2 um-1', 'erg s-1 cm-2 A-1', 'Jy', 'mJy', 'uJy',
                    'W m-2 Hz-1', 'erg s-1 cm-2 Hz-1'}
        The desired units of 'fluxdensity'

    Returns
    -------
    ndarray

    The 'fluxdensity' array converted to units of 'output_funit'.

    """

    #
    # Check parameters
    #

    types = ['W m-2 um-1', 'erg s-1 cm-2 A-1', 'Jy', 'mJy', 'uJy',
             'W m-2 Hz-1', 'erg s-1 cm-2 Hz-1']
    
    check_parameter('convert_fluxdensity', 'wavelength', wavelength, 'ndarray')

    check_parameter('convert_fluxdensity', 'flux_density', flux_density,
                    'ndarray')

    check_parameter('convert_fluxdensity', 'input_wunit', input_wunit, 'str',
                    possible_values=['A', 'nm', 'um'])

    check_parameter('convert_fluxdensity', 'input_funit', input_funit, 'str',
                    possible_values=types)

    check_parameter('convert_fluxdensity', 'output_funit', output_funit, 'str',
                    possible_values=types)    
    

    #
    # Convert input wavelengths to microns
    #

    if input_wunit == 'A':

        twavelength =wavelegth /1e4

    elif input_wunit == 'nm':

        twavelength = wavelength / 1e3

    else:

        twavelength = wavelength


    #
    # Now convert the flux density to erg s-1 cm-2 A-1
    #

    if input_funit == 'W m-2 um-1':

        from_unit = u.W / u.m**2 / u.um
        to_unit = u.erg / u.cm**2 / u.s / u.AA

        tflux_density = flux_density*(from_unit).to(to_unit)        

    elif input_funit == 'W m-2 Hz-1':

        from_unit = u.W / u.m**2 / u.Hz
        to_unit = u.erg / u.cm**2 / u.s / u.AA
        
        tflux_density = flux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))
        
    elif input_funit == 'erg s-1 cm-2 A-':

        tflux_density = flux_density

    elif input_funit == 'erg s-1 cm-2 Hz-1':

        from_unit = u.erg / u.s / u.cm**2 / u.Hz
        to_unit = u.erg / u.cm**2 / u.s / u.AA
        
        tflux_density = flux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))
        
    elif input_funit == 'Jy':

        from_unit = u.Jy
        to_unit = u.erg / u.cm**2 / u.s / u.AA

        tflux_density = flux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))
        

    elif input_funit == 'mJy':

        from_unit = u.Jy
        to_unit = u.erg / u.cm**2 / u.s / u.AA

        tflux_density = 1e3*flux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))

        
    elif input_funit == 'uJy':

        from_unit = u.Jy
        to_unit = u.erg / u.cm**2 / u.s / u.AA

        tflux_density = 1e6*flux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))        
    else:

        tflux_density = flux_density

        
        
    #
    # Now convert to the desired units
    #

    if output_funit == 'W m-2 um-1':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.W / u.m**2 / u.um

        tflux_density = tflux_density*(from_unit).to(to_unit)

    elif output_funit == 'W m-2 Hz-1':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.W / u.m**2 / u.Hz
        
        tflux_density = tflux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))

    elif output_funit == 'erg s-1 cm-2 A-1':

        tflux_density = tflux_density

    elif output_funit == 'erg s-1 cm-2 Hz-1':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.erg / u.cm**2 / u.s / u.Hz
                
        tflux_density = tflux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))
        
    elif output_funit == 'Jy':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.Jy
                
        tflux_density = tflux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))
        
    elif output_funit == 'mJy':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.Jy
                
        tflux_density = tflux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))

        tflux_density *= 1e3
        
    elif output_funit == 'uJy':

        from_unit = u.erg / u.cm**2 / u.s / u.AA
        to_unit = u.Jy
                
        tflux_density = tflux_density*(from_unit).to(to_unit,
                        equivalencies=u.spectral_density(twavelength*u.micron))

        tflux_density *= 1e6
        
    return tflux_density
        

    
def convert_wavelength(wavelength, input_unit, output_unit):

    """
    To convert between wavelength units

    Parameters
    ----------
    wavelength : ndarray
        A (nwave,) array of wavelengths

    input_unit : {'um', 'A', 'nm'}
        The units of `wavelength`
        'um'=microns
        'A'=Angstroms
        'nm'=nanometers

    output_unit : {'um', 'A', 'nm'}
        The desired units of `wavelength`.
        'um'=microns
        'A'=Angstroms
        'nm'=nanometers

    Returns
    -------
    ndarray

    The `wavelength` array converted to units of `output_unit`.

    """

    #
    # Check parameters
    #

    check_parameter('convert_wavelength', 'wavelength', wavelength, 'ndarray')

    check_parameter('convert_wavelength', 'input_unit', input_unit, 'str',
                    possible_values=['A', 'nm', 'um'])

    check_parameter('convert_wavelength', 'output_unit', output_unit, 'str',
                    possible_values=['A', 'nm', 'um'])

    #
    # Convert to Angstroms
    #

    if input_unit == 'nm':

        twavelength =wavelegth * 10

    elif input_unit == 'um':

        twavelength = wavelength * 1e4

    else:

        twavelength = wavelength

    #
    # Now convert to the requested unit
    #

    if output_unit == 'nm':

        twavelength /= 10

    elif output_unit == 'um':

        twavelength /= 1e4

                
    return twavelength

def get_latex_fluxdensity(unit:str):

    """
    To return latex appropriate things

    Parameters
    ----------
    unit : {'W m-2 um-1', 'erg s-1 cm-1 A-1', 'Jy', 'mJy', 'uJy',
            'W m-2 Hz-1', 'erg s-1 cm-2 Hz-1', 'reflectence'}
        A str ASCII unit

    
    Returns
    -------
    str, str, str
    latex_unit : str
        Latex version of the ASCII units

    latex_flux : str
        Either $f_\lambda$ or $f\_nu$ and (`latex_unit`)

    latex_unc : str
        Either $\sigma$ and (`latex_unit`)
    
        
    """

    #
    # Check parameter
    #

    values = ['W m-2 um-1',
              'erg s-1 cm-2 A-1',
              'Jy',
              'mJy',
              'uJy',
              'W m-2 Hz-1',
              'erg s-1 cm-2 Hz-1',
              'reflectance']
    
    check_parameter('get_latex_fluxdensity', 'unit', unit, 'str',
                    possible_values=values)

    units = np.array(['W m-2 um-1',
                      'erg s-1 cm-2 A-1',
                      'Jy',
                      'mJy',
                      'uJy',
                      'W m-2 Hz-1',
                      'erg s-1 cm-2 Hz-1',
                      'reflectance'])

    lunits = np.array([r'(W m$^{-2}$ $\mu$m$^{-1}$)',
                       r'(erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)',
                       '(Jy)',
                       '(mJy)',
                       r'($\mu$Jy)',
                       r'(W m$^{-2}$ Hz$^{-1}$)',
                       r'(erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)',
                       ''])

    ftype = np.array([r'$f_\lambda$',
                      r'$f_\lambda$',
                      r'$f_\\nu$',
                      r'$f_\\nu$',
                      r'$f_\\nu$',
                      r'$f_\\nu$',                      
                      r'$f_\\nu$',
                      'Reflectance'])

    #
    # Find the match, and return to user
    #
    
    z = units == unit
    
    return lunits[z][0], ftype[z][0]+lunits[z][0], r'$\sigma$'+lunits[z][0] 
    



    
