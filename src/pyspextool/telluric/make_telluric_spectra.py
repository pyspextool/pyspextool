import copy
import logging

from pyspextool import config as setup
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric import config as telluric
from pyspextool.io.check import check_parameter
from pyspextool.utils.units import convert_fluxdensity
from pyspextool.telluric.core import make_telluric_spectrum

def make_telluric_spectra(intensity_unit:str='W m-2 um-1'):

    """
    To create telluric correction spectra.

    Parameters
    ----------
    fluxdensity_unit : str, default='W m-2 um-1'
        The requested flux density units of the final output spectrum.  

    Returns
    -------
    None
    Loads data into memory.

        telluric.load['intensity_unit']
        telluric.state['telluric_spectra']
        telluric.state['vega_spectra']
        telluric.state['make_done']    

    """

    #
    # Check the variable
    #

    if telluric.state['method'] == 'deconvolution':
             
        if telluric.state['kernel_done'] is False:

            message = "The kernels have not been created.  Please run '\
            'make_kernels.py."
            raise pySpextoolError(message)

    else:
            
        if telluric.state['load_done'] is False:

            message = "The spectra have not been loaded.  Please run '\
            'load_spectra.py."
            raise pySpextoolError(message)
           
    #
    # Check parameters
    #

    check_parameter('make_telluric_spectra', 'intensity_unit',
                    intensity_unit, 'str',
                    possible_values=setup.state['units'])    
        
    logging.info(f' Making telluric correction spectra...')

    telluric.load['intensity_unit'] = intensity_unit
    
    telluric_spectra = copy.deepcopy(telluric.state['standard_spectra'])
    
    if telluric.load['reflectance'] is False:
        
        telluric.load['intensity_unit'] = intensity_unit
        
        vega_spectra = copy.deepcopy(telluric.state['standard_spectra'])
        vega_spectra[:,2,:] = 1.0
        vega_spectra[:,3,:] = 0
        
        for i in range(telluric.state['standard_norders']):    
            
            standard_wavelength = telluric.state['standard_spectra'][i,0,:]
            standard_fluxdensity = telluric.state['standard_spectra'][i,1,:]
            standard_uncertainty = telluric.state['standard_spectra'][i,2,:]
            
            standard_bmag = telluric.state['standard_bmag']
            standard_vmag = telluric.state['standard_vmag']
            standard_rv = telluric.state['standard_rv']                
            
            vega_wavelength = telluric.state['vega_wavelength']
            vega_fluxdensity = telluric.state['vega_fluxdensity']
            vega_continuum = telluric.state['vega_continuum']
            vega_fitted_continuum = telluric.state['vega_fitted_continuum']
            kernel = telluric.state['kernels'][i]
            scale = float(telluric.state['ew_scale'])
            
            result = make_telluric_spectrum(standard_wavelength,
                                            standard_fluxdensity,
                                            standard_uncertainty,
                                            standard_rv,
                                            standard_bmag,
                                            standard_vmag,
                                            vega_wavelength,
                                            vega_fluxdensity,
                                            vega_continuum,
                                            vega_fitted_continuum,
                                            kernel,
                                            scale)
            
            #
            # Change units to those requested by the user
            #

            flux = result[0]
            unc = result[1]
            vega = result[2]
            
            flux = convert_fluxdensity(standard_wavelength,
                                       result[0],'um','erg s-1 cm-2 A-1',
                                       intensity_unit)
            
            unc = convert_fluxdensity(standard_wavelength,
                                      result[1],'um','erg s-1 cm-2 A-1',
                                      intensity_unit)
            
            vega = convert_fluxdensity(standard_wavelength,
                                       result[2],'um','erg s-1 cm-2 A-1',
                                       intensity_unit)
            
            # Updates labels
            
            telluric_spectra[i,1,:] = flux
            telluric_spectra[i,2,:] = unc
            vega_spectra[i,1,:] = vega
            
        #
        # Store the results
        #

        telluric.state['telluric_spectra'] = telluric_spectra
        telluric.state['vega_spectra'] = vega_spectra    

    else:

        telluric.load['intensity_unit'] = 'reflectance'
        
        for i in range(telluric.state['object_norders']):    
            
            standard_fluxdensity = telluric.state['standard_spectra'][i,1,:]
            standard_uncertainty = telluric.state['standard_spectra'][i,2,:]

        
            telluric_spectra[i,1,:] = 1/standard_fluxdensity
            telluric_spectra[i,2,:] = 1/standard_fluxdensity**2 * \
                standard_uncertainty
            
    #
    # Store the results
    #

    telluric.state['telluric_spectra'] = telluric_spectra

    telluric.state['make_done'] = True
    
