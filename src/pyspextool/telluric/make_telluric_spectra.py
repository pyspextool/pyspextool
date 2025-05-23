import numpy as np
import copy
import logging

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.io.check import check_parameter
from pyspextool.utils.units import convert_fluxdensity
from pyspextool.utils.irplanck import irplanck
from pyspextool.telluric.core import make_telluric_spectrum
from pyspextool.utils.interpolate import linear_interp1d
from pyspextool.utils.interpolate import linear_bitmask_interp1d


def make_telluric_spectra(intensity_unit:str='W m-2 um-1',
                          verbose:bool=True,
                          new=False):

    """
    To create telluric correction spectra.

    Parameters
    ----------
    intensity_unit : str, default='W m-2 um-1'
        The requested intensity units of the final output spectrum.  

    Returns
    -------
    None
    Loads data into memory.

        telluric.load['intensity_unit']
        telluric.state['telluric_spectra']
        telluric.state['model_spectra']
        telluric.state['make_done']    

    """

    #
    # Check the variable
    #

    if tc.state['method'] == 'deconvolution':
             
        if tc.state['kernel_done'] is False:

            message = "The kernels have not been created.  Please run '\
            'make_kernels.py."
            raise pySpextoolError(message)

    else:
            
        if tc.state['load_done'] is False:

            message = "The spectra have not been loaded.  Please run '\
            'load_spectra.py."
            raise pySpextoolError(message)
           
    #
    # Check parameters
    #

    check_parameter('make_telluric_spectra', 'intensity_unit',
                    intensity_unit, 'str',
                    possible_values=setup.state['units'])    
        
    logging.info(' Making telluric correction spectra.')

    tc.state['intensity_unit'] = intensity_unit
    
    telluric_spectra = copy.deepcopy(tc.state['standard_spectra'])

    #
    # Go by correction type
    #

    model_spectra = copy.deepcopy(tc.state['standard_spectra'])
    model_spectra[:,1,:] = 1.0
    model_spectra[:,2,:] = 1.0
    model_spectra[:,3,:] = 0
        
    if tc.state['type'] == 'A0V':

        tc.state['intensity_unit'] = intensity_unit
                
        for i in range(tc.state['standard_norders']):    
            
            standard_wavelength = tc.state['standard_spectra'][i,0,:]
            standard_fluxdensity = tc.state['standard_spectra'][i,1,:]
            standard_uncertainty = tc.state['standard_spectra'][i,2,:]
            
            standard_bmag = tc.state['standard_bmag']
            standard_vmag = tc.state['standard_vmag']
            standard_rv = tc.state['standard_rv']                
            
            vega_wavelength = tc.state['vega_wavelength']
            vega_fluxdensity = tc.state['vega_fluxdensity']
            vega_continuum = tc.state['vega_continuum']
            vega_fitted_continuum = tc.state['vega_fitted_continuum']
            kernel = tc.state['kernels'][i]

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
                                            tc.state['control_points'][i][0,:],
                                            tc.state['control_points'][i][1,:],
                                            int(tc.state['standard_orders'][i]),
                                            intensity_unit,
                                            new=new)
            
            # Store results
            
            telluric_spectra[i,1,:] = result[0]
            telluric_spectra[i,2,:] = result[1]
            model_spectra[i,1,:] = result[2]
            
    if tc.state['type'] == 'basic':

        tc.state['intensity_unit'] = intensity_unit
        
        for i in range(tc.state['object_norders']):    

            standard_wavelength = tc.state['standard_spectra'][i,0,:]
            standard_fluxdensity = tc.state['standard_spectra'][i,1,:]
            standard_uncertainty = tc.state['standard_spectra'][i,2,:]

        
            telluric_spectrum = 1/standard_fluxdensity
            telluric_unc = 1/standard_fluxdensity**2 * standard_uncertainty

            # Flux calibrate using a Planck function

            planck_atzlambda = irplanck(setup.state['vega_zlambda'],
                                        tc.state['standard_teff'])

            dmag = tc.state['standard_vmag']+setup.state['vega_zmag']
            scale = setup.state['vega_zfd']*10**(-0.4*(dmag)) / planck_atzlambda

            planck = irplanck(standard_wavelength,
                              tc.state['standard_teff'])*scale
            
            telluric_spectrum *= planck 
            telluric_unc *= planck

            #
            # Change units to those requested by the user
            #
            
            flux = convert_fluxdensity(standard_wavelength,
                                       telluric_spectrum,
                                       'um','erg s-1 cm-2 A-1',
                                       intensity_unit)
            
            unc = convert_fluxdensity(standard_wavelength,
                                      telluric_unc,
                                      'um','erg s-1 cm-2 A-1',
                                      intensity_unit)
            
            planck = convert_fluxdensity(standard_wavelength,
                                         planck,
                                         'um','erg s-1 cm-2 A-1',
                                         intensity_unit)
            
            # Updates labels
            
            telluric_spectra[i,1,:] = flux
            telluric_spectra[i,2,:] = unc
            model_spectra[i,1,:] = planck
                               
    if tc.state['type'] == 'reflectance':

        tc.state['intensity_unit'] = 'reflectance'

        vega_spectra = copy.deepcopy(tc.state['standard_spectra'])
        vega_spectra[:,1,:] = 1.0
        vega_spectra[:,2,:] = 1.0
        vega_spectra[:,3,:] = 0
        
        for i in range(tc.state['object_norders']):    
            
            standard_fluxdensity = tc.state['standard_spectra'][i,1,:]
            standard_uncertainty = tc.state['standard_spectra'][i,2,:]

        
            telluric_spectra[i,1,:] = 1/standard_fluxdensity
            telluric_spectra[i,2,:] = 1/standard_fluxdensity**2 * \
                standard_uncertainty

    #
    # Interpolate the telluric spectra onto the object wavelengths
    #

    object_spectra = tc.state['object_spectra']
    tmp = copy.deepcopy(tc.state['object_spectra'])    

    for i in range(tc.state['object_norders']):

        for j in range(tc.state['object_napertures']):

            std_idx = tc.state['standard_orders'] == \
                tc.state['object_orders'][i]
            obj_idx = i + j*tc.state['object_napertures']
            
            ri, ru = linear_interp1d(telluric_spectra[std_idx,0,:],
                                     telluric_spectra[std_idx,1,:],
                                     object_spectra[obj_idx,0,:],
                                     input_u=telluric_spectra[std_idx,2,:])

            tmp[obj_idx,1,:] = ri
            tmp[obj_idx,2,:] = ru

            rm = linear_bitmask_interp1d(telluric_spectra[std_idx,0,:],
                    telluric_spectra[std_idx,3,:].astype(np.uint8),
                                         object_spectra[obj_idx,0,:])
            tmp[obj_idx,3,:] = rm
            
    #
    # Store the results.  Pre fill ewcorrected and shifted in case the
    # user chooses to not correct the H lines or shift the standard.
    #

    tc.state['model_spectra'] = model_spectra
    tc.state['rawtc_spectra'] = np.array(tmp)
    tc.state['ewcorrectedtc_spectra'] = np.array(tmp)
    tc.state['shiftedtc_spectra'] = np.array(tmp)        
    tc.state['shifts'] = np.zeros((tc.state['object_norders'],
                                   tc.state['object_napertures']))
    
    tc.state['make_done'] = True


    
