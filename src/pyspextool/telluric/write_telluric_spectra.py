from astropy.io import fits
import os
import logging

from pyspextool import config as setup
from pyspextool.telluric import config 
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.utils.units import get_latex_fluxdensity


def write_telluric_spectra(
    write_model_spectra:bool=False,
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float |int =None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    Write telluric-corrected spectra in pySpextool FITS files to disk.

    Will write the telluric-corrected spectra to disk, and optionally the
    telluric correction spectra and the (Vega) model spectra.  

    Parameters
    ----------
    write_model_spectra : {False, True}
        Set to True to write the modified (Vega) model to disk.
        Set to False to not write the modified (Vega) model to disk.    

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
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].
        
    Returns
    -------
    None
        Writes pyspextool FITS files and QA plots to disk.
    
    """

    #
    # Check parameters and keywords
    #

    check_parameter('write_telluric_spectra', 'write_model_spectra', \
                    write_model_spectra, 'bool')

    check_parameter('write_telluric_spectra', 'verbose', verbose, \
                    ['NoneType','bool'])
    
    check_parameter('write_telluric_spectra', 'qa_show', qa_show, \
                    ['NoneType','bool'])

    check_parameter('write_telluric_spectra', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])
    
    check_parameter('write_telluric_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
        
    check_parameter('write_telluric_spectra', 'qa_write', qa_write, \
                    ['NoneType','bool'])

    check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)

    #
    # Grab the standard's hdrinfo dictionary
    #

    std_hdrinfo = get_headerinfo(config.state['standard_astropyheader'])

    #
    # Write the (Vega) model to disk if requested.
    #

    if write_model_spectra is True:

        # Create the header. Make a new hdrinfo dictionary using the standard
        # hdrinfo dictionary
        
        hdrinfo = {}
        
        hdrinfo['NORDERS'] = std_hdrinfo['NORDERS']

        hdrinfo['ORDERS'] = std_hdrinfo['ORDERS']

        hdrinfo['NAPS'] = std_hdrinfo['NAPS']                

        hdrinfo['MODULE']=std_hdrinfo['MODULE']

        hdrinfo['MODULE'][0] = 'telluric'

        hdrinfo['VERSION'] = std_hdrinfo['VERSION']

        filename = [config.state['telluric_output_filename']+'_model.fits', 
                    ' File name']
        hdrinfo['FILENAME'] = filename

        mag = float('{:.3f}'.format(config.state['standard_bmag']))
        hdrinfo['STD_BMAG'] = (mag, 'Telluric standard B mag')

        mag = float('{:.3f}'.format(config.state['standard_vmag']))
        hdrinfo['STD_VMAG'] = (mag, 'Telluric standard V mag')

        # Now deal with the units

        latex = get_latex_fluxdensity(config.state['intensity_unit']) 
        
        hdrinfo['XUNITS'] = std_hdrinfo['XUNITS']
        hdrinfo['YUNITS'] = std_hdrinfo['YUNITS']
        hdrinfo['YUNITS'][0] = config.state['intensity_unit']
        hdrinfo['LXUNITS'] = std_hdrinfo['LXUNITS']
        hdrinfo['LYUNITS'] = std_hdrinfo['LYUNITS']        
        hdrinfo['LYUNITS'][0] = latex[0]
        hdrinfo['LXLABEL'] = std_hdrinfo['LXLABEL']
        hdrinfo['LYLABEL'] = std_hdrinfo['LYLABEL']        
        hdrinfo['LYLABEL'][0] = latex[1]

        # Create the basic header

        phdu = fits.PrimaryHDU()
        newhdr = phdu.header

        # Add our keywords
    
        keys = list(hdrinfo.keys())
    
        for i in range(len(keys)):
                            
                newhdr[keys[i]] = (hdrinfo[keys[i]][0],
                                   hdrinfo[keys[i]][1])
        
        #
        # Write the file out
        #

        full_path = os.path.join(setup.state['proc_path'],
                                 config.state['telluric_output_filename']+\
                                 '_model.fits')    

        fits.writeto(full_path, config.state['model_spectra'], newhdr,
                     overwrite=True)
        
        logging.info(' Wrote file '+os.path.basename(full_path) + \
                     ' to the proc directory.')

    #
    # Now create a header for the output filefile.
    #
        
    hdrinfo = {}

    hdrinfo['RA'] = std_hdrinfo['RA']

    hdrinfo['DEC'] = std_hdrinfo['DEC']

    hdrinfo['AVE_AM'] = std_hdrinfo['AVE_AM']

    hdrinfo['MODE'] = std_hdrinfo['MODE']

    hdrinfo['SLTW_ARC'] = std_hdrinfo['SLTW_ARC']
        
    hdrinfo['NORDERS'] = std_hdrinfo['NORDERS']

    hdrinfo['ORDERS'] = std_hdrinfo['ORDERS']

    hdrinfo['NAPS'] = std_hdrinfo['NAPS']                
    
    hdrinfo['MODULE'] = std_hdrinfo['MODULE']
    hdrinfo['MODULE'][0] = 'telluric'
    
    hdrinfo['VERSION'] = std_hdrinfo['VERSION']

    filename = [config.state['telluric_output_filename']+'.fits', ' File name']
    hdrinfo['FILENAME'] = filename

    mag = float('{:.3f}'.format(config.state['standard_bmag']))
    hdrinfo['STD_BMAG'] = (mag, 'Telluric standard B mag')
    
    mag = float('{:.3f}'.format(config.state['standard_vmag']))
    hdrinfo['STD_VMAG'] = (mag, 'Telluric standard V mag')
        
    value = [config.state['standard_fullfilename'], 
             'Telluric standard input file']
    hdrinfo['TC_SFILE'] = value
    
    value = [config.state['standard_name'], 'Telluric standard']
    hdrinfo['TC_STDID'] = value
    
    value = [config.state['standard_sptype'], 'Telluric standard spectral type']
    hdrinfo['TC_STDST'] = value
    
    value = [int(config.state['standard_teff']), 
             'Telluric standard effective temperature (K)']
    hdrinfo['TC_STDT'] = value
    
    value = [float('{:.3f}'.format(config.state['standard_bmag'])),
             'Telluric standard B mag']
    hdrinfo['TC_STDB'] = value
    
    value = [float('{:.3f}'.format(config.state['standard_vmag'])),
             'Telluric standard V mag']
    hdrinfo['TC_STDV'] = value
    
    if 'standard_rv' in list(config.state.keys()):
        value = [float('{:.2f}'.format(config.state['standard_rv'])),
                 'Telluric standard radial velocity (km s-1)']
        hdrinfo['TC_STDRV'] =  value
    else: hdrinfo['TC_STDRV'] =  [0.,'Telluric standard radial velocity (km s-1)']
            
    value = [config.state['correction_type'],'Telluric correction type']
    hdrinfo['TC_TYPE'] = value
    
    value = [config.state['kernel_method'], 'Kernel creation method']
    hdrinfo['TC_METH'] = value
    
    if config.state['kernel_method'] == 'deconvolution' and config.state['correction_type'] == 'A0 V':
        
        value = [float('{:.5f}'.format(config.state['max_deviation'])),
                 'Telluric maxmimum % deviation of Vega-data']
        hdrinfo['TC_MXDEV'] = value
        
        value = [float('{:.5f}'.format(config.state['rms_deviation'])),
                 'Telluric RMS deviation of Vega-data']
        hdrinfo['TC_RMS'] = value
            
    # Add units to the hdrinfo

    result = get_latex_fluxdensity(config.state['intensity_unit'])[0]
    
    hdrinfo['XUNITS'] = std_hdrinfo['XUNITS']
    hdrinfo['YUNITS'] = std_hdrinfo['YUNITS']
    hdrinfo['YUNITS'][0] = config.state['intensity_unit']+' / DN s-1'
    
    hdrinfo['LXUNITS'] = std_hdrinfo['LXUNITS']
    hdrinfo['LYUNITS'] = std_hdrinfo['LYUNITS']        
    hdrinfo['LYUNITS'][0] = result+' / DN s$^{-1}$'
    
    hdrinfo['LXLABEL'] = std_hdrinfo['LXLABEL']
    hdrinfo['LYLABEL'] = std_hdrinfo['LYLABEL']        
    hdrinfo['LYLABEL'][0] = 'Ratio ('+result+' / DN s$^{-1}$)'

    # Create the FITS header

    # Create the basic headers
        
    phdu = fits.PrimaryHDU()
    hdr = phdu.header


    
    # Add our keywords
    
    keys = list(hdrinfo.keys())
    
    for i in range(len(keys)):
        
        if keys[i] == 'COMMENT':
            
            pass
        
        else:
            
            hdr[keys[i]] = (hdrinfo[keys[i]][0], 
                            hdrinfo[keys[i]][1])

    #
    # Write the file out
    #
    
    telluric_fullpath = os.path.join(
        setup.state['proc_path'],
        config.state['telluric_output_filename']+\
        '.fits')    

    
    fits.writeto(
        telluric_fullpath,
        config.state['rawtc_spectra'],
        hdr,
        overwrite=True)
    
    logging.info(' Wrote file '+os.path.basename(telluric_fullpath) + \
                 " to the proc directory.")
    

    



    
