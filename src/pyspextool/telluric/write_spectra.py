from astropy.io import fits
import os
import logging
import numpy as np

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils.add_entry import add_entry
from pyspextool.utils.units import get_latex_fluxdensity
from pyspextool.plot.plot_spectra import plot_spectra

def write_spectra(write_model_spectra:bool=False,
                  write_telluric_spectra:bool=True,
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

    write_telluric_spectra : {True, False}
        Set to True to write the telluric correction spectra to disk.
        Set to False to not write the telluric correction spectra to disk.

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

    check_parameter('write_spectra', 'write_model_spectra', write_model_spectra,
                    'bool')

    check_parameter('write_spectra', 'write_telluric_spectra',
                    write_telluric_spectra, 'bool')    
    
    check_parameter('write_spectra', 'verbose', verbose, ['NoneType','bool'])
    
    check_parameter('write_spectra', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('write_spectra', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])
    
    check_parameter('write_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
        
    check_parameter('write_spectra', 'qa_write', qa_write, ['NoneType','bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    #
    # Write the (Vega) model to disk if requested.
    #

    if write_model_spectra is True:

        # Create the header. Make a new hdrinfo dictionary using the standard
        # hdrinfo dictionary
        
        std_hdrinfo = tc.state['standard_hdrinfo']

        hdrinfo = {}
        
        hdrinfo['NORDERS'] = std_hdrinfo['NORDERS']

        hdrinfo['ORDERS'] = std_hdrinfo['ORDERS']

        hdrinfo['NAPS'] = std_hdrinfo['NAPS']                

        hdrinfo['MODULE']=std_hdrinfo['MODULE']

        hdrinfo['MODULE'][0] = 'telluric'

        hdrinfo['VERSION'] = std_hdrinfo['VERSION']

        filename = [tc.state['output_filename']+'_model.fits', ' File name']
        hdrinfo['FILENAME'] = filename

        mag = float('{:.3f}'.format(tc.state['standard_bmag']))
        hdrinfo['STD_BMAG'] = (mag, 'Telluric standard B mag')

        mag = float('{:.3f}'.format(tc.state['standard_vmag']))
        hdrinfo['STD_VMAG'] = (mag, 'Telluric standard V mag')

        # Now deal with the units

        latex = get_latex_fluxdensity(tc.state['intensity_unit']) 
        
        hdrinfo['XUNITS'] = std_hdrinfo['XUNITS']
        hdrinfo['YUNITS'] = std_hdrinfo['YUNITS']
        hdrinfo['YUNITS'][0] = tc.state['intensity_unit']
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
                                 tc.state['output_filename']+\
                                 '_vega.fits')    

        fits.writeto(full_path, tc.state['model_spectra'], newhdr,
                     overwrite=True)
        
        logging.info(' Wrote file '+os.path.basename(full_path) + ' to proc/.')

    #
    # Now create a headerinfo that is the basis for the _telluric.fits file and 
    # the telluric corrected spectrum file
    #

    obj_hdrinfo = tc.state['object_hdrinfo']

    hdrinfo = {}
        
    hdrinfo['NORDERS'] = obj_hdrinfo['NORDERS']

    hdrinfo['ORDERS'] = obj_hdrinfo['ORDERS']

    hdrinfo['NAPS'] = obj_hdrinfo['NAPS']                
    
    hdrinfo['MODULE'] = obj_hdrinfo['MODULE']
    hdrinfo['MODULE'][0] = 'telluric'
    
    hdrinfo['VERSION'] = obj_hdrinfo['VERSION']

    filename = [tc.state['output_filename']+'_telluric.fits', ' File name']
    hdrinfo['FILENAME'] = filename

    mag = float('{:.3f}'.format(tc.state['standard_bmag']))
    hdrinfo['STD_BMAG'] = (mag, 'Telluric standard B mag')
    
    mag = float('{:.3f}'.format(tc.state['standard_vmag']))
    hdrinfo['STD_VMAG'] = (mag, 'Telluric standard V mag')
        
    value = [tc.state['standard_file'], 'Telluric standard input file']
    hdrinfo['TC_SFILE'] = value
    
    value = [tc.state['standard_name'], 'Telluric standard']
    hdrinfo['TC_STDID'] = value
    
    value = [tc.state['standard_sptype'], 'Telluric standard spectral type']
    hdrinfo['TC_STDST'] = value
    
    value = [int(tc.state['standard_teff']), 
             'Telluric standard effective temperature (K)']
    hdrinfo['TC_STDT'] = value
    
    value = [float('{:.3f}'.format(tc.state['standard_bmag'])),
             'Telluric standard B mag']
    hdrinfo['TC_STDB'] = value
    
    value = [float('{:.3f}'.format(tc.state['standard_vmag'])),
             'Telluric standard V mag']
    hdrinfo['TC_STDV'] = value
    
    value = [float('{:.2f}'.format(tc.state['standard_rv'])),
             'Telluric standard radial velocity (km s-1)']
    hdrinfo['TC_STDRV'] =  value
            
    value = [tc.state['type'],'Telluric correction type']
    hdrinfo['TC_TYPE'] = value
    
    value = [tc.state['method'], 'Kernel creation method']
    hdrinfo['TC_METH'] = value
    
    if tc.state['method'] == 'deconvolution':
        
        value = [float('{:.5f}'.format(tc.state['max_deviation'])),
                 'Telluric maxmimum % deviation of Vega-data']
        hdrinfo['TC_MXDEV'] = value
        
        value = [float('{:.5f}'.format(tc.state['rms_deviation'])),
                 'Telluric RMS deviation of Vega-data']
        hdrinfo['TC_RMS'] = value
            
    # Deal with the shifts
            
    orders = tc.state['object_orders']
    for i in range(tc.state['object_norders']):
                
        name = 'TC_SO' + str(orders[i]).zfill(3)
        value = ', '.join(tc.state['shifts'][i,:].astype(str))
        comment = 'Standard shift (pixels) for order ' + \
            str(orders[i]).zfill(3)
        hdrinfo[name] = [value, comment]
        
    #
    # Write the telluric correction spectra to disk.
    #
    
    if write_telluric_spectra is True:
        
        # Add units to the hdrinfo

        result = get_latex_fluxdensity(tc.state['intensity_unit'])[0]
        
        hdrinfo['XUNITS'] = obj_hdrinfo['XUNITS']
        hdrinfo['YUNITS'] = obj_hdrinfo['YUNITS']
        hdrinfo['YUNITS'][0] = tc.state['intensity_unit']+' / DN s-1'

        hdrinfo['LXUNITS'] = obj_hdrinfo['LXUNITS']
        hdrinfo['LYUNITS'] = obj_hdrinfo['LYUNITS']        
        hdrinfo['LYUNITS'][0] = result

        hdrinfo['LXLABEL'] = obj_hdrinfo['LXLABEL']
        hdrinfo['LYLABEL'] = obj_hdrinfo['LYLABEL']        
        hdrinfo['LYLABEL'][0] = 'Ratio ('+result+' / DN s$^{-1}$)'

        #
        # Create the FITS header
        #
        
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
        
        telluric_fullpath = os.path.join(setup.state['proc_path'],
                                         tc.state['output_filename']+\
                                         '_telluric.fits')    
        
        fits.writeto(telluric_fullpath,
                     tc.state['shiftedtc_spectra'],
                     hdr,
                     overwrite=True)
        
        logging.info(' Wrote file '+os.path.basename(telluric_fullpath) + \
                     " to the proc directory.")
        
    
    #
    # Write the corrected spectra to disk
    #

    #  Add two additional keywords to the hdrinfo dictionary


    value = [float('{:.2f}'.format(tc.state['delta_airmass'])),
             'Telluric Average airmass difference']
    hdrinfo = add_entry(hdrinfo, 'TC_STDV', 'after', 'TC_dAM', value)

    value = [float('{:.2f}'.format(tc.state['delta_angle'])),\
             'Telluric angular separation (deg) of obj and std']
    hdrinfo = add_entry(hdrinfo, 'TC_dAM', 'after', 'TC_dAN', value)

    
    # Rename header for ease of reading
            
    obj_hdrinfo = tc.state['object_hdrinfo']
    
    # Store the history

    old_history = obj_hdrinfo['HISTORY']

    # remove it from the avehdr

    obj_hdrinfo.pop('HISTORY')

    #
    # Now fill in the values from the generic hdrinfo

    for key, value in hdrinfo.items():
        
        obj_hdrinfo[key] = value

                                
    # Deal with the units and plot labels
    
    obj_hdrinfo['YUNITS'][0] = tc.state['intensity_unit']

    latex = get_latex_fluxdensity(tc.state['intensity_unit']) 

    obj_hdrinfo['LYUNITS'][0] = latex[0]    
    obj_hdrinfo['LYLABEL'][0] = latex[1]
    obj_hdrinfo['LULABEL'][0] = latex[2]    

    # Add the S/N ratio 

    for i in range(tc.state['object_norders']):

        name = 'SNRO' + str(tc.state['object_orders'][i]).zfill(3)
        comment = ' Median S/N values for order ' + \
            str(tc.state['object_orders'][i]).zfill(3)

        values = []
        for j in range(tc.state['object_napertures']):        

            idx = i*tc.state['object_napertures']+j
        
            signal = tc.state['corrected_spectra'][idx,1,:]
            noise = tc.state['corrected_spectra'][idx,2,:]
            
            values.append(str(int(np.round(np.nanmedian(signal/noise)))))

        obj_hdrinfo[name] = [", ".join(values), comment]

    #
    # Create the FITS header
    #

    # Create the basic headers

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add our keywords
    
    keys = list(obj_hdrinfo.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            pass

        else:

            hdr[keys[i]] = (obj_hdrinfo[keys[i]][0], obj_hdrinfo[keys[i]][1])
    
    # Add the history

    for hist in old_history:

        hdr['HISTORY'] = hist

    #
    # Write the file out
    #

    correct_fullpath = os.path.join(setup.state['proc_path'],
                                    tc.state['output_filename']+'.fits')
    
    fits.writeto(correct_fullpath,
                 tc.state['corrected_spectra'],
                 hdr,
                 overwrite=True)

    logging.info(' Wrote file '+os.path.basename(correct_fullpath) + \
                 ' to the proc directory.')
    
    #
    # Do the QA plotting
    #

    if qa['show'] is True:

        figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                       setup.plots['landscape_size'][1]*qa['showscale'])

        font_size = setup.plots['font_size']*qa['showscale']
        
        plot_spectra(correct_fullpath,
                     ytype='flux and uncertainty',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(correct_fullpath),
                     showblock=qa['showblock'],
                     plot_number=setup.plots['abeam_spectra'],
                     figure_size=figure_size,
                     font_size=font_size,
                     colors=['green','black'])

        plot_spectra(correct_fullpath,
                     ytype='snr',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(correct_fullpath),
                     showblock=qa['showblock'],
                     plot_number=setup.plots['abeam_snr'],
                     figure_size=figure_size,
                     font_size=font_size,
                     colors=['green','black'])

        
    if qa['write'] is True:

        file_fullpath = os.path.join(setup.state['qa_path'],
                                     tc.state['output_filename']+\
                                     setup.state['qa_extension'])

        plot_spectra(correct_fullpath,
                     ytype='flux and uncertainty',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(correct_fullpath),
                     showblock=qa['showblock'],
                     output_fullpath=file_fullpath,
                     figure_size=setup.plots['landscape_size'],
                     font_size=setup.plots['font_size'],
                     colors=['green','black'])

        file_fullpath = os.path.join(setup.state['qa_path'],
                                     tc.state['output_filename']+\
                                     '_snr'+\
                                     setup.state['qa_extension'])

        plot_spectra(correct_fullpath,
                     ytype='snr',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(correct_fullpath),
                     showblock=qa['showblock'],
                     output_fullpath=file_fullpath,
                     figure_size=setup.plots['landscape_size'],
                     font_size=setup.plots['font_size'],
                     colors=['green','black'])



        if write_telluric_spectra is True:
        
            file_fullpath = os.path.join(setup.state['qa_path'],
                                         tc.state['output_filename']+\
                                         '_telluric'+\
                                         setup.state['qa_extension'])


            plot_spectra(telluric_fullpath,
                         ytype='flux and uncertainty',
                         spectrum_linewidth=setup.plots['spectrum_linewidth'],
                         spine_linewidth=setup.plots['spine_linewidth'],
                         title=os.path.basename(telluric_fullpath),
                         showblock=qa['showblock'],
                         output_fullpath=file_fullpath,
                         figure_size=setup.plots['landscape_size'],
                         font_size=setup.plots['font_size'],
                         colors=['green','black'])
