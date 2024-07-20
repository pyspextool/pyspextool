from astropy.io import fits
import os
import logging

from pyspextool import config as setup
from pyspextool.telluric import config as telluric
from pyspextool.io.check import check_parameter, check_keywords
from pyspextool.utils.units import get_latex_fluxdensity
from pyspextool.plot.plot_spectra import plot_spectra

def write_spectra(write_model_spectra:bool=False,
                  write_telluric_spectra:bool=True,
                  verbose:bool=None,
                  qa_show:bool=None,
                  qa_scale:float |int =None,
                  qa_block:bool=None,
                  qa_write:bool=None,
                  overwrite:bool=True):

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
    
    qa_block : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_scale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window which is (9,6).  This does affect plots written to disk.
        Set to None to default to setup.state['qa_scale'].

    overwrite : {True, False, None}
        Set to True to overwrite FITS file on disk.
        Set to False to not overwrite FITS file on disk.
        Set to None to default to setup.state['overwrite'].
        
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

    check_parameter('write_spectra', 'qa_scale', qa_scale,
                    ['NoneType','float','int'])
    
    check_parameter('write_spectra', 'qa_block', qa_block, ['NoneType','bool'])
        
    check_parameter('write_spectra', 'qa_write', qa_write, ['NoneType','bool'])

    check_parameter('write_spectra', 'overwrite', overwrite, 'bool')    

    keywords = check_keywords(verbose=verbose, qa_show=qa_show,
                              qa_scale=qa_scale, qa_block=qa_block,
                              qa_write=qa_write, overwrite=overwrite)    

    #
    # Write the (Vega) model to disk if requested.
    #

    if write_model_spectra is True:

        # Create the header. Make a new hdrinfo dictionary using the standard
        # hdrinfo dictionary
        
        std_hdrinfo = telluric.state['standard_hdrinfo']

        new_hdrinfo = {}
        
        new_hdrinfo['NORDERS'] = std_hdrinfo['NORDERS']
        new_hdrinfo['ORDERS'] = std_hdrinfo['ORDERS']
        new_hdrinfo['NAPS'] = std_hdrinfo['NAPS']                
        new_hdrinfo['MODULE']=std_hdrinfo['MODULE']
        new_hdrinfo['MODULE'][0] = 'telluric'
        new_hdrinfo['VERSION'] = std_hdrinfo['VERSION']
        filename = [telluric.load['output_filename']+'_vega.fits', ' File name']
        new_hdrinfo['FILENAME'] = filename

        f = '{:.3f}'
        mag = float(f.format(telluric.state['standard_bmag']))
        new_hdrinfo['STD_BMAG'] = (mag, 'Telluric standard B mag')
        mag = float(f.format(telluric.state['standard_vmag']))
        new_hdrinfo['STD_VMAG'] = (mag, 'Telluric standard V mag')

        # Now deal with the units

        latex = get_latex_fluxdensity(telluric.load['intensity_unit']) 
        
        new_hdrinfo['XUNITS'] = std_hdrinfo['XUNITS']
        new_hdrinfo['YUNITS'] = std_hdrinfo['YUNITS']
        new_hdrinfo['YUNITS'][0] = telluric.load['intensity_unit']
        new_hdrinfo['LXUNITS'] = std_hdrinfo['LXUNITS']
        new_hdrinfo['LYUNITS'] = std_hdrinfo['LYUNITS']        
        new_hdrinfo['LYUNITS'][0] = latex[0]
        new_hdrinfo['LXLABEL'] = std_hdrinfo['LXLABEL']
        new_hdrinfo['LYLABEL'] = std_hdrinfo['LYLABEL']        
        new_hdrinfo['LYLABEL'][0] = latex[1]

        # Create the basic header

        phdu = fits.PrimaryHDU()
        newhdr = phdu.header

        # Add our keywords
    
        keys = list(new_hdrinfo.keys())
    
        for i in range(len(keys)):
                            
                newhdr[keys[i]] = (new_hdrinfo[keys[i]][0],
                                   new_hdrinfo[keys[i]][1])
        
        #
        # Write the file out
        #

        full_path = os.path.join(setup.state['proc_path'],
                                 telluric.load['output_filename']+\
                                 '_vega.fits')    

        fits.writeto(full_path, telluric.state['vega_spectra'], newhdr,
                 overwrite=keywords['overwrite'])
        
        logging.info(' Wrote file '+os.path.basename(full_path) + ' to disk.')


    
    if write_telluric_spectra is True:

        hdr = telluric.state['standard_hdrinfo']

        # Store the history

        old_history = hdr['HISTORY']

        # remove it from the avehdr

        hdr.pop('HISTORY')

        # Add things to it
        
        hdr['MODULE'][0] = 'telluric'
        
        hdr['FILENAME'][0] = telluric.load['output_filename']+'_telluric.fits'
        
        hdr['STDFILE'] = [telluric.load['standard_file'],
                          'Telluric standard input file']
        
        hdr['STD_ID'] = [telluric.state['standard_name'],'Telluric standard']

        f = '{:.3f}'
        hdr['STD_BMAG'] = [float(f.format(telluric.state['standard_bmag'])),
                           'Telluric standard B mag']

        hdr['STD_VMAG'] = [float(f.format(telluric.state['standard_vmag'])),
                       'Telluric standard V mag']

        # Deal with the units and plot labels
    
        hdr['YUNITS'][0] = telluric.load['intensity_unit']+' / DN s-1'

        lylabel = get_latex_fluxdensity(telluric.load['intensity_unit'])[0]

        hdr['LYLABEL'][0] = 'Ratio ('+lylabel+' / DN s$^{-1}$)'
        
        # Create the basic headers

        phdu = fits.PrimaryHDU()
        newhdr = phdu.header

        # Add our keywords
    
        keys = list(hdr.keys())
    
        for i in range(len(keys)):
            
            if keys[i] == 'COMMENT':
                
                junk = 1
                
            else:
                
                newhdr[keys[i]] = (hdr[keys[i]][0], hdr[keys[i]][1])
                
        # Do the history
        
        for hist in old_history:

            newhdr['HISTORY'] = hist

        #
        # Write the file out
        #

        full_path = os.path.join(setup.state['proc_path'],
                                 telluric.load['output_filename']+\
                                 '_telluric.fits')    

        fits.writeto(full_path, telluric.state['telluric_spectra'], newhdr,
                 overwrite=keywords['overwrite'])
        
        logging.info(' Wrote file '+os.path.basename(full_path) + ' to disk.')
        
    
    #
    # Write the corrected spectra to disk
    #
    
    # Rename header for ease of reading
            
    hdr = telluric.state['object_hdrinfo']
                
    # Store the history

    old_history = hdr['HISTORY']

    # remove it from the avehdr

    hdr.pop('HISTORY')

    #
    # Add things to the header
    #

    hdr['MODULE'][0] = 'telluric'

    hdr['FILENAME'][0] = telluric.load['output_filename']+'.fits'

    hdr['TC_OFILE'] = [telluric.load['object_file'],
                      'Telluric object input file']

    hdr['TC_SFILE'] = [telluric.load['standard_file'],
                      'Telluric standard input file']

    hdr['TC_STDID'] = [telluric.state['standard_name'],'Telluric standard']


    hdr['TC_STDB'] = [float('{:.3f}'.format(telluric.state['standard_bmag'])),
                       'Telluric standard B mag']

    hdr['TC_STDV'] = [float('{:.3f}'.format(telluric.state['standard_vmag'])),
                       'Telluric standard V mag']

    hdr['TC_STDRV'] =  [float('{:.2f}'.format(telluric.state['standard_rv'])),
                      'Telluric standard radial velocity (km s-1)']

    hdr['TC_dAM'] = [float('{:.2f}'.format(telluric.state['delta_airmass'])),\
                       'Telluric Average airmass difference']    

    hdr['TC_METH'] = [telluric.state['method'],'Telluric method']

    if telluric.state['method'] == 'deconvolution':
    
        hdr['TC_MXDEV'] = [float('{:.5f}'.format(telluric.state['max_deviation'])),
                           'Telluric maxmimum % deviation of Vega-data']

        hdr['TC_RMS'] = [float('{:.5f}'.format(telluric.state['rms_deviation'])),
                         'Telluric RMS deviation of Vega-data']

    # Deal with the units and plot labels
    
    hdr['YUNITS'][0] = telluric.load['intensity_unit']

    latex = get_latex_fluxdensity(telluric.load['intensity_unit']) 

    hdr['LYUNITS'][0] = latex[0]    
    hdr['LYLABEL'][0] = latex[1]
    
    #
    # Create the header
    #

    # Create the basic headers

    phdu = fits.PrimaryHDU()
    newhdr = phdu.header

    # Add our keywords
    
    keys = list(hdr.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            newhdr[keys[i]] = (hdr[keys[i]][0], hdr[keys[i]][1])
    
    # Do the history

    for hist in old_history:

        newhdr['HISTORY'] = hist

    #
    # Write the file out
    #

    full_path = os.path.join(setup.state['proc_path'],
                             telluric.load['output_filename']+'.fits')
    
    fits.writeto(full_path, telluric.state['corrected_spectra'], newhdr,
                 overwrite=keywords['overwrite'])

    logging.info(' Wrote file '+os.path.basename(full_path) + ' to disk.')
        
    #
    # Do the QA plotting
    #

    if keywords['qa_show'] is True:

        plot_spectra(full_path, ytype='flux and uncertainty',
                     line_width=0.5,
                     title=os.path.basename(full_path),
                     block=keywords['qa_block'],
                     colors=['green','black'])
        
    if keywords['qa_write'] is True:

        qafileinfo = {'figsize': (6,4),
                      'filepath': setup.state['qa_path'],
                      'filename': telluric.load['output_filename'],
                      'extension': setup.state['qa_extension']}

        plot_spectra(full_path, ytype='flux and uncertainty',
                     order_numbers=False,
                     line_width=0.5,colors=['green','black'],
                     title=os.path.basename(full_path),
                     file_info=qafileinfo)

    
