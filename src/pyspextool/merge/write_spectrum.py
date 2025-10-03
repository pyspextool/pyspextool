from astropy.io import fits
import os
import logging

from pyspextool import config as setup
from pyspextool.merge import config
from pyspextool.io.check import check_qakeywords
from pyspextool.io.apertures_to_array import apertures_to_array
from pyspextool.plot.plot_spectra import plot_spectra

def write_spectrum(
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float=None,
    qa_showblock:bool=None,
    qa_write:bool=None):
    
    """
    To write merged spectra to disk.

    Parameters
    ----------
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
    
    qa_showscale : float or int, default None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    Returns
    -------
    None

    """


    #
    # Check qa keywords
    #

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Conver the list of apertures to a single array suitable for FITS
    #

    result = apertures_to_array(config.state['merged_spectrum'],
                                1,
                                config.state['napertures'])

    #
    # Get a header going
    #

    hdrinfo = config.state["hdrinfo"]
    hdrinfo['NORDERS'][0] = 1
    hdrinfo['ORDERS'][0] = '1'
    hdrinfo['NAPS'][0] = config.state['napertures']

    old_history = hdrinfo['HISTORY']

    hdrinfo.pop('HISTORY')

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add our keywords
    
    keys = list(hdrinfo.keys())

    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            pass

        else:

            hdr[keys[i]] = (hdrinfo[keys[i]][0], hdrinfo[keys[i]][1])
    
    # Add the history

    for hist in old_history:

        hdr['HISTORY'] = hist

    #
    # Write the file out
    #

    fullpath = os.path.join(setup.state['proc_path'],
                            config.state["outputfile_root"]+'.fits')
    
    fits.writeto(fullpath,
                 result,
                 hdr,
                 overwrite=True)

    logging.info(' Wrote file '+os.path.basename(fullpath) + \
                 ' to the proc directory.')

    #
    # Do the QA plotting
    #

    if qa['show'] is True:

        figure_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                       setup.plots['landscape_size'][1]*qa['showscale'])

        font_size = setup.plots['font_size']*qa['showscale']
        
        plot_spectra(fullpath,
                     ytype='flux and uncertainty',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(fullpath),
                     showblock=qa['showblock'],
                     plot_number=setup.plots['merged'],
                     figure_size=figure_size,
                     font_size=font_size,
                     colors=['green','black'])


    if qa['write'] is True:

        file_fullpath = os.path.join(setup.state['qa_path'],
                                     config.state["outputfile_root"]+
                                     setup.state['qa_extension'])

        plot_spectra(fullpath,
                     ytype='flux and uncertainty',
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],            
                     title=os.path.basename(fullpath),
                     showblock=qa['showblock'],
                     output_fullpath=file_fullpath,
                     figure_size=setup.plots['landscape_size'],
                     font_size=setup.plots['font_size'],
                     colors=['green','black'])


