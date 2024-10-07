import numpy as np
import logging
from astropy.io import fits
import matplotlib.pyplot as pl
from os.path import join


from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.io.check import check_parameter, check_qakeywords, \
    check_file, check_sansfits

from pyspextool.io.files import files_to_fullpath
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io.fitsheader import get_headerinfo
from pyspextool.combine.core import plot_allorders

from pyspextool.pyspextoolerror import pySpextoolError


def load_spectra(files:str | list,
                 output_name:str,
                 verbose:bool=None,
                 qa_show:bool=None,
                 qa_showscale:float | int=None,
                 qa_showblock:bool=None,
                 qa_write:bool=None):


    """
    To load spectra to later be combined.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.

        e.g. ['spc', '1-2']

    output_name : str
        The output file name sans the suffix.
    
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
        the plot window which is (9,6).  This does NOT affect plots written
        to disk.  Set to None to default to setup.state['qa_scale'].    
    

    Returns
    -------
    None
        Loads data into memory

    
    """

    #
    # Check the parameters and QA keywords
    #

    check_parameter('load_spectra', 'files', files, ['str', 'list'])

    check_parameter('load_spectra', 'output_name', output_name, 'str')

    check_parameter('load_spectra', 'verbose', verbose, ['NoneType', 'bool'])
    
    check_parameter('load_spectra', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('load_spectra', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('load_spectra', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('load_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    logging.info(' Combining Spectra\n-----------------------\n')
    logging.info(' Loading the spectra.')

    
    #
    # Create the file names
    #

    results = files_to_fullpath(setup.state['proc_path'],
                                files,
                                setup.state['nint'],
                                '',
                                '.fits')

    input_files = results[0]
    file_readmode = results[1]
    filenames = results[2]

    check_file(input_files)

    combine.state['input_files'] = input_files
    combine.state['filenames'] = filenames
    combine.state['nfiles'] = len(input_files)

    check_sansfits(output_name,'output_name')
    combine.state['output_name'] = output_name
    
    #
    # Read the first file and store useful things
    #

    first_spectra, info = read_spectra_fits(combine.state['input_files'][0])

    combine.state['npixels'] = np.size(first_spectra[0,0,:])
    combine.state['module'] = info['module']    
    combine.state['orders'] = info['orders']
    combine.state['norders'] = len(info['orders'])
    combine.state['napertures'] = info['napertures']
    combine.state['xlabel'] = info['lxlabel']
    combine.state['ylabel'] = info['lylabel']    
    combine.state['xunits'] = info['xunits']


    #
    # Determine the combination parameters
    #

    # Start with the module that created the data

    if combine.state['module'] == 'extract':

        # Now check whether we are combining apertures.
        # 1 file means yes, >1 files means no.

        if combine.state['nfiles'] == 1:

            # Single file means combining all apertures

            combine.state['combine_apertures'] = True
            combine.state['combine_type'] = 'twoaperture'

        else:

            # Multiple file means combining on aperture by aperture basis

            combine.state['combine_apertures'] = False
            combine.state["combine_type"] = "standard"

    elif combine.state['module'] == 'telluric':

        combine.state['combine_apertures'] = False
        combine.state['combine_type'] = 'telluric'

    else:

        message = 'Unknown module.'
        raise pySpextoolError(message)

    #
    #  Compute various values and build the various arrays and lists.
    #

    if combine.state['combine_apertures'] is True:

        combine.state['final_napertures'] = 1
        combine.state['nspectra'] = 2*combine.state['nfiles']
        combine.state['scales'] = np.full(combine.state['nspectra'], 1)

    else:

        combine.state['final_napertures'] = combine.state['napertures']
        combine.state['nspectra'] = combine.state['nfiles']
        shape = (combine.state['nspectra'], combine.state['final_napertures'])
        combine.state['scales'] = np.full(shape,1)

    #
    # Create output arrays, wavelength first, then the others
    #

    shape = (combine.state['final_napertures'], combine.state['norders'],
             combine.state['npixels'])

    wavelengths = np.full(shape, np.nan)

    shape = (combine.state['final_napertures'], combine.state['norders'],
             combine.state['nspectra'], combine.state['npixels'])

    intensities = np.full(shape, np.nan)
    uncertainties = np.full(shape, np.nan)
    bitmasks = np.full(shape, 0, dtype=np.int8)

    #
    # Now start the loop over each file
    #

    headers = []

    keywords  = setup.state['combine_keywords']

    for i in range(combine.state['nfiles']):

        for j in range(combine.state['norders']):

            for k in range(combine.state['napertures']):

                # read the file

                hdul = fits.open(combine.state['input_files'][i])
                hdul[0].verify('silentfix')
                # this was needed to correct hdr problems
                spectra = hdul[0].data
                header = hdul[0].header
                hdul.close()

                # Grab header keywords and store

                info = get_headerinfo(header,keywords=keywords)
                headers.append(info)

                # store the data

                if combine.state['combine_apertures'] is False:

                    idx = j*combine.state['napertures']+k
                    wavelengths[k,j,:] = spectra[idx,0,:]
                    intensities[k,j,i,:] = spectra[idx,1,:]
                    uncertainties[k,j,i,:] = spectra[idx,2,:]
                    bitmasks[k,j,i,:] = np.nan_to_num(spectra[idx,3,:]).astype(np.int8)

                else:

                    idx = j*combine.state['napertures']+k
                    wavelengths[0,j,:] = spectra[idx,0,:]
                    intensities[0,j,k,:] = spectra[idx,1,:]
                    uncertainties[0,j,k,:] = spectra[idx,2,:]
                    bitmasks[0,j,k,:] = np.nan_to_num(spectra[idx,3,:]).astype(np.int8)

    #
    # Load into memory
    #
    
    combine.state['headers'] = headers
    combine.state['wavelengths'] = wavelengths
    combine.state['intensities'] = intensities
    combine.state['uncertainties'] = uncertainties
    combine.state['bitmasks'] = bitmasks
    combine.state['spectra_scaled'] = False
    
    
    #
    # Do the QA
    #

    if qa['write'] is True:

        plot_allorders(setup.plots['combine_spectra'],
                       setup.plots['landscape_size'],
                       setup.plots['font_size'],
                       setup.plots['spectrum_linewidth'],
                       setup.plots['spine_linewidth'],
                       filenames,
                       scalerange=None,
                       title='Raw Spectra')

        
        pl.savefig(join(setup.state['qa_path'],output_name+'_raw'+\
                   setup.state['qa_extension']))
        pl.close()

    if qa['show'] is True:

        scaled_size = (setup.plots['landscape_size'][0]*qa['showscale'],
                       setup.plots['landscape_size'][1]*qa['showscale'])

        scaled_font = setup.plots['font_size']*qa['showscale']
        
        plot_allorders(setup.plots['combine_spectra'],
                       scaled_size,
                       scaled_font,
                       setup.plots['spectrum_linewidth'],
                       setup.plots['spine_linewidth'],
                       filenames,
                       scalerange=None,                       
                       title='Raw Spectra')
        
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False: pl.pause(1)


    #
    # Set the done variables
    #

    combine.state['load_done'] = True
    combine.state['scale_done'] = False
        

    
    
    
    
    
