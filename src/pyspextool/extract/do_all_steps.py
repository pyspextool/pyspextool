import numpy as np
import os

from pyspextool.extract import config
from pyspextool.extract.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.io.files import extract_filestring
from pyspextool.extract.load_image import load_image
from pyspextool.extract.make_spatial_profiles import make_spatial_profiles
from pyspextool.extract.locate_aperture_positions import locate_aperture_positions
from pyspextool.extract.select_orders import select_orders
from pyspextool.extract.trace_apertures import trace_apertures
from pyspextool.extract.define_aperture_parameters import define_aperture_parameters
from pyspextool.extract.extract_apertures import extract_apertures

def do_all_steps(files):

    """
    To extract spectra from images after parameters are set.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.



    Returns
    -------
    None.  Writes FITS files to disk.

    """

    #
    # Check parameter
    #

    check_parameter('do_all_steps', 'files', files, ['str', 'list'])
    
    #
    # Figure out how many files we are talking about
    #

    if config.state['filereadmode'] == 'filename':

        files = extract_filestring(files, 'filename')
        
    else:

        files = extract_filestring(files[1], 'index')

    nfiles = len(files)
    
    # Test for even-ity if the mode is A-B
    
    if nfiles %2 != 0 and config.state['reductionmode'] == 'A-B':

        message = 'The number of images must be even.'
        raise ValueError(message)

    if config.state['reductionmode'] == 'A-B':

        nloop = int(nfiles/2)

    else:
        
        nloop = int(nfiles)

    # Start the loop

    for i in range(nloop):

        if config.state['reductionmode'] == 'A-B':

            subset = files[i*2:i*2+2]

        if config.state['reductionmode'] == 'A':

            subset = files[i]

        # Load the data
            
        load_image([config.state['prefix'],subset],
                   config.user['load']['flatfile'],
                   config.user['load']['wavecalfile'],
                   flat_field=config.user['load']['doflat'],
                   linearity_correction=config.user['load']['doflat'],
                   qaplot=config.user['load']['qaplot'],
                   qafile=config.user['load']['qafile'],
                   reduction_mode=config.state['reductionmode'],
                   do_all_steps=True, verbose=config.user['load']['verbose'])

        # Set continue to 2 to pass the set_extraction_type step
        config.state['continue'] = 2
        
        # Make the Profiles

        make_spatial_profiles(verbose=config.user['profiles']['verbose'],
                              qaplot=config.user['profiles']['qaplot'],
                              qafile=config.user['profiles']['qafile'])

        # Locate the Aperture Positions

        locate_aperture_positions(config.user['locateaps']['apertures'],
                                  method=config.user['locateaps']['method'],
                                  qaplot=config.user['locateaps']['qaplot'],
                                  qafile=config.user['locateaps']['qafile'],
                                  verbose=config.user['locateaps']['verbose'])

        # Select orders

        select_orders(include=config.user['orders']['include'],
                      exclude=config.user['orders']['exclude'],
                      include_all=config.user['orders']['include_all'],
                      verbose=config.user['orders']['verbose'],
                      qaplot=config.user['orders']['qaplot'],
                      qafile=config.user['orders']['qafile'])
        

        # Trace apertures

        trace_apertures(fit_degree=config.user['trace']['fitdeg'],
                        step_size=config.user['trace']['step'],
                        summation_width=config.user['trace']['sumwidth'],
                        centroid_threshold=config.user['trace']['centhresh'],
                        fwhm=config.user['trace']['fwhm'],
                        verbose=config.user['trace']['verbose'],
                        qaplot=config.user['trace']['qaplot'],
                        qafile=config.user['trace']['qafile'])

        # Define aperture parameters

        define_aperture_parameters(config.user['apparms']['apradii'],
                                psf_radius=config.user['apparms']['psfradius'],
                                bg_radius=config.user['apparms']['bgradius'],
                                bg_width=config.user['apparms']['bgwidth'],
                                bg_regions=config.user['apparms']['bgregions'],
                                bg_fit_degree=config.user['apparms']['bgdeg'],
                                qaplot=config.user['apparms']['qaplot'],
                                qafile=config.user['apparms']['qafile'])

        # Extract apertures

        extract_apertures(verbose=config.user['extract']['verbose'])

        
    #
    # Things proceed differently depending on whether you are extracting a
    # point source or an extended source
    #



        
