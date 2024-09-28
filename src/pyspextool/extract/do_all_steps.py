import logging

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter
from pyspextool.io.files import extract_filestring
from pyspextool.extract.load_image import load_image
from pyspextool.extract.make_profiles import make_profiles
from pyspextool.extract.identify_apertures import identify_apertures
from pyspextool.extract.select_orders import select_orders
from pyspextool.extract.trace_apertures import trace_apertures
from pyspextool.extract.define_aperture_parameters import define_aperture_parameters
from pyspextool.extract.extract_apertures import extract_apertures


def do_all_steps(files, verbose=None):
    """
    To extract spectra from images in a loop after parameters are set.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the prefix.
        files[1] is a string giving the index numbers of the files.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file.  

    Returns 
    -------
    list
    A list of str giving the names of the files successfully written to disk.

    """

    #
    # Check if we can proceed
    #

    if extract.state['extract_done'] is False:

        
        message = "extract.state['extract_done'] is False.  "+\
                  "Previous steps not completed."
        raise ValueError(message)

    #
    # Check parameter
    #

    check_parameter('do_all_steps', 'files', files, ['str', 'list'])

    check_parameter('do_all_steps', 'verbose', verbose, ['NoneType','bool'])    


    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if verbose is None:
        verbose = setup.state['verbose']

    if verbose is True:
        logging.getLogger().setLevel(logging.INFO)
        setup.state["verbose"] = True
        
    elif verbose is False:
        logging.getLogger().setLevel(logging.ERROR)
        setup.state["verbose"] = False
    
    #
    # Figure out how many files we are talking about
    #

    if extract.state['filereadmode'] == 'filename':

        files = extract_filestring(files, 'filename')

    else:

        files = extract_filestring(files[1], 'index')

    nfiles = len(files)

    # Test for even-ity if the mode is A-B

    if nfiles % 2 != 0 and extract.state['reductionmode'] == 'A-B':
        message = 'The number of images must be even.'
        raise ValueError(message)

    if extract.state['reductionmode'] == 'A-B':

        nloop = int(nfiles / 2)

    else:

        nloop = int(nfiles)

    #
    # Start the loop
    #

    good_files = []
        
    for i in range(nloop):

        if extract.state['reductionmode'] == 'A-B':
            subset = files[i * 2:i * 2 + 2]

        if extract.state['reductionmode'] == 'A':
            subset = files[i]

        #
        # Load the data
        #

        load_image([extract.state['prefix'], subset], extract.load['flatfile'],
                   extract.load['wavecalfile'],
                   flat_field=extract.load['doflat'],
                   linearity_correction=extract.load['doflat'],
                   qa_show=extract.load['qaplot'],
                   qa_write=extract.load['qafile'],
                   qa_showsize=extract.load['qaplotsize'],
                   reduction_mode=extract.state['reductionmode'],
                   do_all_steps=True, verbose=extract.load['verbose'])

        #
        # Set the extraction type
        #

#        set_extraction_type(extract.type['type'])

        #
        # Make the Profiles
        #

#        make_profiles(verbose=extract.profiles['verbose'],
#                      qa_show=extract.profiles['qaplot'],
#                      qa_write=extract.profiles['qafile'],
#                      qa_showsize=extract.profiles['qaplotsize'])

        #
        # Locate the Aperture Positions
        #

        locate_aperture_positions(extract.apertures['apertures'],
                                  method=extract.apertures['method'],
                                  qa_show=extract.apertures['qaplot'],
                                  qa_showsize=extract.apertures['qaplotsize'],
                                  qa_write=extract.apertures['qafile'],
                                  verbose=extract.apertures['verbose'])
            
        #
        # Select orders
        #

        select_orders(include=extract.orders['include'],
                      exclude=extract.orders['exclude'],
                      include_all=extract.orders['include_all'],
                      verbose=extract.orders['verbose'],
                      qa_show=extract.orders['qaplot'],
                      qa_showsize=extract.orders['qaplotsize'],
                      qa_write=extract.orders['qafile'])


        #
        # Trace apertures
        #

        trace_apertures(fit_degree=extract.trace['fitdeg'],
                        step_size=extract.trace['step'],
                        summation_width=extract.trace['sumwidth'],
                        centroid_threshold=extract.trace['centhresh'],
                        fwhm=extract.trace['fwhm'],
                        verbose=extract.trace['verbose'],
                        qa_show=extract.trace['qaplot'],
                        qa_showsize=extract.trace['qaplotsize'],
                        qa_write=extract.trace['qafile'])

        #
        # Define aperture parameters
        #


        try:
            
            define_aperture_parameters(extract.parameters['apradii'],
                                   psf_radius=extract.parameters['psfradius'],
                                   bg_radius=extract.parameters['bgradius'],
                                   bg_width=extract.parameters['bgwidth'],
                                   bg_regions=extract.parameters['bgregions'],
                                   bg_fit_degree=extract.parameters['bgdeg'],
                                   qa_show=extract.parameters['qaplot'],
                                   qa_showsize=extract.parameters['qaplotsize'],
                                   qa_write=extract.parameters['qafile'])

        except(ValueError) as e:

            message = f"\n\n\nEncountered error `{e}` in "+\
              "define_aperture_parameters.  Moving on to next extraction.\n\n\n"
            
            logging.info(message)
            
            continue
            
        #
        # Extract apertures
        #

        filenames = extract_apertures(verbose=extract.extract['verbose'])
        
        #
        # Store the successful file names
        #

        for item in filenames: good_files.append(item)
            
                
    logging.info(f" Do All Steps Complete.")

    return good_files

        
    #
    # Things proceed differently depending on whether you are extracting a
    # point source or an extended source
    #
