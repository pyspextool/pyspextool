import os


from pyspextool.extract import load_image
from pyspextool.extract import make_profiles
from pyspextool.extract import identify_apertures
from pyspextool.extract import select_orders
from pyspextool.extract import override_aperturesigns
from pyspextool.extract import trace_apertures
from pyspextool.extract import define_aperture_parameters
from pyspextool.extract import extract_apertures

from pyspextool import config as setup
from pyspextool.extract import config as extract


from pyspextool.io.files import files_to_fullpath
from pyspextool.pyspextoolerror import pySpextoolError



def extract(reduction_mode:str,
            files:str | list,
            flat_name:str,
            wavecal_name:str,
            aperture_findinfo:list,
            aperture_radii:int | float | list,
            output_filenames:str=None,
            output_prefix:str='spectra',
            input_extension:str='.fits*',
            load_directory='raw',
            flat_field=True,
            linearity_correction=True,
            detector_info:dict=None,
            write_rectified_orders:bool=False,
            seeing_fwhm:int | float=0.8,
            profile_ybuffer:int=3,
            aperture_signs:list=None,
            include_orders:int | str | list=None,
            exclude_orders:int | str | list=None,
            trace_fitdegree:int=2,
            trace_stepsize:int=5,
            trace_summationwidth:int=5,
            trace_centroidthreshold:int=2,
            bg_annulus:list=None,
            bg_regions:str=None,
            bg_fitdegree:int=0,
            psf_radius:int | float=None,
            fix_badpixels:bool=True,
            use_meanprofile:bool=False,
            badpixel_thresh:int | float=7,                    
            verbose:bool=None,
            qa_show:bool=None,
            qa_showscale:float | int=None,
            qa_showblock:bool=None,
            qa_write:bool=None):

    """



    """



    #
    # Let's get set up for the for loop.  Build the input and output file names.
    #
    
    # Determine the path in which `files` is located.  
    
    if load_directory == 'raw':

        load_path = setup.state['raw_path']
#        image_type = 'raw'
        
    elif load_directory == 'proc':

        load_path = setup.state['proc_path']
#        image_type = 'combined'
    

    # Create the input file names

    results = files_to_fullpath(load_path,
                                files,
                                setup.state['nint'],
                                setup.state['suffix'],
                                input_extension)

    input_fullpaths = results[0]
    file_readmode = results[1]
    n_inputfiles = len(input_fullpaths)

    # Create the full output file names

    if output_filenames is not None:

        # The user wants to use their own file names

        # Convert to a list if a str
        
        if isinstance(output_filenames,str) is True:
            
            output_filenames = [output_filenames]
            
        # Check to make sure the number of files giveb equals the number of
        # files we are loading
            
        if n_inputfiles != len(output_filenames):

            message = 'The number of files in the keyword '+ \
                '`output_filenames` does not match the number of files in '+ \
                'the parameter `files`.'
            raise pySpextoolError(message)

    else:

        # The user wants them created by us.  Do it differently depending
        # on what readmode was selected.
        
        if file_readmode == 'index':
        
            result = files_to_fullpath('',
                                       [output_prefix,files[1]],
                                       setup.state['nint'],
                                       '',
                                       '',
                                       exist=False)

            output_filenames = result[0]        

        else:

            output_filenames = []
            for file in files:

                root = os.path.splitext(os.path.basename(file))
                if root[1] == '.gz':
                    root = os.path.splitext(root[0])
                    
                output_filenames.append(output_prefix+'_'+root[0]+'.fits')

    # If the reduction mode is A-B, is there an even number of files?

    if n_inputfiles % 2 != 0 and reduction_mode == 'A-B':

        message = "The number of images must be even when "+ \
            "`reduction_mode`='A-B'."
        raise pySpextoolError(message)

    if reduction_mode.upper() == 'A-B':

        nloop = int(n_inputfiles / 2)

    else:

        nloop = int(n_inputfiles)

    #
    # start the loop
    #
        
    for i in range(nloop):

        # Load the image

        do_all_steps = False if i == 0 else True


        if reduction_mode.upper() == 'A':

            input_subset = [input_fullpaths[i]]
            output_subset = [output_filenames[i]]            
        
        if reduction_mode.upper() == 'A-B':

            input_subset = input_fullpaths[i * 2:i * 2 + 2]
            output_subset = output_filenames[i * 2:i * 2 + 2]
                    
        input_subset = ','.join(input_subset)
        
        load_image(input_subset,
                   flat_name,
                   wavecal_name,
                   output_filenames=output_subset,
                   output_prefix=output_prefix,
                   input_extension=input_extension,
                   load_directory=load_directory,
                   flat_field=flat_field,
                   linearity_correction=linearity_correction,
                   detector_info=detector_info,
                   write_rectified_orders=write_rectified_orders,
                   do_all_steps=do_all_steps,
                   verbose=verbose,
                   qa_show=qa_show,
                   qa_showscale=qa_showscale,
                   qa_showblock=qa_showblock,
                   qa_write=qa_write)

        # Make the profile

        make_profiles(verbose=verbose,
                      qa_show=qa_show,
                      qa_showscale=qa_showscale,
                      qa_showblock=qa_showblock,
                      qa_write=qa_write)

        # Identify apertures

        identify_apertures(aperture_findinfo,
                           seeing_fwhm=seeing_fwhm,
                           ybuffer=profile_ybuffer,
                           verbose=verbose,
                           qa_show=qa_show,
                           qa_showscale=qa_showscale,
                           qa_showblock=qa_showblock,
                           qa_write=qa_write)

        # Select orders

        if include_orders is not None or exclude_orders is not None:

            select_orders(include=include_orders,
                          exclude=exclude_orders,
                          verbose=verbose,
                          qa_show=qa_show,
                          qa_showscale=qa_showscale,
                          qa_showblock=qa_showblock,
                          qa_write=qa_write)

        # Update aperture positions if requested

        if aperture_signs is not None:

            override_aperturesigns(aperture_signs,
                                   verbose=verbose)

        # Trace the orders

        trace_apertures(fit_degree=trace_fitdegree,
                        step_size=trace_stepsize,
                        summation_width=trace_summationwidth,
                        centroid_threshold=trace_centroidthreshold,
                        seeing_fwhm=seeing_fwhm,
                        verbose=verbose,
                        qa_show=qa_show,
                        qa_showscale=qa_showscale,
                        qa_showblock=qa_showblock,
                        qa_write=qa_write)

        # Define the aperture parameters

        define_aperture_parameters(aperture_radii,
                                   bg_annulus=bg_annulus,
                                   bg_regions=bg_regions,
                                   bg_fit_degree=bg_fitdegree,
                                   psf_radius=psf_radius,
                                   verbose=verbose,
                                   qa_show=qa_show,
                                   qa_showscale=qa_showscale,
                                   qa_showblock=qa_showblock,
                                   qa_write=qa_write)

        # Extract the apertures

        extract_apertures(fix_badpixels=fix_badpixels,
                          use_meanprofile=use_meanprofile,
                          badpixel_thresh=badpixel_thresh,
                          verbose=verbose,
                          qa_show=qa_show,
                          qa_showscale=qa_showscale,
                          qa_showblock=qa_showblock,
                          qa_write=qa_write)





        

        
    

    
