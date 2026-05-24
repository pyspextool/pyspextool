import os

from pyspextool import config as setup

from pyspextool.extract import load_image
from pyspextool.extract import make_profiles
from pyspextool.extract import identify_apertures
from pyspextool.extract import select_orders
from pyspextool.extract import override_aperturesigns
from pyspextool.extract import trace_apertures
from pyspextool.extract import define_aperture_parameters
from pyspextool.extract import extract_apertures
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.io.files import files_to_fullpath
from pyspextool.pyspextoolerror import pySpextoolError

def extract(
    reduction_mode:str,
    filenames:str | list,
    flat_name:str,
    wavecal_name:str,
    aperture_findinfo:list,
    aperture_radii:int | float | list,
    output_filenames:str=None,
    output_prefix:str='spectra',
    input_suffix:str=setup.state['input_suffixes'][0],
    load_directory='raw',
    flat_field=True,
    linearity_correction=True,
    detector_info:dict=None,
    rectification_method:str='cubic',
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
    badpixel_threshold:int | float=7,                    
    verbose:bool=None,
    qa_show:bool=None,
    qa_showscale:float | int=None,
    qa_showblock:bool=None,
    qa_write:bool=None):

    """
    To perform all steps necessary to extract spectra.

    Parameters
    ----------
    reduction_mode : setup.state['reduction_modes']
        The reduction mode to be used.
        If 'A', then no pair subtraction is done.
        If 'A-B', then each pair of images is subtracted.

    filenames : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
    
    flat_name : str
        The full name of a pySpextool flat file.

    wavecal_name : str or None
        The full name of a pySpextool wavecal file.

    aperture_findinfo : list
        A (2,) list if

        list[0] = 'auto' then list[1] is the number of apertures to
        search for.

        list[0] = 'guess' then list[1] is a (2,) list giving guess positions

        list[0] = 'fixed' then list[1] is a (2,) list giving positions

    aperture_radii : int, float, list
        If int or float, then the aperture radius to use for all apertures.
        If list, then a (naps,) list of aperture radii.  

    output_filenames : str, list str, optional
        A str or list of str of output files names.  Only required if 
        `filenames` gives file names instead of a prefix and index numbers.
 
    output_prefix : str, default = 'spectra'
        The prefix for output spectral files if `filenames` is a (2,) list 
        with a prefix and index numbers.

    input_suffix : setup.state['reduction_modes'], default='.fits*'
        The suffix for raw data files.

    load_directory : {'raw', 'proc'}
        The directory in which to look for raw data.

    flat_field : {True, False}
        Set to True to flat field the data.
        Set to False to not flat field the data.

    linearity_correction : {True, False}
        Set to True to correct raw images for non-linearity.
        Set to False to not correct raw images for non-linearity.

    detector_info : dict, optional
        A dictionary with any information that needs to be passed to the
        instrument specific readfits program.  

    'rectification_method' : setup.state['rectification_methods']
        If 'linear', then bi-linear interpolation is used.  See 
        scipy.interpolate.RegularGridInterpolator.

        If 'cubic', k=3 tensor-product spline is used.  See
        scipy.interpolate.RegularGridInterpolator.

    write_rectified_orders: {False, True}
        If False, then the rectified orders are not written to disk.
        If True, then the rectified orders are written to disk.

    seeing_fwhm : float, default = 0.8 (arcseconds).
        The approximate FWHM of the peak to be identified.  Only used 
        if `aperture_findinfo`[0] is 'auto' or 'guess'.

    profile_ybuffer : int, default=3
        The number of pixels on the edge of the orders to ignore.  Useful
        as sometimes there is a downturn that can mess with the finding
        routine.

    aperture_signs : list, optional
        If given, an (naps,) list of aperture signs, e.g. '-','+' that will
        override results computed by the software.

    include_orders : int, list, str, optional
        If the type is int, the single order to include.
        If the type is list, a list of integer orders to include.
        If the type is str, a str giving the orders, e.g. '1-3,4,5'.

    exclude_orders : int, list, str, optional
        If the type is int, the single order to exclude.  
        If the type is list, a list of integer orders to include.
        If the type is str, a str giving the orders, e.g. '1-3,4,5'.

    trace_fitdegree : int, default=2
        The polynomial degree for the fit.

    trace_stepsize : int, default=5
        The step size as the function moves across the array identifying
        the peaks in the spatial dimension.

    trace_summationwidth : int, default=5
        The number of columns to combine in order to increase the S/N
        of the data fit to identify the peaks in the spatial dimension.

    trace_centroidthreshold : int, default=2
        If (fit-guess) > centroid_threshold to peak identification is found
        to fail.
    
    bg_annulus : list, optional
        A (2,) list giving the bg start radius and bg width (arcseconds).

    bg_regions : str, optional
        A string giving bg regions, e.g. '1-2,14,15'.

    bg_fitdegree: int, default=0
        The polynomial degree of the background fit.

    psf_radius:int, float, optional
        The PSF radius (in arcsconds) used in optimal extraction.  

    fix_badpixels : {True, False}
        Set to True to fix bad pixels.  Only used when not optimal extraction.
        Set to False to not fix bad pixels.  Only used when not optimal 
        extraction.

    use_meanprofile {False, True}
        Set to False to use the wavelength-dependent spatial profile for 
        bad pixel finding and/or optimal extraction.

        Set to True to use the mean of the wavelength-dependent spatial 
        profile for bad pixel finding and/or optimal extraction.

    badpixel_threshold : int, float, default=7
        The sigma threshold used to identify bad pixels.
       
    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, optional
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    """

    #
    # Check the parameters and QA keywords
    #

    check_parameter('extract', 'reduction_mode', 
                    reduction_mode, 'str', 
                    possible_values=setup.state['reduction_modes'])

    check_parameter('extract', 'filenames', 
                    filenames, ['str', 'list'], list_types=['str','str'])

    check_parameter('extract', 'flat_name', 
                    flat_name, 'str')

    check_parameter('extract', 'wavecal_name', 
                    wavecal_name, ['NoneType','str'])

    check_parameter('extract', 'aperture_findinfo', 
                    aperture_findinfo, 'list')

    check_parameter('extract', 'aperture_radii', 
                    aperture_radii, ['int', 'float','list'])

    check_parameter('extract', 'output_filenames', 
                    output_filenames, ['NoneType','str','list'])

    check_parameter('extract', 'output_prefix', 
                    output_prefix, 'str')

    check_parameter('extract', 'input_suffix', 
                    input_suffix, 'str', 
                    possible_values=setup.state['input_suffixes'])

    check_parameter('extract', 'load_directory', 
                    load_directory, 'str', possible_values=['raw','proc'])

    check_parameter('extract', 'flat_field', 
                    flat_field, 'bool')

    check_parameter('extract', 'linearity_correction', 
                    linearity_correction, 'bool')

    check_parameter('extract', 'detector_info', 
                    detector_info, ['dict','NoneType'])

    check_parameter('extract', 'rectification_method', 
                    rectification_method, 'str', 
                    possible_values=setup.state['rectification_methods'])

    check_parameter('extract', 'write_rectified_orders', 
                    write_rectified_orders, 'bool')

    check_parameter('extract', 'seeing_fwhm', 
                    seeing_fwhm, ['int','float'])

    check_parameter('extract', 'profile_ybuffer', 
                    profile_ybuffer, 'int')

    check_parameter('extract', 'aperture_signs', 
                    aperture_signs, ['list', 'NoneType'])

    check_parameter('extract', 'include_orders', 
                    include_orders, ['NoneType', 'int', 'list', 'str'])

    check_parameter('extract', 'exclude_orders', 
                    exclude_orders, ['NoneType', 'int', 'list', 'str'])    

    check_parameter("extract", "trace_fitdegree", 
                    trace_fitdegree, "int")
    
    check_parameter("extract", "trace_stepsize", 
                    trace_stepsize, "int")
    
    check_parameter("extract", "trace_summationwidth", 
                    trace_summationwidth, "int")

    check_parameter("extract", "trace_centroidthreshold", 
                    trace_centroidthreshold, "int")

    check_parameter('extract', 'bg_annulus', 
                    bg_annulus, ['ndarray', 'list','NoneType'])    

    check_parameter('extract', 'bg_regions', 
                    bg_regions, ['str','NoneType'])    

    check_parameter('extract', 'bg_fitdegree', 
                    bg_fitdegree, ['int','NoneType'])    
    
    check_parameter('extract', 'psf_radius', 
                    psf_radius, ['int','float', 'NoneType'])

    check_parameter('extract', 'fix_badpixels', 
                    fix_badpixels, 'bool')

    check_parameter('extract', 'use_meanprofile', 
                    use_meanprofile, 'bool')

    check_parameter('extract', 'badpixel_threshold', 
                    badpixel_threshold, ['int','float'])

    check_parameter('extract', 'verbose', 
                    verbose, ['NoneType', 'bool'])
    
    check_parameter('extract', 'qa_write', 
                    qa_write, ['NoneType', 'bool'])

    check_parameter('extract', 'qa_show', 
                    qa_show, ['NoneType', 'bool'])

    check_parameter('extract', 'qa_showscale', 
                    qa_showscale, ['int', 'float', 'NoneType'])

    check_parameter('extract', 'qa_showblock', 
                    qa_showblock, ['NoneType', 'bool'])

    qa = check_qakeywords(
        verbose=verbose,
        show=qa_show,
        showscale=qa_showscale,
        showblock=qa_showblock,
        write=qa_write)



    #
    # Let's get set up for the for loop.  Build the input and output file names.
    #
    
    # Determine the path in which `files` is located.  
    
    if load_directory == 'raw':

        load_path = setup.state['raw_path']
        
    elif load_directory == 'proc':

        load_path = setup.state['proc_path']

    # Create the input file names

    results = files_to_fullpath(
        load_path,
        filenames,
        setup.state['nint'],
        setup.state['suffix'],
        input_suffix)

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
                        [output_prefix,filenames[1]],
                        setup.state['nint'],
                        '',
                        '',
                        exist=False)

            output_filenames = result[0]        

        else:

            output_filenames = []
            for file in filenames:

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
        
        load_image(
            input_subset,
            flat_name,
            wavecal_name,
            output_filenames=output_subset,
            output_prefix=output_prefix,
            input_suffix=input_suffix,
            load_directory=load_directory,
            flat_field=flat_field,
            linearity_correction=linearity_correction,
            detector_info=detector_info,
            rectification_method=rectification_method,
            write_rectified_orders=write_rectified_orders,
            do_all_steps=do_all_steps,
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])


        # Make the profile

        make_profiles(
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])
        
        # Identify apertures

        identify_apertures(
            aperture_findinfo,
            seeing_fwhm=seeing_fwhm,
            ybuffer=profile_ybuffer,
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])


        # Select orders

        if include_orders is not None or exclude_orders is not None:

            select_orders(
                include=include_orders,
                exclude=exclude_orders,
                verbose=qa['verbose'],
                qa_show=qa['show'],
                qa_showscale=qa['showscale'],
                qa_showblock=qa['showblock'],
                qa_write=qa['write'])

        # Update aperture positions if requested

        if aperture_signs is not None:

            override_aperturesigns(
                aperture_signs,
                verbose=qa['verbose'])

        # Trace the orders

        trace_apertures(
            fit_degree=trace_fitdegree,
            step_size=trace_stepsize,
            summation_width=trace_summationwidth,
            centroid_threshold=trace_centroidthreshold,
            seeing_fwhm=seeing_fwhm,
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])
        
        # Define the aperture parameters

        define_aperture_parameters(
            aperture_radii,
            bg_annulus=bg_annulus,
            bg_regions=bg_regions,
            bg_fit_degree=bg_fitdegree,
            psf_radius=psf_radius,
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])

        # Extract the apertures

        extract_apertures(
            fix_badpixels=fix_badpixels,
            use_meanprofile=use_meanprofile,
            badpixel_thresh=badpixel_threshold,
            verbose=qa['verbose'],
            qa_show=qa['show'],
            qa_showscale=qa['showscale'],
            qa_showblock=qa['showblock'],
            qa_write=qa['write'])




        

        
    

    
