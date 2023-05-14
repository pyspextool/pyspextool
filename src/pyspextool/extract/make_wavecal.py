import importlib
import os
import numpy as np

from pyspextool import config as setup

from pyspextool.io.check import check_parameter
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.files import make_full_path
from pyspextool.io.wavecal import read_wavecal_file
from pyspextool.io.wavecal import read_line_list

from pyspextool.extract.find_lines_1dxd import find_lines_1dxd
from pyspextool.extract.wavecal_solution_1d import wavecal_solution_1d
from pyspextool.extract.get_line_guess_position import get_line_guess_position
from pyspextool.extract.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.extract.make_interp_indices_1d import make_interp_indices_1d
from pyspextool.extract.get_spectral_pixelshift import get_spectral_pixelshift
from pyspextool.utils.math import scale_data_stack
from pyspextool.utils.math import median_data_stack

from pyspextool.io.wavecal import write_wavecal_1d


def make_wavecal(files, flat_file, output_name, extension='.fits*',
                 use_stored_solution=False, verbose=None, qa_file=None,
                 qa_plot=None, qa_plotsize=(8, 10), overwrite=True):
    """
    To create a spextool wavecal file.

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a str giving the perfix, files[1] is a str giving the 
        index numbers of the files, e.g. ['spc', '1-2,5-10,13,14'].
        
    flat_file : str
        The full path of a spextool flat field file.

    output_name : str
        The filename of the wavecal file to written to disk.      

    extension : str, default='.fits*', optional
        The file extension

    use_stored_solution : {False, True}
        Set to True to use  the solution stored in the pyspextool file.

    verbose : {None, True, False}, optional
        Set to True/False to override config.setup['verbose']

    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.setup['qa_file'].  If set to True, 
        quality assurance plots will be written to disk.

    qa_plot : {None, True, False}, optional
        Set to True/False to override config.setup['qa_plot'].  If set to True,
        quality assurance plots will be interactively generated.

    qa_plotsize : tuple, default=(8,10)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    Returns
    -------
     None
        Writes a FITS file to disk.
        
    """

    #
    # Check parameters
    #

    check_parameter('make_flat', 'files', files, ['str', 'list'])

    if isinstance(files, list):
        check_parameter('make_flat', 'files[0]', files[0], 'str')
        check_parameter('make_flat', 'files[1]', files[1],
                        ['str', 'list', 'int'])

    check_parameter('make_wavecal', 'flat_file', flat_file, 'str')

    check_parameter('make_wavecal', 'output_name', output_name, 'str')

    check_parameter('make_wavecal', 'use_stored_solution',
                    use_stored_solution, 'bool')

    check_parameter('make_wavecal', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('make_wavecal', 'overwrite', overwrite, 'bool')

    check_parameter('make_wavecal', 'qa_file', qa_file, ['NoneType', 'bool'])

    check_parameter('make_wavecal', 'qa_plot', qa_plot, ['NoneType', 'bool'])

    check_parameter('make_wavecal', 'qa_plotsize', qa_plotsize,
                    ['NoneType', 'tuple'])

    #
    # Check the qa variables and set to system default if need be.
    #

    if qa_file is None:
        qa_file = setup.state['qa_file']

    if qa_plot is None:
        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

    #
    # Let the user know what you are doing.
    #
    
    if verbose is True:
        print('Generating Wavelength Solution')
        print('------------------------------')        
        
    #
    # Load the instrument module for the read_fits program
    #

    module = 'pyspextool.instrument_data.' + setup.state['instrument'] + \
             '_dir.' + setup.state['instrument']

    instr = importlib.import_module(module)

    #
    # Read the flat and wavecal files
    #

    flat_file = os.path.join(setup.state['cal_path'], flat_file)
    flatinfo = read_flat_fits(flat_file)

    wavecalfile = os.path.join(setup.state['instrument_path'],
                               flatinfo['mode'] + '_wavecalinfo.fits')

    wavecalinfo = read_wavecal_file(wavecalfile)

    #
    # At some point, add a check to ensure the arcs have the same mode
    #

    #
    # Load the files into memory
    #

    # Create the file names

    if isinstance(files, str):

        # You are in FILENAME mode

        files = files.replace(" ", "").split(',')
        files = make_full_path(setup.state['raw_path'], files, exist=True)

    else:

        # You are in INDEX mode

        prefix = files[0]
        nums = files[1]

        files = make_full_path(setup.state['raw_path'], nums,
                               indexinfo={'nint': setup.state['nint'],
                                          'prefix': prefix,
                                          'suffix': setup.state['suffix'],
                                          'extension': extension},
                               exist=True)

    # Get the mode name in order to get the array rotation value.

    img, var, hdr, mask = instr.read_fits(files,
                                          setup.state['linearity_info'],
                                    keywords=setup.state['extract_keywords'],
                                          rotate=flatinfo['rotation'],
                                          verbose=verbose)

    #
    # Combine images as necessary and don't worry about the variance
    #

    # Now scale their intensities to a common flux level

    if len(files) > 1:

        if verbose is True:
            print('Scaling images...')

        simgs, svars, scales = scale_data_stack(img, None)

    else:
        simgs = img

    # Now median the scaled images

    if len(files) > 1:

        if verbose is True:
            print('Medianing the images...')

        med, munc = median_data_stack(simgs)

    else:
        med = simgs

    #
    # Let's do the extraction of the "arc" spectra
    #

    # Create wavecal and spatcal images

    wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                             flatinfo['nrows'],
                                             flatinfo['edgecoeffs'],
                                             flatinfo['xranges'],
                                             flatinfo['slith_arc'])

    # Extract the "arc"                                

    appos = np.full((np.size(flatinfo['orders']), 1), flatinfo['slith_arc'] / 2)

    spectra = extract_extendedsource_1dxd(med, med, flatinfo['ordermask'],
                                          flatinfo['orders'], wavecal,
                                          spatcal, appos,
                                          wavecalinfo['apradius'],
                                          linmax_bitmask=None,
                                          badpixel_mask=None, bginfo=None,
                                          verbose=verbose)

    #
    # Find the pixel offset between these spectra and the disk spectra
    #

    # Get the anchor order and spectra.

    z = np.sum(np.where(flatinfo['orders'] == wavecalinfo['xcororder']))

    xanchor = np.arange(int(wavecalinfo['xranges'][z, 0]),
                        int(wavecalinfo['xranges'][z, 1] + 1), dtype=int)
    fanchor = np.squeeze(wavecalinfo['spectra'][z, 1, :])

    # Get the source order

    xsource = np.squeeze(spectra['spectra'][z][0, :])
    fsource = np.squeeze(spectra['spectra'][z][1, :])

    if qa_file is True:

        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': setup.state['qa_path'],
                      'filename': output_name, 'extension': '.pdf'}
    else:
        qafileinfo = None

    offset = get_spectral_pixelshift(xanchor, fanchor, xsource, fsource,
                                     qafileinfo=qafileinfo, qa_plot=qa_plot,
                                     qa_plotsize=qa_plotsize)

    #
    # Are we using the stored solution?
    #

    if use_stored_solution is False:

        #
        # Locate the line positions
        #

        # Get the line list to search for lines

        filename = os.path.join(setup.state['instrument_path'],
                                wavecalinfo['linelist'])

        lineinfo = read_line_list(filename, delta_to_microns=True)

        #  Determine the guess position and search range for each

        lineinfo = get_line_guess_position(wavecalinfo['spectra'],
                                           wavecalinfo['orders'],
                                           flatinfo['xranges'], lineinfo)

        # Add the shift offset to the results

        lineinfo['xguess'] = lineinfo['xguess'] + offset
        lineinfo['range_min_xguess'] = lineinfo['range_min_xguess'] + offset
        lineinfo['range_max_xguess'] = lineinfo['range_max_xguess'] + offset

        # Now find the lines

        if verbose:
            print('Finding the lines...')

        if qa_file is True:

            qafileinfo = {'figsize': (8.5, 11),
                          'filepath': setup.state['qa_path'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

            # Find the lines

        lineinfo = find_lines_1dxd(spectra['spectra'], wavecalinfo['orders'],
                                   lineinfo, flatinfo['slitw_pix'],
                                   qafileinfo=qafileinfo, verbose=verbose)

        #
        # Let's do the actual calibration
        #

        if verbose:
            print('Determining the wavelength solution...')

        # Get set up for either 1d of 1dxd

        if wavecalinfo['wcaltype'] == '1d':

            qafile_figsize = (6, 4)
            xdinfo = None

        else:

            qafile_figsize = (8.5, 11)
            xdinfo = {'homeorder': wavecalinfo['homeorder'],
                      'orderdeg': wavecalinfo['ordrdeg']}

        if qa_file is True:

            qafileinfo = {'figsize': qafile_figsize,
                          'filepath': setup.state['qa_path'],
                          'filename': output_name, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the solution

        solution = wavecal_solution_1d(wavecalinfo['orders'], lineinfo,
                                       wavecalinfo['dispdeg'], xdinfo=xdinfo,
                                       qa_fileinfo=qafileinfo,
                                       qa_plot=qa_plot,
                                       qa_plotsize=qa_plotsize,
                                       verbose=verbose)

    else:

        if verbose:
            print('Using stored solution...')

        solution = {'coeffs': wavecalinfo['coeffs'],
                    'covar': wavecalinfo['covar'],
                    'rms': wavecalinfo['rms'],
                    'nlines': wavecalinfo['nlines'],
                    'ngood': wavecalinfo['ngood'],
                    'nbad': wavecalinfo['nbad']}

    #
    # Creating rectification indices
    #

    indices = []
    for i in range(flatinfo['norders']):
        idxs = make_interp_indices_1d(flatinfo['edgecoeffs'][i, :, :],
                                      flatinfo['xranges'][i, :],
                                      flatinfo['slith_arc'],
                                      array_output=True)

        indices.append(idxs)

        #
    # Write the wavecal file to disk.
    #

    if verbose:
        print('Writing wavecal to disk...')

    if wavecalinfo['wcaltype'] == '1d':

        xdinfo = None

    else:

        xdinfo = {'orderdeg': wavecalinfo['ordrdeg'],
                  'homeorder': wavecalinfo['homeorder']}

    write_wavecal_1d(flatinfo['ncols'], flatinfo['nrows'],
                     flatinfo['orders'], flatinfo['edgecoeffs'],
                     flatinfo['xranges'], solution['coeffs'],
                     solution['covar'], wavecalinfo['dispdeg'],
                     solution['rms'] * 1e4, solution['nlines'],
                     solution['ngood'], solution['nbad'],
                     wavecal, spatcal, indices, flatinfo['rotation'],
                     flat_file, os.path.join(setup.state['cal_path'],
                                             output_name + '.fits'),
                     setup.state['version'],
                     xdinfo=xdinfo,
                     stored_solution=use_stored_solution,
                     overwrite=overwrite)

    if verbose:
        print('Wavecal ' + output_name + '.fits written to disk.\n')
