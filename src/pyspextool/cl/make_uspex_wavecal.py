import os
import numpy as np
import pickle

from pyspextool.calibration.simulate_wavecal_1dxd import simulate_wavecal_1dxd
from pyspextool.calibration.get_line_guess_position import get_line_guess_position
from pyspextool.calibration.find_lines_1dxd import find_lines_1dxd
from pyspextool.calibration.wavecal_solution_1dxd import wavecal_solution_1dxd
from pyspextool.cl import config
from pyspextool.io.files import make_full_path
from pyspextool.io.read_uspex_fits import read_uspex_fits
from pyspextool.io.flat import read_flat_fits
from pyspextool.io.wavecal import read_line_list
from pyspextool.io.wavecal import read_wavecal_file
from pyspextool.io.wavecal import write_wavecal_1dxd
from pyspextool.spectroscopy.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.spectroscopy.get_spectral_pixelshift import get_spectral_pixelshift
from pyspextool.utils.math import median_data_stack
from pyspextool.utils.math import scale_data_stack


def make_uspex_wavecal(files, flatfile, oname, prefix='arc-',
                       suffix='.[ab]', extension='.fits*',
                       input_method='index', use_stored_solution=False,
                       clupdate=True, qafile_shift=True,
                       qafile_findlines=True, qafile_fitlines=True,
                       overwrite=True):
    #
    # Load the files into memory
    #

    # Get the user requested keywords

    keywords = config.state['xspextool_keywords']

    # Add GRAT and DIT so that you can determine the mode.  Users will just have
    # to live with them in their files

    if 'GRAT' not in keywords:
        keywords.append('GRAT')

    if 'DIT' not in keywords:
        keywords.append('DIT')

    # Create the file names

    if input_method == 'index':

        files = make_full_path(config.state['rawpath'], files,
                               indexinfo={'nint': config.state['nint'],
                                          'prefix': prefix,
                                          'suffix': suffix,
                                          'extension': extension},
                               exist=True)

    elif input_method == 'filename':

        files = make_full_path(config.state['rawpath'], files, exist=True)

    else:

        raise ValueError('Unknown input_method.')

    # Load the FITS files into memory

    if clupdate is True:
        print('Loading FITS images...')

    img, var, hdr, mask = read_uspex_fits(files,
                                          config.state['linearity_info'],
                                          keywords=keywords,
                                          clupdate=clupdate)

    #
    # Combine images as necessary
    #

    # Now scale their intensities to a common flux level

    if len(files) > 1:

        if clupdate is True:
            print('Scaling images...')

        simgs, svars, scales = scale_data_stack(img, None)

    else:
        simgs = img

    # Now median the scaled images

    if len(files) > 1:

        if clupdate is True:
            print('Medianing the images...')

        med, munc = median_data_stack(simgs)

    else:
        med = simgs

    #
    # Let's do the extraction of the "arc" spectra
    #

    # Read the flat and wavecal files

    flatfile = os.path.join(config.state['calpath'], flatfile)
    flatinfo = read_flat_fits(flatfile)

    wavecalfile = os.path.join(config.state['packagepath'], 'instruments',
                               config.state['instrument'], 'data',
                               flatinfo['mode'] + '_wavecalinfo.fits')

    wavecalinfo = read_wavecal_file(wavecalfile)

    # Create wavecal and spatcal images

    wavecal, spatcal = simulate_wavecal_1dxd(flatinfo['ncols'],
                                             flatinfo['nrows'],
                                             flatinfo['edgecoeffs'],
                                             flatinfo['xranges'],
                                             flatinfo['slith_arc'])

    # Extract the "arc"                                

    spectra = extract_extendedsource_1dxd(med, var, flatinfo['ordermask'],
                                          flatinfo['orders'], wavecal,
                                          spatcal, flatinfo['slith_arc'] / 2,
                                          wavecalinfo['apradius'],
                                          linmax_bitmask=None,
                                          badpixel_mask=None, bginfo=None,
                                          clupdate=clupdate)

    #
    # Find the pixel offset between these spectra and the disk spectra
    #

    # Get the anchor order and spectra.

    z = flatinfo['orders'] == wavecalinfo['xcororder']

    xanchor = np.arange(int(wavecalinfo['xranges'][z, 0]),
                        int(wavecalinfo['xranges'][z, 1] + 1), dtype=int)
    fanchor = np.squeeze(wavecalinfo['spectra'][z, 1, :])

    # Get the source order

    key = 'OR' + str(int(flatinfo['orders'][z])).zfill(3) + '_AP01'
    xsource = np.squeeze(spectra[key][0, :])
    fsource = np.squeeze(spectra[key][1, :])

    if qafile_shift is True:

        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
                      'filename': oname, 'extension': '.pdf'}
    else:
        qafileinfo = None

    offset = get_spectral_pixelshift(xanchor, fanchor, xsource, fsource,
                                     qafileinfo=qafileinfo)

    #
    # Are we using the stored solution?
    #

    if use_stored_solution is False:

       #
       # Locate the line positions
       #

        # Get the line list to search for lines

        filename = os.path.join(config.state['packagepath'], 'instruments',
                                config.state['instrument'], 'data',
                                wavecalinfo['linelist'])

        lineinfo = read_line_list(filename, delta_to_microns=True)

#        with open('data.sav', 'wb') as f:
#                pickle.dump([spectra, wavecalinfo, lineinfo, flatinfo, offset, wavecal, spatcal], f)
        
#        return
        
#        with open('data.sav', 'rb') as f:
#          spectra, wavecalinfo, lineinfo, flatinfo, offset, wavecal, spatcal,   = pickle.load(f)

        #     Determine the guess position and search range for each

        lineinfo = get_line_guess_position(wavecalinfo['spectra'],
                                           wavecalinfo['orders'],
                                           flatinfo['xranges'], lineinfo)

        # Add the shift offset to the results

        lineinfo['xguess'] = lineinfo['xguess'] + offset
        lineinfo['range_min_xguess'] = lineinfo['range_min_xguess'] + offset
        lineinfo['range_max_xguess'] = lineinfo['range_max_xguess'] + offset

        # Now find the lines

        if clupdate:
            print('Finding the lines...')

        if qafile_findlines is True:

            qafileinfo = {'figsize': (8.5, 11),
                          'filepath': config.state['qapath'],
                          'filename': oname, 'extension': '.pdf'}

        else:
            qafileinfo = None

        # Find the lines

        lineinfo = find_lines_1dxd(spectra, wavecalinfo['orders'], lineinfo,
                                   flatinfo['slitw_pix'], qafileinfo=qafileinfo,
                                   clupdate=clupdate)
        return
        #
        # Let's do the actual calibration
        #

        if clupdate:
            print('Determining the wavelength solution...')

        if qafile_fitlines is True:

            qafileinfo = {'figsize': (8.5, 11),
                          'filepath': config.state['qapath'],
                          'filename': oname, 'extension': '.pdf'}

        else:
            qafileinfo = None

        if wavecalinfo['wcaltype'] == '1D':
            print('hi')

        if wavecalinfo['wcaltype'] == '1dxd':

            solution = wavecal_solution_1dxd(wavecalinfo['orders'],lineinfo,
                                             wavecalinfo['homeorder'],
                                             wavecalinfo['dispdeg'],
                                             wavecalinfo['ordrdeg'],
                                             qafileinfo=qafileinfo,
                                             clupdate=clupdate)

            

            
    else:

#        with open('data.sav', 'rb') as f:
#           spectra, wavecalinfo, lineinfo, flatinfo, offset, wavecal, spatcal,   = pickle.load(f)

        if clupdate:
            print('Using stored solution...')
        
        solution = {'coeffs':wavecalinfo['coeffs'],
                    'covar':wavecalinfo['covar'],
                    'rms':wavecalinfo['rms'],
                    'nlines':wavecalinfo['nlines'],
                    'ngood':wavecalinfo['ngood'],
                    'nbad':wavecalinfo['nbad']}


    #
    # Write the wavecal file to disk.
    #

    if clupdate:
        print('Writing wavecal to disk...')

    if wavecalinfo['wcaltype'] == '1D':
        print('later')

    elif wavecalinfo['wcaltype'] == '1dxd':

        write_wavecal_1dxd(flatinfo['ncols'], flatinfo['nrows'],
                           flatinfo['orders'], flatinfo['edgecoeffs'],
                           flatinfo['xranges'], solution['coeffs'],
                           solution['covar'], wavecalinfo['dispdeg'],
                           wavecalinfo['ordrdeg'], wavecalinfo['homeorder'],
                           solution['rms'] * 1e4, solution['nlines'],
                           solution['ngood'], solution['nbad'],
                           wavecal, spatcal, flatinfo['rotation'], flatfile,
                           os.path.join(config.state['calpath'], oname + '.fits'),
                           config.state['version'],
                           stored_solution=use_stored_solution,
                           overwrite=overwrite)

    else:
        print('unknown wcaltype.')
