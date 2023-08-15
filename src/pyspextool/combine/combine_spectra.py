import os
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
from astropy.io import fits

from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.io.check import *
from pyspextool.io.files import *
from pyspextool.io.read_spectra_fits import read_spectra_fits
from pyspextool.io.fitsheader import average_header_info
from pyspextool.io.fitsheader import get_header_info
from pyspextool.plot.limits import get_stack_range
from pyspextool.utils import math
from pyspextool.utils.split_text import split_text
from pyspextool.plot.plot_spectra import plot_spectra



def combine_spectra(files, output_name, input_path=None, output_path=None,
                    scale_spectra=True, scale_order=None, scale_range=None,
                    scale_range_fraction=0.7, correct_spectral_shape=False,
                    statistic='robust weighted mean', robust_sigma=8,
                    qa_show=None, qa_showsize=(10, 6), qa_write=None,
                    line_width=0.5, verbose=None, overwrite=True):

    """
    To combine raw extracted spectra

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

    input_path : str or None
        An optional path with the data.  Otherwise the default is the proc/
        directory.

    output_path : str or None
        An optional output path.  Otherwise the default is the proc/
        directory.

    scale_spectra : {True, False}
        Set to scale the spectra to a common flux level.  See below.

    scale_order : int or None
        The order number to compute the scale factors.  If None, defaults
        to middle-most order.

    scale_range : list or None
        A (2,) list giving the wavelength range, e.g. [1.1,1.15].  If None, 
        defaults to the central `scale_range_fraction` of `scale_order`.

    scale_range_fraction : float, default=0.7
        A float giving the central fraction of the wavelength range to be
        used for scaling.  

    correct_spectral_shape : {False, True}
        Set to True to correct the spectral shape of the spectra.

    statistic : {'robust weighted mean', 'robust mean', 'weighted mean', 
                 'mean', 'median'}
        For any of the means, the uncertainty is the standard error on the 
        mean, e.g. s/np.sqrt(n), where s is the sample standard deviation and
        n is the number of data points.

        For the median, the uncertainty is given by 1.482*MAD/np.sqrt(n) where
        MAD is the median absolute deviation given by,

        1.482*median(|x_i-x_med|).

    robust_sigma : float or int, default=8
        The sigma threshold of a robust statistic.   Values are identified as
        outliers if

        |x_i - x_med|/MAD > robust_sigma,

        where x_i is the ith data point, x_med is the median of the data, and 
        MAD is the median absolute deviation given by,

        1.482*median(|x_i-x_med|).
    
    qa_show : {None, True, False}
        Set to True/False to override config.state['qa_show'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be interactively generated.

    qa_showsize : tuple, default=(10, 6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}
        Set to True/False to override config.state['qa_write'] in the 
        pyspextool config file.  If set to True, quality assurance 
        plots will be written to disk.

    line_width : float or int, default=0.5
        The line width passed to matplotlib for the qa plots

    verbose : {None, True, False}
        Set to True/False to override config.state['verbose'] in the 
        pyspextool config file.  

    overwrite: {True, False}, optional
        Set to False to to not overwrite a file already on disk.

    Returns
    -------
    None
        Writes a pyspextool FITS file to disk.

    """

    #
    # Check the parameters
    #

    check_parameter('combine_spectra', 'files', files, ['str', 'list'])

    check_parameter('combine_spectra', 'output_name', output_name, 'str')    

    check_parameter('combine_spectra', 'input_path', input_path,
                    ['NoneType', 'str'])

    check_parameter('combine_spectra', 'output_path', output_path,
                    ['NoneType', 'str'])    
    
    check_parameter('combine_spectra', 'scale_spectra', scale_spectra, 'bool')

    check_parameter('combine_spectra', 'scale_order', scale_order,
                    ['NoneType', 'int'])

    check_parameter('combine_spectra', 'scale_range', scale_range,
                    ['NoneType', 'list'])

    check_parameter('combine_spectra', 'scale_range_fraction',
                    scale_range_fraction, 'float')                                            
    check_parameter('combine_spectra', 'correct_spectral_shape',
                    correct_spectral_shape, 'bool')

    check_parameter('combine_spectra', 'statistic', statistic, 'str')

    check_parameter('combine_spectra', 'robust_sigma', robust_sigma,
                    ['int', 'float'])            

    check_parameter('combine_spectra', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'qa_showsize', qa_showsize, 'tuple')    

    check_parameter('combine_spectra', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'line_width', line_width,
                    ['float', 'int'])    

    check_parameter('combine_spectra', 'verbose', verbose, ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'overwrite', overwrite, 'bool')    

    #
    # Check the qa and verbose variables and set to system default if need be.
    #
        
    if qa_write is None:

        qa_write = setup.state['qa_write']

    if qa_show is None:

        qa_show = setup.state['qa_show']

    if verbose is None:
        verbose = setup.state['verbose']

    # Get user paths if need be.
        
    if input_path is None:

        input_path = setup.state['proc_path']

    if output_path is None:

        output_path = setup.state['proc_path']        

    check_path(input_path)
    check_path(output_path)    
        
    #
    # Store user inputs
    #

    combine.load['files'] = files
    combine.load['output_name'] = output_name
    combine.load['input_path'] = input_path
    combine.load['output_path'] = output_path    
    combine.load['scale_spectra'] = scale_spectra
    combine.load['scale_order'] = scale_order
    combine.load['scale_range'] = scale_range
    combine.load['scale_range_fraction'] = scale_range_fraction
    combine.load['correct_spectral_shape'] = correct_spectral_shape
    combine.load['statistic'] = statistic.lower()
    combine.load['robust_sigma'] = robust_sigma    
    combine.load['qa_show'] = qa_show
    combine.load['qa_showsize'] = qa_showsize
    combine.load['qa_write'] = qa_write
    combine.load['line_width'] = line_width
    combine.load['verbose'] = verbose
    combine.load['overwrite'] = overwrite

    #
    # Create the file names
    #

    if isinstance(files, str):

        # You are in FILENAME mode

        files = files.replace(" ", "").split(',')
        input_files = make_full_path(input_path, files, exist=True)

    else:

        # You are in INDEX mode

        prefix = files[0]
        nums = files[1]

        # Create the files to read into memory

        indexinfo = {'nint': setup.state['nint'], 'prefix': prefix,
                     'suffix':'', 'extension': '.fits'}

        input_files = make_full_path(input_path, nums, indexinfo=indexinfo,
                                     exist=True)

    combine.state['input_files'] = input_files
    combine.state['nfiles'] = len(input_files)
    
    #
    # Read the files
    #

    if verbose is True:

        print('Combining Spectra')
        print('-----------------')
        print('Loading spectra...')
        
    load_allorders()

    # Do the qa plotting

    if combine.load['qa_show']:

        plot_allorders(figure_size=qa_showsize, title='raw')

    if combine.load['qa_write']:
        
        qafileinfo = {'figsize': (11,8.5),
                      'filepath': setup.state['qa_path'],
                      'filename': combine.load['output_name'],
                      'extension': setup.state['qa_extension']}

        plot_allorders(file_info=qafileinfo, suffix='_raw', title='raw')
            
    #
    # Scale the spectra
    #

    if scale_spectra is True:

        if verbose is True:

            print('Scaling the spectra...')

        scale_allorders(scale_order, scale_range, scale_range_fraction)

    #
    # Do the qa plotting
    #
    
    if combine.load['qa_show']:

        plot_allorders(figure_size=qa_showsize, title='Scaled',
                       plot_scale_range=True)

    if combine.load['qa_write']:
        
        qafileinfo = {'figsize': (11,8.5),
                      'filepath': setup.state['qa_path'],
                      'filename': combine.load['output_name'],
                      'extension': setup.state['qa_extension']}

        plot_allorders(file_info=qafileinfo, suffix='_scaled',
                       plot_scale_range=True, title='Scaled')

    #
    # Combine them 
    #

    if verbose is True:

        print('Combining spectra...')
        
    combine_allorders()

    #
    # Write the results to disk
    #

    write_file()


def combine_allorders():

    """
    To combine the spectra using the statistic of the user's choice.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    
    #
    # Create output arrays
    #
    
    intensity = np.full_like(combine.state['wavelengths'], np.nan)
    uncertainty = np.full_like(combine.state['wavelengths'], np.nan)
    bitmask = np.full_like(combine.state['wavelengths'], 0)

    #
    #  Loop over apertures and orders
    #

    sigma = combine.load['robust_sigma']
    statistic = combine.load['statistic'].replace(" ", "")
    
    for i in range(combine.state['final_napertures']):


        for j in range(combine.state['norders']):

            data = combine.state['intensities'][i,j,:,:]
            var = combine.state['uncertainties'][i,j,:,:]**2

            
            # First do the data

            if statistic == 'robustweightedmean':

                result = math.mean_data_stack(data, robust=sigma,
                                              weights=1/var, stderr=True)
                
            elif statistic == 'robustmean':

                result = math.mean_data_stack(data, robust=sigma, stderr=True)
            
            elif statistic == 'weightedmean':
                
                result = math.mean_data_stack(data,weights=1/var, stderr=True)
        
            elif statistic == 'mean':

                result = math.mean_data_stack(data,stderr=True)
                
            elif statistic == 'median':                        
                
                result = math.median_data_stack(data,stderr=True)
                
            else:

                message=combine.load['statistic']+ ' is an unknown `statitics`.'
                raise ValueError(message)


            mean = result[0]
            unc = result[1]

            # Now do the mask
            
            array = combine.state['bitmasks'][i,j,:,:]
            mask = math.combine_flag_stack(array)
        
            # Store the results
                    
            intensity[i,j,:] = mean
            uncertainty[i,j,:] = unc
            bitmask[i,j,:] = mask

    #
    # Store all the data
    #
            
    combine.state['intensity'] = intensity
    combine.state['uncertainty'] = uncertainty
    combine.state['bitmask'] = bitmask

        
def load_allorders():

    """
    To combine the spectra using the statistic of the user's choice.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

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
            combine.state['combine_type'] = 'standard'
        

    elif combine.state['module'] == 'telluric':

        combine.state['combine_apertures'] = False
        combine.state['combine_type'] = 'telluric'
                
    else:

        message = 'Unknown module.'
        raise ValueError(message)


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
    bitmasks = np.full(shape, 0)
        
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

                info = get_header_info(header,keywords=keywords)
                headers.append(info)

                # store the data

                if combine.state['combine_apertures'] is False:

                    idx = j*combine.state['napertures']+k
                    wavelengths[k,j,:] = spectra[idx,0,:]
                    intensities[k,j,i,:] = spectra[idx,1,:]
                    uncertainties[k,j,i,:] = spectra[idx,2,:]
                    bitmasks[k,j,i,:] = spectra[idx,3,:]

                else:

                    idx = j*combine.state['napertures']+k
                    wavelengths[0,j,:] = spectra[idx,0,:]
                    intensities[0,j,k,:] = spectra[idx,1,:]
                    uncertainties[0,j,k,:] = spectra[idx,2,:]
                    bitmasks[0,j,k,:] = spectra[idx,3,:]

    # Load into memory
                    
    combine.state['headers'] = headers
    combine.state['wavelengths'] = wavelengths
    combine.state['intensities'] = intensities
    combine.state['uncertainties'] = uncertainties
    combine.state['bitmasks'] = bitmasks

    
def plot_allorders(file_info=None, figure_size=None, plot_scale_range=False,
                   suffix=None, title=None):

    """
    To plot all the spectra

    Parameters
    ----------
    None

    Returns
    -------
    None
        

    """

    if file_info is not None:

        figure_size = file_info['figsize']
    
    #
    # Copy for ease of use
    #

    wavelength = combine.state['wavelengths']
    intensity = combine.state['intensities']
    uncertainty = combine.state['uncertainties']
    bitmask = combine.state['bitmasks']
    
    #
    # Get the ranges
    #

    wavelength_range = [np.nanmin(wavelength),np.nanmax(wavelength)]

    wavelength_ranges = np.empty((combine.state['norders'], 2))    
    intensity_ranges = np.empty((combine.state['norders'], 2))
    uncertainty_ranges = np.empty((combine.state['norders'], 2))        

    # Do the loop
    
    for i in range(combine.state['norders']):

        wavelength_ranges[i,:] = np.array([np.nanmin(wavelength[0,i]),
                                           np.nanmax(wavelength[0,i])])
        
        intensity_ranges[i,:] = get_stack_range(intensity[0,i], savgol=True,
                                                frac=0.05)

        uncertainty_ranges[i,:] = get_stack_range(uncertainty[0,i], savgol=True,
                                                frac=0.05)        

    # Find the minimum and maximum values
        
    array = [intensity_ranges, uncertainty_ranges]
    intensity_range = [np.min(array), np.max(array)]

    # Create the figure
    
    figure = pl.figure(figsize=figure_size)
    ax = figure.add_axes([0.125, 0.11, 0.775, 0.77])
    ax.plot([np.nan], [np.nan])
    ax.set_xlim(wavelength_range)
    ax.set_ylim(intensity_range)
    ax.set_xlabel(combine.state['xlabel'])
#    ax.set_ylabel(ylabel)

    axis_to_data = ax.transAxes + ax.transData.inverted()
    data_to_axis = axis_to_data.inverted()
    

    for i in range(combine.state['norders']):

        color = iter(cm.rainbow(np.linspace(0, 1, combine.state['nspectra'])))
        for j in range(combine.state['nspectra']):

            c = next(color)
            pl.plot(wavelength[0,i,:], intensity[0,i,j,:],c=c,
                    lw=combine.load['line_width'])                              
            pl.plot(wavelength[0,i,:], uncertainty[0,i,j,:],c='gray',
                    lw=combine.load['line_width'])            

        # Get the order number position

        ave_wave = np.mean(wavelength_ranges[i, :])
        xnorm = data_to_axis.transform((ave_wave, 0))[0]

        # Print the text

        ax.text(xnorm, 1.01 + 0.01 * (i % 2), str(combine.state['orders'][i]),
                color='black', transform=ax.transAxes)

    # Add the title

    if title is not None:

        ax.set_title(title, pad=20.0)
    
    if plot_scale_range is True:

        ax.axvline(x=combine.state['scale_range'][0],color='0',ls='--')
        ax.axvline(x=combine.state['scale_range'][1],color='0',ls='--')        

    
    if file_info is not None:


        qafileinfo = {'figsize': (11,8.5),
                      'filepath': setup.state['qa_path'],
                      'filename': combine.load['output_name'],
                      'extension': setup.state['qa_extension']}

        
        
        pl.savefig(os.path.join(file_info['filepath'],
                                file_info['filename']+suffix+\
                                file_info['extension']))
        pl.close()
                                
    else:

        pl.show()
        pl.close()
                    

def scale_allorders(scale_order, scale_range, scale_range_fraction):

    """
    To scale the spectra and uncertainties to a common level

    Parameters
    ----------

    Returns
    -------
    None

    """

    #
    # First determine which order we are using to scale
    #
    
    if scale_order is None:

        # Do it ourselves.  Pick the middle-most order.

        scale_order = int(np.median(combine.state['orders']))

    else:

        # Let the user choose, but make sure it is an order we have.
        
        z = combine.state['orders'] == scale_order
        if np.sum(z) == 0:

            message = '`scale_order`='+str(scale_order)+' not in the files.'
            raise ValueError(message)

    # Store the results
        
    combine.state['scale_order'] = scale_order

    #
    # Now let's figure out the wavelength range over which to do the scaling.
    
    if scale_range is None:

        # Determine the wavelenth range ourselves.  

        z_order = combine.state['orders'] == combine.state['scale_order']    

        # Grab the zeroth aperture wavelength array.
        
        wave = np.squeeze(combine.state['wavelengths'][0, z_order, :])

        min_wave = np.nanmin(wave)
        max_wave = np.nanmax(wave)        
        median_wave = np.nanmedian(wave)

        delta_wave = (max_wave-min_wave)*scale_range_fraction/2

        scale_range = np.array([median_wave-delta_wave,median_wave+delta_wave])

        combine.state['scale_range'] = scale_range

    else:

        # Let the user choose, but make sure it falls in range.

        combine.state['scale_range'] = combine.load['scale_range']
        
    #                    
    # Determine which order we are using to determine the scale factors
    #
        
    z_order = np.where(combine.state['orders'] == combine.state['scale_order'])

    # 
    # Loop over each aperture and order
    #

    scales = np.empty((combine.state['final_napertures'],
                       combine.state['nspectra']))
    
    for i in range(combine.state['final_napertures']):

        # Determine which pixels we are using for the scaling

        wave = np.squeeze(combine.state['wavelengths'][i, z_order, :])
        intensity = np.squeeze(combine.state['intensities'][i,z_order,:,:])

        z_wave = np.where((wave > combine.state['scale_range'][0]) &
                          (wave < combine.state['scale_range'][1]))
            
        # Get scale factors
                
        junk, junk, scale = math.scale_data_stack(intensity[:, z_wave[0]],
                                                   None)
        scales[i,:] = scale
        
        # Now scale each order

        shape = np.shape(intensity)                
        reshape_info = (shape[0], 1)
        tile_info = (1, shape[1])        

        scale_array = np.tile(np.reshape(scale, reshape_info), tile_info)
        #        scale_array = np.absolute(scale_array)
        
        for j in range(combine.state['norders']):

            np.multiply(combine.state['intensities'][i,j,:,:], scale_array,
                        out=combine.state['intensities'][i,j,:,:])

            np.multiply(combine.state['uncertainties'][i,j,:,:],
                        np.sqrt(scale_array),
                        out=combine.state['uncertainties'][i,j,:,:])

    combine.state['scales'] = scales
    

def write_file():

    """
    To write the results to disk.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    #
    # Create an array into which you can store the results
    #
    
    shape = (combine.state['norders']*combine.state['final_napertures'],4, \
             combine.state['npixels'])
    
    array = np.empty(shape)

    #
    # File the array
    #
    
    for i in range(combine.state['norders']):

        for j in range(combine.state['final_napertures']):        

            idx = i+combine.state['final_napertures']*j
            array[idx,0,:]  = combine.state['wavelengths'][j,i]
            array[idx,1,:]  = combine.state['intensity'][j,i]
            array[idx,2,:]  = combine.state['uncertainty'][j,i]
            array[idx,3,:]  = combine.state['bitmask'][j,i]

    #
    # Average the headers together
    #

#    avehdr = average_header_info(combine.state['headers'])

#    avehdr['MODULE'][0] = 'combine'
#    avehdr['FILENAME'][0] = combine.load['output_name']+'.fits'

    
    #
    # Determine useful things depending on combination mode
    #
    
    
    if combine.state['combine_type'] == 'standard':

        napertures_combined = 0
        avehdr = average_header_info(combine.state['headers'])

    elif combine.state['combine_type'] == 'twoaperture':

        napertures_combined = 2
        avehdr = combine.state['headers'][0]
        avehdr['TOTITIME'] = [2*avehdr['TOTITIME'][0],
                              ' Total integration time (sec)']

    elif combine.state['combine_type'] == 'telluric':        

        napertures_combined = 0

    # Store the history

    old_history = avehdr['HISTORY']

    # remove it from the avehdr

    avehdr.pop('HISTORY')

    #
    # Add things to the average header
    #

    avehdr['MODULE'][0] = 'combine'

    avehdr['FILENAME'][0] = combine.load['output_name']+'.fits'


    avehdr['NAPS'][0] = combine.state['final_napertures']

    avehdr['NFLCOMB'] = [combine.state['nfiles'],
                         ' Number of spectra files combined']

    avehdr['NAPCOMB'] = [napertures_combined,
                         ' Number of apertures combined']

    avehdr['COMBSTAT'] = [combine.load['statistic'],
                         ' Combination statistic']

    avehdr['RBTHRESH'] = [combine.load['robust_sigma'],
                         ' Robust threshold (if used)']            

        
    #
    # Create the header
    #
    
    # Create the basic headers

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add our keywords
    
    keys = list(avehdr.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (avehdr[keys[i]][0], avehdr[keys[i]][1])


    # Strip the directories.

    files = [os.path.basename(i) for i in combine.state['input_files']]
                        
    # Deal with the general history

    history = 'The spectra in the file(s) '+', '.join(files)+\
              ' were combined using a '+combine.load['statistic']

    if combine.load['statistic'][0:6] == 'robust':

        history += ' with a sigma threshold of '+\
                   str(combine.load['robust_sigma'])+'.  '

    else:

        history += '.  '

    # Now do the modifications for each aperture

    for i in range(combine.state['final_napertures']):

        history += 'Aperture '+str(i+1)+' modifications: '
        
        if combine.load['scale_spectra'] == True:

            scales = [str(i) for i in combine.state['scales'][i,:]]
            
            history += 'The scale factors were determined using order '+\
              str(combine.state['scale_order'])+' and are: '+\
              ', '.join(scales)+'.  '
        
    # Now split, add the old history, and add to the hdr

    history = split_text(history, length=65)            

    history = old_history + [' '] + history

    for hist in history:

        hdr['HISTORY'] = hist

    #
    # Write the file out
    #

    full_path = os.path.join(combine.load['output_path'],
                             combine.load['output_name']+'.fits')
    
    fits.writeto(full_path, array, hdr, overwrite=combine.load['overwrite'])

    if combine.load['verbose'] is True:
        print('Wrote', os.path.basename(full_path) + ' to disk.')    

            
    #
    # Do the plotting
    #

    if combine.load['qa_show'] is True:

        plot_spectra(full_path, ytype='flux and uncertainty',
                     line_width=combine.load['line_width'],
                     title=os.path.basename(full_path))

        
    if combine.load['qa_write'] is True:

        qafileinfo = {'figsize': combine.load['qa_showsize'],
                      'filepath': setup.state['qa_path'],
                      'filename': combine.load['output_name'],
                      'extension': setup.state['qa_extension']}

        plot_spectra(full_path, ytype='flux and uncertainty',
                     line_width=combine.load['line_width'],
                     title=os.path.basename(full_path), file_info=qafileinfo)            
        if combine.load['verbose'] is True:
            print('Wrote', combine.load['output_name']+\
                  setup.state['qa_extension']+' to disk.')    
