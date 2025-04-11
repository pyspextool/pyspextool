import numpy as np
from astropy.io import fits
from os.path import basename as osbasename, join as osjoin
import logging

from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.utils import math
from pyspextool.io.fitsheader import average_headerinfo
from pyspextool.utils.split_text import split_text
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.plot.plot_spectra import plot_spectra


def combine_spectra(statistic:str='robust weighted mean',
                    robust_sigma:int | float=8,
                    verbose:bool=None,
                    qa_show:bool=None,
                    qa_showscale:float | int=None,
                    qa_showblock:bool=None,
                    qa_write:bool=None):


    """
    To combine the spectra and write the results to disk

    Parameters
    ----------
    statistic : {'robust weighted mean', 'robust mean', 'weighted mean', 
                 'mean', 'median'}
        For the means, the uncertainty is the standard error on the 
        mean, e.g. s/np.sqrt(n) where s is the sample standard deviation and 
        n is the number of data points.

        For the median, the uncertainty is given by 1.482*MAD/np.sqrt(n) where
        MAD is the median absolute deviation and is given by,

        1.482*median(|x_i-x_med|).

    robust_sigma : float or int, default 8
        The sigma threshold of a robust statistic is requested.   Values are 
        identified as outliers if 

        |x_i - x_med|/MAD > robust_sigma,

        where x_i is the ith data point, x_med is the median of the data, and 
        MAD is the median absolute deviation given by,

        1.482*median(|x_i-x_med|).
        
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
        the plot window.  Set to None to default to setup.state['qa_scale'].    
        
    Returns
    -------
    None
        Write a pySpextool FITS file to the proc/ directory.

    """

    #
    # Check to make sure we can proceed.
    #

    if combine.state['scale_done'] is False:

        message = 'Previous steps complete.  Please run combine.scale_spectra.'
        raise pySpextoolError(message)
    
    #
    # Check the parameters and QA keywords
    #

    check_parameter('combine_spectra', 'statistic', statistic, 'str',
                    possible_values=combine.state['statistics'])

    check_parameter('combine_spectra', 'robust_sigma', robust_sigma,
                    ['int', 'float'])            

    check_parameter('combine_spectra', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('combine_spectra', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('combine_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)
    
    #
    # Create output arrays
    #

    intensity = np.full_like(combine.state['wavelengths'], np.nan)
    uncertainty = np.full_like(combine.state['wavelengths'], np.nan)
    bitmask = np.full_like(combine.state['wavelengths'], 0,dtype=np.int8)


    logging.info(' Combining the spectra using a '+statistic+'.')
    
    #
    #  Loop over apertures and orders
    #

    for i in range(combine.state['final_napertures']):

        for j in range(combine.state['norders']):

            data = combine.state['intensities'][i,j,:,:]
            var = combine.state['uncertainties'][i,j,:,:]**2

            # First do the data

            if statistic == 'robust weighted mean':

                result = math.mean_data_stack(data, robust=robust_sigma,
                                              weights=1/var, stderr=True)

            elif statistic == 'robust mean':

                result = math.mean_data_stack(data, robust=sigma, stderr=True)

            elif statistic == 'weighted mean':

                result = math.mean_data_stack(data,weights=1/var, stderr=True)

            elif statistic == 'mean':

                result = math.mean_data_stack(data,stderr=True)

            elif statistic == 'median':                        

                result = math.median_data_stack(data,stderr=True)

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
    # Write the results to disk
    #

    # Create an array into which you can store the results

    
    shape = (combine.state['norders']*combine.state['final_napertures'],4, \
             combine.state['npixels'])
    array = np.empty(shape)

    #
    # File the array
    #
    
    for i in range(combine.state['norders']):

        for j in range(combine.state['final_napertures']):        

            idx = i*combine.state['final_napertures']+j            
            array[idx,0,:]  = combine.state['wavelengths'][j,i]
            array[idx,1,:]  = intensity[j,i]
            array[idx,2,:]  = uncertainty[j,i]
            array[idx,3,:]  = bitmask[j,i]

    #
    # Determine useful things depending on combination mode
    #
        
    if combine.state['combine_type'] == 'standard':

        napertures_combined = 0
        avehdr = average_headerinfo(combine.state['headers'])

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

    avehdr['FILENAME'][0] = combine.state['output_name']+'.fits'


    avehdr['NAPS'][0] = combine.state['final_napertures']

    avehdr['CMB_NFL'] = [combine.state['nfiles'],
                          'Number of spectra files combined']

    avehdr['CMB_NAP'] = [napertures_combined,
                          'Number of apertures combined']

    avehdr['CMB_STAT'] = [statistic,
                          'Combination statistic']

    avehdr['CMB_THSH'] = [robust_sigma,
                         'Combination robust threshold (if used)']   

    # Add the S/N ratio 

    for i in range(combine.state['norders']):

        name = 'SNRO' + str(combine.state['orders'][i]).zfill(3)
        comment = ' Median S/N values for order ' + \
            str(combine.state['orders'][i]).zfill(3)

        values = []
        for j in range(combine.state['final_napertures']):        

            idx = i*combine.state['final_napertures']+j
        
            signal = array[idx,1,:]
            noise = array[idx,2,:]
            
            values.append(str(int(np.round(np.nanmedian(signal/noise)))))

        avehdr[name] = [", ".join(values), comment]

                
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

    files = [osbasename(i) for i in combine.state['input_files']]
                        
    # Deal with the general history

    history = 'The spectra in the file(s) '+', '.join(files)+\
              ' were combined using a '+statistic

    if statistic[0:6] == 'robust':

        history += ' with a sigma threshold of '+\
                   str(robust_sigma)+'.  '

    else:

        history += '.  '

    # Now do the modifications for each aperture

    for i in range(combine.state['final_napertures']):

        history += 'Aperture '+str(i+1)+' modifications: '
        
        if combine.state['spectra_scaled'] == True:

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

    full_path = osjoin(setup.state['proc_path'],
                       combine.state['output_name']+'.fits')

    fits.writeto(full_path, array, hdr, overwrite=True)

    logging.info(' Wrote file '+osbasename(full_path) + ' to proc/.')

    #
    # Do the QA plotting
    #
    
    if qa['show'] is True:

        plot_spectra(full_path,
                     figure_size=setup.plots['landscape_size'],
                     font_size=setup.plots['font_size'], 
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],
                     showblock=qa['showblock'],
                     showscale=qa['showscale'])

    
    if qa['write'] is True:

        qafullpath = osjoin(setup.state['qa_path'],
                            combine.state['output_name']+\
                            setup.state['qa_extension'])
        
        plot_spectra(full_path,
                     figure_size=setup.plots['landscape_size'],
                     font_size=setup.plots['font_size'],
                     spectrum_linewidth=setup.plots['spectrum_linewidth'],
                     spine_linewidth=setup.plots['spine_linewidth'],
                     output_fullpath=qafullpath)
    
            
            
        
            
    


            
