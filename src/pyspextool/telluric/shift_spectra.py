import numpy as np
import logging
from os.path import join as osjoin
import matplotlib.pyplot as pl

from pyspextool import config as setup
from pyspextool.telluric import config as tc
from pyspextool.io.check import check_parameter, check_qakeywords, check_file
from pyspextool.telluric.core import find_shift, plot_shifts
from pyspextool.utils.interpolate import linear_interp1d


def shift_spectra(default_shiftranges:bool=True,
                  user_shiftranges:tuple | list=None,
                  verbose:bool=None,
                  qa_show:bool=None,
                  qa_showscale:float=None,
                  qa_showblock:bool=None,
                  qa_write:bool=None):

    """
    To shift the telluric spectra to minimize residual telluric noise.

    Parameters
    ----------
    default_shiftranges : {True, False}
        Set to True to use the default ranges.
        Set to False to not use the default ranges.

    user_shiftranges : tuple, list, deafult None
        If a tuple, then a (3,) tuple as
        (order number, lower wavelength, upper wavelength).

        If a list, then a (norders,) list of (3,) tuples as
        (order number, lower wavelength, upper wavelength).
            
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
    Loads data into memory:

        tc.state['shift_ranges']
        tc.state['shifts']
        tc.state['shiftedtc_spectra']
    
    """

    #
    # Check the parameters and keywords
    #

    check_parameter('shift_spectra', 'default_shiftranges',
                    default_shiftranges, 'bool')

    check_parameter('shift_spectra', 'user_shiftranges', user_shiftranges,
                    ['tuple', 'list', 'NoneType'])    
                   
    check_parameter('shift_spectra', 'verbose', verbose, ['NoneType','bool'])

    check_parameter('shift_spectra', 'qa_show', qa_show, ['NoneType','bool'])

    check_parameter('shift_spectra', 'qa_showscale', qa_showscale,
                    ['NoneType','float','int'])

    check_parameter('shift_spectra', 'qa_showblock', qa_showblock,
                    ['NoneType','bool'])
    
    check_parameter('shift_spectra', 'qa_write', qa_write, ['NoneType','bool'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    logging.info(" Shifting spectra to minimize telluric noise.")

    #
    # Create an empty shift_ranges array
    #

    shift_orders = tc.state['object_orders']
    shift_ranges = np.full((tc.state['object_norders'],2), np.nan)

    #
    # Now update with the default ranges if requested.
    #
    
    if default_shiftranges is True:

        # Read the telluric_shiftinfo.dat file.
        
        file = osjoin(setup.state['instrument_path'],
                      'telluric_shiftinfo.dat')

        result = check_file(file, raise_error=False)

        if result is None:

            message = ' File telluric_shiftinfo.dat not found.  '+\
                'Skipping default shifts.'
            logging.info( message)

        else:
            
            mode, order, minwave, maxwave = np.loadtxt(file,
                                                       comments='#',
                                                       dtype='str',
                                                       unpack=True)
            
            order = order.astype(int)
            minwave = minwave.astype(float)
            maxwave = maxwave.astype(float)    
            
            #
            # yank out the pertinent information for the mode
            #
            
            z = mode == tc.state['mode']
            if not np.sum(z):

                message = ' Default shift ranges for mode '+tc.state['mode']+\
                    ' not found.  Skipping default shifts.'
                logging.info(message)
            
            default_orders = order[z]
            minwave = minwave[z]
            maxwave = maxwave[z]
        
            for i in range(len(shift_orders)):

                z = np.where(default_orders == shift_orders[i])[0]
                if len(z) == 1:
                    
                    shift_ranges[i,0] = minwave[z][0]
                    shift_ranges[i,1] = maxwave[z][0]                    
                                       
    #
    # Update the shifts_ranges for the user values
    #

    if isinstance(user_shiftranges,tuple):

        z = np.where(shift_orders == user_shiftranges[0])[0]
        if len(z) == 1:

            user_range = check_inrange(user_shiftranges,
                                       verbose=qa['verbose'])
            if user_range is not None:

                shift_ranges[z,:] = user_range

    tc.state['shift_ranges'] = shift_ranges

                
    #
    # Now do the shifts
    #

    shifts = np.zeros((tc.state['object_norders'],
                       tc.state['object_napertures']))

    for i in range(tc.state['object_norders']):

        if np.isnan(shift_ranges[i][1]):

            continue
        
        for j in range(tc.state['object_napertures']):

            idx = i + j*tc.state['object_napertures']

            shifts[i,j] = find_shift(tc.state['object_spectra'][idx,0,:],
                                     tc.state['object_spectra'][idx,1,:],
                                     tc.state['ewcorrectedtc_spectra'][idx,1,:],
                                    list(shift_ranges[i][0:2]))

    tc.state['shifts'] = np.round(shifts, decimals=2)

    #
    # Perform the shifts
    #

    for i in range(tc.state['object_norders']):
        
        for j in range(tc.state['object_napertures']):

            idx = i + j*tc.state['object_napertures']
    
            x = np.arange(len(tc.state['object_spectra'][idx,0,:]))
            tc_f = tc.state['ewcorrectedtc_spectra'][idx,1,:]
            tc_u = tc.state['ewcorrectedtc_spectra'][idx,2,:]            

            rtc_f, rtc_u = linear_interp1d(x+tc.state['shifts'][i,j],
                                           tc_f,
                                           x,
                                           input_u=tc_u)
            
            tc.state['shiftedtc_spectra'][idx,1,:] = rtc_f
            tc.state['shiftedtc_spectra'][idx,2,:] = rtc_u

    #
    # Do the QA plots
    #

    if qa['show'] is True and np.sum(tc.state['shifts']) != 0:
        
        plot_shifts(setup.plots['shifts'],
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    qa['showscale'],
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    tc.state['xlabel'],
                    tc.state['object_orders'],
                    tc.state['object_spectra'],
                    tc.state['ewcorrectedtc_spectra'],
                    tc.state['shiftedtc_spectra'],
                    tc.state['shift_ranges'],
                    tc.state['shifts'])
        
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False:

            pl.pause(1)
        
    if qa['write'] is True and np.sum(tc.state['shifts']) != 0:

        plot_shifts(None,
                    setup.plots['subplot_size'],
                    setup.plots['stack_max'],
                    setup.plots['font_size'],
                    1,
                    setup.plots['spectrum_linewidth'],
                    setup.plots['spine_linewidth'],
                    tc.state['xlabel'],
                    tc.state['object_orders'],
                    tc.state['object_spectra'],
                    tc.state['ewcorrectedtc_spectra'],
                    tc.state['shiftedtc_spectra'],
                    tc.state['shift_ranges'],
                    tc.state['shifts'])
        
        pl.savefig(osjoin(setup.state['qa_path'],
                          tc.state['output_filename']+ \
                          '_shifts' + \
                          setup.state['qa_extension']))
        pl.close()

            
def check_inrange(test:tuple,
                  verbose:bool=True):

    """
    To determine if a requested user shift range falls within the order

    Parameters
    ----------
    test : tuple
        A (3,) tuple where test[0] is the order number, test[1] is the
        lower wavelength limit and test[2] is the upper wavelength limit.

    verbose : {True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
    
    Returns
    -------
    None, tuple
        If user wavelength range falls in the order, a numpy array with the
        wavelength range.

        If user wavelength range does not fall in the order, None.
    
    """

    #
    # Check parameters
    #

    check_parameter('check_inrange', 'test', test, 'tuple')

    check_qakeywords(verbose=verbose)
    
    #
    # Do the check
    #
        
    # Rename the object wavelength ranges for ease.
    
    object_ranges = tc.state['standard_wavelengthranges']
    
    z = np.where(tc.state['object_orders'] == test[0])[0]

    # Does the requested order exist?
    
    if len(z) == 0:

        message = ' Requested telluric shift '+str(tc.state['mode'])+\
            ' order '+str(test[0])+' not extracted.  Skipping this order.'
        logging.info(message)

        return None
           
    # Is the range in order?

    z = int(z)
    lower = np.logical_and(test[1] >= object_ranges[z][0], 
                           test[0] <= object_ranges[z][1])
    
    upper = np.logical_and(test[2] >= object_ranges[z][0], 
                           test[1] <= object_ranges[z][1])
    
    if not lower*upper:
        
        message = ' Requested telluric shift range for '+\
            str(tc.state['mode'])+\
            ' Order '+str(test[0])+' of '+\
            ' and '.join([str(test[1]),str(test[2])])+\
            ' is out of range.  Skipping this order.'

        logging.info(message)

        return None
        
    return np.array(test[1:3])
