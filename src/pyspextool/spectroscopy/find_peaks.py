import numpy as np

from pyspextool.fit.fit_peak1d import *
from pyspextool.io.check_parameter import check_parameter
from pyspextool.utils.arrays import find_index
from pyspextool.utils.arrays import trim_nan

def find_peaks(profiles, method_info, fwhm=0.8):

    """
    To determine the locations of peaks in spatial profiles


    Parameters
    ----------
    profile : list
        A (norders,) list where each element is a dictionary.  
        The 'y' keyword contains the spatial positions in arcseconds.
        The 'p' keyword contains the spatial profile.


    method_info : dict
        The 'method' keyword can be, 'auto', 'guess', or 'fixed'.
        The 'apertures' keyword is the number of apertures to search for 
        if the 'method' is auto and a list of numpy.ndarray of apertures 
        if 'method' is 'guess' or 'fixed'.

    fwhm: float, default 0.8 (arcseconds), optional
        The approximate FWHM of the peak to be identified.  Only used 
        if `method` is 'auto' or 'guess'.


    Returns
    -------
    numpy.ndarray, numpy.ndarray
    
    apertures : An (norders, naps) numpy.ndarray of aperture positions.  
                By Spextool convention, apertures[0,:] are the apertures for
                the order closest to the bottom of the image.

    apsigns : An (norders, naps) numpy.ndarray of aperture signs 
              (1 = positive, -1=negative).  By Spextool convention, 
              apsigns[0,:] are the aperture signs the order closest to the 
              bottom of the image.


    Notes
    -----


    Examples
    --------
    later
    
    """

    #
    # Check parameters
    #
    
    check_parameter('find_apertures', 'profiles', profiles, 'list')
    
    if len(method_info) != 2:

        message = 'method_info must contain two elements.'
        raise ValueError(message)


    # Get set up

    norders = len(profiles)

    

    # Run things by the type of method requested

    if method_info['method'] == 'auto':

        npeaks = method_info['peaks']

        # Create arrays
        
        peaks = np.empty((norders, npeaks),dtype=float)
        signs = np.empty((norders, npeaks),dtype=int)
       
        # loop over the orders
        
        for i in range(norders):

            x = profiles[i]['y']
            y = profiles[i]['p']
            
            # Trim NaNs
            
            good = trim_nan(y)
            x = x[good]
            y = y[good]
           
            # Subtract the background
            
            y = y-np.median(y)        
            absy = np.abs(y)
        
            for j in range(npeaks):

                # Find the largest peak

                idx = int(np.argmax(absy))

                fit = fit_peak1d(x,absy,p0=[absy[idx], x[idx], fwhm/2.354, 0])
                peaks[i,j] = fit['parms'][1]
               
                # Get the aperture sign

                idx = int(find_index(x,fit['parms'][1]))
                
                signs[i,j] = 1 if y[idx] > 0 else -1
               
                # Subtract the results off

                absy = absy-fit['fit']
                
    elif method_info['method'] == 'guess':

        # Get apertures
        
        guess_peaks = np.array(method_info['peaks'])
        npeaks = int(np.size(guess_peaks)/norders)
        
        # Create arrays
        
        peaks = np.empty((norders, npeaks),dtype=float)
        signs = np.empty((norders, npeaks),dtype=int)
                
        # loop over the orders
    
        for i in range(norders):

            x = profiles[i]['y']
            y = profiles[i]['p']
            
            # Trim NaNs
            
            good = trim_nan(y)
            x = x[good]
            y = y[good]
                
            # Subtract the background
            
            y = y-np.median(y)
                        
            # Now loop over each aperture and do the fit
        
            for j in range(npeaks):
                                
                # Get the peak value for the guess
                
                idx = int(find_index(x,guess_peaks[i, j]))
                
                # Do the fit
                
                fit = fit_peak1d(x,y,p0=[y[idx], guess_peaks[i,j],
                                 fwhm/2.354, 0])
                peaks[i,j] = fit['parms'][1]

                # Get the aperture sign

                signs[i,j] = 1 if fit['parms'][0] > 0 else -1

    elif method_info['method'] == 'fixed':

      # Get apertures
        
        guess_peaks = np.array(method_info['peaks'])
        npeaks = int(np.size(guess_peaks)/norders)
        
        # Create arrays
        
        peaks = np.empty((norders, npeaks),dtype=float)
        signs = np.empty((norders, npeaks),dtype=int)
                
        # loop over the orders
    
        for i in range(norders):

            x = profiles[i]['y']
            y = profiles[i]['p']
            
            # Trim NaNs
            
            good = trim_nan(y)
            x = x[good]
            y = y[good]
                
            # Subtract the background
            
            y = y-np.median(y)

            for j in range(npeaks):
                                
                # Get the peak value for the guess
                
                idx = find_index(x,guess_peaks[i, j])
                
                # Store results
                
                peaks[i,j] = guess_peaks[i,j]

                # Get the aperture sign

                signs[i,j] = 1 if y[int(idx)] > 0 else -1

    else:

        message = 'Unknown search method.'
        raise ValueError(message)


    # Now sort everything properly

    for i in range(norders):

        idx = np.argsort(peaks[i,:])
        peaks[i,:] = peaks[i,idx]
        signs[i,:] = signs[i,idx]        

    
    return peaks, signs
