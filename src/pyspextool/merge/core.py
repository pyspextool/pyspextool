import numpy as np
import numpy.typing as npt

from pyspextool.utils.arrays import trim_nan
from pyspextool.io.check import check_parameter
from pyspextool.utils.interpolate import linear_interp1d, linear_bitmask_interp1d
from pyspextool.utils.math import mean_data_stack, combine_flag_stack


def get_spectra_position(
    anchor_wavelength:npt.ArrayLike,
    add_wavelength:npt.ArrayLike):

    """
    Merges two, potentially overlapping, spectra into a single spectrum.

    Parameters
    ----------
    anchor_wavelength : ndarray
        An (nanchor,) array of wavelength values for the anchor spectrum.

    add_wavelength : ndarray
        An (ndd,) array of wavelength values for the add spectrum.

    Returns
    -------
    str {'left', 'left+overlap', 'right', right+overlap', 'inside', 
        'encompasses'}

    """

    #
    # Check parameters and qa values
    # 

    check_parameter("get_spectra_positions", "anchor_wavelength", 
                    anchor_wavelength, "ndarray")

    check_parameter("get_spectra_positions", "add_wavelength", 
                    add_wavelength, "ndarray")


    #
    # Get the minimum and maximum wavelengths
    #

    anchor_minwavelength = np.nanmin(anchor_wavelength)
    anchor_maxwavelength = np.nanmax(anchor_wavelength)

    add_minwavelength = np.nanmin(add_wavelength)
    add_maxwavelength = np.nanmax(add_wavelength)

    #
    # Figure out their relative locations
    #

    if (add_minwavelength < anchor_minwavelength and 
        add_maxwavelength > anchor_minwavelength and 
        add_maxwavelength < anchor_maxwavelength):

        label = 'left+overlap' 
        range = [anchor_minwavelength, add_maxwavelength]

        return label, range


    if (add_minwavelength < anchor_minwavelength and 
        add_maxwavelength < anchor_minwavelength):

        label = 'left'
        range = [add_maxwavelength, anchor_minwavelength]

        return label, range

    if (add_maxwavelength > anchor_maxwavelength and 
        add_minwavelength < anchor_maxwavelength and 
        add_minwavelength > anchor_minwavelength):

        label = 'right+overlap' 
        range = [add_minwavelength, anchor_maxwavelength]
        return label, range


    if (add_maxwavelength > anchor_maxwavelength and 
        add_minwavelength > anchor_maxwavelength):

        label = 'right'
        range = [anchor_maxwavelength, add_minwavelength]

        return label, range

    if (add_minwavelength < anchor_minwavelength and 
        add_maxwavelength > anchor_maxwavelength):

        label = 'encompass'
        range = [anchor_minwavelength, anchor_maxwavelength]

        return label, range

    if (add_minwavelength > anchor_minwavelength and 
        add_maxwavelength < anchor_maxwavelength):
        
        label = 'inside'
        range = [add_minwavelength, add_maxwavelength]

        return label, range
 

def merge_spectra(
    anchor_wavelength:npt.ArrayLike,
    anchor_intensity:npt.ArrayLike,
    add_wavelength:npt.ArrayLike,
    add_intensity:npt.ArrayLike,
    anchor_uncertainty:npt.ArrayLike=None,
    anchor_bitmask:npt.ArrayLike=None,
    add_uncertainty:npt.ArrayLike=None,
    add_bitmask:npt.ArrayLike=None):

    """
    Merges two, potentially overlapping, spectra into a single spectrum.

    Parameters
    ----------
    anchor_wavelength : ndarray
        An (nanchor,) array of wavelength values for the anchor spectrum.

    anchor_intensity : ndarray
        An (nanchor,) array of flux density values for the anchor spectrum.

    add_wavelength : ndarray
        An (nadd,) array of wavelength values for the add spectrum.

    add_intensity : ndarray
        An (nadd,) array of flux density values for the add spectrum.

    anchor_uncertainty : ndarray, default None
        An (nanchor,) array of uncertainty values for the anchor spectrum or 
        None.

    anchor_bitmask : ndarray, default None
        An (nanchor,) array of bit values for the anchor spectrum or None.

    add_uncertainty : ndarray, default None
        An (nadd,) array of uncertainty values for the add spectrum or None.

    add_bitmask : ndarray, default None
        A (nadd,) array of bit values for the add spectrum or None.

    Returns
    -------
    dict
        `"wavelength"` : ndarray
            An (nwaves,) array of merged wavelengths.

        `"intensity"` : ndarray
            An (nwaves,) array of merged flux density values.

        `"uncertainty"` : ndarray or None
            An (nwaves,) array of merged uncertainty values or None.

        `"bitmask"` : ndarray or None
            An (nwaves,) array of merged bitmask values or None.

        `"wavelength"` : ndarray
            An (nwaves,) array of merged wavelengths.
   
    Procedure
    ---------
    If the add spectrum does not overlap in wavelength with the anchor 
    spectrum, the two spectra are concatinated after setting the interior 
    edge pixels of the flux density and uncertainty to NaN and the bitmask 
    to zero.

    If the two spectra do overlap in wavelength, the add spectrum is first 
    linearly interpolated onto the wavelength scale of the anchor spectrum.  
    In the overlap region, the two spectra are (weighted) averaged.  

    There is an implicit assumption that the "add" order cannot cover the 
    entire wavelength range of the anchor_order.  This needs to be fixed.
   
    """

    #
    # Check parameters and qa values
    # 

    check_parameter("merge_spectra", "anchor_wavelength", anchor_wavelength, 
                    "ndarray")

    check_parameter("merge_spectra", "anchor_intensity", anchor_intensity, 
                    "ndarray", ndarray_size=np.size(anchor_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_spectra", "add_wavelength", add_wavelength, 
                    "ndarray")

    check_parameter("merge_spectra", "add_intensity", add_intensity, 
                    "ndarray", ndarray_size=np.size(add_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_spectra", "anchor_uncertainty", anchor_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_spectra", "anchor_bitmask", anchor_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength))

    check_parameter("merge_spectra", "add_uncertainty", add_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_spectra", "add_bitmask", add_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength))

    #
    # Check optional inputs
    #

    if anchor_uncertainty is None or add_uncertainty is None:

        no_uncertainty = True

    else:

        no_uncertainty = False

    if anchor_bitmask is None or add_bitmask is None:

        no_bitmask = True

    else:

        no_bitmask = False

    #
    # Trim NaNs at the edges of the two arrays
    #

    anchor_mask = trim_nan(anchor_intensity, flag=2)
    add_mask = trim_nan(add_intensity, flag=2)

    # Anchor data

    tanchor_wavelength = anchor_wavelength[anchor_mask]
    tanchor_intensity = anchor_intensity[anchor_mask]

    if no_uncertainty is False:

        tanchor_uncertainty = anchor_uncertainty[anchor_mask]

    else:

        tanchor_uncertainty = None

    if no_bitmask is False:

        tanchor_bitmask = anchor_bitmask[anchor_mask]

    else:

        tanchor_bitmask = None

    # Add data 

    tadd_wavelength = add_wavelength[add_mask]
    tadd_intensity = add_intensity[add_mask]

    if no_uncertainty is False:

        tadd_uncertainty = add_uncertainty[add_mask]

    else:

        tadd_uncertainty = None

    if no_bitmask is False:

        tadd_bitmask = add_bitmask[add_mask]

    else:

        tadd_bitmask = None

#    add_ndat = np.size(tadd_wavelength)

    #
    # Get the minimum and maximum wavelengths
    #

    tanchor_minwavelength = np.nanmin(tanchor_wavelength)
    tanchor_maxwavelength = np.nanmax(tanchor_wavelength)

    tadd_minwavelength = np.nanmin(tadd_wavelength)
    tadd_maxwavelength = np.nanmax(tadd_wavelength)

    #
    # Determine the position of the add spectrum relative to the anchor spectrum
    #

    if tadd_minwavelength < tanchor_minwavelength:

        position = 'Left'

    else:

        if tadd_maxwavelength > tanchor_maxwavelength:

            position = 'Right'

        else:

            position = 'Inside'

    #
    # Start the process
    # 

    if position == 'Right':

        result = _merge_onright(tanchor_wavelength,
                                tanchor_intensity,
                                tanchor_uncertainty,
                                tanchor_bitmask,
                                tadd_wavelength,
                                tadd_intensity,
                                tadd_uncertainty,
                                tadd_bitmask)

    if position == 'Left':

        result = _merge_onright(tadd_wavelength,
                               tadd_intensity,
                               tadd_uncertainty,
                               tadd_bitmask,
                               tanchor_wavelength,
                               tanchor_intensity,
                               tanchor_uncertainty,
                               tanchor_bitmask)

    if position == 'Inside':

        result = _merge_inmiddle(tanchor_wavelength,
                                 tanchor_intensity,
                                 tanchor_uncertainty,
                                 tanchor_bitmask,
                                 tadd_wavelength,
                                 tadd_intensity,
                                 tadd_uncertainty,
                                 tadd_bitmask)



    return result


def scale_order():

    """
    


    """




def _merge_inmiddle(anchor_wavelength:npt.ArrayLike,
                   anchor_intensity:npt.ArrayLike,
                   anchor_uncertainty:npt.ArrayLike | None,
                   anchor_bitmask:npt.ArrayLike | None,
                   add_wavelength:npt.ArrayLike,
                   add_intensity:npt.ArrayLike,
                   add_uncertainty:npt.ArrayLike | None,
                   add_bitmask:npt.ArrayLike | None):

    """
    Merges two spectra when the add spectrum is on the right.

    Parameters
    ----------
    anchor_wavelength : ndarray
        An (nanchor,) array of wavelength values for the anchor spectrum.

    anchor_intensity : ndarray
        An (nanchor,) array of flux density values for the anchor spectrum.

    anchor_uncertainty : ndarray or None
        An (nanchor,) array of uncertainty values for the anchor spectrum or 
        None.

    anchor_bitmask : ndarray or None
        An (nanchor,) array of bit values for the anchor spectrum or None.

    add_wavelength : ndarray
        An (nadd,) array of wavelength values for the add spectrum.

    add_intensity : ndarray
        An (nadd,) array of flux density values for the add spectrum.

    add_uncertainty : ndarray, default None
        An (nadd,) array of uncertainty values for the add spectrum or None.

    add_bitmask : ndarray, default None
        A (nadd,) array of bit values for the add spectrum or None.

    Returns
    -------
    dict
        `"wavelength"` : ndarray
            An (nwaves,) array of merged wavelengths.

        `"intensity"` : ndarray
            An (nwaves,) array of merged flux density values.

        `"uncertainty"` : ndarray or None
            An (nwaves,) array of merged uncertainty values or None.

        `"bitmask"` : ndarray or None
            An (nwaves,) array of merged bitmask values or None.

    """

    #
    # Check input parameters
    #

    check_parameter("merge_inmiddle", "anchor_wavelength", anchor_wavelength, 
                    "ndarray")

    check_parameter("merge_inmiddle", "anchor_intensity", anchor_intensity, 
                    "ndarray", ndarray_size=np.size(anchor_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_inmiddle", "add_wavelength", add_wavelength, 
                    "ndarray")

    check_parameter("merge_inmiddle", "add_intensity", add_intensity, 
                    "ndarray", ndarray_size=np.size(add_wavelength))

    check_parameter("merge_inmiddle", "anchor_uncertainty", anchor_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength))

    check_parameter("merge_inmiddle", "anchor_bitmask", anchor_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength))

    check_parameter("merge_inmiddle", "add_uncertainty", add_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength))

    check_parameter("merge_inmiddle", "add_bitmask", add_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength))

    #
    # Set defaults and find min and max wavelength ranges
    #


    ouncertainty = None
    obitmask = None

    add_minwavelength = np.nanmin(add_wavelength)
    add_maxwavelength = np.nanmax(add_wavelength)

    #
    # Start the merging.  
    #
        
    z = np.where(anchor_wavelength < add_minwavelength)[0]
    
    left_wavelength = anchor_wavelength[z]
    left_intensity = anchor_intensity[z]
    
    if anchor_uncertainty is not None:
        
        left_uncertainty = anchor_uncertainty[z]
        
        if anchor_bitmask is not None:
            
            left_bitmask = anchor_bitmask[z]
            
    # Now get the right part of the final spectrum.
            
    z = np.where(anchor_wavelength > add_maxwavelength)[0]

        
    right_wavelength = anchor_wavelength[z]
    right_intensity = anchor_intensity[z]

    if anchor_uncertainty is not None:

        right_uncertainty = anchor_uncertainty[z]

    if anchor_bitmask is not None:

        right_bitmask = anchor_bitmask[z]
        
    # Now do the overlap section.  
    
    if anchor_uncertainty is None:
            
        iflux = linear_interp1d(add_wavelength, 
                                add_intensity,
                                anchor_wavelength, 
                                leave_nans=True)
        
        mean = np.sum(np.vstack((iflux,anchor_intensity)),axis=0)/2
            
        # Figure which pixels are in the overlap by exploiting the fact
        # that the pixels outside of overlap are set to NaNs in the 
        # interpolation.  
        
        idx = trim_nan(iflux, flag=2)

        middle_intensity = mean[idx]
        middle_wavelength = anchor_wavelength[idx]

    else:

        iflux, iunc = linear_interp1d(add_wavelength, 
                                      add_intensity,
                                      anchor_wavelength, 
                                      input_u=add_uncertainty, 
                                      leave_nans=True)
        
        
        weights = 1/np.vstack((anchor_uncertainty,iunc))**2
        result = mean_data_stack(np.vstack((anchor_intensity,iflux)),
                                 weights=weights)
        
        imean = result[0]
        iunc = np.sqrt(result[1])
        
        # Figure which pixels are in the overlap by exploiting the fact
        # that the pixels outside of overlap are set to NaNs in the 
        # interpolation.  
        
        idx = trim_nan(iflux, flag=2)
        
        middle_intensity = imean[idx]
        middle_uncertainty = iunc[idx]
        middle_wavelength = anchor_wavelength[idx]
        
    if anchor_bitmask is not None:
            
        result = linear_bitmask_interp1d(add_wavelength, 
                                         add_bitmask,
                                         anchor_wavelength)
            
        bitmask = combine_flag_stack(np.vstack((result,anchor_bitmask)))

        
        middle_bitmask = bitmask[idx]
            
    #
    # Build the final spectrum
    #
    
    owavelength = np.concatenate((left_wavelength, 
                                  middle_wavelength, 
                                  right_wavelength))
    
    ointensity = np.concatenate((left_intensity, 
                                   middle_intensity, 
                                   right_intensity))
    
    if anchor_uncertainty is not None:
        
        ouncertainty = np.concatenate((left_uncertainty, 
                                       middle_uncertainty, 
                                       right_uncertainty))
        
    else:
        
        ouncertainty = None
        
    if anchor_bitmask is not None:
        
        obitmask = np.concatenate((left_bitmask, 
                                   middle_bitmask, 
                                   right_bitmask))
        
    else:
        
        obitmask = None

    #
    # Build the dictionary of results
    #

    dict = {'wavelength':owavelength,
            'intensity':ointensity,
            'uncertainty':ouncertainty,
            'bitmask':obitmask}

    return dict


            
def _merge_onright(anchor_wavelength:npt.ArrayLike,
                   anchor_intensity:npt.ArrayLike,
                   anchor_uncertainty:npt.ArrayLike | None,
                   anchor_bitmask:npt.ArrayLike | None,
                   add_wavelength:npt.ArrayLike,
                   add_intensity:npt.ArrayLike,
                   add_uncertainty:npt.ArrayLike | None,
                   add_bitmask:npt.ArrayLike | None):

    """
    Merges two spectra when the add spectrum is on the right.

    Parameters
    ----------
    anchor_wavelength : ndarray
        An (nanchor,) array of wavelength values for the anchor spectrum.

    anchor_intensity : ndarray
        An (nanchor,) array of flux density values for the anchor spectrum.

    anchor_uncertainty : ndarray or None
        An (nanchor,) array of uncertainty values for the anchor spectrum or 
        None.

    anchor_bitmask : ndarray or None
        An (nanchor,) array of bit values for the anchor spectrum or None.

    add_wavelength : ndarray
        An (nadd,) array of wavelength values for the add spectrum.

    add_intensity : ndarray
        An (nadd,) array of flux density values for the add spectrum.

    add_uncertainty : ndarray, default None
        An (nadd,) array of uncertainty values for the add spectrum or None.

    add_bitmask : ndarray, default None
        A (nadd,) array of bit values for the add spectrum or None.

    Returns
    -------
    dict
        `"wavelength"` : ndarray
            An (nwaves,) array of merged wavelengths.

        `"intensity"` : ndarray
            An (nwaves,) array of merged flux density values.

        `"uncertainty"` : ndarray or None
            An (nwaves,) array of merged uncertainty values or None.

        `"bitmask"` : ndarray or None
            An (nwaves,) array of merged bitmask values or None.

        `"wavelength"` : ndarray
            An (nwaves,) array of merged wavelengths.

    """

    #
    # Check input parameters
    #

    check_parameter("merge_onright", "anchor_wavelength", anchor_wavelength, 
                    "ndarray")

    check_parameter("merge_onright", "anchor_intensity", anchor_intensity, 
                    "ndarray", ndarray_size=np.size(anchor_wavelength),
                    ndarray_dtype='float')

    check_parameter("merge_onright", "add_wavelength", add_wavelength, 
                    "ndarray")

    check_parameter("merge_onright", "add_intensity", add_intensity, 
                    "ndarray", ndarray_size=np.size(add_wavelength))

    check_parameter("merge_onright", "anchor_uncertainty", anchor_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength))

    check_parameter("merge_onright", "anchor_bitmask", anchor_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(anchor_wavelength))

    check_parameter("merge_onright", "add_uncertainty", add_uncertainty, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength))

    check_parameter("merge_onright", "add_bitmask", add_bitmask, 
                    ["ndarray", "NoneType"], 
                    ndarray_size=np.size(add_wavelength))

    #
    # Set defaults and find min and max wavelength ranges
    #


    ouncertainty = None
    obitmask = None

    anchor_maxwavelength = np.nanmax(anchor_wavelength)
    add_minwavelength = np.nanmin(add_wavelength)

    #
    # Start the merging.  Do they overlap in wavelength?
    #

    if add_minwavelength > anchor_maxwavelength:

        #
        # They do not.  Just merge.  Set inside edge pixels to NaN/0.
        #
        
        # Wavelength

        owavelength = np.concatenate((anchor_wavelength, add_wavelength))
        
        # Flux density
        
        anchor_intensity[-1] = np.nan
        add_intensity[0] = np.nan
        
        ointensity = np.concatenate((anchor_intensity, add_intensity))
        
        # Uncertainty

        if anchor_uncertainty is not None and add_uncertainty is not None:
        
            anchor_uncertainty[-1] = np.nan
            add_uncertainty[0] = np.nan
        
            ouncertainty = np.concatenate((anchor_uncertainty, add_uncertainty))
        
        # Bitmask

        if anchor_bitmask is not None and add_bitmask is not None:

            anchor_bitmask[-1] = 0
            add_bitmask[0] = 0
            
            obitmask = np.concatenate((anchor_bitmask, add_bitmask))
                
    else:
        
        # They do overlap. Start by getting the left part of the final spectrum.
        
        z = np.where(anchor_wavelength < add_minwavelength)[0]
        
        left_wavelength = anchor_wavelength[z]
        left_intensity = anchor_intensity[z]

        if anchor_uncertainty is not None:

            left_uncertainty = anchor_uncertainty[z]

        if anchor_bitmask is not None:

            left_bitmask = anchor_bitmask[z]

        # Now get the right part of the final spectrum.
        
        z = np.where(add_wavelength > anchor_maxwavelength)[0]
        
        right_wavelength = add_wavelength[z]
        right_intensity = add_intensity[z]

        if anchor_uncertainty is not None:

            right_uncertainty = add_uncertainty[z]

        if anchor_bitmask is not None:

            right_bitmask = add_bitmask[z]
        
        # Now do the overlap section.  
    
        if anchor_uncertainty is None:
            
            iflux = linear_interp1d(add_wavelength, 
                                    add_intensity,
                                    anchor_wavelength, 
                                    leave_nans=True)
            
            mean = np.sum(np.vstack((iflux,anchor_intensity)),axis=0)/2
            
            # Figure which pixels are in the overlap by exploiting the fact
            # that the pixels outside of overlap are set to NaNs in the 
            # interpolation.  

            idx = trim_nan(iflux, flag=2)

            middle_intensity = mean[idx]
            middle_wavelength = anchor_wavelength[idx]

        else:

            iflux, iunc = linear_interp1d(add_wavelength, 
                                          add_intensity,
                                          anchor_wavelength, 
                                          input_u=add_uncertainty, 
                                          leave_nans=True)


            weights = 1/np.vstack((anchor_uncertainty,iunc))**2
            result = mean_data_stack(np.vstack((anchor_intensity,iflux)),
                                     weights=weights)

            imean = result[0]
            iunc = np.sqrt(result[1])

            # Figure which pixels are in the overlap by exploiting the fact
            # that the pixels outside of overlap are set to NaNs in the 
            # interpolation.  

            idx = trim_nan(iflux, flag=2)

            middle_intensity = imean[idx]
            middle_uncertainty = iunc[idx]
            middle_wavelength = anchor_wavelength[idx]

        if anchor_bitmask is not None:
            
            result = linear_bitmask_interp1d(add_wavelength, 
                                             add_bitmask,
                                             anchor_wavelength)
            
            bitmask = combine_flag_stack(np.vstack((result,anchor_bitmask)))


            middle_bitmask = bitmask[idx]

        #
        # Build the final spectrum
        #
    
        owavelength = np.concatenate((left_wavelength, 
                                      middle_wavelength, 
                                      right_wavelength))
        
        ointensity = np.concatenate((left_intensity, 
                                       middle_intensity, 
                                       right_intensity))
        
        if anchor_uncertainty is not None:
            
            ouncertainty = np.concatenate((left_uncertainty, 
                                           middle_uncertainty, 
                                           right_uncertainty))

        else:

            ouncertainty = None

        if anchor_bitmask is not None:
            
            obitmask = np.concatenate((left_bitmask, 
                                       middle_bitmask, 
                                       right_bitmask))

        else:

            obitmask = None

    #
    # Build the dictionary of results
    #

    dict = {'wavelength':owavelength,
            'intensity':ointensity,
            'uncertainty':ouncertainty,
            'bitmask':obitmask}

    return dict


