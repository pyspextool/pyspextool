import numpy as np
import pyspextool as ps
from pyspextool.merge.core import merge_spectra
import pytest
from pyspextool.pyspextoolerror import pySpextoolError

def test_merge_spectra():

    #
    # Testing add to the right
    #

    anchor_wavelength  = np.array([1,2,3,4,5,6])
    anchor_intensity = np.array([1,np.nan,1,1,1,1])
    anchor_uncertainty = np.array([1.,1,1,1,1,1])
    anchor_bitmask = np.array([0,0,1,0,0,0])

    add_wavelength  = np.array([4.5, 5.5, 6.5, 7])
    add_intensity = np.array([1. , 1  , 1  , 1])
    add_uncertainty = np.array([2.,2,2,2])
    add_bitmask = np.array([0,1,0,0])


    # Data+bitmask+uncertainty w/ overlap

    result = merge_spectra(anchor_wavelength,
                           anchor_intensity,
                           add_wavelength,
                           add_intensity,
                           anchor_uncertainty=anchor_uncertainty,
                           add_uncertainty=add_uncertainty,
                           anchor_bitmask=anchor_bitmask,
                           add_bitmask=add_bitmask)

    np.testing.assert_array_equal(result['wavelength'], 
                                  np.array([1,2,3,4,5,6,6.5,7]))

    np.testing.assert_array_equal(result['intensity'],
                                  np.array([1,np.nan,1,1,1,1,1,1]))


    np.testing.assert_array_equal(result['bitmask'],
                                  np.array([0,0,1,0,1,1,0,0]))


    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([1,1,1,1,0.9258201,0.9258201,2,2]))

    # Data+bitmask+uncertainty w/o overlap

    anchor_wavelength  = np.array([1,2,3,4,5,6.])
    anchor_intensity = np.array([1,np.nan,1,1,1,1.])
    anchor_uncertainty = np.array([1,1,1,1,1,1.])
    anchor_bitmask = np.array([0,1,0,0,0,0])

    add_wavelength  = np.array([10 ,11,12.])
    add_intensity = np.array([1.,1,1.])
    add_uncertainty = np.array([1,1,1.])
    add_bitmask = np.array([1,1,0])

    result = merge_spectra(anchor_wavelength,
                           anchor_intensity,
                           add_wavelength,
                           add_intensity,
                           anchor_uncertainty=anchor_uncertainty,
                           add_uncertainty=add_uncertainty,
                           anchor_bitmask=anchor_bitmask,
                           add_bitmask=add_bitmask)


    np.testing.assert_array_equal(result['wavelength'], 
                                  np.array([1,2,3,4,5,6,10, 11, 12]))

    np.testing.assert_array_equal(result['intensity'],
                                  np.array([1,np.nan,1,1,1,np.nan,np.nan,1,1]))


    np.testing.assert_array_equal(result['bitmask'],
                                  np.array([0,1,0,0,0,0,0,1,0]))

    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([1,1,1,1,1,np.nan,np.nan,1,1]))


    #
    # Testing add to the left
    #

    # Data+bitmask+uncertainty w/o overlap


    add_wavelength  = np.array([1,2,3,4,5,6])
    add_intensity = np.array([1,np.nan,1,1,1,1])
    add_uncertainty = np.array([1.,1,1,1,1,1])
    add_bitmask = np.array([0,0,1,0,0,0])

    anchor_wavelength  = np.array([4.5, 5.5, 6.5, 7])
    anchor_intensity = np.array([1. , 1  , 1  , 1])
    anchor_uncertainty = np.array([2.,2,2,2])
    anchor_bitmask = np.array([0,1,0,0])


    result = merge_spectra(add_wavelength,
                           add_intensity,
                           anchor_wavelength,
                           anchor_intensity,
                           anchor_uncertainty=add_uncertainty,
                           add_uncertainty=anchor_uncertainty,
                           anchor_bitmask=add_bitmask,
                           add_bitmask=anchor_bitmask)

    
    np.testing.assert_array_equal(result['wavelength'], 
                                  np.array([1,2,3,4,5,6,6.5,7]))

    np.testing.assert_array_equal(result['intensity'],
                                  np.array([1,np.nan,1,1,1,1,1,1]))


    np.testing.assert_array_equal(result['bitmask'],
                                  np.array([0,0,1,0,1,1,0,0]))


    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([1,1,1,1,0.9258201,0.9258201,2,2]))


    # Data+bitmask+uncertainty w/o overlap

    add_wavelength  = np.array([1,2,3,4,5,6.])
    add_intensity = np.array([1,np.nan,1,1,1,1.])
    add_uncertainty = np.array([1,1,1,1,1,1.])
    add_bitmask = np.array([0,1,0,0,0,0])

    anchor_wavelength  = np.array([10 ,11,12.])
    anchor_intensity = np.array([1.,1,1.])
    anchor_uncertainty = np.array([1,1,1.])
    anchor_bitmask = np.array([1,1,0])

    result = merge_spectra(add_wavelength,
                           add_intensity,
                           anchor_wavelength,
                           anchor_intensity,
                           anchor_uncertainty=add_uncertainty,
                           add_uncertainty=anchor_uncertainty,
                           anchor_bitmask=add_bitmask,
                           add_bitmask=anchor_bitmask)


    np.testing.assert_array_equal(result['wavelength'], 
                                  np.array([1,2,3,4,5,6,10, 11, 12]))

    np.testing.assert_array_equal(result['intensity'],
                                  np.array([1,np.nan,1,1,1,np.nan,np.nan,1,1]))


    np.testing.assert_array_equal(result['bitmask'],
                                  np.array([0,1,0,0,0,0,0,1,0]))

    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([1,1,1,1,1,np.nan,np.nan,1,1]))

    #
    # Testing add to middle    
    #

    anchor_wavelength  = np.array([1,2,3,4,5,6])
    anchor_intensity = np.array([1,np.nan,1,1,1,1])
    anchor_uncertainty = np.array([1.,1,1,1,1,1])
    anchor_bitmask = np.array([0,0,1,0,0,0])

    add_wavelength  = np.array([2.5, 3.5, 4.5, 5.5])
    add_intensity = np.array([1. , 1  , 1  , 1])
    add_uncertainty = np.array([2.,2,2,2])
    add_bitmask = np.array([0,1,0,0])

    result = merge_spectra(anchor_wavelength,
                           anchor_intensity,
                           add_wavelength,
                           add_intensity,
                           anchor_uncertainty=anchor_uncertainty,
                           add_uncertainty=add_uncertainty,
                           anchor_bitmask=anchor_bitmask,
                           add_bitmask=add_bitmask)

    np.testing.assert_array_equal(result['wavelength'], 
                                  np.array([1,2,3,4,5,6]))

    np.testing.assert_array_equal(result['intensity'],
                                  np.array([1,np.nan,1,1,1,1]))


    np.testing.assert_array_equal(result['bitmask'],
                                  np.array([0,0,1,1,0,0]))


    np.testing.assert_allclose(result['uncertainty'],
                                  np.array([1,1,0.9258201,0.9258201,0.9258201, 1]))


