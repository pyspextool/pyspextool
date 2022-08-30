import numpy as np
from unittest import TestCase
from pyspextool.utils.math import *

#
# =========================== bit_set test =============================
#
def test_bit_set():


    array = np.array([3,4,1])
    result = bit_set(array,0)
    np.testing.assert_array_equal(result,[1,0,1])
    
    array = np.array([3,4,1])
    result = bit_set(array,[0,3])
    np.testing.assert_array_equal(result,[1,0,1])
    
    array = np.array([3,4,1])
    result = bit_set(array,[2,3])
    np.testing.assert_array_equal(result,[0,1,0])

#
# ============================ round test ==============================
#
def test_round():

    array = np.array([-3.5,-2.5,2.5,3.5])
    result = round(array)
    np.testing.assert_array_equal(result,[-4,-3,3,4])

#    
# ======================== combine_flag_stack ==========================
#
def test_combine_flag_stack():

    # Do a spectral stack
    
    spec1 = np.array([0,4,4,2])
    spec2 = np.array([0,0,3,1])
    stack = np.stack((spec1,spec2))

    result = combine_flag_stack(stack)
    np.testing.assert_array_equal(result,[0,4,7,3])


    # Do an image stack
    
    img1 = np.array([[0,2,0],[3,0,4],[0,0,0]])
    img2 = np.array([[1,0,0],[1,0,0],[0,0,0]])
    stack = np.stack((img1,img2))
    result = combine_flag_stack(stack)

    np.testing.assert_array_equal(result,np.array([[1,2,0],[3,0,4],[0,0,0]]))

#    
# ========================== median_data_stack ==========================
#
def test_median_data_stack():

    # Do a spectral stack with mask

    ss = np.array([[1,2,3,7],[0,5,8,2],[2,9,4,6]])
    msk = np.ones((3,4),dtype=int)
    msk[0,0] = 0
    med, unc = median_data_stack(ss, mask=msk)
    np.testing.assert_array_equal(med,[1,5,4,6])
    np.testing.assert_allclose(unc,[1.04835651,2.56793853,0.85597951,
                                    0.85597951],0.001)

    # Do an image stack with mask

    
    istack = np.array([[[1,2,3],[4,5,6],[7,8,9]],
                       [[6,3,1],[9,2,4],[1,5,0]],
                       [[3,4,9],[5,7,7],[3,9,1]],
                       [[1,6,5],[2,1,9],[5,2,7]]])              
    msk = np.ones((4,3,3),dtype=int)
    msk[0,0,0] = 0
    med, unc = median_data_stack(istack, mask=msk)
    
    np.testing.assert_array_equal(med,np.array([[3.,3.5,4.],[4.5,3.5,6.5],
                                                    [4.,6.5,4.]]))

    np.testing.assert_allclose(unc,np.array([[1.71195902,0.7413,1.4826],
                                             [1.11195,1.4826,1.11195],
                                             [1.4826,1.4826,2.59455]]))

#    
# ========================== scale_data_stack ============================
#
def test_scale_data_stack():

    # Spectral stack data set

    spec1 = np.array([1,1,1,1])
    spec2 = np.array([2,2,2,2])
    spec3 = np.array([4,4,4,4])
    specstack = np.stack((spec1, spec2, spec3))    
    
    # spectral stack, scale to median, no variance

    scaledstack, scaledvar, scales = scale_data_stack(specstack,None)
    
    np.testing.assert_array_equal(scaledstack,np.stack((spec2,spec2,spec2)))
    assert scaledvar == None
    np.testing.assert_array_equal(scales,[2,1,0.5])
    
    # spectral stack, scale to median, with variance

    scaledstack, scaledvar, scales = scale_data_stack(specstack,specstack)

    np.testing.assert_array_equal(scaledstack,np.stack((spec2,spec2,spec2)))
    np.testing.assert_array_equal(scaledvar,np.stack((spec3,spec2,spec1)))
    np.testing.assert_array_equal(scales,[2,1,0.5])    
    
    # spectral stack, scale to first spectrum, no variance

    scaledstack, scaledvar, scales = scale_data_stack(specstack,None,index=0)

    np.testing.assert_array_equal(scaledstack,np.stack((spec1,spec1,spec1)))
    assert scaledvar == None
    np.testing.assert_array_equal(scales,[1,0.5,0.25])

    
    # spectral stack, scale to first spectrum, with variance

    scaledstack, scaledvar, scales = scale_data_stack(specstack,specstack,
                                                      index=0)
    
    np.testing.assert_array_equal(scaledstack,np.stack((spec1,spec1,spec1)))
    np.testing.assert_array_equal(scaledvar,np.stack((spec1,1/spec2,1/spec3)))
    np.testing.assert_array_equal(scales,[1,0.5,0.25])

    # Image stack data

    img1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
    img2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
    img3 = np.array([[4,4,4],[4,4,4],[4,4,4]])
    imgstack = np.stack((img1,img2,img3))    
    
    # image stack, scale to median, no variance

    scaledstack, scaledvar, scales = scale_data_stack(imgstack, None)   

    np.testing.assert_array_equal(scaledstack,np.stack((img2,img2,img2)))
    assert scaledvar == None
    np.testing.assert_array_equal(scales,[2,1,0.5])

    # image stack, scale to median, with variance

    scaledstack, scaledvar, scales = scale_data_stack(imgstack, imgstack)   

    np.testing.assert_array_equal(scaledstack,np.stack((img2,img2,img2)))
    np.testing.assert_array_equal(scaledvar,np.stack((img3,img2,img1)))
    np.testing.assert_array_equal(scales,[2,1,0.5])        
    
    # image stack, scale to first image, no variance

    scaledstack, scaledvar, scales = scale_data_stack(imgstack, None, index=0)
    
    np.testing.assert_array_equal(scaledstack,np.stack((img1,img1,img1)))
    assert scaledvar == None
    np.testing.assert_array_equal(scales,[1,0.5,0.25])

    # image stack, scale to first image, no variance

    scaledstack, scaledvar, scales = scale_data_stack(imgstack, imgstack,
                                                      index=0)   

    np.testing.assert_array_equal(scaledstack,np.stack((img1,img1,img1)))
    np.testing.assert_array_equal(scaledvar,np.stack((img1,1/img2,1/img3)))
    np.testing.assert_array_equal(scales,[1,0.5,0.25])


#
#  ============================== moments =============================
#
def test_moments():

    x = 1
    assert x == 1
    
    
