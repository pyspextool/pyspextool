from pyspextool.io import files
from pyspextool.io.files import *
import pytest
from os.path import exists


def test_check_file():

    result = check_file(files.__file__)
    assert exists(result) == True

    with pytest.raises(ValueError) as info:
        result = check_file(files.__file__+'a')


def test_extract_filestring():

    files = '1-3,5,7,10-12'
    result = extract_filestring(files, 'index')
    assert result == [1, 2, 3, 5, 7, 10, 11, 12]

    files = 'spec1.fits,spec2.fits'
    result = extract_filestring(files, 'filename')
    assert result == ['spec1.fits', 'spec2.fits']

    files = 'spec1.fits,spec2.fits'
    with pytest.raises(ValueError) as info:
        result = extract_filestring(files, 'filenames')
    

def test_make_full_path():

    files = '1-5'
    dir = '../data/'
    result = make_full_path(dir, files, indexinfo={'nint': 5,
                                                   'prefix': 'spc-',
                                                   'suffix': '.[ab].fits',
                                                   'extension': ''})
    assert result == ['../data/spc-00001.[ab].fits',
                      '../data/spc-00002.[ab].fits',
                      '../data/spc-00003.[ab].fits',
                      '../data/spc-00004.[ab].fits',
                      '../data/spc-00005.[ab].fits']
