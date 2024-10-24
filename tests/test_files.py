from os.path import join
from pyspextool.io import files
from pyspextool.io.files import check_file, extract_filestring, make_full_path
#from pyspextool.io.convert_fits import convert_to_fits, spectrum_isplottable
import pytest
from os.path import exists


def test_check_file():

    result = check_file(files.__file__)
    print(result)
    assert exists(result) is True

#    with pytest.raises(ValueError):
#        result = check_file(files.__file__+'a')


def test_extract_filestring():

    files = '1-3,5,7,10-12'
    result = extract_filestring(files, 'index')
    assert result == [1, 2, 3, 5, 7, 10, 11, 12]

    files = 'spec1.fits,spec2.fits'
    result = extract_filestring(files, 'filename')
    assert result == ['spec1.fits', 'spec2.fits']

    files = 'spec1.fits,spec2.fits'
    with pytest.raises(ValueError):
        result = extract_filestring(files, 'filenames')


def test_make_full_path():
    files = '626-628'
    dir = 'tests/test_data/processed/spex-SXD/proc/'
    result = make_full_path(dir, files, indexinfo={'nint': 4,
                                                   'prefix': 'spectra',
                                                   'suffix': '.fits',
                                                   'extension': ''})
    print(result)
    assert result == ['tests/test_data/processed/spex-SXD/proc/spectra0626.fits',
                      'tests/test_data/processed/spex-SXD/proc/spectra0627.fits',
                      'tests/test_data/processed/spex-SXD/proc/spectra0628.fits']


#@pytest.mark.parametrize(
#    "file,output_path,out_file",
#    [
#        (
#            "tests/test_data/processed/spex-SXD/proc/combspec626-635.fits",
#            "tests/test_data/processed/spex-SXD/proc/",
#            "HD_160365_2003Jul07.fits",
#        ),
#        (
#            "tests/test_data/processed/spex-SXD/proc/combspec636-645.fits",
#            "tests/test_data/processed/spex-SXD/proc/",
#            "HD_165029_2003Jul07.fits",
#        ),
#        (
#            "tests/test_data/processed/uspex-prism/proc/combspec1-2.fits",
#            "tests/test_data/processed/uspex-prism/proc/",
#            "2010-1707_2022Oct19.fits",
#        ),
#        (
#            "tests/test_data/processed/uspex-prism/proc/combspec7-8.fits",
#            "tests/test_data/processed/uspex-prism/proc/",
#            "HD_193689_2022Oct19.fits",
#        ),
#        (
#            "tests/test_data/processed/uspex-SXD/proc/combspec1-8.fits",
#            "tests/test_data/processed/uspex-SXD/proc/",
#            "HD100906_G9w+_2015Jun03.fits",
#        ),
#        (
#            "tests/test_data/processed/uspex-SXD/proc/combspec11-18.fits",
#            "tests/test_data/processed/uspex-SXD/proc/",
#            "HD101369_A0V_2015Jun03.fits",
#        ),
#    ],
#)
#
#
#def test_convert_to_fits(file, output_path, out_file):
#    convert_to_fits(file, output_path=output_path)
#    assert exists(join(output_path, out_file)) is True


#@pytest.mark.parametrize(
#    "file",
#    [
#        "tests/test_data/processed/spex-SXD/proc/HD_160365_2003Jul07.fits",
#        "tests/test_data/processed/spex-SXD/proc/HD_165029_2003Jul07.fits",
#        "tests/test_data/processed/uspex-prism/proc/2010-1707_2022Oct19.fits",
#        "tests/test_data/processed/uspex-prism/proc/HD_193689_2022Oct19.fits",
#        "tests/test_data/processed/uspex-SXD/proc/HD100906_G9w+_2015Jun03.fits"#,
#        "tests/test_data/processed/uspex-SXD/proc/HD101369_A0V_2015Jun03.fits",
#    ],
#)

#@pytest.mark.skip(reason='issues with organization of prism data currently prevent this test')
#def test_spectrum_isplottable(file):
#    result = spectrum_isplottable(file, raise_error=True, show_plot=False)
#    assert result is True
