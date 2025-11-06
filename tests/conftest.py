import pytest


@pytest.fixture
def raw_setup():
    raw_setup = {
        "spex_lxd": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-LXD/data/",
            "cal_path": "tests/test_data/raw/spex-LXD/cals/",
            "proc_path": "tests/test_data/raw/spex-LXD/proc/",
            "qa_path": "tests/test_data/raw/spex-LXD/qa/",
        },
        "spex_prism": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-prism/data/",
            "cal_path": "tests/test_data/raw/spex-prism/cals/",
            "proc_path": "tests/test_data/raw/spex-prism/proc/",
            "qa_path": "tests/test_data/raw/spex-prism/qa/",
        },
        "spex_sxd": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-SXD/data/",
            "cal_path": "tests/test_data/raw/spex-SXD/cals/",
            "proc_path": "tests/test_data/raw/spex-SXD/proc/",
            "qa_path": "tests/test_data/raw/spex-SXD/qa/",
        },
        "uspex_lxd": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-LXD/data/",
            "cal_path": "tests/test_data/raw/uspex-LXD/cals/",
            "proc_path": "tests/test_data/raw/uspex-LXD/proc/",
            "qa_path": "tests/test_data/raw/uspex-LXD/qa/",
        },
        "uspex_prism": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-prism/data/",
            "cal_path": "tests/test_data/raw/uspex-prism/cals/",
            "proc_path": "tests/test_data/raw/uspex-prism/proc/",
            "qa_path": "tests/test_data/raw/uspex-prism/qa/",
        },
        "uspex_sxd": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-SXD/data/",
            "cal_path": "tests/test_data/raw/uspex-SXD/cals/",
            "proc_path": "tests/test_data/raw/uspex-SXD/proc/",
            "qa_path": "tests/test_data/raw/uspex-SXD/qa/",
        },
    }
    return raw_setup


@pytest.fixture
def proc_setup():
    proc_setup = {
        "spex_lxd": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-LXD/data/",
            "cal_path": "tests/test_data/raw/spex-LXD/cals/",
            "proc_path": "tests/test_data/processed/spex-LXD/proc/",
            "qa_path": "tests/test_data/processed/spex-LXD/qa/",
            "standard_file": "combspec36-39.fits",
        },
        "spex_prism": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-prism/data/",
            "cal_path": "tests/test_data/raw/spex-prism/cals/",
            "proc_path": "tests/test_data/processed/spex-prism/proc/",
            "qa_path": "tests/test_data/processed/spex-prism/qa/",
            "standard_file": "spectra0009.fits",
        },
        "spex_sxd": {
            "instrument": "spex",
            "raw_path": "tests/test_data/raw/spex-SXD/data/",
            "cal_path": "tests/test_data/raw/spex-SXD/cals/",
            "proc_path": "tests/test_data/processed/spex-SXD/proc/",
            "qa_path": "tests/test_data/processed/spex-SXD/qa/",
            "standard_file": "spectra0640.fits",
        },
        "uspex_lxd": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-LXD/data/",
            "cal_path": "tests/test_data/raw/uspex-LXD/cals/",
            "proc_path": "tests/test_data/processed/uspex-LXD/proc/",
            "qa_path": "tests/test_data/processed/uspex-LXD/qa/",
            "standard_file": "HD223352.fits",
        },
        "uspex_prism": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-prism/data/",
            "cal_path": "tests/test_data/raw/uspex-prism/cals/",
            "proc_path": "tests/test_data/processed/uspex-prism/proc/",
            "qa_path": "tests/test_data/processed/uspex-prism/qa/",
            "standard_file": "spectra00007.fits",
        },
        "uspex_sxd": {
            "instrument": "uspex",
            "raw_path": "tests/test_data/raw/uspex-SXD/data/",
            "cal_path": "tests/test_data/raw/uspex-SXD/cals/",
            "proc_path": "tests/test_data/processed/uspex-SXD/proc/",
            "qa_path": "tests/test_data/processed/uspex-SXD/qa/",
            "standard_file": "spectra00010.fits"
        },
    }
    return proc_setup


@pytest.fixture
def postextraction_setup():
    """
    Fixture for PostExtractionTests data.
    Contains extracted spectra ready for combining and telluric correction.
    """
    postextraction_setup = {
        "spex_prism": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-LowRes15/",
            "object_files": ['spectra', '17-28'],
            "standard_files": ['spectra', '29-40'],
            "standard_name": "HD 101060",
        },
        "spex_sxd": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-ShortXD/",
            "object_files": ['spectra', '626-635'],
            "standard_files": ['spectra', '636-645'],
            "standard_name": "HD 165029",
        },
        "spex_lxd": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-LongXD2.1/",
            "object_files": ['spectra', '585-594'],
            "standard_files": ['spectra', '595-614'],
            "standard_name": "HD 165029",
        },
        "uspex_prism": {
            "instrument": "uspex",
            "proc_path": "tests/test_data/PostExtractionTests/uspex-prism/",
            "object_files": ['spectra', '1-2'],
            "standard_files": ['spectra', '7-8'],
            "standard_name": "HD 223352",
        },
        "uspex_lxd_short": {
            "instrument": "uspex",
            "proc_path": "tests/test_data/PostExtractionTests/uspex-LXD_Short/",
            "object_files": ['spectra', '11-20'],
            "standard_files": ['spectra', '1-10'],
            "standard_name": "HD 17778",
        },
    }
    return postextraction_setup
