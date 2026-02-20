import os
import pytest

_TEST_DATA_PATH = os.path.join(os.path.dirname(__file__), "test_data")
_TEST_DATA_AVAILABLE = os.path.isdir(_TEST_DATA_PATH) and bool(os.listdir(_TEST_DATA_PATH))


def pytest_collection_modifyitems(config, items):
    """Skip tests that require the test_data submodule when it is not available."""
    if not _TEST_DATA_AVAILABLE:
        skip_marker = pytest.mark.skip(reason="test_data submodule not cloned")
        requires_test_data_files = {
            "test_combine.py",
            "test_locate_orders.py",
            "test_make_flat.py",
            "test_make_wavecal.py",
            "test_override_aperturesigns.py",
            "test_telluric.py",
        }
        requires_test_data_tests = {
            ("test_files.py", "test_inoutfiles_to_fullpaths"),
        }
        for item in items:
            basename = os.path.basename(str(item.fspath))
            if basename in requires_test_data_files:
                item.add_marker(skip_marker)
            elif (basename, item.originalname) in requires_test_data_tests:
                item.add_marker(skip_marker)


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
        "spex_lxd_1_9": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-LongXD1.9/",
            "object_files": ['spectra', '200-205'],
            "standard_files": ['spectra', '176-189'],
            "standard_name": "HD 7927",
        },
        "spex_lxd_2_1": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-LongXD2.1/",
            "object_files": ['spectra', '585-594'],
            "standard_files": ['spectra', '595-614'],
            "standard_name": "HD 165029",
        },
        "spex_lxd_2_3": {
            "instrument": "spex",
            "proc_path": "tests/test_data/PostExtractionTests/spex-LongXD2.3/",
            "object_files": ['spectra', '555-564'],
            "standard_files": ['spectra', '565-574'],
            "standard_name": "BS3314",
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
        "uspex_lxd_long": {
            "instrument": "uspex",
            "proc_path": "tests/test_data/PostExtractionTests/uspex-LXD_Long/",
            "object_files": ['spectra', '23-32'],
            "standard_files": ['spectra', '33-42'],
            "standard_name": "HD 223352",
        },
        "uspex_sxd": {
            "instrument": "uspex",
            "proc_path": "tests/test_data/PostExtractionTests/uspex-ShortXD/",
            "object_files": ['spectra', '1-6'],
            "standard_files": ['spectra', '9-16'],
            "standard_name": "HD 222332",
        },
    }
    return postextraction_setup
