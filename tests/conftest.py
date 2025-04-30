import pytest


@pytest.fixture
def raw_paths():
    raw_paths = {
        "spex_prism": {
            "raw_path": "tests/test_data/raw/spex-prism/data/",
            "cal_path": "tests/test_data/raw/spex-prism/cals/",
            "proc_path": "tests/test_data/raw/spex-prism/proc/",
            "qa_path": "tests/test_data/raw/spex-prism/qa/",
        },
        "spex_sxd": {
            "raw_path": "tests/test_data/raw/spex-SXD/data/",
            "cal_path": "tests/test_data/raw/spex-SXD/cals/",
            "proc_path": "tests/test_data/raw/spex-SXD/proc/",
            "qa_path": "tests/test_data/raw/spex-SXD/qa/",
        },
        "uspex_prism": {
            "raw_path": "tests/test_data/raw/uspex-prism/data/",
            "cal_path": "tests/test_data/raw/uspex-prism/cals/",
            "proc_path": "tests/test_data/raw/uspex-prism/proc/",
            "qa_path": "tests/test_data/raw/uspex-prism/qa/",
        },
        "uspex_sxd": {
            "raw_path": "tests/test_data/raw/uspex-SXD/data/",
            "cal_path": "tests/test_data/raw/uspex-SXD/cals/",
            "proc_path": "tests/test_data/raw/uspex-SXD/proc/",
            "qa_path": "tests/test_data/raw/uspex-SXD/qa/",
        },
    }
    return raw_paths


@pytest.fixture
def proc_paths():
    proc_paths = {
        "spex_prism": {
            "raw_path": "tests/test_data/raw/spex-prism/data/",
            "cal_path": "tests/test_data/raw/spex-prism/cals/",
            "proc_path": "tests/test_data/processed/spex-prism/proc/",
            "qa_path": "tests/test_data/processed/spex-prism/qa/",
        },
        "spex_sxd": {
            "raw_path": "tests/test_data/raw/spex-SXD/data/",
            "cal_path": "tests/test_data/raw/spex-SXD/cals/",
            "proc_path": "tests/test_data/processed/spex-SXD/proc/",
            "qa_path": "tests/test_data/processed/spex-SXD/qa/",
        },
        "uspex_prism": {
            "raw_path": "tests/test_data/raw/uspex-prism/data/",
            "cal_path": "tests/test_data/raw/uspex-prism/cals/",
            "proc_path": "tests/test_data/processed/uspex-prism/proc/",
            "qa_path": "tests/test_data/processed/uspex-prism/qa/",
        },
        "uspex_sxd": {
            "raw_path": "tests/test_data/raw/uspex-SXD/data/",
            "cal_path": "tests/test_data/raw/uspex-SXD/cals/",
            "proc_path": "tests/test_data/processed/uspex-SXD/proc/",
            "qa_path": "tests/test_data/processed/uspex-SXD/qa/",
        },
    }
    return proc_paths
