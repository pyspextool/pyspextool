import os
import pytest
from pyspextool import config as setup
from pyspextool.setup_utils import (
    pyspextool_setup,
    set_paths,
    set_instrument,
    set_qa_state,
)


cwd = os.path.abspath(os.getcwd())

spex_prism_paths = {
    "raw_path": "tests/test_data/raw/spex-prism/data/",
    "cal_path": "tests/test_data/raw/spex-prism/cals/",
    "proc_path": "tests/test_data/raw/spex-prism/proc/",
    "qa_path": "tests/test_data/raw/spex-prism/qa/",
}

spex_sxd_paths = {
    "raw_path": "tests/test_data/raw/spex-SXD/data/",
    "cal_path": "tests/test_data/raw/spex-SXD/cals/",
    "proc_path": "tests/test_data/raw/spex-SXD/proc/",
    "qa_path": "tests/test_data/raw/spex-SXD/qa/",
}

uspex_sxd_paths = {
    "raw_path": "tests/test_data/raw/uspex-SXD/data/",
    "cal_path": "tests/test_data/raw/uspex-SXD/cals/",
    "proc_path": "tests/test_data/raw/uspex-SXD/proc/",
    "qa_path": "tests/test_data/raw/uspex-SXD/qa/",
}

uspex_prism_paths = {
    "raw_path": "tests/test_data/raw/uspex-prism/data/",
    "cal_path": "tests/test_data/raw/uspex-prism/cals/",
    "proc_path": "tests/test_data/raw/uspex-prism/proc/",
    "qa_path": "tests/test_data/raw/uspex-prism/qa/",
}


def test_pyspextool_setup_defaults():
    pyspextool_setup()

    assert setup.state["verbose"] == False

    assert setup.state["raw_path"] == cwd
    assert setup.state["cal_path"] == cwd
    assert setup.state["proc_path"] == cwd

    assert setup.state["instrument"] == "uspex"

    assert setup.state["qa_path"] == cwd
    assert setup.state["qa_extension"] == ".pdf"
    assert setup.state["qa_show"] == False
    assert setup.state["qa_write"] == True


@pytest.mark.parametrize(
    "paths", [spex_prism_paths, spex_sxd_paths, uspex_prism_paths, uspex_sxd_paths]
)
def test_set_paths(paths):
    set_paths(paths)

    assert setup.state["raw_path"] == paths["raw_path"]
    assert setup.state["cal_path"] == paths["cal_path"]
    assert setup.state["proc_path"] == paths["proc_path"]


#  # TODO: add tests for bad paths
# def test_set_paths_bad():


@pytest.mark.parametrize("instrument", ["spex", "uspex"])
def test_set_instrument(instrument):
    set_instrument(instrument)

    if instrument == "uspex":
        assert setup.state["instrument"] == "uspex"
        assert setup.state["lincormax"] == 35000
    elif instrument == "spex":
        assert setup.state["instrument"] == "spex"
        assert setup.state["lincormax"] == 4000


def test_set_instrument_state_bad():
    with pytest.raises(ValueError):
        set_instrument("not_an_instrument")

    with pytest.raises(TypeError):
        set_instrument(5)


def test_set_qa_state():
    set_qa_state()

    assert setup.state["qa_show"] == False
    assert setup.state["qa_write"] == True
    assert setup.state["qa_path"] == cwd
    assert setup.state["qa_extension"] == ".pdf"
