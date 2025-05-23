import os
import pytest
from pyspextool import config as setup
from pyspextool.setup_utils import (
    pyspextool_setup,
    set_paths,
    set_instrument,
    set_qa_state,
    mishu,
)


def test_pyspextool_setup_defaults(raw_setup):
    spex_prism_paths = raw_setup["spex_prism"]
    pyspextool_setup(
        raw_path=spex_prism_paths["raw_path"],
        cal_path=spex_prism_paths["cal_path"],
        proc_path=spex_prism_paths["proc_path"],
        qa_path=spex_prism_paths["qa_path"],
    )

    assert setup.state["verbose"] is True

    assert setup.state["raw_path"] == os.path.abspath(spex_prism_paths["raw_path"])
    assert setup.state["cal_path"] == os.path.abspath(spex_prism_paths["cal_path"])
    assert setup.state["proc_path"] == os.path.abspath(spex_prism_paths["proc_path"])

    assert setup.state["instrument"] == "uspex"

    assert setup.state["qa_path"] == os.path.abspath(spex_prism_paths["qa_path"])
    assert setup.state["qa_extension"] == ".png"
    assert setup.state["qa_show"] is False
    assert setup.state["qa_write"] is False


@pytest.mark.parametrize(
    "setup_name", ["spex_lxd","spex_prism", "spex_sxd", "uspex_lxd","uspex_prism", "uspex_sxd"]
)
def test_set_paths(setup_name, raw_setup):
    setup_dict = raw_setup[setup_name]
    set_paths(
        setup_dict["raw_path"],
        setup_dict["cal_path"],
        setup_dict["proc_path"],
        setup_dict["qa_path"],
    )

    assert setup.state["raw_path"] == os.path.abspath(setup_dict["raw_path"])
    assert setup.state["cal_path"] == os.path.abspath(setup_dict["cal_path"])
    assert setup.state["proc_path"] == os.path.abspath(setup_dict["proc_path"])


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


# def test_set_instrument_state_bad():
#    with pytest.raises(ValueError):
#        set_instrument("SpeX")
#
#    with pytest.raises(ValueError):
#        set_instrument("not_an_instrument")
#
#    with pytest.raises(TypeError):
#        set_instrument(5)


def test_set_qa_state():
    set_qa_state(False, 1, False, False, ".pdf")

    assert setup.state["qa_show"] is False
    assert setup.state["qa_write"] is False
    assert setup.state["qa_extension"] == ".pdf"


@pytest.mark.parametrize(
    "file",
    [
        "uspex_lincorr.fits",
        "uspex_bias.fits",
        "spex_lincorr.fits",
        "Vega50000.fits",
    ],
)
def test_pooch_cache(file):
    assert mishu.is_available(file)
