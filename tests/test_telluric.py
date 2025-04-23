from pyspextool.telluric.load_spectra import load_data, load_modeinfo, load_vegamodel
from pyspextool.telluric import config as tc
from pyspextool.setup_utils import pyspextool_setup
import pytest


@pytest.mark.parametrize(
    ("instrument_name", "mode", "path_name"),
    [
        ("uspex", "ShortXD", "uspex_sxd"),
        ("uspex", "Prism", "uspex_prism"),
        #("uspex", "LXD_short", ""),
        #("uspex", "LXD_long", ""),
        ("spex", "ShortXD", "spex_sxd"),
        ("spex", "LowRes15", "spex_prism"),
        #("spex", "LongXD1.9", ""),
    ],
)
def test_load_vegamodel(instrument_name, mode, path_name, paths):
    path = paths[path_name]
    pyspextool_setup(
        instrument=instrument_name,
        raw_path=path["raw_path"],
        cal_path=path["cal_path"],
        proc_path=path["proc_path"],
        qa_path=path["qa_path"],
    )
    tc.state["mode"] = mode
    load_data()
    load_modeinfo()
    load_vegamodel()
    print(f"model:{tc.state["model"]}")
    print(f"method:{tc.state["method"]}")
