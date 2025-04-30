from pyspextool.telluric.load_spectra import _load_vegamodel
from pyspextool.setup_utils import pyspextool_setup
import pytest


@pytest.mark.parametrize(
    ("instrument_name", "mode", "path_name", "standard_file", "new", "vega_file"),
    [
        ("uspex", "ShortXD", "uspex_sxd","spectra00011.fits",False,"Vega50000.fits"),
        ("uspex", "ShortXD", "uspex_sxd","spectra00011.fits",True,"Vega50000_new.fits"),
        #("uspex", "Prism", "uspex_prism",""),
        #("spex", "ShortXD", "spex_sxd",""),
        #("spex", "LowRes15", "spex_prism",""),
    ],
)
def test_load_vegamodel(instrument_name, mode, path_name, standard_file, new, vega_file, proc_paths):
    path = proc_paths[path_name]
    pyspextool_setup(
        instrument=instrument_name,
        raw_path=path["raw_path"],
        cal_path=path["cal_path"],
        proc_path=path["proc_path"],
        qa_path=path["qa_path"],
    )
    result = _load_vegamodel(standard_file, new=new)
    assert result["vega_file"] == vega_file