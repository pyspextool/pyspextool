import importlib
import numpy as np
import pyspextool as ps
from pyspextool import config as setup
from pyspextool.extract.flat import read_flatcal_file, locate_orders


def test_locate_orders():

    # Get set up

    ps.pyspextool_setup(
        raw_path="tests/test_data/raw/uspex-SXD/data/",
        qa_path="tests/test_data/raw/uspex-SXD/qa/",
        cal_path="tests/test_data/raw/uspex-SXD/cals/",
        proc_path="tests/test_data/raw/uspex-SXD/proc/",
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    # Get flat info

    # Load a single flat frame
    # This flat is missing the TCS-AM keyword and didn't work with previous version
    # of locate_orders. https://github.com/pyspextool/pyspextool/issues/221
    file = setup.state["raw_path"] + "/sbd.2023B006.231021.flat.00293.a.fits"

    module = (
        "pyspextool.instruments."
        + setup.state["instrument"]
        + "."
        + setup.state["instrument"]
    )

    instr = importlib.import_module(module)

    data = instr.read_fits([file], setup.state["linearity_info"])

    hdr = data[2]

    mode = hdr[0]["MODE"][0]
    modefile = setup.state["instrument_path"] + "/" + mode + "_flatinfo.fits"

    modeinfo = read_flatcal_file(modefile)

    result = locate_orders(
        data[0],
        modeinfo["guesspos"],
        modeinfo["xranges"],
        modeinfo["step"],
        modeinfo["slith_range"],
        modeinfo["edgedeg"],
        modeinfo["ybuffer"],
        modeinfo["flatfrac"],
        modeinfo["comwidth"],
        qa_show=False,
    )

    coeffs = result[0]
    assert np.isclose(coeffs[0][0][0], 368.719)
    assert np.isclose(coeffs[1][0][0], 616.042)
    assert np.isclose(coeffs[2][0][0], 861.734)
    assert np.isclose(coeffs[3][0][0], 1146.96056)
    assert np.isclose(coeffs[4][0][0], 1495.087)
    assert np.isclose(coeffs[5][0][0], 1928.40968)
    assert np.isclose(coeffs[6][0][0], 2468.47831)

    xranges = result[1]
    assert (xranges[0] == (4, 2043)).all()
    assert (xranges[1] == (4, 2043)).all()
    assert (xranges[2] == (4, 2043)).all()
    assert (xranges[3] == (4, 2043)).all()
    assert (xranges[4] == (4, 2043)).all()
    assert (xranges[5] == (200, 2043)).all()
    assert (xranges[6] == (660, 2043)).all()
