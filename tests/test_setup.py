from pyspextool.setup_utils import *
import pytest

spex_prism_paths = {
    'raw_path':  'tests/test_data/spex-prism/data/',
    'cal_path':  'tests/test_data/spex-prism/cals/',
    'proc_path': 'tests/test_data/spex-prism/proc/',
    'qa_path': 'tests/test_data/spex-prism/qa/',
    }


@pytest.mark.parametrize("instrument, paths", [('spex', spex_prism_paths)])
def test_pyspextool_setup(instrument, paths):
    pyspextool_setup(instrument, paths=paths, verbose=True)


def test_set_instrument_state():
    set_instrument('uspex')

    assert setup.state['instrument'] == 'uspex'
    assert setup.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument


