from pyspextool.setup_utils import *
import pytest


cwd = os.path.abspath(os.getcwd())

spex_prism_paths = {
    'raw_path':  'tests/test_data/spex-prism/data/',
    'cal_path':  'tests/test_data/spex-prism/cals/',
    'proc_path': 'tests/test_data/spex-prism/proc/',
    'qa_path': 'tests/test_data/spex-prism/qa/',
    }

spex_sxd_paths = {
    'raw_path':  'tests/test_data/spex-sxd/data/',
    'cal_path':  'tests/test_data/spex-sxd/cals/',
    'proc_path': 'tests/test_data/spex-sxd/proc/',
    'qa_path': 'tests/test_data/spex-sxd/qa/',
    }

uspex_sxd_paths = {
    'raw_path':  'tests/test_data/uspex-sxd/data/',
    'cal_path':  'tests/test_data/uspex-sxd/cals/',
    'proc_path': 'tests/test_data/uspex-sxd/proc/',
    'qa_path': 'tests/test_data/uspex-sxd/qa/',
    }

uspex_prism_paths = {
    'raw_path':  'tests/test_data/uspex-prism/data/',
    'cal_path':  'tests/test_data/uspex-prism/cals/',
    'proc_path': 'tests/test_data/uspex-prism/proc/',
    'qa_path': 'tests/test_data/uspex-prism/qa/',
    }


@pytest.mark.parametrize("instrument, paths", [('spex', spex_prism_paths)])
def test_pyspextool_setup(instrument, paths):
    pyspextool_setup(instrument, paths=paths, verbose=True)


@pytest.mark.parametrize("paths", [(spex_prism_paths)]) #, spex_sxd_paths, uspex_prism_paths, uspex_sxd_paths)])
def test_set_paths(paths):
    # TODO: add tests for bad paths
    state = set_paths(paths)

    #assert state['paths'] == paths


@pytest.mark.parametrize("instrument", [('spex')]) #, uspex
def test_set_instrument_state(instrument):
    set_instrument(instrument)

    if instrument == 'uspex':
        assert setup.state['instrument'] == 'uspex'
        assert setup.state['lincormax'] == 35000
    elif instrument == 'spex':
        assert setup.state['instrument'] == 'spex'
        assert setup.state['lincormax'] == 4000

    # TODO: add tests for not_supported instrument


def test_set_qa_state():
    state = set_qa_state()

    assert state['qa_path'] == cwd
    assert state['qa_plot'] == False
    assert state['qa_extension'] == '.pdf'
    assert state['qa_file'] == True