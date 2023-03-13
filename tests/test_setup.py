from pyspextool.extract import config
from pyspextool.setup_utils import *


def test_pyspextool_setup():
    set_paths(verbose=True)


def test_set_instrument_state():
    set_instrument('uspex')

#    assert config.state['instrument_name'] == 'uspex' # this is not consistent with config structure
    assert config.user['setup']['instrument'] == 'uspex'
    assert config.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument


