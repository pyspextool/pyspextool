from pyspextool import config as setup
from pyspextool.setup_utils import *


def test_pyspextool_setup():
    set_parameters(verbose=True)


def test_set_instrument_state():
    set_instrument('uspex')

    assert setup.state['instrument'] == 'uspex'
    assert setup.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument


