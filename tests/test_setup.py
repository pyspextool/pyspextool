from pyspextool.extract import config
from pyspextool.extract.setup_utils import *


def test_spex_setup():
    pyspextool_setup(instrument='uspex', verbose=True)


def test_set_instrument_state():
    set_instrument('uspex')

    assert config.state['instrument_name'] == 'uspex'
    #assert config.state['readfits']        == 'read_uspex_fits' # DEPRECATED
    assert config.state['lincormax']       == 35000

    # TODO: add tests for not_supported instrument
    


