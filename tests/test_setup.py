from pyspextool.cl import config
from pyspextool.cl.setup import *


def test_set_instrument_state():
    set_instrument_state('uspex')

    assert config.state['instrument_name'] == 'uspex'
    assert config.state['readfits'] == 'read_uspex_fits'
    assert config.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument


def test_spex_setup():
    setup(instrument_name='uspex', clupdate=True)