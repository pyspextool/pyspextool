from pyspextool.cl.setup import *
from pyspextool.cl import config

def test_set_instrument_state():
    instrument = 'uspex'
    set_instrument_state(instrument)

    assert config.state['instrument'] == 'uspex'
    assert config.state['readfits'] == 'read_uspex_fits'
    assert config.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument

