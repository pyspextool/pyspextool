from pyspextool.setup_utils import *

# This should be a dictionary

base_folder = 'tests/test_data/spex-prism'
# what instrument and mode we are using
instrument = 'spex'
mode = 'prism'

# set output pathways
# set output pathways
raw_path = base_folder+'/data/'
cal_path = base_folder+'/cals/'
proc_path = base_folder+'/proc/'
qa_path = base_folder+'/qa/'


def test_pyspextool_setup():
    pyspextool_setup(instrument,raw_path=raw_path, cal_path=cal_path, proc_path=proc_path, qa_path=qa_path,verbose=True)



def test_set_instrument_state():
    set_instrument('uspex')

    assert setup.state['instrument'] == 'uspex'
    assert setup.state['lincormax'] == 35000

    # TODO: add tests for not_supported instrument


