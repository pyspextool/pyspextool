import numpy as np
from pyspextool.extract.override_aperturesigns import override_aperturesigns
from pyspextool.setup_utils import pyspextool_setup
from pyspextool.extract import config as extract

def test_override_aperturesigns(raw_setup):

    setup_dict = raw_setup['spex_prism']

    pyspextool_setup(
        instrument=setup_dict["instrument"],
        raw_path=setup_dict["raw_path"],
        cal_path=setup_dict["cal_path"],
        proc_path=setup_dict["proc_path"],
        qa_path=setup_dict["qa_path"],
    )

    extract.state['average_aperturesigns'] = [1,-1]
    
    result = override_aperturesigns('+,+')

    np.testing.assert_array_equal(result,np.array([1,1]))

