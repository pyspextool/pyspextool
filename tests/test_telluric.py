#from pyspextool.telluric.load_spectra import _load_vegamodel
from pyspextool.setup_utils import pyspextool_setup
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric import telluric
import pytest
import os
import glob

@pytest.mark.parametrize("setup_name", ["uspex_sxd"])
def test_telluric(setup_name, proc_setup):
    setup_dict = proc_setup[setup_name]
    proc_path = setup_dict["proc_path"]
    qa_extension = '.png'
    qa_showblock = True
    pyspextool_setup(
                    proc_path=proc_path, 
                    qa_extension=qa_extension,
                    qa_showblock=qa_showblock,
                    verbose=False)
    telluric(['spectra','1-6'] ,'cspectra9-16.fits', 'HD17778', \
                     'telluric', 'tcspectra')

    telluric_path = os.path.join(proc_path, 'telluric.fits')
    assert os.path.exists(telluric_path)

    # CLEANUP
    # remove generated files
    os.remove(telluric_path)
    for files in glob.glob(os.path.join(proc_path, 'tcspectra*.fits')):
        os.remove(files)


#    [
#        ("spex_lxd", False, "Vega50000.fits"),
#        ("spex_lxd", True, "Vega50000_new.fits"),
#        ("spex_prism", False, "Vega50000.fits"),
#        ("spex_prism", True, "Vega50000_new.fits"),
#        ("spex_sxd", False, "Vega50000.fits"),
#        ("spex_sxd", True, "Vega50000_new.fits"),
#        ("uspex_lxd", False, "Vega50000.fits"),
#        ("uspex_lxd", True, "Vega50000_new.fits"),
#        ("uspex_prism", False, "Vega50000.fits"),
#        ("uspex_prism", True, "Vega50000_new.fits"),
#        ("uspex_sxd",False,"Vega50000.fits"),
#        ("uspex_sxd",True,"Vega50000_new.fits"),
#    ],
#)
#def test_load_vegamodel(setup_name, new, vega_file, proc_setup):
#    setup_dict = proc_setup[setup_name]
#    pyspextool_setup(
#        instrument=setup_dict["instrument"],
#        raw_path=setup_dict["raw_path"],
#        cal_path=setup_dict["cal_path"],
#        proc_path=setup_dict["proc_path"],
#        qa_path=setup_dict["qa_path"],
#    )
#    result = _load_vegamodel(setup_dict["standard_file"], new=new)
#    
#    
#    assert result["vega_file"] == vega_file

