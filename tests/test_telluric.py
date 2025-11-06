#from pyspextool.telluric.load_spectra import _load_vegamodel
from pyspextool.setup_utils import pyspextool_setup
from pyspextool.pyspextoolerror import pySpextoolError
from pyspextool.telluric import telluric
from pyspextool.combine import combine
import pytest
import os
import glob

@pytest.mark.parametrize("setup_name", ["uspex_sxd"])
def test_telluric(setup_name, proc_setup):
    """Test telluric correction with uspex SXD data (existing test)"""
    setup_dict = proc_setup[setup_name]
    proc_path = setup_dict["proc_path"]
    qa_extension = '.png'
    qa_showblock = False  # Changed to False to avoid blocking
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


@pytest.mark.parametrize("setup_name", [
    "spex_prism",      # Minimum requirement - crucial
    "spex_sxd",        # Optional
    "spex_lxd",        # Preferred
    "uspex_prism",     # Optional
    "uspex_lxd_short", # Optional
])
def test_telluric_postextraction(setup_name, postextraction_setup):
    """
    Test telluric correction for various instrument modes using PostExtractionTests data.
    
    This test:
    1. Sets up the processing path
    2. Combines object spectra
    3. Combines standard spectra  
    4. Runs telluric correction
    5. Verifies telluric correction file is created
    6. Cleans up generated files
    """
    setup_dict = postextraction_setup[setup_name]
    proc_path = setup_dict["proc_path"]
    
    # Setup pyspextool
    pyspextool_setup(
        proc_path=proc_path,
        qa_extension='.png',
        qa_showblock=False,
        verbose=False
    )
    
    # Combine object spectra
    object_combined = f'cspectra_object_{setup_name}.fits'
    combine(
        files=setup_dict["object_files"],
        output_name=object_combined.replace('.fits', ''),
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Combine standard spectra
    standard_combined = f'cspectra_standard_{setup_name}.fits'
    combine(
        files=setup_dict["standard_files"],
        output_name=standard_combined.replace('.fits', ''),
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Run telluric correction
    telluric(
        object_filenames=object_combined,
        standard_filename=standard_combined,
        standard_info=setup_dict["standard_name"],
        telluric_filename='telluric',
        corrected_filenames='tcspectra_object',
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Verify telluric correction file was created
    telluric_path = os.path.join(proc_path, 'telluric.fits')
    assert os.path.exists(telluric_path), f"Telluric correction file not created for {setup_name}"
    
    # CLEANUP - remove generated files
    os.remove(telluric_path)
    
    # Remove combined spectra
    for pattern in [f'cspectra_object_{setup_name}.fits', 
                     f'cspectra_standard_{setup_name}.fits',
                     'tcspectra_object*.fits']:
        for filepath in glob.glob(os.path.join(proc_path, pattern)):
            os.remove(filepath)


@pytest.mark.parametrize("setup_name", ["spex_prism"])
def test_telluric_precomputed(setup_name, postextraction_setup):
    """
    Test applying a pre-computed telluric correction file to individual spectra.
    
    This test:
    1. Creates a telluric correction file from combined spectra
    2. Applies that correction to individual uncorrected spectra
    3. Verifies the corrected individual spectra are created
    4. Cleans up generated files
    
    This tests the workflow where a telluric correction is computed once
    and then applied to multiple individual observations.
    """
    setup_dict = postextraction_setup[setup_name]
    proc_path = setup_dict["proc_path"]
    
    # Setup pyspextool
    pyspextool_setup(
        proc_path=proc_path,
        qa_extension='.png',
        qa_showblock=False,
        verbose=False
    )
    
    # Step 1: Create telluric correction file from combined spectra
    # Combine object spectra
    object_combined = f'cspectra_object_{setup_name}_precomp.fits'
    combine(
        files=setup_dict["object_files"],
        output_name=object_combined.replace('.fits', ''),
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Combine standard spectra
    standard_combined = f'cspectra_standard_{setup_name}_precomp.fits'
    combine(
        files=setup_dict["standard_files"],
        output_name=standard_combined.replace('.fits', ''),
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Create telluric correction file
    telluric(
        object_filenames=object_combined,
        standard_filename=standard_combined,
        standard_info=setup_dict["standard_name"],
        telluric_filename='telluric_precomp',
        corrected_filenames='tcspectra_combined',
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    telluric_path = os.path.join(proc_path, 'telluric_precomp.fits')
    assert os.path.exists(telluric_path), "Pre-computed telluric correction file not created"
    
    # Step 2: Apply pre-computed telluric correction to individual spectra
    # Use the first few individual object spectra
    telluric(
        object_filenames=setup_dict["object_files"],
        standard_filename=standard_combined,
        standard_info=setup_dict["standard_name"],
        telluric_filename='telluric_precomp',  # Reuse the pre-computed correction
        corrected_filenames='tcspectra_individual',
        verbose=False,
        qa_show=False,
        qa_write=False
    )
    
    # Verify corrected individual spectra were created
    corrected_files = glob.glob(os.path.join(proc_path, 'tcspectra_individual*.fits'))
    assert len(corrected_files) > 0, "No corrected individual spectra were created"
    
    # CLEANUP - remove generated files
    os.remove(telluric_path)
    os.remove(os.path.join(proc_path, object_combined))
    os.remove(os.path.join(proc_path, standard_combined))
    
    # Remove corrected spectra
    for filepath in glob.glob(os.path.join(proc_path, 'tcspectra_combined*.fits')):
        os.remove(filepath)
    for filepath in corrected_files:
        os.remove(filepath)


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

