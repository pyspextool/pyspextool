import glob
import os
import pyspextool as ps


def test_telluric():

    raw_path = "tests/test_data/raw/spex-prism/data"
    proc_path = "tests/test_data/processed/spex-prism/proc"
    qa_path = "tests/test_data/processed/spex-prism/qa"
    cal_path = "tests/test_data/processed/spex-prism/cals"

    ps.pyspextool_setup(
        raw_path=raw_path,
        qa_path=qa_path,
        cal_path=cal_path,
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.telluric.telluric_correction(
        "combspec1-8.fits",
        "HD 101060",
        "combspec9-14.fits",
        "test",
        fluxdensity_units="Jy",
        input_path=proc_path,
        output_path=proc_path,
    )

    fits_files = glob.glob(os.path.join(proc_path, "test.fits"))
    png_files = glob.glob(os.path.join(qa_path, "test.png"))

    assert len(fits_files) == 1
    assert len(png_files) == 1

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)
