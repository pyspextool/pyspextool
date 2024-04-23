import glob
import os
import pyspextool as ps


def test_combine():

    proc_path = "tests/test_data/processed/spex-prism/proc"
    qa_path = "tests/test_data/processed/spex-prism/qa"

    ps.pyspextool_setup(
        instrument="spex",
        raw_path="tests/test_data/raw/spex-prism/data",
        qa_path=qa_path,
        cal_path="tests/test_data/processed/spex-prism/cals",
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.combine.combine_spectra(["spectra", "1-7"], "spectra1-7", statistic="mean")

    fits_files = glob.glob(os.path.join(proc_path, "spectra1-7.fits"))
    png_files = glob.glob(os.path.join(qa_path, "*.png"))

    assert len(fits_files) == 1
    assert len(png_files) == 3

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)
