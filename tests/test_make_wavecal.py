import glob
import os
import pyspextool as ps


def test_make_wavecal():

    ps.pyspextool_setup(
        raw_path="tests/test_data/raw/uspex-prism/data/",
        qa_path="tests/test_data/raw/uspex-prism/qa/",
        cal_path="tests/test_data/raw/uspex-prism/cals/",
        proc_path="tests/test_data/raw/uspex-prism/proc/",
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.extract.make_flat(
        ["sbd.2022B046.221019.flat.", "15-19"], "flat15-19", verbose=True
    )

    ps.extract.make_wavecal(
        ["sbd.2022B046.221019.arc.", "20"],
        "flat15-19.fits",
        "wavecal20",
        use_stored_solution=False,
    )

    fits_files = glob.glob(
        os.path.join("tests/test_data/raw/uspex-prism/cals/", "*.fits")
    )

    png_files = glob.glob(os.path.join("tests/test_data/raw/uspex-prism/qa/", "*.png"))

    assert len(fits_files) == 2
    assert len(png_files) == 4

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)
