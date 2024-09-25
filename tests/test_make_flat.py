import glob
import os
import pyspextool as ps

#
# Running the test gives:
#
#  /opt/homebrew/lib/python3.11/site-packages/scipy/stats/_stats_py.py:1069: DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)

# which have tracked down to the moments.py function and the call to stats.  But I cannot figure out what is going on.  Deal with later.
#


def test_make_flat():

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

    ps.extract.make_flat(['sbd.2022B046.221019.flat.','15-19'],'flat15-19',
                          verbose=True)

    fits_files = glob.glob(
        os.path.join("tests/test_data/raw/uspex-prism/cals/", "flat15-19.fits")
    )

    png_files = glob.glob(os.path.join("tests/test_data/raw/uspex-prism/qa/", "flat*.png"))

    assert len(fits_files) == 1
    assert len(png_files) == 2

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)
