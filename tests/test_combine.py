import glob
import os
import pyspextool as ps


def test_combine():

    #
    # uSpeX SXD
    #
    
    proc_path = "tests/test_data/processed/uspex-SXD/proc"
    qa_path = "tests/test_data/processed/uspex-SXD/qa"

    ps.pyspextool_setup("uspex",
        qa_path=qa_path,
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.combine.combine(["spectra", "1-8"], "cspectra1-8", statistic="mean")
    ps.combine.combine(["spectra", "11-16"], "cspectra11-16", statistic="mean")
    

    fits_files = glob.glob(os.path.join(proc_path, "cspectra*.fits"))
    png_files = glob.glob(os.path.join(qa_path, "*.png"))

    assert len(fits_files) == 2
    assert len(png_files) == 6

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)

    #
    # uSpeX prism
    #
    
    proc_path = "tests/test_data/processed/uspex-prism/proc"
    qa_path = "tests/test_data/processed/uspex-prism/qa"

    ps.pyspextool_setup("uspex",
        qa_path=qa_path,
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.combine.combine(["spectra", "1-2"], "cspectra1-2")
    ps.combine.combine(["spectra", "7-8"], "cspectra7-8")
    

    fits_files = glob.glob(os.path.join(proc_path, "cspectra*.fits"))
    png_files = glob.glob(os.path.join(qa_path, "*.png"))

    assert len(fits_files) == 2
    assert len(png_files) == 6

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)



    #
    # SpeX SXD
    #
    
    proc_path = "tests/test_data/processed/spex-SXD/proc"
    qa_path = "tests/test_data/processed/spex-SXD/qa"

    ps.pyspextool_setup("spex",
        qa_path=qa_path,
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.combine.combine(["spectra", "626-635"], "cspectra626-635")
    ps.combine.combine(["spectra", "636-645"], "cspectra636-645")
    

    fits_files = glob.glob(os.path.join(proc_path, "cspectra*.fits"))
    png_files = glob.glob(os.path.join(qa_path, "*.png"))

    assert len(fits_files) == 2
    assert len(png_files) == 6

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)


    #
    # SpeX prism
    #
    
    proc_path = "tests/test_data/processed/spex-prism/proc"
    qa_path = "tests/test_data/processed/spex-prism/qa"

    ps.pyspextool_setup("spex",
        qa_path=qa_path,
        proc_path=proc_path,
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    ps.combine.combine(["spectra", "1-7"], "cspectra1-7")
    ps.combine.combine(["spectra", "8-14"], "cspectra8-14")
    

    fits_files = glob.glob(os.path.join(proc_path, "cspectra*.fits"))
    png_files = glob.glob(os.path.join(qa_path, "*.png"))

    assert len(fits_files) == 2
    assert len(png_files) == 6

    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)

        
