import glob
import os
import pyspextool as ps

basefold = "tests/test_data/"
rawfold = os.path.join(basefold, "raw/")
procfold = os.path.join(basefold, "processed/")
test_instruments = ["spex-prism", "spex-SXD","uspex-prism", "uspex-SXD"]

def test_make_wavecal():

    inst = 'uspex-prism'
    ps.pyspextool_setup(
        raw_path=os.path.join(rawfold,inst,'data'),
        cal_path=os.path.join(procfold,inst,'cals'),
        proc_path=os.path.join(procfold,inst,'proc'),
        qa_path=os.path.join(procfold,inst,'qa'),
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    # expecting to find flat15-19.fits in the processed cals directory

    ps.extract.make_wavecal(
        ["sbd.2022B046.221019.arc.", "20"],
        "flat15-19.fits",
        "wavecal20",
        use_stored_solution=False,
    )

    fits_files = glob.glob(os.path.join(procfold,inst,'cals',"wavecal*.fits"))

    png_files = glob.glob(os.path.join(procfold,inst,'qa',"wavecal*.png"))    

    assert len(fits_files) == 1
    assert len(png_files) == 2
    
    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)
