import glob
import os
import pyspextool as ps


def test_combine():

    ps.pyspextool_setup(raw_path='./',
                        qa_path='tests/test_data/misc_data/',
                        cal_path='tests/test_data/misc_data/',
                        proc_path='tests/test_data/misc_data/',
                        verbose=True, qa_show=False, qa_write=True,
                        qa_extension='.png')

    ps.combine.combine_spectra(['spectra','1-7'], 'spectra1-7',
                                statistic='mean')
    

    fits_files = glob.glob(os.path.join("tests/test_data/misc_data/","spectra1-7.fits"))
    png_files = glob.glob(os.path.join("tests/test_data/misc_data/","*.png"))

    assert len(fits_files) == 1
    assert len(png_files) == 3
    
    # CLEANUP
    # remove generated files
    for files in fits_files + png_files:
        os.remove(files)

