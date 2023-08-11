import glob
import os
import pyspextool as ps


ps.pyspextool_setup(raw_path='./',
                    qa_path='tests/test_data/misc_data/',
                    cal_path='tests/test_data/misc_data/',
                    proc_path='tests/test_data/misc_data/',
                    verbose=True, qa_show=False, qa_write=True,
                    qa_extension='.png')


ps.telluric.telluric_correction('spectra1-8.fits', 'HD 101369',
                                'spectra11-18.fits',
                                'test', fluxdensity_units='Jy',
                                input_path='tests/test_data/misc_data/',
                                output_path='tests/test_data/misc_data/')

fits_files = glob.glob(os.path.join("tests/test_data/misc_data/","test.fits"))
png_files = glob.glob(os.path.join("tests/test_data/misc_data/","test.png"))

assert len(fits_files) == 1
assert len(png_files) == 1

# CLEANUP
    # remove generated files
for files in fits_files + png_files:
    os.remove(files)

