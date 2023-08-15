import os
import glob
import pandas
from pyspextool.batch.batch import readDriver, batchReduce, processFolder
from pyspextool.batch.batch import HEADER_DATA, REDUCTION_FOLDERS


basefold = "tests/test_data/"

test_instruments = ["spex-prism", "spex-SXD"]
tfold = os.path.join(basefold, "spex-prism")


def test_inputs():
    # make sure code can find test data
    assert os.path.exists(basefold)
    for inst in test_instruments:
        assert os.path.exists(os.path.join(basefold, inst))
        for fold in REDUCTION_FOLDERS:
            assert os.path.exists(os.path.join(basefold, inst, fold))
    dfiles = glob.glob(os.path.join(tfold, "data/*.fits"))
    assert len(dfiles) > 0


def test_processFolder():
    result = processFolder(os.path.join(tfold, "data/"))
    assert isinstance(result, pandas.core.frame.DataFrame)
    assert len(result) > 0
    for k in list(HEADER_DATA.keys()):
        assert k in list(result.columns)


# def test_organizeLegacy():
#     pass


# def test_writeLog():
#     # write a file and then delete it
#     pass


# def test_readDriver():
#     # read in existing driver in test folder
#     pass


# def test_writeDriver():
#     # write a file and then delete it
#     pass


# def test_makeQApage():
#     # write a file and then delete it
#     pass


def test_batchReduce_spex_prism():
    base_folder = "./tests/test_data/spex-prism/"

    # first, make sure proc and cal folders _do not have_ FITS files
    proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
    cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
    qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
    assert len(proc_fits_files) == 0
    assert len(cal_fits_files) == 0
    assert len(qa_pdf_files) == 0

    driver_file = os.path.join(base_folder, "proc/driver_test.txt")
    par = readDriver(driver_file)
    batchReduce(par)

    # now, make sure proc and qa folders _have_ FITS files and qa has no PDF files
    proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
    cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
    qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
    assert len(proc_fits_files) == 17
    assert len(cal_fits_files) == 2
    assert len(qa_pdf_files) == 60

    # CLEANUP
    # remove generated files
    for files in proc_fits_files + cal_fits_files + qa_pdf_files:
        os.remove(files)
