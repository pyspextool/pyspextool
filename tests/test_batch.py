import os
import glob
import pandas
import pytest
from pyspextool.batch.batch import (
    readDriver,
    writeDriver,
    writeLog,
    batchReduce,
    processFolder,
)
from pyspextool.batch.batch import (
    HEADER_DATA,
    REDUCTION_FOLDERS,
    BATCH_PARAMETERS,
    LOG_PARAMETERS,
    # OBSERVATION_SET_KEYWORD,
)


basefold = "tests/test_data/raw"

test_instruments = ["spex-prism", "spex-SXD"]
tfold = os.path.join(basefold, "spex-prism")
rawfold = os.path.join(basefold, "raw/")
procfold = os.path.join(basefold, "processed/")


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


def test_writeLog():
    logfile = os.path.join(tfold, "qa", "log_test.csv")

    # process files and generate a test log
    dp = processFolder(os.path.join(tfold, "data/"))
    writeLog(dp, logfile)
    assert os.path.exists(logfile)

    # read in log and test contents
    dp = pandas.read_csv(logfile)
    assert isinstance(dp, pandas.core.frame.DataFrame)
    assert len(dp) > 0
    for x in LOG_PARAMETERS["COLUMNS"]:
        assert x in list(dp.columns)

    # CLEANUP
    # remove generated files
    os.remove(logfile)


def test_writeReadDriver():
	logfile = os.path.join(tfold, "qa", "log_test.csv")
	driver_file = os.path.join(tfold, "proc/driver_test.txt")

# process files and generate a driver file
	dp = processFolder(os.path.join(tfold, "data/"))
	writeLog(dp,logfile)
	dp = pandas.read_csv(logfile)
    writeDriver(
        dp,
        driver_file,
        data_folder=os.path.join(tfold, "data"),
        check=False,
        create_folders=False,
        verbose=False,
    )
    assert os.path.exists(driver_file)

    # read in existing driver in test folder and check
    par = readDriver(driver_file)
    assert isinstance(par, dict)
    for x in BATCH_PARAMETERS:
        assert x in list(par.keys())

# CLEANUP
# remove generated files
	os.remove(logfile)
	os.remove(driver_file)


# def test_makeQApage():
# 	 # write a file and then delete it
# 	 pass


@pytest.mark.skip(reason='issues with organization of prism data currently prevent this test')
def test_batchReduce_spex_prism():
    base_folder = "./tests/test_data/spex-prism/"

    # first, make sure proc and cal folders _do not have_ FITS files
    proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
    cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
    qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
    assert len(proc_fits_files) == 0
    assert len(cal_fits_files) == 0
    assert len(qa_pdf_files) == 0

    # generate driver file
    driver_file = os.path.join(base_folder, "proc/driver_test.txt")
    # dp = processFolder(os.path.join(base_folder, "data"),verbose=False)
    # writeDriver(dp,driver_file,data_folder=os.path.join(base_folder, "data"),check=False,create_folders=False,verbose=False)

    par = readDriver(driver_file)
    batchReduce(par,verbose=True)

    # now, make sure proc and qa folders _have_ FITS files and qa has no PDF files
    proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
    cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
    qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
    assert len(proc_fits_files) == 17
    assert len(cal_fits_files) == 2
    assert len(qa_pdf_files) == 55

    # CLEANUP
    # remove generated files
    for files in proc_fits_files + cal_fits_files + qa_pdf_files:
        os.remove(files)

# def test_batchReduce_uspex_prism():
# 	base_folder = "./tests/test_data/uspex-prism/"

# 	# first, make sure proc and cal folders _do not have_ FITS files
# 	proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
# 	cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
# 	qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
# 	assert len(proc_fits_files) == 0
# 	assert len(cal_fits_files) == 0
# 	assert len(qa_pdf_files) == 0

# # generate driver file
# 	driver_file = os.path.join(base_folder, "proc/driver_test.txt")
# 	# dp = processFolder(os.path.join(base_folder, "data"),verbose=False)
# 	# writeDriver(dp,driver_file,data_folder=os.path.join(base_folder, "data"),check=False,create_folders=False,verbose=False)

# 	par = readDriver(driver_file)
# 	batchReduce(par,verbose=True)

# 	# now, make sure proc and qa folders _have_ FITS files and qa has no PDF files
# 	proc_fits_files = glob.glob(os.path.join(base_folder, "proc/*.fits"))
# 	cal_fits_files = glob.glob(os.path.join(base_folder, "cals/*.fits"))
# 	qa_pdf_files = glob.glob(os.path.join(base_folder, "qa/*.pdf"))
# 	assert len(proc_fits_files) == 7
# 	assert len(cal_fits_files) == 2
# 	assert len(qa_pdf_files) == 25

# 	# CLEANUP
# 	# remove generated files
# 	for files in proc_fits_files + cal_fits_files + qa_pdf_files:
# 		os.remove(files)
