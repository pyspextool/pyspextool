import pyspextool as ps
from pyspextool.batch import batch
import os
import pandas as pd
import logging


testdatafold = os.path.join(batch.CODEDIR,"../../tests/test_data/")
if os.path.exists(testdatafold)==False: testdatafold = os.path.join(batch.CODEDIR,"tests","test_data/")
if os.path.exists(testdatafold)==False: 
    raise ValueError('Cannot find test data in {} or {}'.format(os.path.join(batch.CODEDIR,"../../tests/test_data/"),os.path.join(batch.CODEDIR,"tests","test_data/")))

rawfold = os.path.join(testdatafold, "raw/")
procfold = os.path.join(testdatafold, "processed/")


def test_batch_log(testdatafold=testdatafold,rawfold=rawfold,procfold=procfold,verbose=batch.ERROR_CHECKING):
    '''
    Purpose
    -------
    Checks that logs can be made and read

    ''' 
# process files and generate a test log
    test_instruments = ['spex-prism','spex-LXD','spex-SXD','uspex-prism','uspex-LXD','uspex-SXD']
    log_test_file = 'log_test.csv'
    driver_test_file = 'driver_test.txt'
    if verbose==True: logging.info('...checking that code can generate log')
    for inst in test_instruments:
        rfold = os.path.join(rawfold, inst)
        pfold = os.path.join(procfold, inst)
        logfile = os.path.join(pfold, "qa", log_test_file)
        print(inst,rfold,pfold,logfile)
        dp = batch.processFolder(os.path.join(rfold, "data/"),verbose=False)
        batch.writeLog(dp,logfile,verbose=False)
        assert os.path.exists(logfile), 'Could not find log file {} after creation'.format(logfile)
        dp = pd.read_csv(logfile)
        assert isinstance(dp, pd.core.frame.DataFrame), 'Log file parameter error: not a pandas dataframe for log file {}'.format(logfile)
        assert len(dp)>0
        for x in batch.LOG_PARAMETERS['COLUMNS']:
            assert x in list(dp.columns), 'Could not find required column {} in log file {}'.format(x,logfile)
        for x in batch.SIMBAD_COLS:
            assert x in list(dp.columns), 'Could not find required SIMBAD column {} in log file {}'.format(x,logfile)
        if verbose==True: logging.info('\t{}: PASS'.format(inst))

# CLEANUP
# remove generated files
        os.remove(logfile)

def test_batch_driver(testdatafold=testdatafold,rawfold=rawfold,procfold=procfold,verbose=batch.ERROR_CHECKING):
    '''
    Purpose
    -------
    Checks that driver file can be made and read

    ''' 
# process files and generate a test log
    test_instruments = ['spex-prism','spex-LXD','spex-SXD','uspex-prism','uspex-LXD','uspex-SXD']
    log_test_file = 'log_test.csv'
    driver_test_file = 'driver_test.txt'
    if verbose==True: logging.info('...checking that code can generate log')
    for inst in test_instruments:
        rfold = os.path.join(rawfold, inst)
        pfold = os.path.join(procfold, inst)
# make log firs
        logfile = os.path.join(pfold, "qa", log_test_file)
        driverfile = os.path.join(pfold, "proc", driver_test_file)
        dp = batch.processFolder(os.path.join(rfold, "data/"),verbose=False)
        batch.writeLog(dp,logfile,verbose=False)
# read log and make driver file
        dp = pd.read_csv(logfile)
        batch.writeDriver(dp,driverfile,data_folder=os.path.join(rfold, "data"),check=False,create_folders=False,verbose=False)
        assert os.path.exists(driverfile), 'Could not find driver file {} after creation'.format(driverfile)
# check driver file parameters
        par = batch.readDriver(driverfile,verbose=False)
        assert isinstance(par, dict), 'Driver file parameter error: not a dictionary for driver file {}'.format(driverfile)
        for x in batch.BATCH_PARAMETERS:
            assert x in list(par.keys()), 'Could not find required parameter {} in driver file {}'.format(x,driverfile)
        if verbose==True: logging.info('\t{}: PASS'.format(inst))

# CLEANUP
# remove generated files
        os.remove(logfile)
        os.remove(driverfile)

    return

def test_batch_fullpass(testdatafold=testdatafold,rawfold=rawfold,procfold=procfold,verbose=batch.ERROR_CHECKING):
    '''
    Purpose
    -------
    runs a full pass reduction

    ''' 
    test_instruments = ['spex-prism']
    log_test_file = 'log_test.csv'
    driver_test_file = 'driver_test.txt'
    pass

