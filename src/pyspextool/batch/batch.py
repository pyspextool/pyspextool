# Functions for batch reduction using pyspextool
# General procedure:
# - read in files in a given folder and process into a log file and batch driving file
#  -- subprocesses: files --> database, database --> log, database --> batch driver
# - user checks logs and batch driving file to ensure they are correct
# - run through reduction steps based on batch driving file
#  -- subprocesses: read driver --> do all calibrations --> extract all spectra 
#		--> combine spectral groups --> telluric correction --> stitch (SXD/LXD)
#
# REMAINING TASKS
# - xtellcor integration
# - stitching integration
# - finish docstrings
# - integrate testing scripts
# - TESTING

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table
from astroquery.xmatch import XMatch
from astroquery.simbad import Simbad
from astropy.utils.exceptions import AstropyWarning
import astropy.units as u
import copy
import os
import glob
import pandas as pd
pd.set_option('mode.chained_assignment',None) # suppress pandas warnings
import numpy as np
import re
import shutil
import sys
import yaml
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

import pyspextool as ps
from pyspextool import config as setup
from pyspextool.io.files import extract_filestring,make_full_path
from pyspextool.utils.arrays import numberList

VERSION = '2024 Aug 5'

ERROR_CHECKING = True
DIR = os.path.dirname(os.path.abspath(__file__))
CODEDIR = ps.__path__[0]

XMatch.TIMEOUT = 180 # time out for XMatch search
Simbad.TIMEOUT = 60
MOVING_MAXSEP = 15 # max separation for fixed target in arcsec
MOVING_MAXRATE = 10 # max moving rate for fixed target in arcsec/hr
INSTRUMENT_DATE_SHIFT = Time('2014-07-01').mjd # date spex --> uspex
ARC_NAME = 'arc lamp'
FLAT_NAME = 'flat lamp'
OBSERVATION_SET_KEYWORD = 'OBS_SET'
DRIVER_DEFAULT_FILENAME = 'driver.txt'
BACKUP_INSTRUMENT_DATA_URL = 'https://drive.google.com/drive/folders/1hxWS9a4G41qh9d_tk_eAx2aOmJrKvPe6?usp=sharing'
BACKUP_REDUCTION_DATA_URL = 'https://drive.google.com/drive/folders/1fK345gOCKHtffluiFjPLUuRYbS7tK_yl?usp=sharing'

# modes to ignore
IGNORE_MODES = ['LongXD2.3','LongXD2.1','LongXD1.9','LongSO','LXD2.3','LXD2.1','LXD1.9','LXD_short','LXD_long','LowRes60']
IGNORE_SLITS = ['0.3x60','0.5x60','0.8x60','1.6x60','3.0x60']

# this is the data to extract from header, with options for header keywords depending on when file was written
HEADER_DATA={
	'UT_DATE': ['DATE_OBS'],
	'UT_TIME': ['TIME_OBS'],
#	'MJD': ['MJD_OBS'],
	'TARGET_NAME': ['OBJECT','TCS_OBJ'],
	'RA': ['RA','TCS_RA'],
	'DEC': ['DEC','TCS_DEC'],
	'HA': ['HA','TCS_HA'],
	'PA': ['POSANGLE'],
	'AIRMASS': ['AIRMASS','TCS_AM'],
	'INTEGRATION': ['ITIME'],
	'COADDS': ['CO_ADDS'],
	'SLIT': ['SLIT'],
	'MODE': ['GRAT'],
	'BEAM': ['BEAM'],
	'PARALLACTIC': ['TCS_PA'],
	'PROGRAM': ['PROG_ID'],
	'OBSERVER': ['OBSERVER'],
	'TARGET_TYPE': ['DATATYPE'],
}

# these are default parameters for generating the log
LOG_PARAMETERS = {
	'FILENAME': 'log.csv',
	'HTML_TEMPLATE_FILE' : os.path.join(DIR,'log_template.txt'),
	'BASE_HTML': r'<!DOCTYPE html>\n<html lang="en">\n\n<head>\n<title>SpeX Log</title>\n</head>\n\n<body>\n<h1>SpeX Log</h1>\n\n[TABLE]\n\n</body>\n</html>',	
	'COLUMNS' : ['FILE','BEAM','TARGET_NAME','TARGET_TYPE','FIXED-MOVING','RA','DEC','UT_DATE','UT_TIME','MJD','AIRMASS','INTEGRATION','COADDS','INSTRUMENT','SLIT','MODE','PROGRAM','OBSERVER'],
}

# these are default batch reduction parameters
BATCH_PARAMETERS = {
	'INSTRUMENT': 'spex',
	'DATA_FOLDER': '',
	'CALS_FOLDER': '',
	'PROC_FOLDER': '',
	'QA_FOLDER': '',
	'FLAT_FILE_PREFIX': 'flat',
	'ARC_FILE_PREFIX': 'arc',
	'SCIENCE_FILE_PREFIX': 'spc',
	'SPECTRA_FILE_PREFIX': 'spectra',
	'COMBINED_FILE_PREFIX': 'combspec',
	'CALIBRATED_FILE_PREFIX': 'calspec',
	'STITCHED_FILE_PREFIX': 'stitched',
	'PROGRAM': '',
	'OBSERVER': '',
	'qa_write': True,
	'qa_show': False, 
	'PLOT_TYPE': '.pdf', 
	'CALIBRATIONS': True,
	'EXTRACT': True,
	'COMBINE': True,
	'FLUXTELL': True,
	'REREDUCE': False,
	'STITCH': True,
	'OVERWRITE': False,
}

# these are required parameters for target-calibration sequences
OBSERVATION_PARAMETERS_REQUIRED = {
	'MODE': 'SXD',
	'TARGET_TYPE': 'fixed ps',
	'TARGET_REDUCTION_MODE': 'A-B',
	'TARGET_NAME': 'target',
	'TARGET_PREFIX': 'spc',
	'TARGET_FILES': '1-2',
	'TARGET_FLAT_FILE': 'flat.fits',
	'TARGET_WAVECAL_FILE': 'wavecal.fits',
	'STD_REDUCTION_MODE': 'A-B',
	'STD_NAME': 'standard',
	'STD_SPT': '',
	'STD_B': -99,
	'STD_V': -99.,
	'STD_PREFIX': 'spc',
	'STD_FILES': '1-2',
	'STD_FLAT_FILE': 'flat.fits',
	'STD_WAVECAL_FILE': 'wavecal.fits',
	}

# these are optional parameters for target-calibration sequences
# note that all with have a TARGET and STD prefix attached when actual code is run
OBSERVATION_PARAMETERS_OPTIONAL = {
#	'USE_STORED_SOLUTION': False,
	'ORDERS': '3-9',
#	'TARGET_TYPE': 'ps',
#	'REDUCTION_MODE': 'A-B',
	'NPOSITIONS': 2,
	'APERTURE_POSITIONS': [4.,11.5],
	'APERTURE_METHOD': 'auto',
	'APERTURE': 1.5,
	'PSF_RADIUS': 1.5,
	'BACKGROUND_RADIUS': 2.5,
	'BACKGROUND_WIDTH': 4,
	'SCALE_RANGE': [1.0,1.5],
}

# these are default parameters for QA page creation
QA_PARAMETERS = {
	'FILENAME': 'index.html',
	'TEMPLATE_FILE': os.path.join(DIR,'qa_template.txt'),
	'SOURCE_TEMPLATE_FILE': os.path.join(DIR,'qa_source_template.txt'),
	'SINGLE_PLOT_TEMPLATE_FILE': os.path.join(DIR,'qa_singleplot_template.txt'),
	'CSS_FILE' : os.path.join(DIR,'qa.css'),
	'PLOT_TYPE': '.pdf',
	'NIMAGES': 3,
	'IMAGE_WIDTH': 300,
	'MKWC_ARCHIVE_URL': 'http://mkwc.ifa.hawaii.edu/forecast/mko/archive/index.cgi',
	'HTML_TABLE_HEAD': '<table width=100%>\n <tr>\n',
	'HTML_TABLE_TAIL': ' </tr>\n</table>\n',
}

# columns to grab from SIMBAD
# SIMBAD_COLS = {
# 	'SIMBAD_SEP': 'angDist',
# 	'SIMBAD_NAME': 'main_id',
# 	'SIMBAD_TYPE': 'sp_type',
# 	'SIMBAD_BMAG': 'Bmag',
# 	'SIMBAD_VMAG': 'Vmag',
# 	'SIMBAD_GMAG': 'Gmag',
# 	'SIMBAD_JMAG': 'Jmag',
# 	'SIMBAD_HMAG': 'Hmag',
# 	'SIMBAD_KMAG': 'Kmag',
# }

# SIMBAD COLS - first value is the votable name, second value is the returned name
SIMBAD_COLS = {
#	'SIMBAD_SEP': 'SEP',
	'SIMBAD_NAME': ['main_id','MAIN_ID'],
	'SIMBAD_TYPE': ['sptype','SP_TYPE'],
	'SIMBAD_CLASS': ['otype','OTYPE'],
	'SIMBAD_BMAG': ['fluxdata(B)','FLUX_B'],
	'SIMBAD_VMAG': ['fluxdata(V)','FLUX_V'],
	'SIMBAD_GMAG': ['fluxdata(G)','FLUX_G'],
	'SIMBAD_JMAG': ['fluxdata(J)','FLUX_J'],
	'SIMBAD_HMAG': ['fluxdata(H)','FLUX_H'],
	'SIMBAD_KMAG': ['fluxdata(K)','FLUX_K'],
}

SIMBAD_EXCLUDE = ['Planet','Galaxy','QSO','Cluster*','Association','Stream','MouvGroup','GlobCluster','OpenCluster','BrightestCG','RadioG','LowSurfBrghtG','BlueCompactG','StarburstG','StarburstG','EmissionG','AGN','Seyfert','Seyfert1','Seyfert2','LINER','InteractingG','PairG','GroupG','Compact_Gr_G','ClG','protoClG','SuperClG','Void','Transient','HI','Maser','gammaBurst','PartofG','Unknown','Region','**']
SIMBAD_RADIUS = 30*u.arcsecond

# legacy download constants
LEGACY_FOLDER = 'irtfdata.ifa.hawaii.edu'
IRSA_FOLDER = 'IRTF'
IRSA_PREFIX = 'sbd.'
REDUCTION_FOLDERS = ['data','cals','proc','qa']

# observing schedule
OBSERVER_SCHEDULE_FILE = os.path.join(DIR,'observer_schedule.csv')

##################################
####### BASIC TEST OF CODE #######
##################################

def test(verbose=True):
	'''
	Purpose
	-------
	Checks to make sure installation is correct and necessary data files are present

	Parameters
	----------
	None

	Outputs
	-------
	Will raise warning errors if any aspects of code aren't working, otherwise returns True

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> batch.test(verbose=True)
	...checking that code can find instrument data
		PASS
	...checking that code can find reduction data
		PASS
	...checking that code can find test data
		PASS
	...checking that code can process data files
		PASS
	...checking that code can generate log
		PASS
	...checking that code can generate driver
		PASS
	*** Batch reduction tests have all passed - happy reducing! ***

	Dependencies
	------------
	processfolder()
	writeLog()
	writeDriver()
	readDriver()
	glob
	os
	pandas
	'''	
	testdatafold = os.path.join(CODEDIR,"../../tests/test_data/")
	instrumentdatafold = os.path.join(CODEDIR,"instrument_data/")
	reductiondatafold = os.path.join(CODEDIR,"data/")
	test_instruments = ['spex-prism', 'spex-SXD','uspex-prism','uspex-SXD']
	rawfold = os.path.join(testdatafold, "raw/")
	procfold = os.path.join(testdatafold, "processed/")
	log_test_file = 'log_test.csv'
	driver_test_file = 'driver_test.txt'

# make sure code can find main instrument data files
	if verbose==True: print('...checking that code can find instrument data')
	assert os.path.exists(instrumentdatafold), 'could not find instrument data folder {}, try downloading from {}'.format(instrumentdatafold,BACKUP_INSTRUMENT_DATA_URL)
	for x in ['spex','uspex']:
		testfold = os.path.join(instrumentdatafold,'{}_dir'.format(x))
		assert os.path.exists(testfold), 'could not find instrument data folder {}, try downloading from {}'.format(testfold,BACKUP_INSTRUMENT_DATA_URL)
		testfile = os.path.join(testfold,'{}_bdpxmk.fits'.format(x))
		assert os.path.exists(testfile), 'could not find bad pixel mask file {}, try downloading from {}'.format(testfile,BACKUP_INSTRUMENT_DATA_URL)
		fstat = os.stat(testfile)
		assert fstat.st_size > 1000000, 'file {} is too small, try downloading from {}'.format(testfile,BACKUP_INSTRUMENT_DATA_URL)
		if verbose==True: print('\t{}: PASS'.format(x))

# make sure code can find main reduction data files
	if verbose==True: print('...checking that code can find reduction data')
	assert os.path.exists(reductiondatafold), 'could not find reduction data folder {}, try downloading from {}'.format(reductiondatafold,BACKUP_REDUCTION_DATA_URL)
	for x in [100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,60000,75000]:
		atmfile = os.path.join(reductiondatafold,'atran{:.0f}.fits'.format(x))
		assert os.path.exists(atmfile), 'could not find atmosphere transmission file {}, try downloading from {}'.format(atmfile,BACKUP_REDUCTION_DATA_URL)
		fstat = os.stat(atmfile)
		assert fstat.st_size > 1000000, 'atmosphere transmission file {} is too small, try downloading from {}'.format(atmfile,BACKUP_REDUCTION_DATA_URL)
	for x in [5000,50000]:
		vfile = os.path.join(reductiondatafold,'Vega{:.0f}.fits'.format(x))
		assert os.path.exists(vfile), 'Vega spectrum file {} not found in data folder {}, try downloading from {}'.format(vfile,BACKUP_REDUCTION_DATA_URL)
		fstat = os.stat(vfile)
		assert fstat.st_size > 1000000, 'Vega spectrum file {} is too small, try downloading from {}'.format(vfile,BACKUP_REDUCTION_DATA_URL)
	if verbose==True: print('\tPASS')

# make sure code can find raw test data
	if verbose==True: print('...checking that code can find and process test data')
	assert os.path.exists(testdatafold), 'Could not find test data folder {}'.format(testdatafold)
	for inst in test_instruments:
		rfold = os.path.join(rawfold, inst)
		print(rawfold,rfold)
		assert os.path.exists(rfold), 'Could not find test data folder {}'.format(os.path.join(rawfold, inst))
		for fold in REDUCTION_FOLDERS:
			assert os.path.exists(os.path.join(rfold, fold)), 'Could not find test data folder {}'.format(os.path.join(rawfold, inst, fold))
		dfiles = glob.glob(os.path.join(rfold, "data/*.fits"))
		assert len(dfiles) > 0, 'no data files found in {}'.format(os.path.join(rfold, "data/*.fits"))
		result = processFolder(os.path.join(rfold, "data/"),verbose=False)
		assert isinstance(result, pd.core.frame.DataFrame), 'Could not process data folder {}'.format(os.path.join(rfold, "data/"))
		assert len(result) > 0, 'no data files found in {}'.format(os.path.join(rfold, "data/"))
		for x in list(HEADER_DATA.keys()):
			assert x in list(result.columns), 'Could not find required header data {} in data files in {}'.format(x,os.path.join(rfold, "data/"))
		if verbose==True: print('\t{}: PASS'.format(inst))
#	if verbose==True: print('\tPASS')

# process files and generate a test log
	if verbose==True: print('...checking that code can generate log')
	for inst in test_instruments:
		rfold = os.path.join(rawfold, inst)
		pfold = os.path.join(procfold, inst)
		logfile = os.path.join(pfold, "qa", log_test_file)
		dp = processFolder(os.path.join(rfold, "data/"),verbose=False)
		writeLog(dp,logfile,verbose=False)
		assert os.path.exists(logfile), 'Could not find log file {} after creation'.format(driver_file)
		dp = pd.read_csv(logfile)
		assert isinstance(dp, pd.core.frame.DataFrame), 'Log file parameter error: not a pandas dataframe for log file {}'.format(log_file)
		assert len(dp)>0
		for x in LOG_PARAMETERS['COLUMNS']:
			assert x in list(dp.columns), 'Could not find required column {} in log file {}'.format(x,logfile)
		for x in SIMBAD_COLS:
			assert x in list(dp.columns), 'Could not find required SIMBAD column {} in log file {}'.format(x,logfile)
		if verbose==True: print('\t{}: PASS'.format(inst))

# CLEANUP
# remove generated files

# process files and generate a driver file
	if verbose==True: print('...checking that code can generate driver')
	for inst in test_instruments:
		rfold = os.path.join(rawfold, inst)
		pfold = os.path.join(procfold, inst)
		logfile = os.path.join(pfold, "qa", log_test_file)
		driver_file = os.path.join(pfold, "proc", driver_test_file)
#	dp = processFolder(os.path.join(tfold, "data/"),verbose=False)
		dp = pd.read_csv(logfile)
		writeDriver(dp,driver_file,data_folder=os.path.join(rfold, "data"),check=False,create_folders=False,verbose=False)
		assert os.path.exists(driver_file), 'Could not find driver file {} after creation'.format(driver_file)

# read in existing driver in test folder and check
		par = readDriver(driver_file,verbose=False)
		assert isinstance(par, dict), 'Driver file parameter error: not a dictionary for driver file {}'.format(driver_file)
		for x in BATCH_PARAMETERS:
			assert x in list(par.keys()), 'Could not find required parameter {} in driver file {}'.format(x,driver_file)
		if verbose==True: print('\t{}: PASS'.format(inst))

# CLEANUP
# remove generated files
	if verbose==True: print('...cleaning up')
	for inst in test_instruments:
		rfold = os.path.join(rawfold, inst)
		pfold = os.path.join(procfold, inst)
		logfile = os.path.join(pfold, "qa", log_test_file)
		os.remove(logfile)
		driver_file = os.path.join(pfold, "proc", driver_test_file)
		os.remove(driver_file)
	if verbose==True: print('\n*** Batch reduction tests have all passed - happy reducing! ***\n')

	return

##################################
######### PROCESS FOLDERS ########
####### AND ORGANIZED DATA #######
##################################

def processFolder(folder,verbose=False):
	'''
	Purpose
	-------
	Reads in fits files from a data folder and organizes into pandas dataframe 
	with select keywords from headers defined in HEADER_DATA

	Parameters
	----------
	folder : str
		full path to the data folder containing the fits files

	verbose : bool, default=False
		Set to True to return verbose output

	Outputs
	-------
	pandas DataFrame containing the main header elements of the fits data files

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/adam/data/spex/20200101/data/'
	>>> dp = batch.processFolder(dpath)
	>>> print(dp[:4])

					   FILE	 UT_DATE		  UT_TIME TARGET_NAME		   RA  \
		0	spc0001.a.fits  2003-05-21  05:46:06.300712   1104+1959  11:04:01.12   
		1	spc0002.b.fits  2003-05-21  05:48:25.591305   1104+1959  11:04:01.08   
		2	spc0003.b.fits  2003-05-21  05:50:42.924982   1104+1959  11:04:01.06   
		3	spc0004.a.fits  2003-05-21  05:53:02.194063   1104+1959  11:04:01.08   

	Dependencies
	------------
	astropy.coordinates
	astropy.io
	astropy.time
	astropy.units
	glob
	os.path
	pandas
	'''
# error checking folder 
	if os.path.isdir(folder) == False: 
		raise ValueError('Cannot folder find {}'.format(folder))
	if folder[-1] != '/': folder=folder+'/'

# read in files, with accounting for .gz files
	files = glob.glob(folder+'/*.fits')
	gzflag = 0
	if len(files) == 0: 
		files = glob.glob(folder+'/*.fits.gz')
		gzflag = 1
	if len(files) == 0: 
		raise ValueError('Cannot find any .fits data files in {}'.format(folder))

# generate pandas data frame
	dp = pd.DataFrame()
	dp['FILE'] = [os.path.basename(f) for f in files]
	for k in list(HEADER_DATA.keys()): dp[k] = ['']*len(dp)

# run through each file header and populate the dataframe
	for i,f in enumerate(files):
		hdu = fits.open(f,ignore_missing_end=True)
		hdu[0].verify('silentfix')
		header = hdu[0].header
		hdu.close()
		for k in list(HEADER_DATA.keys()): 
			ref = ''
			for ii in HEADER_DATA[k]:
				if ii in list(header.keys()) and ref=='': ref=ii
			if ref!='': dp.loc[i,k] = header[ref]
			if ref=='' and verbose==True and i==0:
				print('Could not find keywords {} in file {}'.format(HEADER_DATA[k],f))
# update some of the mode names
		if 'SHORTXD' in dp.loc[i,'MODE'].upper(): dp.loc[i,'MODE'] = 'SXD'
		if 'lowres' in dp.loc[i,'MODE'].lower(): dp.loc[i,'MODE'] = 'Prism'		
		if 'arc' in dp.loc[i,'FILE']: 
			dp.loc[i,'TARGET_NAME'] = ARC_NAME
			dp.loc[i,'TARGET_TYPE'] = 'calibration'
		if 'flat' in dp.loc[i,'FILE']: 
			dp.loc[i,'TARGET_NAME'] = FLAT_NAME
			dp.loc[i,'TARGET_TYPE'] = 'calibration'
# set standards based on integration time
# this is not a great way to do this!
		if dp.loc[i,'TARGET_TYPE']=='': 
			dp.loc[i,'TARGET_TYPE']='target'
			if 'SXD' in dp.loc[i,'MODE'] and float(dp.loc[i,'INTEGRATION'])<=30.: dp.loc[i,'TARGET_TYPE']='standard'
			if 'Prism' in dp.loc[i,'MODE'] and float(dp.loc[i,'INTEGRATION'])<=10.: dp.loc[i,'TARGET_TYPE']='standard'

# if guess above was wrong, reset all standards back to targets
	dpt = dp[dp['TARGET_TYPE']=='target']
	dpt.reset_index(inplace=True)
	if len(dpt)==0: 
		dp.loc[dp['TARGET_TYPE']=='standard','TARGET_TYPE'] = 'target'
	dpt = dp[dp['TARGET_TYPE']=='standard']
	dpt.reset_index(inplace=True)
	if len(dpt)==0: 
		if verbose==True: print('Warning: no standards present in list; be sure to check driver file carefully')

# ignore certain modes
	dp.loc[dp['MODE'].isin(IGNORE_MODES),'TARGET_TYPE'] = 'ignore'
	# for igmode in IGNORE_MODES:
	# 	dp.loc[dp['MODE']==igmode,'TARGET_TYPE'] = 'ignore'

# ignore cases where the slit is 60"
	# for i in range(len(dp)):
	# 	if 'x60' in dp['SLIT'].iloc[i]: dp['TARGET_TYPE'].iloc[i] = 'ignore'
	dp.loc[dp['SLIT'].isin(IGNORE_SLITS),'TARGET_TYPE'] = 'ignore'

# fix to coordinates
	for k in ['RA','DEC']:
		dp[k] = [x.strip() for x in dp[k]]

# generate ISO 8601 time string and sort
# --> might need to do some work on this
	dp['DATETIME'] = [dp.loc[i,'UT_DATE']+'T'+dp.loc[i,'UT_TIME'] for i in range(len(dp))]
	dp['MJD'] = [Time(d,format='isot', scale='utc').mjd for d in dp['DATETIME']]
	dp['INSTRUMENT'] = ['spex']*len(dp)
	if dp.loc[0,'MJD']>INSTRUMENT_DATE_SHIFT: dp['INSTRUMENT'] = ['uspex']*len(dp)
	dp.sort_values('DATETIME',inplace=True)
	dp.reset_index(inplace=True,drop=True)

# fixed/moving - determined by range of position and motion of soruce
# NOTE: THIS FAILS IF USER DOESN'T CHANGE NAME OF SOURCE
	dp['FIXED-MOVING'] = ['fixed']*len(dp)
	names = list(set(list(dp['TARGET_NAME'])))
	if ARC_NAME in names: names.remove(ARC_NAME)
	if FLAT_NAME in names: names.remove(FLAT_NAME)
	if len(names)==0:
		if verbose==True: print('Warning: no science files identified in {}'.format(folder))
		return dp
	else:
		for n in names:
			dps = dp[dp['TARGET_NAME']==n]
			dpsa = dps[dps['BEAM']=='A']
			dpsa.reset_index(inplace=True)
			if len(dpsa)>1:	
				try:
					pos1 = SkyCoord(dpsa.loc[0,'RA']+' '+dpsa.loc[0,'DEC'],unit=(u.hourangle, u.deg))
					pos2 = SkyCoord(dpsa.loc[len(dpsa)-1,'RA']+' '+dpsa.loc[len(dpsa)-1,'DEC'],unit=(u.hourangle, u.deg))
					dx = pos1.separation(pos2).to(u.arcsec)
					time1 = Time(dpsa.loc[0,'DATETIME'],format='isot', scale='utc')
					time2 = Time(dpsa.loc[len(dpsa)-1,'DATETIME'],format='isot', scale='utc')
					dt = (time2-time1).to(u.hour)
					if dx.value > MOVING_MAXSEP or ((dx/dt).value>MOVING_MAXRATE and dt.value > 0.1):
						dp.loc[dp['TARGET_NAME']==n,'FIXED-MOVING'] = 'moving'
					if verbose==True: print('{}: dx={:.1f} arc, dt={:.1f} hr, pm={:.2f} arc/hr = {}'.format(n,dx,dt,dx/dt,dp.loc[dp['TARGET_NAME']==n,'FIXED-MOVING']))
				except Exception as e:
					if verbose==True: print('\tWarning: exception {} encountered for file {}; assuming fixed source'.format(e,dpsa['FILE'],iloc[0]))

# program and PI info (based on code by Evan Watson)
# NOTE: in this realization we're assuming its the same program and PI for the entire folder
# this overrules what is header
	obs_sch = pd.read_csv(OBSERVER_SCHEDULE_FILE, usecols=['MJD START','MJD END','PROGRAM','PI'])
	mjd = dp.loc[1,'MJD']
	obs_sch = obs_sch[obs_sch['MJD START']<mjd]
	obs_sch = obs_sch[obs_sch['MJD END']>mjd]
	obs_sch.reset_index(inplace=True)
	if len(obs_sch)>0:
		dp['PROGRAM'] = [obs_sch.loc[0,'PROGRAM']]*len(dp)
		dp['OBSERVER'] = [obs_sch.loc[0,'PI']]*len(dp)

# some additional cleanup
	dp['TARGET_NAME'] = [x.replace('.','_').replace(',','_') for x in dp['TARGET_NAME']]

	return dp




def organize(folder,outfolder='',verbose=ERROR_CHECKING,makecopy=False,overwrite=False):
	'''
	Purpose
	-------
	Organizes set of fits data files into distinct reduction datasets
	Specifically moves files from infolder into a data folder within a top level folder set by date
	If data are from IRTF Legacy archive or IRSA Archive, will automatically go to these settings
	Also creates proc, cals, and qa folders

	Parameters
	----------
	folder : str
		full path to the folder containing the data

	outfolder : str, default=''
		string indicating top-level organization folder under which data, cals, proc, and qa folders will be created
		if empty string, current directory is assumed 

	makecopy : bool, default=True
		If True, then copies data files; otherwise moves them

	overwrite : bool, default=False
		If True, overwrites files that are already in target directory, otherwise skips

	verbose : bool, default=ERROR_CHECKING
		Set to True to return verbose output

	Outputs
	-------
	None

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/me/spex/reductions/data/'
	>>> batch.organize(dpath)

	From command line:

	> python batch.cl --organize

	Dependencies
	------------
	astropy.io.fits
	glob
	os.path
	shutil
	batch.organizeLegacy()
	batch.organizeIRSA()
	'''	

# First look to see if LEGACY_FOLDER or IRSA_FOLDER is in the infolder path
	if LEGACY_FOLDER in folder: 
		if verbose==True: print('Detected IRTF Legacy archive path {} in input path; using organizeLegacy'.format(LEGACY_FOLDER))
		organizeLegacy(folder,outfolder=outfolder,verbose=verbose,makecopy=makecopy,overwrite=overwrite)
	if IRSA_FOLDER in folder: 
		if verbose==True: print('Detected IRTF IRSA archive path {} in input path; using organizeLegacy'.format(IRSA_FOLDER))
		organizeIRSA(folder,outfolder=outfolder,verbose=verbose,makecopy=makecopy,overwrite=overwrite)

# Generate a list of all fits files in folder
	fits_paths = glob.glob(os.path.join(folder, '*.fits'))
	fits_paths.sort()
	fits_files = [os.path.basename(f) for f in fits_paths]

# define folders based on file date (requires file read - can take a while)
	fits_outfolders = []
	for f in fits_paths:
		ofld = ''
		hdu = fits.open(f,ignore_missing_end=True)
		hd = hdu[0].header
		hdu.close()
		for k in HEADER_DATA['UT_DATE']:
			if k in list(hd.keys()): ofld+=hd[k].replace('/','').replace('-','')
		for k in HEADER_DATA['PROGRAM']:
			if k in list(hd.keys()): ofld+='-'+hd[k].replace('/','').replace('-','').replace(' ','')
		fits_outfolders.append(ofld)
	all_outfolders = list(set(fits_outfolders))

# make the necessary output folders
	for f in all_outfolders:
		output_folder = os.path.join(outfolder,f)
		if os.path.exists(output_folder)== False: 
			os.makedirs(output_folder)
			if verbose==True: print('Creating folder {}'.format(output_folder))

# create other folders
		for rf in REDUCTION_FOLDERS:
			rfolder = os.path.join(output_folder,rf)
			if os.path.exists(rfolder) == False:
				os.makedirs(rfolder)
				if verbose==True: print('Creating folder {}'.format(rfolder))
		data_folder = os.path.join(output_folder,REDUCTION_FOLDERS[0])

# run through fits files and copy relevant files over
		for i,fl in enumerate(fits_paths):
			if fits_outfolders[i]==f:
				if os.path.exists(os.path.join(data_folder,os.path.basename(fl)))==False or overwrite==True: 
					if makecopy==False: shutil.move(fl,data_folder)
					else: shutil.copy2(fl,data_folder)
				else: print('WARNING: fits file {} already exists in {}; move this file or set overwrite to True'.format(os.path.basename(fl),data_folder))

		if verbose==True: print('\nOrganization complete; data and reduction folders can be found in {}\n\n'.format(output_folder))
	return



def organizeLegacy(folder,outfolder='',verbose=ERROR_CHECKING,makecopy=False,overwrite=False):
	'''
	Purpose
	-------
	Organizes files downloaded from the IRTF Legacy database into distinct datasets
	Specifically moves files into a data folder within a top level folder set by date and program number
	Also creates proc, cals, and qa folders

	Parameters
	----------
	folder : str
		full path to the folder containing the downloaded LEGACY_FOLDER, into which other folders will be organized

	outfolder : str, default='dataset'
		string indicating top-level organization folder under which data, cals, proc, and qa folders will be created 

	makecopy : bool, default=True
		If True, then copies data files; otherwise moves them

	overwrite : bool, default=False
		If True, overwrites files that are already in target directory, otherwise skips

	verbose : bool, default=ERROR_CHECKING
		Set to True to return verbose output

	Outputs
	-------
	None

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/me/spex/reductions/'
	>>> batch.organizeLegacy(dpath)

	From command line:

	> python batch.cl --organize-legacy

	Dependencies
	------------
	astropy.io.fits
	glob
	os.path
	shutil
	'''	

# Look for 'irtfdata.ifa.hawaii.edu' folder in the current directory
	base_folder = os.path.join(folder, LEGACY_FOLDER)
	if not os.path.exists(base_folder):
		print("Error: {} folder not found in the current directory, no action".format(LEGACY_FOLDER))
		return

# Generate a list of all unique folders containing fits files
	fits_folders = glob.glob(os.path.join(base_folder, '**', '*.fits'), recursive=True)
	fits_folders = [os.path.dirname(fits_file) for fits_file in fits_folders]
	fits_folders = list(set(fits_folders))
	fits_folders.sort()

# now go through each of the folders, generate data folders, and move the data 	
	cntr,prefix0 = 1,''
	for i,f in enumerate(fits_folders):
# generate file name from UT date of first fits file if not provided
		if verbose==True: print('\nOrganizing data folder {}'.format(f))
		if outfolder=='':
			hdu = fits.open(glob.glob(os.path.join(f,'*.fits'))[0],ignore_missing_end=True)
			hd = hdu[0].header
			hdu.close()
			for k in HEADER_DATA['UT_DATE']:
				if k in list(hd.keys()): prefix = hd[k].replace('/','').replace('-','')
			prefix = os.path.join(folder,prefix)
		else: prefix = os.path.join(folder,prefix)
		if prefix!=prefix0:
			prefix0 = copy.deepcopy(prefix)
			cntr=1
		output_folder = os.path.abspath(prefix+'-{:.0f}'.format(cntr))
		cntr+=1
		if os.path.exists(output_folder)== False: 
			os.makedirs(output_folder)
			if verbose==True: print('Creating folder {}'.format(output_folder))

# if path exists, just move and rename fits folder, otherwise move files, being careful of overwriting 
		data_folder = os.path.join(output_folder,'data')
		if os.path.exists(data_folder) == False: 
			if makecopy==False: shutil.move(f, data_folder)
			else: shutil.copytree(f,data_folder)
			if verbose==True: print('Moved legacy data folder {} to {}'.format(f,data_folder))
		else: 
			if len(glob.glob(os.path.join(f,'*.fits')))==0 or overwrite==True:
				for df in glob.glob(os.path.join(f,'*.fits')): 
					if makecopy==False: shutil.move(df,data_folder)
					else: shutil.copy2(df,data_folder)
				if verbose==True: print('Moved fits files in legacy data folder {} to {}'.format(f,data_folder))
			else: print('WARNING: fits files already exist in {}; move these files or set overwrite to True'.format(data_folder))

# create other folders
		for rf in REDUCTION_FOLDERS[1:]:
			rfolder = os.path.join(output_folder,rf)
			if os.path.exists(rfolder) == False:
				os.makedirs(rfolder)
				if verbose==True: print('Creating folder {}'.format(rfolder))

		if verbose==True: print('\n\nOrganization complete; data and reduction folders can be found in {}\n\n'.format(output_folder))
	return



def organizeIRSA(folder,outfolder='',verbose=ERROR_CHECKING,makecopy=True,overwrite=False):
	'''
	Purpose
	-------
	Organizes files downloaded from the IRTF IRSA archive database into distinct datasets based on date and program
	Specifically moves files into a data folder within a top level folder set by date and program number
	Also creates proc, cals, and qa folders

	Parameters
	----------
	folder : str
		full path to the folder containing the downloaded IRSA_FOLDER, into which other folders will be organized

	outfolder : str, default=''
		string indicating top-level organization folder under which data, cals, proc, and qa folders will be created 
		if empty string, assumes current directory

	makecopy : bool, default=True
		If True, then copies data files; otherwise moves them

	overwrite : bool, default=False
		If True, overwrites files that are already in target directory, otherwise skips

	verbose : bool, default=ERROR_CHECKING
		Set to True to return verbose output

	Outputs
	-------
	None

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/me/spex/reductions/'
	>>> batch.organizeIRSA(dpath)

	From command line:
	> python batch.cl --organize-irsa


	Dependencies
	------------
	glob
	os.path
	shutil
	'''	

# Look for 'irtfdata.ifa.hawaii.edu' folder in the current directory
	base_folder = os.path.join(folder, IRSA_FOLDER)
	if not os.path.exists(base_folder):
		raise ValueError("ERROR: IRSA archive output folder {} not found in the current directory".format(IRSA_FOLDER))

# Generate a list of all fits files
	fits_paths = glob.glob(os.path.join(base_folder, '**', '*.fits'), recursive=True)
	fits_paths.sort()
	fits_files = [os.path.basename(f) for f in fits_paths]
# first check
	if fits_files[0][:len(IRSA_PREFIX)]!=IRSA_PREFIX:
		raise ValueError("ERROR: IRSA archive file prefix {} not found with first file".format(IRSA_PREFIX))

# define folders based on program number and date based on typical fits file name
	padstr = ''
	fits_programs = [x.split('.')[1] for x in fits_files]
	fits_dates = [x.split('.')[2] for x in fits_files]
	if len(fits_dates[0])==6: padstr='20'
	fits_folders = [padstr+x.split('.')[2]+'-'+x.split('.')[1] for x in fits_files]
	fits_folders = list(set(list(fits_folders)))
	fits_folders.sort()

# make the necessary output folders
	for f in fits_folders:
		output_folder = os.path.join(folder,f)
		if os.path.exists(output_folder)== False: 
			os.makedirs(output_folder)
			if verbose==True: print('Creating folder {}'.format(output_folder))

# create other folders
		for rf in REDUCTION_FOLDERS:
			rfolder = os.path.join(output_folder,rf)
			if os.path.exists(rfolder) == False:
				os.makedirs(rfolder)
				if verbose==True: print('Creating folder {}'.format(rfolder))
		data_folder = os.path.join(output_folder,REDUCTION_FOLDERS[0])

# run through fits files and copy relevant file over
		mtchstr = f.split('-')[1]+'.'+f.split('-')[0][len(padstr):]
		mtch = [s for s in fits_paths if mtchstr in s]
		if verbose==True: print('Copying {} files in {}'.format(len(mtch),data_folder))
		for i,fl in enumerate(mtch): 
			if os.path.exists(os.path.join(data_folder,os.path.basename(fl)))==False or overwrite==True: 
				if makecopy==False: shutil.move(fl,data_folder)
				else: shutil.copy2(fl,data_folder)
			else: print('WARNING: fits file {} already exists in {}; move this file or set overwrite to True'.format(os.path.basename(fl),data_folder))

		if verbose==True: print('\n\nOrganization complete; data and reduction folders can be found in {}\n\n'.format(output_folder))
	return


##############################
####### LOG GENERATION #######
##############################

def logTablestyle(val):
	'''
	Plotting style for pandas dataframe for html display
	NOT ACTUALLY SURE THIS IS DOING ANYTHING
	'''
	color='white'
	if isinstance(val,str) == True:
		if 'LowRes' in val: color='blue'
		if 'SXD' in val: color='green'
		if 'LXD' in val: color='red'
	return 'background-color: {}'.format(color)


def writeLog(dp,log_file='',options={},verbose=ERROR_CHECKING):
	'''
	Purpose
	-------
	Takes in dataFrame from processFolder and generates a log file based on 
	select header values defined in LOG_COLUMNS parameter in LOG_PARAMETERS

	Parameters
	----------
	dp : pandas DataFrame
		output of processFolder

	log_file : str, default=LOG_DEFAULT_FILENAME
		full path to output log file; 
		allowed file formats are .csv, .htm, .html, .json, .psv, .tab, .tex, .tsv, .txt, .xls, .xlsx 

	options : dict, default = {}
		options for log generation parameters that override defaults given in LOG_PARAMETERS,

	verbose : bool, default=ERROR_CHECKING
		Set to True to return verbose output

	Outputs
	-------
	No explicit returns; log file is output as specified by file extension

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/adam/data/spex/20200101/data/'
	>>> dp = batch.processFolder(dpath)
	>>> batch.writeLog(dp,dpath+'log.html',verbose=True)
		
		log written to /Users/adam/projects/spex_archive/testing/spex-prism/log.html
		
	Dependencies
	------------
	copy
	os.path
	pandas
	'''
# update default parameters with options
	parameters = copy.deepcopy(LOG_PARAMETERS)
	parameters.update(options)
	if log_file != '': parameters['FILENAME'] = log_file

# downselect columns to save out; by default
	if len(parameters['COLUMNS']) > 0:
		try: dpout = dp[parameters['COLUMNS']]
		except: 
			dpout = copy.deepcopy(dp)
			if verbose==True: print('Warning: could not select subset of columns {} from data frame, saving all columns'.format(parameters['COLUMNS']))
	else: dpout = copy.deepcopy(dp)

# fix for RA and DEC for proper formatting
# catch is for rare cases where telescope data is missed
	for x in ['RA','DEC']:
		if x in list(dpout.columns):
			for i,y in enumerate(dpout[x]):
				try: 
					if str(y)[0] not in ['-','+']: dpout.loc[i,x]='+{}'.format(str(y))
				except Exception as e: 
					if verbose==True: print('\n\tWARNING: {} error: {} for file {} is missing\n'.format(e,x,dpout.loc[i,'FILE']))

# FIX FOR PATHOLOGICAL CASE WHERE NO TARGET NAMES ARE GIVEN
	if 'TARGET_NAME' in list(dpout.columns) and 'TARGET_TYPE' in list(dpout.columns):
		dpouts = dpout[dpout['TARGET_TYPE'] != 'calibration']
		dpouts.reset_index(inplace=True)
		tnames = list(set(list(dpouts['TARGET_NAME'])))
		tnames = [str(x) for x in tnames]
		if len(tnames)==1 and len(dpouts)>50:
			if verbose==True: print('WARNING: only one target name {} suggesting an issue with labeling'.format(tnames[0]))
			if 'RA' in list(dpout.columns) and 'DEC' in list(dpout.columns):
				for i in range(len(dpout)):
					if tnames[0] in dpout.loc[i,'TARGET_NAME']:
						dpout.loc[i,'TARGET_NAME'] = ('J{}{}'.format((dpout.loc[i,'RA'])[:6],(dpout.loc[i,'DEC'])[:6])).replace('J+','J').replace(':','')
				if verbose==True: print('Replaced target names with coordinate short names')

# add in SIMBAD information using astroquery.simbad
# try statement to catch cases where simbad is not accessible
	for x in list(SIMBAD_COLS.keys()): dpout.loc[:,x] = ['']*len(dpout)
	try:
		sb = Simbad()
		for x in SIMBAD_COLS.keys(): sb.add_votable_fields(SIMBAD_COLS[x][0])
		if 'RA' in list(dpout.columns) and 'DEC' in list(dpout.columns):
			dpouts = dpout[dpout['TARGET_TYPE'] != 'calibration']
			dpouts.reset_index(inplace=True)
			tnames = list(set(list(dpouts['TARGET_NAME'])))
			tnames = [str(x) for x in tnames]
			for i,tnm in enumerate(tnames):
				dpsel = dpout[dpout['TARGET_NAME']==tnm]
				dpsel.reset_index(inplace=True)
				if len(dpsel)>0:
					src_coord = SkyCoord(dpsel.loc[0,'RA']+' '+dpsel.loc[0,'DEC'],unit=(u.hourangle, u.deg))
					t_match = sb.query_region(src_coord,radius=SIMBAD_RADIUS)
					if isinstance(t_match,type(None))==False:
						if len(t_match)>0:
							for x in SIMBAD_EXCLUDE: t_match.remove_rows(np.where(t_match['OTYPE']==x))
						if len(t_match)>0:
							t_match['SIMBAD_SEP'] = [src_coord.separation(SkyCoord(str(t_match['RA'][lp]),str(t_match['DEC'][lp]),unit=(u.hourangle,u.degree))).arcsecond for lp in np.arange(len(t_match))]
							t_match.sort(['SIMBAD_SEP'])
							for x in list(SIMBAD_COLS.keys()):
								dpout.loc[dpout['TARGET_NAME']==tnm,x] = t_match[SIMBAD_COLS[x][1]][0]
				# t = Table([[src_coord.ra.deg],[src_coord.dec.deg]],names=('ra','dec'))
			# t_match = XMatch.query(t,u'simbad',SIMBAD_RADIUS,colRA1='ra',colDec1='dec',columns=["**", "+_r"])
			# t_match.sort(['angDist'])
			# if len(t_match)>0:
			# 	for x in SIMBAD_EXCLUDE: t_match.remove_rows(np.where(t_match['main_type']==x))
			# print(t_match[0])
			# print(t_match.columns)
			# if len(t_match)>0:
			# 	for x in list(SIMBAD_COLS.keys()):
				# 		dpout.loc[dpout['TARGET_NAME']==tnm,x] = t_match[SIMBAD_COLS[x]][0]
	except: print('WARNING: There was a problem in trying to match targets to SIMBAD via astroquery.simbad; check internet connection')

# add on a NOTES column
	dpout['NOTES'] = ['']*len(dpout)
	
# save depends on file name
	ftype = parameters['FILENAME'].split('.')[-1]
	if ftype in ['xls','xlsx']: dpout.to_excel(parameters['FILENAME'],index=False)	
	elif ftype in ['csv']: dpout.to_csv(parameters['FILENAME'],sep=',',index=False)	
	elif ftype in ['tsv','txt','tab']: dpout.to_csv(parameters['FILENAME'],sep='\t',index=False)	
	elif ftype in ['psv']: dpout.to_csv(parameters['FILENAME'],sep='|',index=False)	
	elif ftype in ['tex']: dpout.to_latex(parameters['FILENAME'],longtable=True,index=False)	
	elif ftype in ['json']: dpout.to_json(parameters['FILENAME'])	
	elif ftype in ['htm','html']:
		dpout.style.apply(logTablestyle)
		if parameters['HTML_TEMPLATE_FILE'] != '' and os.path.exists(parameters['HTML_TEMPLATE_FILE']):
			with open(parameters['HTML_TEMPLATE_FILE'],'r') as f: dphtml = f.read()
		else: dphtml = copy.deepcopy(parameters['BASE_HTML'])
		for x in ['INSTRUMENT','UT_DATE']:
			if x in list(dp.columns): 
				dphtml = dphtml.replace('[{}]'.format(x),dp.loc[0,x])
		dphtml = dphtml.replace('[TABLE]',dpout.to_html(classes='table table-striped text-center',index=False,bold_rows=True))
		with open(parameters['FILENAME'],'w') as f: f.write(dphtml)
	else: raise ValueError('Could not write out to {}; unknown file format'.format(parameters['FILENAME']))

	if verbose==True: print('log written to {}'.format(parameters['FILENAME']))
	return


#############################
##### BATCH DRIVER FILE #####
#############################

def readDriver(driver_file,options={},verbose=ERROR_CHECKING):
	'''
	Purpose
	-------
	Reads in batch driver file into a dictionary for processing

	Parameters
	----------
	driver_file : str
		full path to output file for batch driver

	options : dict, default = {}
		options for batch reduction parameters that override defaults given in BATCH_PARAMETERS,
		but are themselves overridden by what is in batch driver file

	verbose : bool, default=ERROR_CHECKING
		set to True to return verbose output

	Outputs
	-------
	Returns dictionary of batch reduction parameters

	Example
	-------
	>>> from pyspextool.batch import batch
	>>> import yaml
	>>> driver_file = '/Users/adam/projects/spex_archive/testing/150512/proc/driver.txt'
	>>> params = batch.readDriver(driver_file)
	>>> print(yaml.dump(params))
		
		ARC_FILE_PREFIX: arc-
		CALIBRATED_FILE_PREFIX: calspec
		CAL_FOLDER: /Users/adam/projects/spex_archive/testing/150512/cals/
		CAL_SETS: 34-39,70-75,92-97,114-119,142-147
		COMBINED_FILE_PREFIX: combspec
		DATA_FOLDER: /Users/adam/projects/spex_archive/testing/150512/data/
		...
		
	Dependencies
	------------
	copy
	os.path
	pandas
	yaml
	'''
# update default parameters with passed options
	parameters = copy.deepcopy(BATCH_PARAMETERS)
	parameters.update(options)

# check presence of file
	if os.path.exists(driver_file)==False:
		raise ValueError('Cannot find batch driver file {}; check path'.format(driver_file))

# read in and clean out comments and newlines
	f = open(driver_file,'r')
	lines = f.readlines()
	f.close()
	lines = list(filter(lambda x: (x != ''), lines))
	lines = list(filter(lambda x: (x != '\n'), lines))
	lines = list(filter(lambda x: (x[0] != '#'), lines))

# extract keyword = value lines outside target observations and use to update parameters
	klines = list(filter(lambda x: ('=' in x), lines))
	klines = list(filter(lambda x: OBSERVATION_SET_KEYWORD not in x, klines))
	klines = [l.replace(' ','').replace('\n','') for l in klines]
	new_kpar = dict([l.strip().split('=') for l in klines])
	parameters.update(new_kpar)

# insert the individual science sets
	slines = list(filter(lambda x: OBSERVATION_SET_KEYWORD in x, lines))
	for i,sline in enumerate(slines):
		spar = copy.deepcopy(OBSERVATION_PARAMETERS_REQUIRED)
		for k in list(OBSERVATION_PARAMETERS_OPTIONAL.keys()):
			for x in ['STD','TARGET']: spar['{}_{}'.format(x,k)] = OBSERVATION_PARAMETERS_OPTIONAL[k]
# update with global batch parameters			
		for k in list(spar.keys()):
			if k in list(parameters.keys()):
				if parameters[k]!='': spar[k] = parameters[k]
# update with passed required parameters			
		rpars = sline.replace('\n','').split('\t')
		if len(rpars) < len(list(OBSERVATION_PARAMETERS_REQUIRED.keys())):
			if verbose==True: print('Warning: line\n{}\ncontains fewer than the required parameters\n{}; skipping'.format(sline,OBSERVATION_PARAMETERS_REQUIRED))
		else:
			for ii,k in enumerate(list(OBSERVATION_PARAMETERS_REQUIRED.keys())): spar[k] = rpars[ii+1]
# add optional parameters			
			opars = list(filter(lambda x: ('=' in x), rpars))
			new_opar = dict([l.strip().split('=') for l in opars])
			for k in list(new_opar.keys()): new_opar[k.upper()] = new_opar.pop(k)
# process special cases
			for k in ['GUESS','FIXED']:
				if k in list(new_opar.keys()):
					new_opar['TARGET_APERTURE_METHOD'] = k.lower()
					new_opar['TARGET_APERTURE_POSITIONS'] = [float(x) for x in new_opar[k].split(',')]
			if 'SUB' in list(new_opar.keys()):
				new_opar['TARGET_REDUCTION_MODE'] = new_opar['SUB'].strip()
			for k in ['ORDER','ORDERS']:
				if k in list(new_opar.keys()):
					for x in ['TARGET','STD']: new_opar['{}_ORDERS'.format(x)] = new_opar[k]
			spar.update(new_opar)
# some instrument specific adjustments
# THIS MIGHT NEED TO BE MANAGED THROUGH OTHER PYSPEXTOOL FUNCTIONS
			if spar['MODE'] in ['Prism','LowRes15']: 
				for x in ['TARGET','STD']: spar['{}_ORDERS'.format(x)] = '1'
#				if verbose==True: print('Updated prism mode orders to 1')
			for x in ['STD','TARGET']:
				if parameters['INSTRUMENT'] == 'spex' and 'SXD' in spar['MODE'] and spar['{}_ORDERS'.format(x)][-1]=='9':
					spar['{}_ORDERS'.format(x)] = spar['{}_ORDERS'.format(x)][:-1]+'8'
#				if verbose==True: print('Updated spex SXD mode orders to {}'.format(spar['ORDERS']))

		parameters['{}-{}'.format(OBSERVATION_SET_KEYWORD,str(i+1).zfill(4))] = spar

# some specific parameters that need format conversation
	for k in ['qa_write','qa_show']:
		if k in list(parameters.keys()): parameters[k] = bool(parameters[k])


# report out parameters
	# if verbose==True:
	# 	print('\nDriver parameters:')
	# 	print(yaml.dump(parameters))

	return parameters


def writeDriver(dp,driver_file='driver.txt',data_folder='',options={},create_folders=False,comment='',check=ERROR_CHECKING,verbose=ERROR_CHECKING):
	'''
	Purpose
	-------
	Writes the batch driver file based on pandas DataFrame input from processFolder and 
	additional user-supplied options

	Parameters
	----------
	dp : pandas DataFrame
		output of processFolder

	driver_file : str, default = 'driver.txt'
		full path to output file for batch driver

	data_folder : str, default = ''
		full path to data folder containing raw science files

	options : dict, default = {}
		options for driver parameter that override defaults given in BATCH_PARAMETERS

	create_folders: bool, default = False
		if set to True, create proc, cals, and qa folders if they don't exist

	comment: str, default = ''
		comment string to add to head of driver file for information purposes

	check : bool, default=ERROR_CHECKING
		set to True to check contents of driver file by reading it in and returnin parameters

	verbose : bool, default=ERROR_CHECKING
		set to True to return verbose output

	Outputs
	-------
	If check==True, returns dictionary of batch reduction parameters (output of readDriver)
	otherwise returns nothing; in either case, writes contents to file specified by driver_file


	Example
	-------
	(1) Create from processFolder output:

	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/adam/projects/spex_archive/testing/spex-prism/data/'
	>>> driver_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt'
	>>> dp = batch.processFolder(dpath,verbose=False)
	>>> pars = batch.writeDriver(dp,driver_file=driver_file,data_folder=dpath,create_folders=True,check=True,verbose=True)
			
			Created CALS_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/cals/
			Created PROC_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/proc/
			Created QA_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/qa/

			Batch instructions written to /Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt, please this check over before proceeding

			Driver parameters:
			ARC_FILE_PREFIX: arc
			CALIBRATED_FILE_PREFIX: calspec
			CALS_FOLDER: /Users/adam/projects/spex_archive/testing/spex-prism/cals
			CAL_SETS: 15-20
			...

	(2) Create from previously created log file:
	>>> from pyspextool.batch import batch
	>>> import pandas
	>>> dpath = '/Users/adam/projects/spex_archive/testing/spex-prism/data/'
	>>> log_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/logs.csv'
	>>> driver_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt'
	>>> dp = pandas.read_csv(log_file,sep=',')
	>>> pars = batch.writeDriver(dp,driver_file=driver_file,data_folder=dpath,create_folders=True, check=True, verbose=True)

			Batch instructions written to /Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt, please this check over before proceeding

			Driver parameters:
			ARC_FILE_PREFIX: arc
			CALIBRATED_FILE_PREFIX: calspec
			CALS_FOLDER: /Users/adam/projects/spex_archive/testing/spex-prism/cals
			CAL_SETS: 15-20
			...
		
	Dependencies
	------------
	copy
	numpy
	os.path
	pandas
	readDriver()
	'''
# manage default options
	driver_param = copy.deepcopy(BATCH_PARAMETERS)
	driver_param.update(options)
	if data_folder != '': driver_param['DATA_FOLDER'] = data_folder

# check folders
	if os.path.exists(driver_param['DATA_FOLDER'])==False:
		raise ValueError('Cannot find data folder {}: you must include this in batch driver file'.format(driver_param['DATA_FOLDER']))
	for x in ['CALS_FOLDER','PROC_FOLDER','QA_FOLDER']:
# if blank, create cals, proc, and qa folders parallel to data folder
		if driver_param[x]=='': 
			driver_param[x] = os.path.join(os.path.split(driver_param['DATA_FOLDER'])[0],(x.split('_')[0]).lower())
		if os.path.exists(driver_param[x])==False:
			if create_folders==True: 
				os.mkdir(driver_param[x])
				if verbose==True: print('Created {} folder {}'.format(x,driver_param[x]))
			else: raise ValueError('Folder {} does not exist; either create or set keyword create_folders to True'.format(driver_param[x]))

# get information from first line of log
	dpc = copy.deepcopy(dp)
	for x in ['INSTRUMENT','UT_DATE','OBSERVER','PROGRAM']:
		driver_param[x] = dpc.loc[0,x]

# get rid of all ignored modes
#	for ignm in IGNORE_MODES:
#		dpc = dpc[dpc['MODE']!=ignm]

# some file name work
# THIS NEEDS TO BE UPDATED WITH INSTRUMENT SPECIFIC INFORMATION
	dpc['PREFIX'] = [x.replace('.a.fits','').replace('.b.fits','') for x in dpc['FILE']]
	n=4
	if driver_param['INSTRUMENT']=='uspex': n=5
	dpc['FILE NUMBER'] = [int(x[-n:]) for x in dpc['PREFIX']]
	dpc['PREFIX'] = [x[:-n] for x in dpc['PREFIX']]

# set up default prefixes for flats and arcs
	dps = dpc[dpc['TARGET_NAME']==FLAT_NAME]
	dps.reset_index(inplace=True)
	driver_param['FLAT_FILE_PREFIX'] = dps.loc[0,'PREFIX']
	dps = dpc[dpc['TARGET_NAME']==ARC_NAME]
	dps.reset_index(inplace=True)
	driver_param['ARC_FILE_PREFIX'] = dps.loc[0,'PREFIX']
	dps = dpc[dpc['TARGET_NAME']!=FLAT_NAME]
	dps = dps[dps['TARGET_NAME']!=ARC_NAME]
	dps.reset_index(inplace=True)
	driver_param['SCIENCE_FILE_PREFIX'] = dps.loc[0,'PREFIX']

# write out instructions
	try: f = open(driver_file,'w')
	except: raise ValueError('Cannot write to {}, check file path and permissions'.format(driver_file))
	f.write('# pyspextool batch reduction driver file for {} observations on {}\n'.format(dpc.loc[0,'INSTRUMENT'],dpc.loc[0,'UT_DATE']))
	f.write('# Code version {}\n'.format(VERSION))
	if comment!='': f.write('# {}\n'.format(comment.replace('\n','\n# ')))

# observational information
	f.write('\n# Observational information\n')
	for x in ['INSTRUMENT','UT_DATE','PROGRAM','OBSERVER']:
		f.write('{} = {}\n'.format(x,str(driver_param[x])))

# files & folders
	f.write('\n# Folders and files\n')
	for x in ['DATA_FOLDER','CALS_FOLDER','PROC_FOLDER','QA_FOLDER']:
		f.write('{} = {}\n'.format(x,str(driver_param[x])))
	for x in ['SCIENCE_FILE_PREFIX','SPECTRA_FILE_PREFIX','COMBINED_FILE_PREFIX','CALIBRATED_FILE_PREFIX','STITCHED_FILE_PREFIX']:
		f.write('{} = {}\n'.format(x,str(driver_param[x])))

# cals - determine here but print below
	# f.write('\n# Lamp calibrations\n')
	# for x in ['FLAT_FILE_PREFIX','ARC_FILE_PREFIX']:
	# 	f.write('{} = {}\n'.format(x,driver_param[x]))
	cal_sets=''
	dpcal = dpc[dpc['TARGET_TYPE']=='calibration']
	dpcal.reset_index(inplace=True)
	fnum = np.array(dpcal['FILE NUMBER'])
	calf1 = fnum[np.where(np.abs(fnum-np.roll(fnum,1))>1)]
	calf2 = fnum[np.where(np.abs(fnum-np.roll(fnum,-1))>1)]
	for i in range(len(calf1)): cal_sets+='{:.0f}-{:.0f},'.format(calf1[i],calf2[i])
#	f.write('CAL_SETS = {}\n'.format(cal_sets[:-1]))

# observations
	dps = dpc[dpc['TARGET_TYPE']=='target']
	dps.reset_index(inplace=True)
	tnames = list(set(list(dps['TARGET_NAME'])))
	tnames = [str(x) for x in tnames]
	fnums = []
	for n in tnames:
		dps = dpc[dpc['TARGET_NAME']==n]
		dps.reset_index(inplace=True)
		fnums.append(int(dps.loc[0,'FILE NUMBER']))
	tnames = [x for _,x in sorted(zip(fnums,tnames))]

# no names - close it up
	if len(tnames)==0:
		if verbose==True: print('Warning: no science files identified; you may need to update the observational logs')
		f.close()
		return

# loop over names
	f.write('\n# Science observations')
	f.write('\n#\tMode\tTarget Type\tTarget Sub\tTarget Name\tTarget Prefix\tTarget Files\tTarget Flat\tTarget Wavecal\tStd Sub\tStd Name\tStd SpT\tStd B\tStd V\tStd Prefix\tStd Files\tStd Flat\tStd Wavecal\tOptions (key=value)\n'.zfill(len(OBSERVATION_SET_KEYWORD)))
	for n in tnames:
		dps = dpc[dpc['TARGET_NAME']==n]
		dps = dps[dps['TARGET_TYPE']=='target']
		dps.reset_index(inplace=True)
		src_coord = SkyCoord(dps.loc[0,'RA']+' '+dps.loc[0,'DEC'],unit=(u.hourangle, u.deg))

# loop over modes, dropping ignored modes
		srcmodes = list(set(list(dps['MODE'])))
		for ignm in IGNORE_MODES: srcmodes = [i for i in srcmodes if i != ignm]
		srcmodes.sort()
		for m in srcmodes:
			dpsrc = dps[dps['MODE']==m]
			dpsrc.reset_index(inplace=True)
			line='{}\t{}\t{}-ps\tA-B\t{}\t{}'.format(OBSERVATION_SET_KEYWORD,str(m),str(dpsrc.loc[0,'FIXED-MOVING']),n,str(dpsrc.loc[0,'PREFIX']))
#			print(n,m,np.array(dpsrc['FILE NUMBER']),numberList(np.array(dpsrc['FILE NUMBER'])))
			line+='\t{}'.format(numberList(np.array(dpsrc['FILE NUMBER'])))

# assign calibrations
			dpcals = dpcal[dpcal['MODE']==m]
			dpcals.reset_index(inplace=True)
			if len(dpcals)==0: 
				if verbose==True: print('WARNING: no calibration files associated with mode {} for source {}'.format(dpsrc.loc[0,'MODE'],n))
				cs = cal_sets.split(',')[0]
			else:
				cnum = np.array(dpcals['FILE NUMBER'])
				cnum1 = cnum[np.where(np.abs(cnum-np.roll(cnum,1))>1)]
				cnum2 = cnum[np.where(np.abs(cnum-np.roll(cnum,-1))>1)]
				ical = np.argmin(np.abs(cnum1-np.median(dpsrc['FILE NUMBER'])))
				cs = '{:.0f}-{:.0f}'.format(cnum1[ical],cnum2[ical])
			line+='\tflat{}.fits\twavecal{}.fits'.format(cs,cs)

# assign flux cals based on closest in airmass (0.2), time (2 hr) and position (10")
			dpflux = dpc[dpc['TARGET_TYPE']=='standard']
			dpflux = dpflux[dpflux['MODE']==dpsrc.loc[0,'MODE']]
			dpflux.reset_index(inplace=True)
			ftxt = '\tA-B\tUNKNOWN\t{}\t0-0\tflat{}.fits\twavecal{}.fits'.format(str(driver_param['SCIENCE_FILE_PREFIX']),cal_sets.split(',')[ical],cal_sets.split(',')[ical])
			if len(dpflux)==0: 
				if verbose==True: print('WARNING: no calibration files associated with mode {} for source {}'.format(dpsrc.loc[0,'MODE'],n))
				line+=ftxt
			else:
				dpflux['DIFF1'] = [np.abs(float(x)-np.nanmedian(dpsrc['AIRMASS']))/0.2 for x in dpflux['AIRMASS']]
				dpflux['DIFF2'] = [np.abs(float(x)-np.nanmedian(dpsrc['MJD']))*12. for x in dpflux['MJD']]
				dpflux['DIFF3'] = [src_coord.separation(SkyCoord(dpflux.loc[i,'RA']+' '+dpflux.loc[i,'DEC'],unit=(u.hourangle, u.deg))).to(u.arcsecond).value/10. for i in range(len(dpflux))]
				dpflux['DIFF'] = dpflux['DIFF1']+dpflux['DIFF2']+dpflux['DIFF3']
				fref = dpflux.loc[np.argmin(dpflux['DIFF']),'FILE NUMBER']
				tname = str(dpflux.loc[np.argmin(dpflux['DIFF']),'TARGET_NAME'])
# temporary fix while sorting out why this doesn't work in testing				
				if 'SIMBAD_TYPE' in list(dpflux.keys()): tspt =  str(dpflux.loc[np.argmin(dpflux['DIFF']),'SIMBAD_TYPE'])
				else: tspt = ''
				if 'SIMBAD_BMAG' in list(dpflux.keys()): tbmag = str(dpflux.loc[np.argmin(dpflux['DIFF']),'SIMBAD_BMAG'])
				else: tbmag = 8.
				if 'SIMBAD_VMAG' in list(dpflux.keys()): tvmag = str(dpflux.loc[np.argmin(dpflux['DIFF']),'SIMBAD_VMAG'])
				else: tvmag = 8.
				dpfluxs = dpflux[dpflux['TARGET_NAME']==tname]
				dpfluxs.reset_index(inplace=True)
				fnum = np.array(dpfluxs['FILE NUMBER'])
				if len(fnum)<2: 
					if verbose==True: print('Fewer than 2 flux calibrator files for source {}'.format(tname))
					line+=ftxt
				else:				
					if len(dpcals)>0: 
						ical = np.argmin(np.abs(cnum1-np.median(dpfluxs['FILE NUMBER'])))
#						ical = np.argmin(np.abs(np.array(dpcals['FILE NUMBER'])-np.median(dpfluxs['FILE NUMBER'])))
						cs = '{:.0f}-{:.0f}'.format(cnum1[ical],cnum2[ical])
					line+='\tA-B\t{}\t{}\t{}\t{}\t{}\t{}\tflat{}.fits\twavecal{}.fits'.format(tname,tspt,tbmag,tvmag,str(dpfluxs.loc[0,'PREFIX']),numberList(fnum),cs,cs)

					# if len(fnum)==2: 
					# 	line+='\t{}\t{}\t{:.0f}-{:.0f}\tflat{}.fits\twavecal{}.fits'.format(tname,str(dpfluxs['PREFIX'].iloc[0]),fnum[0],fnum[1],cal_sets.split(',')[ical],cal_sets.split(',')[ical])
					# else:
					# 	telf1 = fnum[np.where(np.abs(fnum-np.roll(fnum,1))>1)]
					# 	telf2 = fnum[np.where(np.abs(fnum-np.roll(fnum,-1))>1)]
					# 	if len(telf1)==0: line+=ftxt
					# 	elif len(telf1)==1: line+='\t{}\t{}\t{:.0f}-{:.0f}\tflat{}.fits\twavecal{}.fits'.format(tname,str(dpfluxs['PREFIX'].iloc[0]),telf1[0],telf2[0],cal_sets.split(',')[ical],cal_sets.split(',')[ical])
					# 	else:
					# 		i = np.argmin(np.abs(fref-telf1))
					# 		line+='\t{}\t{}\t{:.0f}-{:.0f}\tflat{}.fits\twavecal{}.fits'.format(tname,str(dpfluxs['PREFIX'].iloc[i]),telf1[i],telf2[i],cal_sets.split(',')[ical],cal_sets.split(',')[ical])
			f.write(line+'\n')

# print out cal sets
	f.write('\n# Lamp calibrations\n')
	for x in ['FLAT_FILE_PREFIX','ARC_FILE_PREFIX']:
		f.write('{} = {}\n'.format(x,driver_param[x]))
	f.write('CAL_SETS = {}\n'.format(cal_sets[:-1]))

# all other options
		# if len(extraction_options)>0:
		# 	for k in list(extraction_options.keys()): line+='\t{}={}'.format(k,extraction_options[k])

# space for notes

	f.write('\n# Notes\n')
	for i in range(5): f.write('# \n')

# done!
	f.close()

# report out location
	if verbose==True: 
		print('\nBatch instructions written to {}, please this check over before proceeding'.format(driver_file))

# if check=True, read back in and return parameters
	if check==True:
		return readDriver(driver_file,verbose=verbose)

	return


#############################
##### QUALITY ASSURANCE #####
#############################


def makeQApage(driver_input,log_input,image_folder='images',output_folder='',log_html_name='',options={},show_options=False,verbose=ERROR_CHECKING):
	'''
	Purpose
	-------
	Generates the quality assurance page to evaluate a night of reduced data

	Parameters
	----------
	driver_input : str or dict
		either the full path to the batch driver file or the dict output of readDriver()

	log_input : str or pandas DataFrame
		either the full path to the log csv file or the pandas DataFrame output of processFolder()

	output_folder : str, default=''
		by default the page is saved in QA folder; this option specifies the full path 
		to an alternate output folder (e.g. public www folder); where the entire contents 
		of the QA folder is copied; use carefully!

	log_html_name : str, default=''
		full path to the log html file

	options : dict, default = {}
		options for driver parameter that override defaults given in QA_PARAMETERS

	verbose : bool, default=ERROR_CHECKING
		set to True to return verbose output

	Outputs
	-------
	No explicit output, generates an html page for QA review


	Example
	-------
	TBD
		
	Dependencies
	------------
	copy
	numpy
	os.path
	pandas
	readDriver()
	'''
# update default parameters with options
	qa_parameters = copy.deepcopy(QA_PARAMETERS)
	qa_parameters.update(options)

# if necessary read in driver
	if isinstance(driver_input,str)==True:
		if os.path.exists(driver_input)==False: raise ValueError('Cannot find driver file {}'.format(driver_input))
		try: driver = readDriver(driver_input)
		except: raise ValueError('Unable to read in driver file {}'.format(driver_input))
	elif isinstance(driver_input,dict)==True:
		driver = copy.deepcopy(driver_input)
	else: raise ValueError('Do not recognize format of driver input {}'.format(driver_input))
	for x in ['DATA_FOLDER','CALS_FOLDER','PROC_FOLDER','QA_FOLDER']: qa_parameters[x] = driver[x]

# set up image folder
	if image_folder!='':
		image_folder+='/'
		if os.path.exists(os.path.join(qa_parameters['QA_FOLDER'],image_folder))==False:
			try: os.mkdir(os.path.join(qa_parameters['QA_FOLDER'],image_folder))
			except:
				print('WARNING: could not generate image folder {}; reverting to qa folder {}'.format(image_folder,qa_parameters['QA_FOLDER']))
				image_folder = ''


# if necessary read in log
	if isinstance(log_input,str)==True:
		if os.path.exists(log_input)==False: raise ValueError('Cannot find log file {}'.format(log_input))
		if log_input.split('.')[-1] != 'csv': raise ValueError('Log input must be a csv file, you passed {}'.format(log_input))
		obslog = pd.read_csv(log_input)
	elif isinstance(log_input,pd.DataFrame)==True: 
		obslog = copy.deepcopy(log_input)
	else: raise ValueError('Do not recognize format of observing log input {}'.format(log_input))

# show the QA parameters if desired
	if show_options==True: print('\nQA Parameters:\n{}\n'.format(yaml.dump(qa_parameters)))

# set up pyspextool parameters
	ps.pyspextool_setup(driver['INSTRUMENT'])

# read in templates and split index
	for x in ['TEMPLATE_FILE','SOURCE_TEMPLATE_FILE']:
		if os.path.exists(qa_parameters[x])==False: raise ValueError('Cannot find QA html template file {}'.format(qa_parameters[x]))
	with open(qa_parameters['TEMPLATE_FILE'],'r') as f: index_txt = f.read()
	output_text,index_bottom = index_txt.split('[SOURCES]')
	with open(qa_parameters['SOURCE_TEMPLATE_FILE'],'r') as f: source_txt = f.read()
	with open(qa_parameters['SINGLE_PLOT_TEMPLATE_FILE'],'r') as f: single_txt = f.read()

# replace relevant parameters in top of index
	output_text = output_text.replace('[CSS_FILE]',os.path.basename(qa_parameters['CSS_FILE']))
	if log_html_name=='': log_html_name = os.path.basename(log_input).replace('.csv','.html')
	output_text = output_text.replace('[LOG_HTML]',log_html_name)
# NEED TO UPDATE WEATHER LINK
	output_text = output_text.replace('[WEATHER_HTML]',qa_parameters['MKWC_ARCHIVE_URL'])

# replace parameters from driver
#	page_parameters = {}
	for x in list(driver.keys()):
		output_text = output_text.replace('[{}]'.format(x),str(driver[x]))

# now cycle through each of the science sets
	scikeys = list(filter(lambda x: OBSERVATION_SET_KEYWORD in x,list(driver.keys())))
	scikeys.sort()
	imodes = []
	for sci in scikeys:
		stxt = copy.deepcopy(source_txt)
		sci_param = driver[sci]
		imodes.append(sci_param['MODE'])
# data from the driver file		
		for x in list(OBSERVATION_PARAMETERS_REQUIRED.keys()):
			stxt = stxt.replace('['+x+']',str(sci_param[x]))
# target data from the observing log
		sci_log = obslog[obslog['TARGET_NAME']==sci_param['TARGET_NAME']]
		sci_log.reset_index(inplace=True)
		if len(sci_log)>0:
			for x in ['RA','DEC','AIRMASS','SLIT']:
				stxt = stxt.replace('['+x+']',str(sci_log.loc[0,x]))
# nicer coordinates
			s = SkyCoord(sci_log.loc[0,'RA']+' '+sci_log.loc[0,'DEC'], unit=(u.hourangle, u.deg))
			coord = s.to_string('hmsdms',sep=':')
			desig = 'J'+s.to_string('hmsdms',sep='',pad=True).replace('.','').replace(' ','')
			stxt = stxt.replace('[DESIGNATION]',desig)
			if 'fixed' in sci_param['TARGET_TYPE']:
				stxt = stxt.replace('[COORDINATE]',coord)
				stxt = stxt.replace('[TARGET_ALADIN_URL]',\
					'[<a href="https://aladin.cds.unistra.fr/AladinLite/?target={}&fov=0.12&survey=CDS%2FP%2F2MASS%2Fcolor" target="_blank">Aladin</a>]'.format(desig.replace('+','%2B').replace(' ','++')))
				stxt = stxt.replace('[TARGET_SIMBAD_URL]',\
					'[<a href="https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=2&Radius.unit=arcmin&submit=submit+query&CoordList=" target="_blank">SIMBAD</a>]'.format(desig.replace('+','%2B').replace(' ','++')))
				stxt = stxt.replace('[SBDB]','')
			else:
				stxt = stxt.replace('[COORDINATE]','{} at MJD {}'.format(coord,str(sci_log.loc[0,'MJD'])))
				stxt = stxt.replace('[TARGET_ALADIN_URL]','')
				stxt = stxt.replace('[TARGET_SIMBAD_URL]','')
				stxt = stxt.replace('[SBDB]','[<a href="https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr={}" target="_blank">JPL SBDB</a>]'.format(sci_param['TARGET_NAME'].replace(' ','')))
			stxt = stxt.replace('[STD_SIMBAD_URL]',\
				'[<a href="https://simbad.u-strasbg.fr/simbad/sim-id?Ident={}&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id" target="_blank">SIMBAD</a>]'.format(sci_param['STD_NAME'].replace('+','%2B').replace(' ','++')))
			stxt = stxt.replace('[DESIGNATION_URL]',desig.replace('+','%2B'))
			stxt = stxt.replace('[UT_START]',str(sci_log.loc[0,'UT_TIME']).split('.')[0])
			stxt = stxt.replace('[UT_END]',str(sci_log.loc[len(sci_log)-1,'UT_TIME']).split('.')[0])
#			sci_log['TINT'] = [sci_log.loc['INTEGRATION',i]*sci_log.loc['COADDS',i] for i in range(len(sci_log))]
# NOTE: THIS CURRENTLY DOESN'T TAKE INTO ACCOUNT COADDS AS THE COMMAND ABOVE HAD AN ERROR
			stxt = stxt.replace('[INTEGRATION]','{:.1f}'.format(np.nansum(sci_log['INTEGRATION'])))

# std data from the observing log
		sci_log = obslog[obslog['TARGET_NAME']==sci_param['STD_NAME']]
		sci_log.reset_index(inplace=True)
		if len(sci_log)>0:
			stxt = stxt.replace('[STD_AIRMASS]',str(sci_log.loc[0,'AIRMASS']))
			stxt = stxt.replace('[STD_TYPE]','({})'.format(str(sci_log.loc[0,'SIMBAD_TYPE'])))

# final calibrated file (check for presence, otherwise combined file)
		ptxt = copy.deepcopy(qa_parameters['HTML_TABLE_HEAD'])
# fix for distributed file list
		tmp = re.split('[,-]',sci_param['TARGET_FILES'])
		fsuf = '{}-{}'.format(tmp[0],tmp[-1])

#		imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],driver['CALIBRATED_FILE_PREFIX'])+'{}*{}'.format(sci_param['TARGET_FILES'],qa_parameters['PLOT_TYPE']))
		imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}{}'.format(driver['CALIBRATED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
		if len(imfile)==0:
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}{}'.format(driver['CALIBRATED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
#			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
		if len(imfile)>0:
			ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
			fitsfile = os.path.join(qa_parameters['PROC_FOLDER'],os.path.basename(imfile[0]).replace('.pdf','.fits'))
			if os.path.exists(fitsfile)==True:
				ptxt = ptxt.replace('[FITS]','[<a href="{}" download="{}">{}</a>]'.format(fitsfile,os.path.basename(fitsfile),os.path.basename(fitsfile)))
			else: ptxt.replace('[FITS]','')
		else: ptxt=''
#		ptxt+=copy.deepcopy(qa_parameters['HTML_TABLE_TAIL'])
		stxt = stxt.replace('[CALIBRATED_FILE]',ptxt)

# insert other target and calibrator files
		for x in ['TARGET','STD']:	

# combined files
			ptxt = copy.deepcopy(qa_parameters['HTML_TABLE_HEAD'])
# fix for distributed file list
			tmp = re.split('[,-]',sci_param['{}_FILES'.format(x)])
			fsuf = '{}-{}'.format(tmp[0],tmp[-1])
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}*_raw{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)==0:
				imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}*_raw{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)>0:
				ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}*_scaled{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)==0:
				imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}*_scaled{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)>0:
				ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)==0:
				imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}{}'.format(driver['COMBINED_FILE_PREFIX'],fsuf,qa_parameters['PLOT_TYPE'])))
			if len(imfile)>0: 
				ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
			ptxt+=copy.deepcopy(qa_parameters['HTML_TABLE_TAIL'])
			stxt = stxt.replace('[{}_COMBINED_FILES]'.format(x),ptxt)

# individual files
			indexinfo = {'nint': setup.state['nint'], 'prefix': driver['SPECTRA_FILE_PREFIX'],'suffix': '', 'extension': '.pdf'}
			files = make_full_path(qa_parameters['QA_FOLDER'],sci_param['{}_FILES'.format(x)], indexinfo=indexinfo)
			files = list(filter(lambda x: os.path.exists(x),files))
			if len(files)==0:
				files = make_full_path(os.path.join(qa_parameters['QA_FOLDER'],image_folder),sci_param['{}_FILES'.format(x)], indexinfo=indexinfo)
				files = list(filter(lambda x: os.path.exists(x),files))
			files.sort()
			ptxt = copy.deepcopy(qa_parameters['HTML_TABLE_HEAD'])
			for i,f in enumerate(files):
				if i>0 and np.mod(i,qa_parameters['NIMAGES'])==0: ptxt+=' </tr>\n <tr> \n'
				ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(f))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
			ptxt+=copy.deepcopy(qa_parameters['HTML_TABLE_TAIL'])
			stxt = stxt.replace('[{}_SPECTRA_FILES]'.format(x),ptxt)

# traces and apertures
			fnums = extract_filestring(sci_param['{}_FILES'.format(x)],'index')
			ptxt = copy.deepcopy(qa_parameters['HTML_TABLE_HEAD'])
			cnt = 0
			for f in fnums:
				if cnt>0 and np.mod(cnt,qa_parameters['NIMAGES'])==0: ptxt+=' </tr>\n <tr> \n'
				imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}.a_*_trace{}'.format(driver['SCIENCE_FILE_PREFIX'],str(f).zfill(setup.state['nint']),qa_parameters['PLOT_TYPE'])))
				if len(imfile)==0: 
					imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}.a_*_trace{}'.format(driver['SCIENCE_FILE_PREFIX'],str(f).zfill(setup.state['nint']),qa_parameters['PLOT_TYPE'])))
				if len(imfile)>0: 
					imfile.sort()
					ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
					cnt+=1
				if cnt>0 and np.mod(cnt,qa_parameters['NIMAGES'])==0: ptxt+=' </tr>\n <tr> \n'
				imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'{}{}.a_*_apertureparms{}'.format(driver['SCIENCE_FILE_PREFIX'],str(f).zfill(setup.state['nint']),qa_parameters['PLOT_TYPE'])))
				if len(imfile)==0: 
					imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'{}{}.a_*_apertureparms{}'.format(driver['SCIENCE_FILE_PREFIX'],str(f).zfill(setup.state['nint']),qa_parameters['PLOT_TYPE'])))
				if len(imfile)>0: 
					imfile.sort()
					ptxt+=copy.deepcopy(single_txt).replace('[IMAGE]',os.path.join(image_folder,os.path.basename(imfile[0]))).replace('[IMAGE_WIDTH]',str(qa_parameters['IMAGE_WIDTH']))
					cnt+=1
			ptxt+=copy.deepcopy(qa_parameters['HTML_TABLE_TAIL'])
			stxt = stxt.replace('[{}_TRACE_APERTURE_FILES]'.format(x.split('_')[0]),ptxt)

# add it into the page
		output_text+=stxt

# list of modes
	imodes = list(set(imodes))
	s=''
	if len(imodes)>0:
		imodes.sort()
		s = str(imodes[0])
		if len(imodes)>1: 
			for x in imodes[1:]: s+=', {}'.format(str(x))
	output_text = output_text.replace('[INSTRUMENT_MODES]',s)

# insert calibrations
	cal_sets = str(driver['CAL_SETS']).split(',')

# flats
	loopstr = '<table>\n <tr>\n'
	cnt=0
	for cs in cal_sets:
		if cnt>0 and np.mod(cnt,qa_parameters['NIMAGES'])==0: loopstr+=' </tr>\n <tr> \n'
		imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'flat{}_locateorders{}'.format(cs,str(qa_parameters['PLOT_TYPE']))))
		if len(imfile)==0: 
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'flat{}_locateorders{}'.format(cs,str(qa_parameters['PLOT_TYPE']))))
		if len(imfile)>0: 
			imfile.sort()
			if qa_parameters['PLOT_TYPE']=='.pdf': loopstr+='  <td align="center">\n   <embed src="{}" width={} height={}>\n   <br>{}\n   </td>\n'.format(os.path.join(image_folder,os.path.basename(imfile[0])),str(qa_parameters['IMAGE_WIDTH']),str(qa_parameters['IMAGE_WIDTH']),os.path.join(image_folder,os.path.basename(imfile[0])))
			else: loopstr+='  <td align="center">\n   <img src="{}" width={}>\n   <br>{}\n   </td>\n'.format(os.path.join(image_folder,os.path.basename(imfile[0])),qa_parameters['IMAGE_WIDTH'],os.path.basename(os.path.join(image_folder,os.path.basename(imfile[0]))))
			cnt+=1
	loopstr+=' </tr>\n</table>\n\n'
	index_bottom = index_bottom.replace('[FLATS]',loopstr)

# wavecals
# note: only PDF			
	loopstr = '<table>\n <tr>\n'
	cnt=0
	for cs in cal_sets:
		if cnt>0 and np.mod(cnt,qa_parameters['NIMAGES'])==0: loopstr+=' </tr>\n <tr> \n'
		imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'wavecal{}_shift.pdf'.format(cs)))
		if len(imfile)==0: 
			imfile = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],image_folder,'wavecal{}_shift.pdf'.format(cs)))
		if len(imfile)>0: 
			imfile.sort()
			loopstr+='  <td align="center">\n   <embed src="{}" width={} height={}>\n   <br>{}\n   </td>\n'.format(os.path.join(image_folder,os.path.basename(imfile[0])),str(qa_parameters['IMAGE_WIDTH']),str(qa_parameters['IMAGE_WIDTH']),os.path.join(image_folder,os.path.basename(imfile[0])))
			cnt+=1
	loopstr+=' </tr>\n</table>\n'
	index_bottom = index_bottom.replace('[WAVECALS]',loopstr)

# WRITE OUT HTML FILE
	output_text+=index_bottom
	outfile = os.path.join(qa_parameters['QA_FOLDER'],qa_parameters['FILENAME'])
	try: 
		with open(outfile,'w') as f: f.write(output_text)
	except: raise ValueError('Unable to write QA page to {}; check path and permissions'.format(outfile))

# move all the image files into image folder
	if image_folder != '':
		imfiles = glob.glob(os.path.join(qa_parameters['QA_FOLDER'],'*{}'.format(qa_parameters['PLOT_TYPE'])))
		for f in imfiles: 
#			print('\n',qa_parameters['QA_FOLDER'],image_folder,f,os.path.join(qa_parameters['QA_FOLDER'],image_folder,f))
			shutil.move(f,os.path.join(qa_parameters['QA_FOLDER'],image_folder,os.path.basename(f)))

# add in CSS file
	shutil.copy2(qa_parameters['CSS_FILE'],qa_parameters['QA_FOLDER'])

# copy entire tree to a separate output folder if specified
# RIGHT NOW JUST COPIES ENTIRE FOLDER; COULD BE DONE MORE SURGICALLY
	if output_folder!='' and output_folder!=qa_parameters['QA_FOLDER']:
		if os.path.exists(output_folder)==True: 
			try: shutil.copytree(qa_parameters['QA_FOLDER'],output_folder)
			except: 
				print('WARNING: could not copy contents of {} into {}; check path or permissions'.format(qa_parameters['QA_FOLDER'],output_folder))
				output_folder = qa_parameters['QA_FOLDER']
		else: output_folder = qa_parameters['QA_FOLDER']

	if verbose==True: 
		print('\nAcesss QA page at {}\n'.format(os.path.join(qa_parameters['QA_FOLDER'],qa_parameters['FILENAME'])))

	return


#############################
###### BATCH REDUCTION ######
#############################


def batchReduce(parameters,verbose=ERROR_CHECKING):
	'''
	THIS FUNCTION AND DOCSTRING NEED TO BE UPDATED
	Purpose
	-------
	Primary script for conducting a batch reduction of a data folder

	Parameters
	----------
	parameters : dict or str
		either the dict output of readDriver() or the full path to the batch driver file

	verbose : bool, default=ERROR_CHECKING
		set to True to return verbose output

	Outputs
	-------
	If check==True, returns dictionary of batch reduction parameters (output of readDriver)
	otherwise returns nothing; in either case, writes contents to file specified by driver_file


	Example
	-------
	(1) Create from processFolder output:

	>>> from pyspextool.batch import batch
	>>> dpath = '/Users/adam/projects/spex_archive/testing/spex-prism/data/'
	>>> driver_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt'
	>>> dp = batch.processFolder(dpath,verbose=False)
	>>> pars = batch.writeDriver(dp,driver_file=driver_file,data_folder=dpath,create_folders=True,check=True,verbose=True)
			
			Created CALS_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/cals/
			Created PROC_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/proc/
			Created QA_FOLDER folder /Users/adam/projects/spex_archive/testing/spex-prism/qa/

			Batch instructions written to /Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt, please this check over before proceeding

			Driver parameters:
			ARC_FILE_PREFIX: arc
			CALIBRATED_FILE_PREFIX: calspec
			CALS_FOLDER: /Users/adam/projects/spex_archive/testing/spex-prism/cals
			CAL_SETS: 15-20
			...

	(2) Create from previously created log file:
	>>> from pyspextool.batch import batch
	>>> import pandas
	>>> dpath = '/Users/adam/projects/spex_archive/testing/spex-prism/data/'
	>>> log_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/logs.csv'
	>>> driver_file = '/Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt'
	>>> dp = pandas.read_csv(log_file,sep=',')
	>>> pars = batch.writeDriver(dp,driver_file=driver_file,data_folder=dpath,create_folders=True, check=True, verbose=True)

			Batch instructions written to /Users/adam/projects/spex_archive/testing/spex-prism/proc/driver.txt, please this check over before proceeding

			Driver parameters:
			ARC_FILE_PREFIX: arc
			CALIBRATED_FILE_PREFIX: calspec
			CALS_FOLDER: /Users/adam/projects/spex_archive/testing/spex-prism/cals
			CAL_SETS: 15-20
			...
		
	Dependencies
	------------
	copy
	numpy
	os.path
	pandas
	readDriver()
	'''
# default parameters
	if 'VERBOSE' not in list(parameters.keys()): parameters['VERBOSE']=verbose
	parameters['VERBOSE'] = (parameters['VERBOSE'] or verbose)

# set up instrument
	ps.pyspextool_setup(parameters['INSTRUMENT'],raw_path=parameters['DATA_FOLDER'], cal_path=parameters['CALS_FOLDER'], \
		proc_path=parameters['PROC_FOLDER'], qa_path=parameters['QA_FOLDER'],qa_extension=parameters['PLOT_TYPE'],verbose=parameters['VERBOSE'])

# reduce all calibrations
	cal_sets = parameters['CAL_SETS'].split(',')
	if parameters['CALIBRATIONS']==True:
		for cs in cal_sets:

# flats
			indexinfo = {'nint': setup.state['nint'], 'prefix': parameters['FLAT_FILE_PREFIX'],\
				'suffix': '.a', 'extension': '.fits'}
			temp_files = make_full_path(parameters['DATA_FOLDER'],cs, indexinfo=indexinfo)
			input_files = list(filter(lambda x: os.path.exists(x),temp_files))
			if len(input_files)==0: raise ValueError('Cannot find any of the flat files {}'.format(temp_files))
			fnum = [int(os.path.basename(f).replace(parameters['FLAT_FILE_PREFIX'],'').replace('.a.fits','').replace('.b.fits','')) for f in input_files]
			fstr = '{}'.format(numberList(fnum))
# if not present or overwrite, make the file
			if os.path.exists(os.path.join(parameters['CALS_FOLDER'],'flat{}.fits'.format(cs)))==False or parameters['OVERWRITE']==True:
				ps.extract.make_flat([parameters['FLAT_FILE_PREFIX'],fstr],'flat{}'.format(cs),qa_show=parameters['qa_show'],qa_write=parameters['qa_write'],verbose=parameters['VERBOSE'])
			else:
				if parameters['VERBOSE']==True: print('\nflat{}.fits already created, skipping (use --overwrite option to remake)'.format(cs))

# wavecal
			indexinfo = {'nint': setup.state['nint'], 'prefix': parameters['ARC_FILE_PREFIX'],\
				'suffix': '.a', 'extension': '.fits'}
			temp_files = make_full_path(parameters['DATA_FOLDER'],cs, indexinfo=indexinfo)
			input_files = list(filter(lambda x: os.path.exists(x),temp_files))
			if len(input_files)==0: raise ValueError('Cannot find any of the arc files {}'.format(temp_files))
			fnum = [int(os.path.basename(f).replace(parameters['ARC_FILE_PREFIX'],'').replace('.a.fits','').replace('.b.fits','')) for f in input_files]
			fstr = '{}'.format(numberList(fnum))
# if not present or overwrite, make the file
			if os.path.exists(os.path.join(parameters['CALS_FOLDER'],'wavecal{}.fits'.format(cs)))==False or parameters['OVERWRITE']==True:
				ps.extract.make_wavecal([parameters['ARC_FILE_PREFIX'],fstr],'flat{}.fits'.format(cs),'wavecal{}'.format(cs),\
					qa_show=parameters['qa_show'],qa_write=parameters['qa_write'],verbose=parameters['VERBOSE'])
			else:
				if parameters['VERBOSE']==True: print('\nwavecal{}.fits already created, skipping (use --overwrite option to remake)'.format(cs))

# extract all sources and standards
	bkeys = list(filter(lambda x: OBSERVATION_SET_KEYWORD not in x,list(parameters.keys())))
	scikeys = list(filter(lambda x: OBSERVATION_SET_KEYWORD in x,list(parameters.keys())))
	if len(scikeys)==0: 
		if parameters['VERBOSE']==True: print('No science files to reduce')
		return

### REDUCTION ###
	if parameters['REREDUCE']==True: 
		if parameters['VERBOSE']==True: print('NOTE: rereduce is set; will overwrite prior reductions ') 
		parameters['OVERWRITE']=True
	if parameters['EXTRACT']==True:
		for k in scikeys:
			spar = parameters[k]
			for kk in bkeys:
				if kk not in list(spar.keys()): spar[kk] = parameters[kk]


# SCIENCE TARGET
# first check if file is present - skip if not overwriting
# FUTURE: could also add file size: os.path.getsize(file1)>100:
			indexinfo = {'nint': setup.state['nint'], 'prefix': spar['SPECTRA_FILE_PREFIX'], 'suffix': '', 'extension': '.fits'}
			file1 = make_full_path(setup.state['proc_path'], extract_filestring(spar['TARGET_FILES'],'index')[0], indexinfo=indexinfo)[0]
			if os.path.exists(file1)==True and parameters['OVERWRITE']==False:
				if parameters['VERBOSE']==True: print('\n{} already created; skipping set {} (use --overwrite option to reextract)'.format(file1,spar['TARGET_FILES']))
				continue

# reduce the first pair of targets
			ps.extract.load_image([spar['TARGET_PREFIX'],extract_filestring(spar['TARGET_FILES'],'index')[:2]],\
				spar['TARGET_FLAT_FILE'], spar['TARGET_WAVECAL_FILE'],reduction_mode=spar['TARGET_REDUCTION_MODE'], \
				flat_field=True, linearity_correction=True,qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# set extraction method
			ps.extract.set_extraction_type(spar['TARGET_TYPE'].split('-')[-1])

# make spatial profiles
			ps.extract.make_spatial_profiles(qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# NEED A CHECK IN HERE TO STOP EXTRACTION OF THE SPATIAL PROFILE IS MOSTLY FLAT (PEAKS < 3 * BACK STDEV)
# HOW DO WE ACCESS THE TRACE?

# identify aperture positions
			if spar['TARGET_APERTURE_METHOD']=='auto':
				ps.extract.locate_aperture_positions(spar['TARGET_NPOSITIONS'], method=spar['TARGET_APERTURE_METHOD'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])
			else:
				ps.extract.locate_aperture_positions(spar['TARGET_APERTURE_POSITIONS'], method=spar['TARGET_APERTURE_METHOD'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# select orders to extract (set above)
			ps.extract.select_orders(include=spar['TARGET_ORDERS'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# trace apertures
			ps.extract.trace_apertures(qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# define the aperture - psf, width, background
# NOTE - ONLY WORKS FOR POINT SOURCE, NEED TO FIX FOR XS?
			ps.extract.define_aperture_parameters(spar['TARGET_APERTURE'], psf_radius=spar['TARGET_PSF_RADIUS'],bg_radius=spar['TARGET_BACKGROUND_RADIUS'],\
						bg_width=spar['TARGET_BACKGROUND_WIDTH'],qa_show=spar['qa_show'],qa_write=spar['qa_write'],)

# extract away
			ps.extract.extract_apertures(verbose=spar['VERBOSE'])

# conduct extraction of all remaining files
			fnum = extract_filestring(spar['TARGET_FILES'],'index')[2:]
			if len(fnum)>0:
				ps.extract.do_all_steps([spar['TARGET_PREFIX'],'{}'.format(numberList(fnum))],verbose=spar['VERBOSE'])

# FLUX STANDARD
# first check if file is present - skip if not overwriting
# FUTURE: could also add file size: os.path.getsize(file1)>100:
			indexinfo = {'nint': setup.state['nint'], 'prefix': spar['SPECTRA_FILE_PREFIX'], 'suffix': '', 'extension': '.fits'}
			file1 = make_full_path(setup.state['proc_path'], extract_filestring(spar['STD_FILES'],'index')[0], indexinfo=indexinfo)[0]
			if os.path.exists(file1)==True and parameters['OVERWRITE']==False:
				if parameters['VERBOSE']==True: print('\n{} already created; skipping set {} (use --overwrite option to reextract)'.format(file1,spar['STD_FILES']))
				continue

# now reduce the first pair of standards
# NOTE: ASSUMING THESE ARE BRIGHT POINT SOURCES WITH AUTO APERTURE FINDING
			ps.extract.load_image([spar['STD_PREFIX'],extract_filestring(spar['STD_FILES'],'index')[:2]],\
				spar['STD_FLAT_FILE'], spar['STD_WAVECAL_FILE'],reduction_mode=spar['STD_REDUCTION_MODE'], \
				flat_field=True, linearity_correction=True,qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# set extraction method
			ps.extract.set_extraction_type('ps')

# make spatial profiles
			ps.extract.make_spatial_profiles(qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# identify aperture positions
			# if spar['STD_APERTURE_METHOD']=='auto':
			# 	ps.extract.locate_aperture_positions(spar['STD_NPOSITIONS'], method=spar['STD_APERTURE_METHOD'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])
			# else:
			# 	ps.extract.locate_aperture_positions(spar['STD_APERTURE_POSITIONS'], method=spar['STD_APERTURE_METHOD'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])
			ps.extract.locate_aperture_positions(spar['STD_NPOSITIONS'], method='auto', qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=parameters['VERBOSE'])

# select orders to extract (set above)
			ps.extract.select_orders(include=spar['STD_ORDERS'], qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# trace apertures
			ps.extract.trace_apertures(qa_show=spar['qa_show'],qa_write=spar['qa_write'], verbose=spar['VERBOSE'])

# define the aperture - psf, width, background
# NOTE - ONLY WORKS FOR POINT SOURCE, NEED TO FIX FOR XS?
			ps.extract.define_aperture_parameters(spar['STD_APERTURE'], psf_radius=spar['STD_PSF_RADIUS'],bg_radius=spar['STD_BACKGROUND_RADIUS'],\
						bg_width=spar['STD_BACKGROUND_WIDTH'], qa_show=spar['qa_show'],qa_write=spar['qa_write'])

# extract away
			ps.extract.extract_apertures(verbose=verbose)

# conduct extraction of all remaining files
			fnum = extract_filestring(spar['STD_FILES'],'index')[2:]
			if len(fnum)>0:
				ps.extract.do_all_steps([spar['STD_PREFIX'],'{}'.format(numberList(fnum))],verbose=spar['VERBOSE'])

## COMBINE SPECTRAL FILES ###
# NOTE: SPECTRA_FILE_PREFIX is currently hardcoded in code to be 'spectra' - how to change?
	if parameters['COMBINE']==True:
		for k in scikeys:
			spar = parameters[k]
			for kk in bkeys:
				if kk not in list(spar.keys()): spar[kk] = parameters[kk]

#			if spar['MODE'] in ['SXD','ShortXD']:

# fix for distributed file list
			tmp = re.split('[,-]',spar['TARGET_FILES'])
			fsuf = '{}-{}'.format(tmp[0],tmp[-1])
			outfile = '{}{}'.format(spar['COMBINED_FILE_PREFIX'],fsuf)

# first check if file is present - skip if not overwriting
			if os.path.exists(os.path.join(parameters['PROC_FOLDER'],outfile))==True and parameters['OVERWRITE']==False:
				if parameters['VERBOSE']==True: print('\n{}{}.fits already created, skipping (use --overwrite option to remake)'.format(spar['COMBINED_FILE_PREFIX'],fsuf))			
# combine
			else:
				ps.combine.combine_spectra(['spectra',spar['TARGET_FILES']],outfile,
					scale_spectra=True,scale_range=spar['TARGET_SCALE_RANGE'],correct_spectral_shape=False,qa_show=spar['qa_show'],qa_write=spar['qa_write'],verbose=spar['VERBOSE'])

# telluric star
#		if spar['MODE']=='SXD':

# fix for distributed file list
			tmp = re.split('[,-]',spar['STD_FILES'])
			fsuf = '{}-{}'.format(tmp[0],tmp[-1])
			outfile = '{}{}'.format(spar['COMBINED_FILE_PREFIX'],fsuf)

# first check if file is present - skip if not overwriting
			if os.path.exists(os.path.join(parameters['PROC_FOLDER'],outfile))==True and parameters['OVERWRITE']==False:
				if parameters['VERBOSE']==True: print('\n{}{}.fits already created, skipping (use --overwrite option to remake)'.format(spar['COMBINED_FILE_PREFIX'],fsuf))
# combine
			else:
				ps.combine.combine_spectra(['spectra',spar['STD_FILES']],outfile,
					scale_spectra=True,scale_range=spar['STD_SCALE_RANGE'],correct_spectral_shape=True,qa_show=spar['qa_show'],qa_write=spar['qa_write'],verbose=spar['VERBOSE'])

## FLUX CALIBRATE ###
	if parameters['FLUXTELL']==True:
		for k in scikeys:
			spar = parameters[k]
			for kk in bkeys:
				if kk not in list(spar.keys()): spar[kk] = parameters[kk]

			standard_data = {'id': spar['STD_NAME'], 'sptype': spar['STD_SPT'],'bmag': float(spar['STD_B']),'vmag': float(spar['STD_V'])}
#			print(standard_data)

# fix for distributed file list
			tmp = re.split('[,-]',spar['TARGET_FILES'])
			tsuf = '{}-{}'.format(tmp[0],tmp[-1])
			tmp = re.split('[,-]',spar['STD_FILES'])
			csuf = '{}-{}'.format(tmp[0],tmp[-1])
			outfile = '{}{}.fits'.format(spar['CALIBRATED_FILE_PREFIX'],tsuf)

# first check if file is present - skip if not overwriting <-- NOT NEEDED
			# if os.path.exists(outfile)==True and parameters['OVERWRITE']==False:
			# 	if parameters['VERBOSE']==True: print('\n{}{}.fits already created, skipping (use --overwrite option to remake)'.format(spar['CALIBRATED_FILE_PREFIX'],tsuf))				
			# else:

# prep file names
			objfile = '{}{}.fits'.format(spar['COMBINED_FILE_PREFIX'],tsuf)
			stdfile = '{}{}.fits'.format(spar['COMBINED_FILE_PREFIX'],csuf)

# depends on fixed or moving source which method we use
			if (spar['TARGET_TYPE'].split('-')[0]).strip()=='fixed': 
				if spar['VERBOSE']==True: print('\nfixed target; doing standard telluric correction and flux calibration')
				reflectance = False
			elif (spar['TARGET_TYPE'].split('-')[0]).strip()=='moving': 
				if spar['VERBOSE']==True: print('\nmoving target; doing reflectance telluric correction and flux calibration assuming G2 V standard')
				reflectance = True
			else: 
				if spar['VERBOSE']==True: print('\nTarget type {} not recognized, skipping flux/telluric correction'.format(spar['TARGET_TYPE'].split('-')[0]))			
				continue

# NOTE: NEED TO PARAMETERIZE DEFAULTS FOR WRITE TELLURIC AND WRITE MODEL
			ps.telluric.telluric_correction(objfile,stdfile,standard_data,outfile,reflectance=reflectance,
				write_telluric_spectra=True,write_model_spectra=True,
				qa_show=spar['qa_show'],qa_write=spar['qa_write'],verbose=spar['VERBOSE'],overwrite=spar['OVERWRITE'])				

## ORDER STITCHING ###
# NOT YET IMPLEMENTED  #
	if parameters['STITCH']==True:
		for k in scikeys:
			spar = parameters[k]
			for kk in bkeys:
				if kk not in list(spar.keys()): spar[kk] = parameters[kk]

			if spar['MODE'] not in ['Prism','LowRes15']:
				if spar['VERBOSE']==True: 
					print('\n** NOTE: STITCHING HAS NOT BEEN IMPLEMENTED **')

	return



# external function call
if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('Call this function from the command line as python batch.py [data_folder] [cals_folder] [proc_folder]')
	else:
		dp = processFolder(sys.argv[1])
		dp = processFolder(sys.argv[1])
		if len(sys.argv) > 4: log_file = sys.argv[4]
		else: log_file = sys.argv[3]+'log.html'
		print('\nWriting log to {}'.format(log_file))
		writeLog(dp,log_file)
		if len(sys.argv) > 5: driver_file = sys.argv[5]
		else: driver_file = sys.argv[3]+'driver.txt'
		print('\nWriting batch driver file to {}'.format(driver_file))
		writeDriver(dp,driver_file,data_folder=sys.argv[1],cals_folder=sys.argv[2],proc_folder=sys.argv[3])
		txt = input('\nCheck the driver file {} and press return when ready to proceed...\n')
		par = readDriver(driver_file)
		print('\nReducing spectra')
		batchReduce(par)
		print('\nReduction complete: processed files are in {}'.format(sys.argv[3]))
