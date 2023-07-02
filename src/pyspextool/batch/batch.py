# Functions for batch reduction using pyspextool
# General procedure:
# - read in files in a given folder and process into a log file and batch driving file
#  -- subprocesses: files --> database, database --> log, database --> batch driver
# - user checks logs and batch driving file to ensure they are correct
# - run through reduction steps based on batch driving file
#  -- subprocesses: read driver --> do all calibrations --> extract all spectra 
#		--> combine spectral groups --> telluric correction --> stitch (SXD/LXD)

import copy
import os
import glob
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

import pyspextool as ps
from pyspextool import config as setup
from pyspextool.io.files import extract_filestring,make_full_path

# this is the data to extract from header, with options for header keywords depending on epoch
HEADER_DATA={
	'UT_DATE': ['DATE_OBS'],
	'UT_TIME': ['TIME_OBS'],
#	'MJD': ['MJD_OBS'],
	'SOURCE_NAME': ['OBJECT','TCS_OBJ'],
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
	'TYPE': ['DATATYPE'],
}
ERROR_CHECKING = True
MOVING_MAXSEP = 15 # max separation for fixed target in arcsec
MOVING_MAXRATE = 10 # max moving rate for fixed target in arcsec/hr
LOG_COLUMNS = ['FILE','BEAM','SOURCE_NAME','UT_DATE','UT_TIME','RA','DEC','TYPE','FIXED-MOVING','AIRMASS','INTEGRATION','COADDS','SLIT','MODE','PROGRAM','OBSERVER']
LOG_TEMPLATE_HTML = r'<!DOCTYPE html><html lang="en"><link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css"><head><title>[TITLE]</title></head><body><h1>[TITLE]</h1>[TABLE]</body></html>'
INSTRUMENT_DATE_SHIFT = Time('2014-07-01').mjd # date spex --> uspex
ARC_SRCNAME = 'arc lamp'
FLAT_SRCNAME = 'flat lamp'
BATCH_PARAMETERS = {
	'INSTRUMENT': '',
	'DATA_FOLDER': './',
	'CAL_FOLDER': './',
	'PROC_FOLDER': './',
	'QA_FOLDER': '',
	'ARC_FILE_PREFIX': 'arc',
	'FLAT_FILE_PREFIX': 'flat',
	'CAL_SETS': '',
}
SCI_PARAMETERS_REQUIRED = ['MODE','TARGET_TYPE','TARGET_NAME','TARGET_FILES','FLAT_FILE','WAVECAL_FILE','STD_NAME','STD_FILES']
SCI_PARAMETERS_OPTIONAL = {
	'SCIENCE_FILE_PREFIX': 'spc',
	'SPECTRAL_FILE_PREFIX': 'spec',
	'COMBINED_FILE_PREFIX': 'combspec',
	'CALIBRATED_FILE_PREFIX': 'calspec',
	'STITCHED_FILE_PREFIX': 'xd1d',	
	'USE_STORED_SOLUTION': False,
	'ORDERS': '3-9',
	'SOURCE_TYPE': 'ps',
	'REDUCTION_MODE': 'A-B',
	'NPOSITIONS': 2,
	'APERTURE_POSITIONS': [3.7,11.2],
	'APERTURE_METHOD': 'auto',
	'PS_APERTURE': 1.5,
	'PSF_RADIUS': 1.5,
	'BACKGROUND_RADIUS': 2.5,
	'BACKGROUND_WIDTH': 4,
	'SCALE_RANGE': [1.0,1.5],
}
QA_PLOT_TYPE = '.pdf'
DIR = os.path.dirname(os.path.abspath(__file__))
QA_CSS = os.path.join(DIR,'qa.css')
QA_INDEX_TEMPLATE_FILE = os.path.join(DIR,'qa_template.txt')
QA_SOURCE_TEMPLATE_FILE = os.path.join(DIR,'qa_source_template.txt')
LOG_TEMPLATE_FILE = os.path.join(DIR,'log_template.txt')
QA_HTML_OUTPUT_NAME = 'index.html'
MKWC_ARCHIVE_URL = 'http://mkwc.ifa.hawaii.edu/forecast/mko/archive/index.cgi'

def tablestyle(val):
	'''
	Plotting style for pandas dataframe for html display
	'''
	color='white'
	if type(val)==str:
		if 'LowRes' in val: color='blue'
		if 'SXD' in val: color='green'
		if 'LXD' in val: color='red'
	return 'background-color: {}'.format(color)


def process_folder(folder,verbose=False):
	'''
	Function reads in fits files from a data folder and organizes into 
	pandas dataframe
	'''
# error checking folder 
	if os.path.isdir(folder) == False: 
		raise ValueError('Cannot folder find {} or {}'.format(inp,folder))
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
		hdu = fits.open(f)
		hdu[0].verify('silentfix')
		header = hdu[0].header
		hdu.close()
		for k in list(HEADER_DATA.keys()): 
			ref = ''
			for ii in HEADER_DATA[k]:
				if ii in list(header.keys()) and ref=='': ref=ii
			if ref!='': dp[k].iloc[i] = header[ref]
			if ref=='' and verbose==True and i==0:
				print('Could not find keywords {} in file {}'.format(HEADER_DATA[k],f))
# update some of the mode names
		if 'SHORTXD' in dp['MODE'].iloc[i].upper(): dp['MODE'].iloc[i] = 'SXD'
		if 'lowres' in dp['MODE'].iloc[i].lower(): dp['MODE'].iloc[i] = 'Prism'		
		if 'arc' in dp['FILE'].iloc[i]: 
			dp['SOURCE_NAME'].iloc[i] = ARC_SRCNAME
			dp['TYPE'].iloc[i] = 'calibration'
		if 'flat' in dp['FILE'].iloc[i]: 
			dp['SOURCE_NAME'].iloc[i] = FLAT_SRCNAME
			dp['TYPE'].iloc[i] = 'calibration'
# set standards based on integration time
# this is not a great way to do this!
		if dp['TYPE'].iloc[i]=='': 
			dp['TYPE'].iloc[i]='target'
			if 'SXD' in dp['MODE'].iloc[i] and float(dp['INTEGRATION'].iloc[i])<=30.: dp['TYPE'].iloc[i]='standard'
			if 'Prism' in dp['MODE'].iloc[i] and float(dp['INTEGRATION'].iloc[i])<=10.: dp['TYPE'].iloc[i]='standard'

# if guess above was wrong, reset all standards back to targets
	dpt = dp[dp['TYPE']=='target']
	if len(dpt)==0: 
		dp.loc[dp['TYPE']=='standard','TYPE'] = 'target'
	dpt = dp[dp['TYPE']=='standard']
	if len(dpt)==0: 
		if verbose==True: print('Warning: no standards present in list; be sure to check driver file carefully')

# generate ISO 8601 time string and sort
# --> might need to do some work on this
	dp['DATETIME'] = [dp['UT_DATE'].iloc[i]+'T'+dp['UT_TIME'].iloc[i] for i in range(len(dp))]
	dp['MJD'] = [Time(d,format='isot', scale='utc').mjd for d in dp['DATETIME']]
	dp['INSTRUMENT'] = ['spex']*len(dp)
	if dp['MJD'].iloc[0]>INSTRUMENT_DATE_SHIFT: dp['INSTRUMENT'] = ['uspex']*len(dp)
	dp.sort_values('DATETIME',inplace=True)
	dp.reset_index(inplace=True,drop=True)

# fixed/moving - determined by range of position and motion of soruce
# NOTE: THIS FAILS IF USER DOESN'T CHANGE NAME OF SOURCE
	dp['FIXED-MOVING'] = ['fixed']*len(dp)
	names = list(set(list(dp['SOURCE_NAME'])))
	names.remove(ARC_SRCNAME)
	names.remove(FLAT_SRCNAME)
	if len(names)==0:
		if verbose==True: print('Warning: no science files identified in {}'.format(folder))
	else:
		for n in names:
			dps = dp[dp['SOURCE_NAME']==n]
			dpsa = dps[dps['BEAM']=='A']
			if len(dpsa)>1:	
				pos1 = SkyCoord(dpsa['RA'].iloc[0]+' '+dpsa['DEC'].iloc[0],unit=(u.hourangle, u.deg))
				pos2 = SkyCoord(dpsa['RA'].iloc[-1]+' '+dpsa['DEC'].iloc[-1],unit=(u.hourangle, u.deg))
				dx = pos1.separation(pos2).to(u.arcsec)
				time1 = Time(dpsa['DATETIME'].iloc[0],format='isot', scale='utc')
				time2 = Time(dpsa['DATETIME'].iloc[-1],format='isot', scale='utc')
				dt = (time2-time1).to(u.hour)
				if dx.value > MOVING_MAXSEP or ((dx/dt).value>MOVING_MAXRATE and dt.value > 0.1):
					dp.loc[dp['SOURCE_NAME']==n,'FIXED-MOVING'] = 'moving'
				if verbose==True: print('{}: dx={:.1f} arc, dt={:.1f} hr, pm={:.2f} arc/hr = {}'.format(n,dx,dt,dx/dt,dp.loc[dp['SOURCE_NAME']==n,'FIXED-MOVING']))

	return dp


def write_log(dp,log_file,columns=LOG_COLUMNS,html_template=LOG_TEMPLATE_HTML,verbose=ERROR_CHECKING):
	'''
	Writes log file based on data Frame from process_folder
	'''
# downselect columns to save out
	try: dpout = dp[columns]
	except: 
		dpout=copy.deepcopy(dp)
		if verbose==True: print('Warning: could not select subset of columns {} from data frame, saving all columns'.format(columns))

# save depends on file name
	ftype = log_file.split('.')[-1]
	if ftype in ['xls','xlsx']: dpout.to_excel(log_file,index=False)	
	elif ftype in ['csv']: dpout.to_csv(log_file,sep=',',index=False)	
	elif ftype in ['tsv','txt']: dpout.to_csv(log_file,sep='\t',index=False)	
	elif ftype in ['tex']: dpout.to_latex(log_file,longtable=True,index=False)	
	elif ftype in ['htm','html']: 
		dpout.style.apply(tablestyle)
		dphtml = copy.deepcopy(html_template)
		dphtml = dphtml.replace('[TITLE]','{} Log for {}'.format(dp['INSTRUMENT'].iloc[0],dp['UT_DATE'].iloc[0]))
		dphtml = dphtml.replace('[TABLE]',dpout.to_html(classes='table table-striped text-center',index=False,bold_rows=True))
		with open(log_file,'w') as f: f.write(dphtml)
	else: raise ValueError('Could not write out to {}; unknown file format'.format(log_file))

	if verbose==True: print('log written to {}'.format(log_file))
	return


def write_driver(dp,driver_file,data_folder=BATCH_PARAMETERS['DATA_FOLDER'],cal_folder=BATCH_PARAMETERS['CAL_FOLDER'],\
	proc_folder=BATCH_PARAMETERS['PROC_FOLDER'],qa_folder=BATCH_PARAMETERS['QA_FOLDER'],\
	science_file_prefix='',arc_file_prefix='',flat_file_prefix='',\
	combined_file_prefix=SCI_PARAMETERS_OPTIONAL['COMBINED_FILE_PREFIX'],\
	spectral_file_prefix=SCI_PARAMETERS_OPTIONAL['SPECTRAL_FILE_PREFIX'],\
	calibrated_file_prefix=SCI_PARAMETERS_OPTIONAL['CALIBRATED_FILE_PREFIX'],\
	stitched_file_prefix=SCI_PARAMETERS_OPTIONAL['STITCHED_FILE_PREFIX'],\
	default_source_type=SCI_PARAMETERS_OPTIONAL['SOURCE_TYPE'],\
	extraction_options={},comment='',create_folders=False,verbose=ERROR_CHECKING):
	''']
	Writes driver file based on database of files
	'''
# check folders
	if qa_folder=='': qa_folder = cal_folder
	for x in [data_folder,cal_folder,proc_folder,qa_folder]:
		if os.path.exists(x)==False:
			if create_folders==True: os.mkdir(x)
			else: raise ValueError('Folder {} does not exist; either create or set keyword create_folders to True'.format(x))

# some file name work
	dpc = copy.deepcopy(dp)
	dpc['PREFIX'] = [x.replace('.a.fits','').replace('.b.fits','') for x in dpc['FILE']]
	n=4
	if dpc['INSTRUMENT'].iloc[0]=='uspex': n=5
	dpc['FILE NUMBER'] = [int(x[-n:]) for x in dpc['PREFIX']]
	dpc['PREFIX'] = [x[:-n] for x in dpc['PREFIX']]

# set up default prefixes
	if arc_file_prefix=='':
		dps = dpc[dpc['SOURCE_NAME']==ARC_SRCNAME]
		arc_file_prefix = dps['PREFIX'].iloc[0]
	if flat_file_prefix=='':
		dps = dpc[dpc['SOURCE_NAME']==FLAT_SRCNAME]
		flat_file_prefix = dps['PREFIX'].iloc[0]
	if science_file_prefix=='':
		dps = dpc[dpc['SOURCE_NAME']!=FLAT_SRCNAME]
		dps = dps[dps['SOURCE_NAME']!=ARC_SRCNAME]
		science_file_prefix = dps['PREFIX'].iloc[0]


# write out instructions
	f = open(driver_file,'w')
	f.write('# Batch reduction driver file for {} observations on {}\n'.format(dpc['INSTRUMENT'].iloc[0],dpc['UT_DATE'].iloc[0]))
	if comment!='': f.write('# {}\n'.format(comment))

# folders
	f.write('\n# File and folder information\n')
	f.write('INSTRUMENT = {}\n'.format(dpc['INSTRUMENT'].iloc[0]))
	f.write('DATA_FOLDER = {}\n'.format(os.path.abspath(data_folder)+'/'))
	f.write('CAL_FOLDER = {}\n'.format(os.path.abspath(cal_folder)+'/'))
	f.write('PROC_FOLDER = {}\n'.format(os.path.abspath(proc_folder)+'/'))
	f.write('QA_FOLDER = {}\n'.format(os.path.abspath(qa_folder)+'/'))

# file prefixes
	f.write('SCIENCE_FILE_PREFIX = {}\n'.format(science_file_prefix))
	f.write('SPECTRAL_FILE_PREFIX = {}\n'.format(spectral_file_prefix))
	f.write('COMBINED_FILE_PREFIX = {}\n'.format(combined_file_prefix))
	f.write('CALIBRATED_FILE_PREFIX = {}\n'.format(calibrated_file_prefix))
	f.write('STITCHED_FILE_PREFIX = {}\n'.format(stitched_file_prefix))

# cals
	f.write('\n# Lamp calibrations\n')
	f.write('ARC_FILE_PREFIX = {}\n'.format(arc_file_prefix))
	f.write('FLAT_FILE_PREFIX = {}\n'.format(flat_file_prefix))
	cal_sets=''
	dpcal = dpc[dpc['TYPE']=='calibration']
	fnum = np.array(dpcal['FILE NUMBER'])
	calf1 = fnum[np.where(np.abs(fnum-np.roll(fnum,1))>1)]
	calf2 = fnum[np.where(np.abs(fnum-np.roll(fnum,-1))>1)]
	for i in range(len(calf1)): cal_sets+='{:.0f}-{:.0f},'.format(calf1[i],calf2[i])
	f.write('CAL_SETS = {}\n'.format(cal_sets[:-1]))

# observations
	dps = dpc[dpc['TYPE']=='target']
	names = list(set(list(dps['SOURCE_NAME'])))
	names.sort()

# no names - close it up
	if len(names)==0:
		if verbose==True: print('Warning: no science files identified')
		f.close()
		return

# loop over names
	f.write('\n# Science observations')
	f.write('\n# \tMode\tTarget Type\tTarget Name\tFiles\tFlat Filename\tWavecal Filename\tStd Name\tStd Files\tOptions (key=value)\n')
	for n in names:
		dpsrc = dpc[dpc['SOURCE_NAME']==n]
		line='SCI\t{}\t{} {}\t{}'.format(dpsrc['MODE'].iloc[0],dpsrc['FIXED-MOVING'].iloc[0],default_source_type,n)
		fnum = np.array(dpsrc['FILE NUMBER'])
		line+='\t{:.0f}-{:.0f}'.format(np.nanmin(fnum),np.nanmax(fnum))

# assign calibrations
		dpcals = dpcal[dpcal['MODE']==dpsrc['MODE'].iloc[0]]
		if len(dpcals)==0: 
			if verbose==True: print('Warning: no calibration files associated with mode {} for source {}'.format(dps['MODE'].iloc[0],n))
			i=0
		else:
			i = np.argmin(np.abs(calf1-np.median(dpsrc['FILE NUMBER'])))
		line+='\tflat{}.fits\twavecal{}.fits'.format(cal_sets.split(',')[i],cal_sets.split(',')[i])

# assign flux cals
		dpflux = dpc[dpc['TYPE']=='standard']
		dpflux = dpflux[dpflux['MODE']==dpsrc['MODE'].iloc[0]]
		if len(dpflux)==0: 
			if verbose==True: print('Warning: no calibration files associated with mode {} for source {}'.format(dps['MODE'].iloc[0],n))
			line+='\tUNKNOWN\t0-0'
		else:
			dpflux['DIFF1'] = [np.abs(x-np.nanmedian(dpsrc['AIRMASS'])) for x in dpflux['AIRMASS']]
			dpflux['DIFF2'] = [np.abs(x-np.nanmedian(dpsrc['MJD'])) for x in dpflux['MJD']]
			dpflux['DIFF'] = dpflux['DIFF1']+dpflux['DIFF2']
			fref = dpflux['FILE NUMBER'].iloc[np.argmin(dpflux['DIFF'])]
			tname = dpflux['SOURCE_NAME'].iloc[np.argmin(dpflux['DIFF'])]
			dpfluxs = dpflux[dpflux['SOURCE_NAME']==tname]
			fnum = np.array(dpfluxs['FILE NUMBER'])
			if len(fnum)==2: line+='\t{}\t{:.0f}-{:.0f}'.format(tname,fnum[0],fnum[1])
			elif len(fnum)<2: 
				if verbose==True: print('Fewer than 2 flux calibrator files for source {}'.format(tname))
				line+='\tUNKNOWN\t0-0'
			else:
				telf1 = fnum[np.where(np.abs(fnum-np.roll(fnum,1))>1)]
				telf2 = fnum[np.where(np.abs(fnum-np.roll(fnum,-1))>1)]
				if len(telf1)==0: line+='\tUNKNOWN\t0-0'
				elif len(telf1)==1: line+='\t{}\t{:.0f}-{:.0f}'.format(tname,telf1[0],telf2[0])
				else:
					i = np.argmin(np.abs(fref-telf1))
					line+='\t{}\t{:.0f}-{:.0f}'.format(tname,telf1[i],telf2[i])

# all other options
		if len(extraction_options)>0:
			for k in list(extraction_options.keys()): line+='\t{}={}'.format(k,extraction_options[k])
		f.write(line+'\n')
	f.close()

	if verbose==True: print('Batch instructions written to {}, please this check over before proceeding'.format(driver_file))
	return


def read_driver(driver_file,default_parameters=BATCH_PARAMETERS,verbose=ERROR_CHECKING):
	'''
	Reads in driver file in dictionary for processing
	'''
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

# extract of keyword = value lines and use to update default parameters
	klines = list(filter(lambda x: ('=' in x), lines))
	klines = [l.replace(' ','').replace('\n','') for l in klines]
	new_parameters = dict([l.strip().split('=') for l in klines])
	parameters = copy.deepcopy(BATCH_PARAMETERS)
	parameters.update(new_parameters)

# set the individual science sets
	slines = list(filter(lambda x: (x[:4]=='SCI\t'), lines))
	for i,sline in enumerate(slines):
		spar = copy.deepcopy(SCI_PARAMETERS_OPTIONAL)
		pars = sline.replace('\n','').split('\t')
		if len(pars) < len(SCI_PARAMETERS_REQUIRED):
			if verbose==True: print('Warning: line {} contains fewer than the required parameters {}; skipping'.format(sline,SCI_REQUIRED))
		else:
# required parameters			
			for ii,k in enumerate(SCI_PARAMETERS_REQUIRED):
				spar[k] = pars[ii+1]
# fill in global parameters			
			for k in list(spar.keys()):
				if k in list(parameters.keys()):
					if parameters[k]!='': spar[k] = parameters[k]
# add optional parameters			
			ksline = list(filter(lambda x: ('=' in x), sline))
			ksline = [l.replace(' ','').replace('\n','') for l in ksline]
			new_spar = dict([l.strip().split('=') for l in ksline])
			spar.update(new_spar)
# some instrument specific adaptations
			if 'Prism' in spar['MODE']:
				spar['ORDERS'] = '1'
#				if verbose==True: print('Updated prism mode orders to 1')
			if parameters['INSTRUMENT'] == 'spex' and 'SXD' in spar['MODE'] and spar['ORDERS'][-1]=='9':
				spar['ORDERS'] = spar['ORDERS'][:-1]+'8'
				if verbose==True: print('Updated spex SXD mode orders to {}'.format(spar['ORDERS']))

			parameters['SCI{}'.format(str(i+1).zfill(3))] = spar

	return parameters


def batch_reduce(parameters,qa_plot=True,qa_file=True,verbose=ERROR_CHECKING):
	'''
	Reduces dataset based on instructions in driver file
	'''
# set up instrument
	ps.pyspextool_setup(parameters['INSTRUMENT'],raw_path=parameters['DATA_FOLDER'], cal_path=parameters['CAL_FOLDER'], proc_path=parameters['PROC_FOLDER'], qa_path=parameters['QA_FOLDER'],qa_extension=QA_PLOT_TYPE,verbose=verbose)

# reduce all calibrations
	cal_sets = parameters['CAL_SETS'].split(',')
	for cs in cal_sets:

# flats
		indexinfo = {'nint': setup.state['nint'], 'prefix': parameters['FLAT_FILE_PREFIX'],\
			'suffix': '.a', 'extension': '.fits'}
		temp_files = make_full_path(parameters['DATA_FOLDER'],cs, indexinfo=indexinfo)
		input_files = list(filter(lambda x: os.path.exists(x),temp_files))
		if len(input_files)==0: raise ValueError('Cannot find any of the flat files {}'.format(temp_files))
		fnum = [int(os.path.basename(f).replace(parameters['FLAT_FILE_PREFIX'],'').replace('.a.fits','').replace('.b.fits','')) for f in input_files]
		fstr = '{:.0f}-{:.0f}'.format(np.nanmin(fnum),np.nanmax(fnum))
# NOTE: qa_file force to False due to ploting error in plot_image
		ps.extract.make_flat([parameters['FLAT_FILE_PREFIX'],fstr],'flat{}'.format(cs),qa_plot=qa_plot,qa_file=True,verbose=verbose)

# wavecal
		indexinfo = {'nint': setup.state['nint'], 'prefix': parameters['ARC_FILE_PREFIX'],\
			'suffix': '.a', 'extension': '.fits'}
		temp_files = make_full_path(parameters['DATA_FOLDER'],cs, indexinfo=indexinfo)
		input_files = list(filter(lambda x: os.path.exists(x),temp_files))
		if len(input_files)==0: raise ValueError('Cannot find any of the arc files {}'.format(temp_files))
		fnum = [int(os.path.basename(f).replace(parameters['ARC_FILE_PREFIX'],'').replace('.a.fits','').replace('.b.fits','')) for f in input_files]
		fstr = '{:.0f}-{:.0f}'.format(np.nanmin(fnum),np.nanmax(fnum))
# HARD CODED choice of use_stored_solution for spex prism mode
# needs to be done on a case-by-case basis due to potentially mixed calibrations
# would be better to move this directly into ps.extrac.make_wavecal function
		hdu = fits.open(input_files[0])
		hdu[0].verify('silentfix')
		header = hdu[0].header
		hdu.close()
		use_stored_solution = False
		if 'lowres' in header['GRAT'].lower(): use_stored_solution = True
		ps.extract.make_wavecal([parameters['ARC_FILE_PREFIX'],fstr],'flat{}.fits'.format(cs),'wavecal{}'.format(cs),\
			qa_plot=qa_plot,qa_file=qa_file,use_stored_solution=use_stored_solution,verbose=verbose)

# extract all sources and standards
	scikeys = list(filter(lambda x: x[:4]=='SCI0',list(parameters.keys())))
	if len(scikeys)==0: 
		if verbose==True: print('No science files to reduce')
		return
	for k in scikeys:
		spar = parameters[k]

# SCIENCE TARGET
# reduce the first pair of targets
# NOTE: WOULD BE HELPFUL TO ADD CHECK THAT FILES HAVEN'T ALREADY BEEN REDUCED
		ps.extract.load_image([spar['SCIENCE_FILE_PREFIX'],extract_filestring(spar['TARGET_FILES'],'index')[:2]],\
			spar['FLAT_FILE'], spar['WAVECAL_FILE'],reduction_mode=spar['REDUCTION_MODE'], \
			flat_field=True, linearity_correction=True,qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# set extraction method
		ps.extract.set_extraction_type(spar['TARGET_TYPE'].split(' ')[-1])

# make spatial profiles
		ps.extract.make_spatial_profiles(qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# identify aperture positions
		ps.extract.locate_aperture_positions(spar['NPOSITIONS'], method=spar['APERTURE_METHOD'], qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# select orders to extract (set above)
		ps.extract.select_orders(include=spar['ORDERS'], qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# trace apertures
		ps.extract.trace_apertures(qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# define the aperture - psf, width, background
# NOTE - ONLY WORKS FOR POINT SOURCE, NEED TO FIX FOR XS?
		ps.extract.define_aperture_parameters(spar['PS_APERTURE'], psf_radius=spar['PSF_RADIUS'],bg_radius=spar['BACKGROUND_RADIUS'],\
						bg_width=spar['BACKGROUND_WIDTH'], qa_plot=qa_plot, qa_file=qa_file)

# extract away
		ps.extract.extract_apertures(verbose=verbose)

# conduct extraction of all remaining files
		fnum = extract_filestring(spar['TARGET_FILES'],'index')[2:]
		if len(fnum)>0:
			ps.extract.do_all_steps([spar['SCIENCE_FILE_PREFIX'],'{:.0f}-{:.0f}'.format(np.nanmin(fnum),np.nanmax(fnum))],verbose=verbose)

# FLUX STANDARD
# now reduce the first pair of standards
# NOTE: WOULD BE HELPFUL TO ADD CHECK THAT FILES HAVEN'T ALREADY BEEN REDUCED
		ps.extract.load_image([spar['SCIENCE_FILE_PREFIX'],extract_filestring(spar['STD_FILES'],'index')[:2]],\
			spar['FLAT_FILE'], spar['WAVECAL_FILE'],reduction_mode=spar['REDUCTION_MODE'], \
			flat_field=True, linearity_correction=True,qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# set extraction method
		ps.extract.set_extraction_type(spar['TARGET_TYPE'].split(' ')[-1])

# make spatial profiles
		ps.extract.make_spatial_profiles(qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# identify aperture positions
		ps.extract.locate_aperture_positions(spar['NPOSITIONS'], method=spar['APERTURE_METHOD'], qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# select orders to extract (set above)
		ps.extract.select_orders(include=spar['ORDERS'], qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# trace apertures
		ps.extract.trace_apertures(qa_plot=qa_plot, qa_file=qa_file, verbose=verbose)

# define the aperture - psf, width, background
# NOTE - ONLY WORKS FOR POINT SOURCE, NEED TO FIX FOR XS?
		ps.extract.define_aperture_parameters(spar['PS_APERTURE'], psf_radius=spar['PSF_RADIUS'],bg_radius=spar['BACKGROUND_RADIUS'],\
						bg_width=spar['BACKGROUND_WIDTH'], qa_plot=qa_plot, qa_file=qa_file)

# extract away
		ps.extract.extract_apertures(verbose=verbose)

# conduct extraction of all remaining files
		fnum = extract_filestring(spar['STD_FILES'],'index')[2:]
		if len(fnum)>0:
			ps.extract.do_all_steps([spar['SCIENCE_FILE_PREFIX'],'{:.0f}-{:.0f}'.format(np.nanmin(fnum),np.nanmax(fnum))],verbose=verbose)


# COMBINE SPECTRAL FILES
# NOTE: SPECTRA_FILE_PREFIX is currently hardcoded in code to be 'spectra' - how to change?
# ALSO: currently not plotting SXD spectra due to an error in combined file code
		if spar['MODE']=='SXD':
			ps.combine.combine_spectra(['spectra',spar['TARGET_FILES']],'{}{}'.format(spar['COMBINED_FILE_PREFIX'],spar['TARGET_FILES']),
				scale_spectra=True,scale_range=spar['SCALE_RANGE'],correct_spectral_shape=False,qa_plot=False,qa_file=False,verbose=verbose)
		else:
			ps.combine.combine_spectra(['spectra',spar['TARGET_FILES']],'{}{}'.format(spar['COMBINED_FILE_PREFIX'],spar['TARGET_FILES']),
				scale_spectra=True,scale_range=spar['SCALE_RANGE'],correct_spectral_shape=False,qa_plot=qa_plot,qa_file=qa_file,verbose=verbose)

# telluric star
		if spar['MODE']=='SXD':
			ps.combine.combine_spectra(['spectra',spar['STD_FILES']],'{}{}'.format(spar['COMBINED_FILE_PREFIX'],spar['STD_FILES']),
				scale_spectra=True,scale_range=spar['SCALE_RANGE'],correct_spectral_shape=False,qa_plot=False,qa_file=False,verbose=verbose)
		else:
			ps.combine.combine_spectra(['spectra',spar['STD_FILES']],'{}{}'.format(spar['COMBINED_FILE_PREFIX'],spar['STD_FILES']),
				scale_spectra=True,scale_range=spar['SCALE_RANGE'],correct_spectral_shape=False,qa_plot=qa_plot,qa_file=qa_file,verbose=verbose)

# FLUX CALIBRATE - NOT CURRENTLY IMPLEMENTED
# depends on fixed or moving source which method we use
		if (spar['TARGET_TYPE'].split(' ')[0]).strip()=='fixed':
			if verbose==True: print('fixed target; doing standard telluric correction and flux calibration')
		elif (spar['TARGET_TYPE'].split(' ')[0]).strip()=='moving':
			if verbose==True: print('moving target; doing reflectance telluric correction and flux calibration')
		else:
			if verbose==True: print('Target type {} not recognized, skipping telluric correction and flux calibration'.format(spar['TARGET_TYPE'].split(' ')[0]))

	return


# QA page maker
def makeQApage(driver_input,log_input,log_html_name='log.html',nimages=4,image_width=400,verbose=ERROR_CHECKING):
	'''
	Generates an html file for the qa folder to review reductions
	'''
# if necessary read in driver
	if isinstance(driver_input,str)==True:
		if os.path.exists(driver_input)==False: raise ValueError('Cannot find driver file {}'.format(driver_input))
		try: driver = read_driver(driver_input)
		except: raise ValueError('Unable to read in driver file {}'.format(driver_input))
	elif isinstance(driver_input,dict)==True:
		driver = copy.deepcopy(driver_input)
	else: raise ValueError('Do not recognize format of driver input {}'.format(driver_input))

# if necessary read in log
	if isinstance(log_input,str)==True:
		if os.path.exists(log_input)==False: raise ValueError('Cannot find log file {}'.format(log_input))
		if log_input.split('.')[-1] != 'csv': raise ValueError('Log file must end in csv, you passed {}'.format(log_input))
		obslog = pd.read_csv(log_input)
	elif isinstance(log_input,pd.DataFrame)==True: 
		obslog = copy.deepcopy(log_input)
	else: raise ValueError('Do not recognize format of observing log input {}'.format(log_input))

# read in templates and split index
	with open(QA_INDEX_TEMPLATE_FILE,'r') as f: index_txt = f.read()
	index_top,index_bottom = index_txt.split('[SOURCES]')
	with open(QA_SOURCE_TEMPLATE_FILE,'r') as f: source_txt = f.read()

# extract key parameters
	page_parameters = {}
	for x in ['INSTRUMENT','CAL_SETS','QA_FOLDER','PROC_FOLDER']:
		page_parameters[x] = driver[x]

# replace relevant parameters in top of index
	index_top = index_top.replace('[CSS]',QA_CSS)
	index_top = index_top.replace('[LOG_HTML]',log_html_name)
# NEED TO UPDATE WEATHER LINK
	index_top = index_top.replace('[WEATHER_HTML]',MKWC_ARCHIVE_URL)
	for x in ['UT_DATE','PROGRAM','OBSERVER']:
		index_top = index_top.replace('['+x+']',obslog[x].iloc[0])
	imodes = list(set(list(obslog['MODE'])))
	imodes.sort()
	s = imodes[0]
	if len(imodes)>1: [s.append(', {}'.format(x)) for x in imodes[1:]]
	index_top = index_top.replace('[INSTRUMENT_MODES]',s)

	output = index_top

# now cycle through each of the science sets
	scikeys = list(filter(lambda x: x[:4]=='SCI0',list(driver.keys())))
	scikeys.sort()
	for sci in scikeys:
# fill in source info		
		sci_param = driver[sci]
		sci_txt = copy.deepcopy(source_txt)
		for x in ['TARGET_NAME','TARGET_TYPE','TARGET_FILES','STD_NAME','STD_FILES','FLAT_FILE','WAVECAL_FILE','SCIENCE_FILE_PREFIX','SPECTRAL_FILE_PREFIX','COMBINED_FILE_PREFIX','CALIBRATED_FILE_PREFIX','STITCHED_FILE_PREFIX']:
			sci_txt = sci_txt.replace('['+x+']',sci_param[x])
		sci_log = obslog[obslog['SOURCE_NAME']==sci_param['TARGET_NAME']]
		if len(sci_log)>0:
			for x in ['RA','DEC','AIRMASS','MODE','SLIT']:
				sci_txt = sci_txt.replace('['+x+']',str(sci_log[x].iloc[0]))
			designation = 'J'+sci_log['RA'].iloc[0].replace(':','').replace(' ','').replace('.','')
			if int(sci_log['DEC'].iloc[0][:2]+'1')>0: designation+='%2B'
			designation+=sci_log['DEC'].iloc[0].replace(':','').replace(' ','').replace('.','')
			sci_txt = sci_txt.replace('[DESIGNATION]',designation)
			sci_txt = sci_txt.replace('[UT_START]',sci_log['UT_TIME'].iloc[0])
			sci_txt = sci_txt.replace('[UT_END]',sci_log['UT_TIME'].iloc[-1])
#			sci_log['TINT'] = [sci_log.loc['INTEGRATION',i]*sci_log.loc['COADDS',i] for i in range(len(sci_log))]
			sci_txt = sci_txt.replace('[INTEGRATION]','{:.1f}'.format(np.nansum(sci_log['INTEGRATION'])))
# final calibrated files
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],sci_param['CALIBRATED_FILE_PREFIX'])+'{}*{}'.format(sci_param['TARGET_FILES'],QA_PLOT_TYPE))
		if len(imfile)>0: sci_txt = sci_txt.replace('[CALIBRATED_FILE]',os.path.basename(imfile[0]))
# combined files
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],sci_param['COMBINED_FILE_PREFIX'])+'{}*_raw{}'.format(sci_param['TARGET_FILES'],QA_PLOT_TYPE))
#		print(os.path.join(page_parameters['QA_FOLDER'],sci_param['COMBINED_FILE_PREFIX'])+'{}*_raw{}'.format(sci_param['TARGET_FILES'],QA_PLOT_TYPE),imfile)
		if len(imfile)>0: sci_txt = sci_txt.replace('[TARGET_COMBINED_RAW_FILE]',os.path.basename(imfile[0]))
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],sci_param['COMBINED_FILE_PREFIX'])+'{}*_scaled{}'.format(sci_param['TARGET_FILES'],QA_PLOT_TYPE))
		if len(imfile)>0: sci_txt = sci_txt.replace('[TARGET_COMBINED_SCALED_FILE]',os.path.basename(imfile[0]))
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],sci_param['COMBINED_FILE_PREFIX'])+'{}*_raw{}'.format(sci_param['STD_FILES'],QA_PLOT_TYPE))
		if len(imfile)>0: sci_txt = sci_txt.replace('[STD_COMBINED_RAW_FILE]',os.path.basename(imfile[0]))
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],sci_param['COMBINED_FILE_PREFIX'])+'{}*_scaled{}'.format(sci_param['STD_FILES'],QA_PLOT_TYPE))
		if len(imfile)>0: sci_txt = sci_txt.replace('[STD_COMBINED_SCALED_FILE]',os.path.basename(imfile[0]))
# individual files
		fnums = extract_filestring(sci_param['TARGET_FILES'],'index')
		loopstr = '<table>\n <tr>\n'
		cnt=0
		for i in fnums:
			if cnt>0 and np.mod(cnt,nimages)==0: loopstr+=' </tr>\n <tr> \n'
			imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],'spectra*{:.0f}{}'.format(i,QA_PLOT_TYPE)))
			if len(imfile)>0: 
				imfile.sort()
				loopstr+='  <td align="center">\n   <embed src="{}" width=300 height=300>\n   <br>{}\n   </td>\n'.format(os.path.basename(imfile[0]),os.path.basename(imfile[0]))
				cnt+=1
		loopstr+=' </tr>\n</table>'
		sci_txt = sci_txt.replace('[TARGET_SPECTRA_FILES]',loopstr)

		fnums = extract_filestring(sci_param['STD_FILES'],'index')
		loopstr = '<table>\n <tr>\n'
		cnt=0
		for i in fnums:
			if cnt>0 and np.mod(cnt,nimages)==0: loopstr+=' </tr>\n <tr> \n'
			imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],'spectra*{:.0f}{}'.format(i,QA_PLOT_TYPE)))
			if len(imfile)>0: 
				imfile.sort()
				loopstr+='  <td align="center">\n   <embed src="{}" width=300 height=300>\n   <br>{}\n   </td>\n'.format(os.path.basename(imfile[0]),os.path.basename(imfile[0]))
				cnt+=1
		loopstr+=' </tr>\n</table>'
		sci_txt = sci_txt.replace('[STD_SPECTRA_FILES]',loopstr)

# traces and apertures - TO BE DONE
		sci_txt = sci_txt.replace('[TRACE_APERTURE_FILES]','')

		output+=sci_txt

# insert calibrations
	cal_sets = driver['CAL_SETS'].split(',')

# flats
	loopstr = '<table>\n <tr>\n'
	cnt=0
	for cs in cal_sets:
		if cnt>0 and np.mod(cnt,nimages)==0: loopstr+=' </tr>\n <tr> \n'
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],'flat{}_locateorders{}'.format(cs,QA_PLOT_TYPE)))
		if len(imfile)>0: 
			imfile.sort()
			loopstr+='  <td align="center">\n   <embed src="{}" width=300 height=300>\n   <br>{}\n   </td>\n'.format(os.path.basename(imfile[0]),os.path.basename(imfile[0]))
			cnt+=1
	loopstr+=' </tr>\n</table>\n\n'
	index_bottom = index_bottom.replace('[FLATS]',loopstr)

# flats
	# loopstr = ''
	# for cs in cal_sets:
	# 	imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],'flat{}_locateorders{}'.format(cs,QA_PLOT_TYPE)))
	# 	if len(imfile)>0: 
	# 		imfile.sort()
	# 		loopstr+='<table><tr><td align="center">\n <embed src="{}" width=300 height=300>\n <br>{}\n</td></tr></table>'.format(os.path.basename(imfile[0]),os.path.basename(imfile[0]))
	# index_bottom = index_bottom.replace('[FLATS]',loopstr)

# wavecals
	loopstr = '<table>\n <tr>\n'
	cnt=0
	for cs in cal_sets:
		if cnt>0 and np.mod(cnt,nimages)==0: loopstr+=' </tr>\n <tr> \n'
		imfile = glob.glob(os.path.join(page_parameters['QA_FOLDER'],'wavecal{}_shift.pdf'.format(cs)))
		if len(imfile)>0: 
			imfile.sort()
			loopstr+='  <td align="center">\n   <embed src="{}" width=300 height=400>\n   <br>{}\n   </td>\n'.format(os.path.basename(imfile[0]),os.path.basename(imfile[0]))
			cnt+=1
	loopstr+=' </tr>\n</table>\n'
	index_bottom = index_bottom.replace('[WAVECALS]',loopstr)

# WRITE OUT HTML FILE
	output+=index_bottom
	with open(os.path.join(page_parameters['QA_FOLDER'],QA_HTML_OUTPUT_NAME),'w') as f:
		f.write(output)
	return


# external function call
if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('Call this function from the command line as python batch.py [data_folder] [cal_folder] [proc_folder]')
	else:
		dp = process_folder(sys.argv[1])
		dp = process_folder(sys.argv[1])
		if len(sys.argv) > 4: log_file = sys.argv[4]
		else: log_file = sys.argv[3]+'log.html'
		print('\nWriting log to {}'.format(log_file))
		write_log(dp,log_file)
		if len(sys.argv) > 5: driver_file = sys.argv[5]
		else: driver_file = sys.argv[3]+'driver.txt'
		print('\nWriting batch driver file to {}'.format(driver_file))
		write_driver(dp,driver_file,data_folder=sys.argv[1],cal_folder=sys.argv[2],proc_folder=sys.argv[3])
		txt = input('\nCheck the driver file {} and press return when ready to proceed...\n')
		par = read_driver(driver_file)
		print('\nReducing spectra')
		batch_reduce(par)
		print('\nReduction complete: processed files are in {}'.format(sys.argv[3]))
