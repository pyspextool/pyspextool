# pyspextool batch driver function
# usage methods
#
# >> python batch_driver.py [base folder] 
# >> python batch_driver.py [data folder] [cal folder] [proc folder] [qa folder]
# >> python batch_driver.py --reduce-only [base folder]
# >> python batch_driver.py --qa-only [base folder]
#

import argparse
import copy
import os
import pandas
import sys
from pyspextool.batch import batch

LOG_FILE_PREFIX_DEFAULT = 'log'
DRIVER_FILE_DEFAULT = 'driver.txt'
VERSION = batch.VERSION
#VERSION = '2024.02.25'
AUTHORS = [
	'Adam Burgasser',
	'Jean Marroquin',
	'Evan Watson']

class runBatch():
	def __init__(self):
		parser = argparse.ArgumentParser(description='Runs the batch mode for pyspextool')
		parser.add_argument('inputs', nargs='*',
			help='main inputs, with the following options:'+\
			'\n\t(1) no inputs assumes user is in currently in a base folder that contains a "data" folder with raw files, in which "cals", "proc", and "qa" folders exist or will be created;'+\
			'\n\t(2) single input giving the full path to a base folder that contains the "data" folder and optionally "cals", "proc", and "qa" folders;'+\
			'\n\t(3) two inputs giving the full path to the data folder and then full path to a base folder in which "cals", "proc", and "qa" folders exist or will be created; or'+\
			'\n\t(4) three or four inputs giving the full paths to the "data", "cals", "proc", and optionally "qa" folders, in that order')
		parser.add_argument('-l', metavar='log-file-prefix', nargs=1, default='',
			required=False, help='set log file prefix (default is "log")')
		parser.add_argument('-d', metavar='driver-filename', nargs=1, default='',
			required=False, help='set driver file name (default is "driver.txt")')
		parser.add_argument('--test', action='store_true',default=False,
			required=False, help='set to test the pyspextool installation')
		parser.add_argument('--organize', action='store_true',default=False,
			required=False, help='organize the legacy archive download data files')
		parser.add_argument('--organize-legacy', action='store_true',default=False,
			required=False, help='organize the legacy archive download data files')
		parser.add_argument('--organize-irsa', action='store_true',default=False,
			required=False, help='organize the IRSA archive download data files')
		parser.add_argument('--log-only', action='store_true',default=False,
			required=False, help='set to just save a log file')
		parser.add_argument('--driver-only', action='store_true',default=False,
			required=False, help='set to just generate driver file')
		parser.add_argument('--rebuild-log', action='store_true',default=False,
			required=False, help='set to rebuild the log file')
		parser.add_argument('--rebuild-driver', action='store_true',default=False,
			required=False, help='set to rebuild the driver file from the log')
		parser.add_argument('--reduce-only', action='store_true',default=False,
			required=False, help='set to just reduce (no log or driver generation)')
		parser.add_argument('--rereduce', action='store_true',default=False,
			required=False, help='set to just re-reduce (no log, driver, or calibrations generation, overwrite prior extractions)')
		parser.add_argument('--qa-only', action='store_true',default=False,
			required=False, help='set to just generate QA page')
		# parser.add_argument('--base', action='store_true',default=False,
		# 	required=False, help='set to establish a base folder with data, proc and cals folder within (cals and proc are created if they do not exist)')
		parser.add_argument('--qaplot-off', action='store_false',default=True,
			required=False, help='set to turn OFF QA plots')
		parser.add_argument('--overwrite', action='store_true',default=False,
			required=False, help='set to automatically overwrite files, including log and driver')
		parser.add_argument('--no-pause', action='store_true',default=False,
			required=False, help='set to remove all pauses for user input')
		# parser.add_argument('--verbose', action='store_true',default=False,
		# 	required=False, help='set to return verbose feedback')
		parser.add_argument('--no-cals', action='store_true',default=False,
			required=False, help='set to skip calibration file creation')
		parser.add_argument('--no-extract', action='store_true',default=False,
			required=False, help='set to skip spectral extraction')
		parser.add_argument('--no-combine', action='store_true',default=False,
			required=False, help='set to skip spectral file combination')
		parser.add_argument('--no-fluxtell', action='store_true',default=False,
			required=False, help='set to skip flux and telluric correction')
		parser.add_argument('--no-stitch', action='store_true',default=False,
			required=False, help='set to skip order stitching')
		parser.add_argument('--quiet', action='store_true',default=False,
			required=False, help='set to return minimal feedback')
		parser.add_argument('--version', action='store_true',default=False,
			required=False, help='report version number of batch reduction code')

		args = vars(parser.parse_args())
#		print(args)
		folders = args['inputs']
		driver_file = args['d']
		log_file_prefix = args['l']
		base_folder = ''
#		print(args)
		# raise()

# give version number (only)
		if args['quiet']==False: print('\npyspextool batch reduction code version {}\n'.format(VERSION))
		if args['version']==True: return

# test
		if args['test']==True: 
			batch.test()
			return

# if nothing passed, assume we are using the local folder as the base folder
		if len(folders)<1: 
			base_folder='./'
			folders = [batch.REDUCTION_FOLDERS[0]]

# one input - assume to be base folder 
		if len(folders) == 1:
			if base_folder=='': base_folder=folders[0]
			if os.path.exists(base_folder)==False: raise ValueError('Cannot find base folder {}'.format(base_folder))
#			elif os.path.exists(os.path.join(folders[0],batch.REDUCTION_FOLDERS[0]))==False: raise ValueError('Cannot find data folder under {}'.format(folders[0]))
#			bfold = copy.deepcopy(folders[0])
			folders[0] = os.path.join(base_folder,batch.REDUCTION_FOLDERS[0])

# organize legacy data?
		if args['organize']==True or args['organize_legacy']==True:
			batch.organizeLegacy(base_folder,verbose=(not args['quiet']),overwrite=args['overwrite'],makecopy=True)
			if args['quiet']==False: print('\n\nFinished IRTF Legacy archive file organization\n\n')
			return

# organize irsa data?
		if args['organize_irsa']==True:
			batch.organizeIRSA(base_folder,verbose=(not args['quiet']),overwrite=args['overwrite'],makecopy=True)
			if args['quiet']==False: print('\n\nFinished IRTF IRSA archive file organization\n\n')
			return

# two inputs - assume to be data folder and base folder 
		if len(folders) == 2:
			if os.path.exists(folders[0])==False: raise ValueError('Cannot find data folder {}'.format(folders[0]))
			if os.path.exists(folders[1])==False: raise ValueError('Cannot find base folder {}'.format(folders[1]))
			if base_folder=='': base_folder = copy.deepcopy(folders[1])
			folders = [folders[0]]

# set or create data, cals, proc, and qa folders 
		for i,nm in enumerate(batch.REDUCTION_FOLDERS):
			nfold = os.path.join(base_folder,nm+'/')
			if len(folders) <= i: folders.append(nfold)
#			folders[i] = os.path.abspath(folders[i])
			if os.path.exists(folders[i])==False:
				os.mkdir(folders[i])
				if args['quiet']==False: print('\nCreated {} folder {}'.format(nm,folders[i]))

# check folders 
#		folders = [os.path.abspath(f) for f in folders]
		for i,f in enumerate(folders):
			if f=='': raise ValueError('Empty path name for {} folder'.format(batch.REDUCTION_FOLDERS[i]))
			if os.path.exists(f)==False: raise ValueError('Cannot find {} folder {}'.format(batch.REDUCTION_FOLDERS[i],f)) 

# generate log csv and html files and put in qa folder
		if log_file_prefix=='': log_file_prefix = copy.deepcopy(LOG_FILE_PREFIX_DEFAULT)
		log_file_prefix = os.path.join(folders[3],log_file_prefix)
		if os.path.exists(log_file_prefix+'.html')==True and os.path.exists(log_file_prefix+'.csv')==True and args['overwrite']==False and args['rebuild_log']==False:
			print('\nWARNING: html log file {} and csv log file {} already exists; use --overwrite if you want to overwrite or --rebuild-log to rebuild'.format(log_file_prefix+'.html',log_file_prefix+'.csv'))
		else:
			if args['driver_only']==False and args['reduce_only']==False and args['qa_only']==False:
				dp = batch.processFolder(folders[0])
				for x in ['.csv','.html']:
					if os.path.exists(log_file_prefix+x) and args['overwrite']==False and args['rebuild_log']==False:
						print('\nWARNING: {} log file {} already exists so not saving; use --overwrite to overwrite'.format(x,log_file_prefix+x))
					else:
#						if args['quiet']==False: print('\nWriting log to {}'.format(log_file_prefix+x))
						batch.writeLog(dp,log_file_prefix+x)

# query to pause and check log
			if args['no_pause']==False and args['log_only']==False: txt = input('\n\nCheck the LOG FILES {} and {} and press return when you are ready to proceed, or type CNTL-C to abort...\n\n'.format(log_file_prefix+'.csv',log_file_prefix+'.html'))

		if args['log_only']==True: 
			print('\n\nLog files {} and {} created.'.format(log_file_prefix+'.csv',log_file_prefix+'.html'))
			return


# generate driver file and put in proc folder
		if driver_file=='': driver_file = copy.deepcopy(DRIVER_FILE_DEFAULT)
		driver_file = os.path.join(folders[2],driver_file)
		if os.path.exists(driver_file)==True and args['overwrite']==False and args['rebuild_driver']==False:
			print('\nWARNING: driver file {} already exists so not saving; use --overwrite to overwrite or --rebuild-driver to rebuild'.format(driver_file))
		else:
			if args['reduce_only']==False and args['qa_only']==False:
				if os.path.exists(log_file_prefix+'.csv')==True:
					dp = pandas.read_csv(log_file_prefix+'.csv')
				else:
					dp = batch.processFolder(folders[0])
					print('\nWARNING: could not find log file {}, this may be a problem later'.format(log_file_prefix+'.csv'))
				if args['quiet']==False: print('\nGenerating driver file and writing to {}'.format(driver_file))
#				print(driver_file,folders[0])
				batch.writeDriver(dp,driver_file,data_folder=folders[0],verbose=(not args['quiet']),check=True,create_folders=True)

# query to pause and check driver
			if args['no_pause']==False and args['driver_only']==False: txt = input('\n\nCheck the DRIVER FILE {} and press return when you are ready to proceed, or type CNTL-C to abort...\n\n'.format(driver_file))

		if args['driver_only']==True: 
			print('\n\nDriver file {} created.'.format(driver_file))
			return

# reduction - only need the driver file for this
##
## NOTE - SOMETHING WEIRD HERE WHERE CALIBRATION FILE CREATION BREAKS IF VERBOSE IS NOT SET
##
		if args['qa_only']==False:
			# if driver_file=='': driver_file = folders[0]
			# if os.path.isdir(driver_file)==True: raise ValueError('Parameter you passed - {} - is a directory; provide path to driver file'.format(driver_file)) 
			if os.path.exists(driver_file)==False: raise ValueError('Cannot find driver file {}'.format(driver_file)) 
			par = batch.readDriver(driver_file)

# stop if there are no science or calibration files
#			print(par)
#			if 'OBS_SET' not in list(par.keys()):
			if 'OBS_SET' not in [x[:7] for x in list(par.keys())]:
				print('No science files listed in the driver file {}; recheck this file before proceeding'.format(driver_file))
				return
			if 'CAL_SETS' not in list(par.keys()):
				print('No calibration sets listed in the driver file {}; recheck this file before proceeding'.format(driver_file))
				return

# add in additional keywords for specific reduction steps:
			if args['no_cals']==True: par['CALIBRATIONS']=False
			if args['no_extract']==True: par['EXTRACT']=False
			if args['no_combine']==True: par['COMBINE']=False
			if args['no_fluxtell']==True: par['FLUXTELL']=False
			if args['no_stitch']==True: par['STITCH']=False
			if args['rereduce']==True: par['REREDUCE']=True
			if args['overwrite']==True: par['OVERWRITE']=False

			if args['quiet']==False: print('\n\nReducing spectra\n\n')
			batch.batchReduce(par,verbose=(not args['quiet']))

		if args['reduce_only']==True: 
			print('\n\nReduction completed!\n\n')
			return

# make qa plots
#		if args['qa_only']==True:
		# if len(folders)<2: raise ValueError('For option --qa-only, you must provide driver file and log file (csv)')
		# if driver_file=='': driver_file = folders[0]
		# if os.path.isdir(driver_file)==True:
		# 	if len(folders)<4: raise ValueError('Could not determine driver file from inputs')
		# 	driver_file = folders[2]+'driver.txt'

		# if log_file_prefix=='': log_file_prefix = 'log'
		# # if os.path.isdir(log_file_prefix)==True:
		# # 	if len(folders)<4: raise ValueError('Could not determine log file from inputs')
		# # 	log_file_prefix = folders[3]+'log'
		# if os.path.exists(log_file_prefix+'.csv')==False: log_file_prefix = folders[3]+'log'

		if os.path.exists(driver_file)==False: raise ValueError('Cannot find driver file {}'.format(driver_file)) 
		if os.path.exists(log_file_prefix+'.csv')==False: raise ValueError('Cannot find log file {}'.format(log_file_prefix+'.csv')) 

		batch.makeQApage(driver_file,log_file_prefix+'.csv',verbose=(not args['quiet']))

		print('\n\nQA page created, please review\n\n')

		return

if __name__ == '__main__':
	app = runBatch()
