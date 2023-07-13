# pyspextool batch driver function
# usage methods
#
# >> python batch_driver.py [base folder] 
# >> python batch_driver.py [data folder] [cal folder] [proc folder] [qa folder]
# >> python batch_driver.py --reduce-only [base folder]
# >> python batch_driver.py --qa-only [base folder]
#

import argparse
import os
import pandas
import sys
from pyspextool.batch.batch import process_folder,write_log,write_driver,read_driver,batch_reduce,makeQApage

DATA_FOLDER = 'data/'
CALS_FOLDER = 'cals/'
PROC_FOLDER = 'proc/'
QA_FOLDER = 'qa/'

class runBatch():
	def __init__(self):
		parser = argparse.ArgumentParser(description='Runs the batch mode for pyspextool')
		parser.add_argument('inputs', nargs='*',
			help='main inputs, with the following options:\n(1) provide four paths to the data, cals, proc, and qa folders in that order (if qa folder not provided, it will be assigned to cals);\n(2) provide one path to a "base" folder inside which a "data" folder already exists (proc, cals, and qa folders will be created if they do not exist);\n(3) if --reduce-only is set, provide the path to the driver file name or include it with the -d keyword')
		parser.add_argument('-l', metavar='log-file-prefix', nargs=1, default='',
			required=False, help='set log file prefix (default is log)')
		parser.add_argument('-d', metavar='driver-filename', nargs=1, default='',
			required=False, help='set driver file prefix (default is "driver")')
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
		parser.add_argument('--qa-only', action='store_true',default=False,
			required=False, help='set to just generate QA page')
		# parser.add_argument('--base', action='store_true',default=False,
		# 	required=False, help='set to establish a base folder with data, proc and cals folder within (cals and proc are created if they do not exist)')
		parser.add_argument('--qaplot-off', action='store_false',default=True,
			required=False, help='set to turn OFF QA plots')
		parser.add_argument('--overwrite', action='store_true',default=False,
			required=False, help='set to automatically overwrite files')
		parser.add_argument('--verbose', action='store_true',default=False,
			required=False, help='set to return verbose feedback')

		args = vars(parser.parse_args())
#		print(args)
		folders = args['inputs']
		driver_file = args['d']
		log_file_prefix = args['l']
#		print(args)
		# raise()

# if nothing passed, assume we are using the local folder
		if len(folders)<1: folders=[os.path.abspath('./')]

# build folders 
		if len(folders) ==1:
			if os.path.exists(folders[0])==False: raise ValueError('Cannot find base folder {}'.format(folders[0]))
			elif os.path.exists(os.path.join(folders[0],DATA_FOLDER))==False: raise ValueError('Cannot find data folder under {}'.format(folders[0]))
			else:
# cals				
				nfold = os.path.join(folders[0],CALS_FOLDER)
				if os.path.exists(nfold)==False: 
					os.mkdir(nfold)
					if args['verbose']==True: print('\nCreated cals folder {}'.format(nfold))
				if len(folders)>=2: folders[1] = nfold
				else: folders.append(nfold)
# proc				
				nfold = os.path.join(folders[0],PROC_FOLDER)
				if os.path.exists(nfold)==False: 
					os.mkdir(nfold)
					if args['verbose']==True: print('\nCreated proc folder {}'.format(nfold))
				if len(folders)>=3: folders[2] = nfold
				else: folders.append(nfold)
# qa				
				nfold = os.path.join(folders[0],QA_FOLDER)
				if os.path.exists(nfold)==False: 
					os.mkdir(nfold)
					if args['verbose']==True: print('\nCreated qa folder {}'.format(nfold))
				if len(folders)>=4: folders[3] = nfold
				else: folders.append(nfold)
# data
				folders[0] = os.path.join(folders[0],DATA_FOLDER)

# check folders 
		if len(folders)<4: raise ValueError('Need to specify three folder paths for data, cals, and proc')
		folders = [os.path.abspath(f) for f in folders]
		for f in folders:
			if f!='' and os.path.exists(f)==False: raise ValueError('Cannot find folder {}'.format(f)) 

# generate log csv and html files and put in qa folder
		if log_file_prefix=='': log_file_prefix = 'log'
		log_file_prefix = os.path.join(folders[3],log_file_prefix)
		if os.path.exists(log_file_prefix+'.html')==True and os.path.exists(log_file_prefix+'.csv')==True and args['overwrite']==False and args['rebuild_log']==False:
			print('\nWARNING: html log file {} and csv log file {} already exists; use --overwrite if you want to overwrite or --rebuild-log to rebuild'.format(log_file_prefix+'.html',log_file_prefix+'.csv'))
		else:
			if args['driver_only']==False and args['reduce_only']==False and args['qa_only']==False:
				dp = process_folder(folders[0])
				for x in ['.csv','.html']:
					if os.path.exists(log_file_prefix+x) and args['overwrite']==False and args['rebuild_log']==False:
						print('\nWARNING: {} log file {} already exists; use --overwrite to overwrite'.format(x,log_file_prefix+x))
					else:
						if args['verbose']==True: print('\nWriting log to {}'.format(log_file_prefix+x))
						write_log(dp,log_file_prefix+x)
		if args['log_only']==True: 
#			print('\n\nLog files {} and {} created\n\n'.format(log_file_prefix+'.csv',log_file_prefix+'.html'))
			return

# query for next step?
#		txt = input('\n\nCheck the log file {} and press return when you are ready to proceed, or type CNTL-C to abort...\n\n'.format(log_file_prefix+'.csv'))

# generate driver file and put in proc folder
		if driver_file=='': driver_file = 'driver.txt'
		driver_file = os.path.join(folders[2],driver_file)
		if os.path.exists(driver_file)==True and args['overwrite']==False and args['rebuild_driver']==False:
			print('\nWARNING: driver file {} already exists; use --overwrite to overwrite or --rebuild-driver to rebuild'.format(driver_file))
		else:
			if args['reduce_only']==False and args['qa_only']==False:
				if os.path.exists(log_file_prefix+'.csv')==True:
					dp = pandas.read_csv(log_file_prefix+'.csv')
				else:
					dp = process_folder(folders[0])
					print('\nWARNING: could not find log file {}, this may be a problem later'.format(log_file_prefix+'.csv'))
				if args['verbose']==True: print('\nGenerating driver file and writing to {}'.format(driver_file))
#				print(driver_file,folders[0])
				write_driver(dp,driver_file,data_folder=folders[0],verbose=args['verbose'],check=True,create_folders=True,exclude_lxd=True)

		if args['driver_only']==True: 
			print('\n\nDriver file {} created'.format(driver_file))
			print('Review this file before reducing\n\n')
			return

# query for next step
		txt = input('\n\nCheck the driver file {} and press return when you are ready to proceed, or type CNTL-C to abort...\n\n'.format(driver_file))

# reduction - only need the driver file for this
##
## NOTE - SOMETHING WEIRD HERE WHERE CALIBRATION FILE CREATION BREAKS IF VERBOSE IS NOT SET
##
		if args['qa_only']==False:
			# if driver_file=='': driver_file = folders[0]
			# if os.path.isdir(driver_file)==True: raise ValueError('Parameter you passed - {} - is a directory; provide path to driver file'.format(driver_file)) 
			if os.path.exists(driver_file)==False: raise ValueError('Cannot find driver file {}'.format(driver_file)) 
			par = read_driver(driver_file)
			if args['verbose']==True: print('\n\nReducing spectra\n\n')
			batch_reduce(par,verbose=args['verbose'])

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

		makeQApage(driver_file,log_file_prefix+'.csv',verbose=args['verbose'])

		print('\n\nQA page created, please review\n\n')

		return

if __name__ == '__main__':
	app = runBatch()
