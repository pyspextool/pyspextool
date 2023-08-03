from pyspextool.batch.batch import *
import os
import glob
import pandas

basefold = os.path.dirname(ps.__file__)+'/test_data/'
test_instruments = ['spex-prism','spex-SXD']
tfold = os.path.join(basefold,'spex-prism')

def test_inputs():
# make sure code can find test data
	assert os.path.exists(basefold)
	for inst in test_instruments:
		assert os.path.exists(os.path.join(basefold,inst))
		for fold in REDUCTION_FOLDERS: assert os.path.exists(os.path.join(basefold,inst,fold))
	dfiles = glob.glob(os.path.join(tfold,'data/*.fits'))
	assert len(dfiles)>0
	
def test_processFolder():
	result = processFolder(os.path.join(tfold,'data/'))
	assert isinstance(result,pandas.core.frame.DataFrame)
	assert len(result)>0
	for k in list(HEADER_DATA.keys()): assert k in list(result.columns)

def test_organizeLegacy():
	pass

def test_writeLog():
# write a file and then delete it
	pass

def test_readDriver():
# read in existing driver in test folder
	pass

def test_writeDriver():
# write a file and then delete it
	pass

def test_makeQApage():
# write a file and then delete it
	pass

def test_batchReduce():
	pass



