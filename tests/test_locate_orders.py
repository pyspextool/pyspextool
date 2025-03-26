import importlib
import glob
import os
import pyspextool as ps
from pyspextool import config as setup

def test_locate_orders():

    
    # Get set up

    ps.pyspextool_setup(
        raw_path="tests/test_data/raw/uspex-SXD/data/",
        qa_path="tests/test_data/raw/uspex-SXD/qa/",
        cal_path="tests/test_data/raw/uspex-SXD/cals/",
        proc_path="tests/test_data/raw/uspex-SXD/proc/",
        verbose=True,
        qa_show=False,
        qa_write=True,
        qa_extension=".png",
    )

    # Get flat info

    


    # Load a single flat frame

    file = setup.state['raw_path']+'/sbd.2023B006.231021.flat.00293.a.fits'

    module = 'pyspextool.instruments.' + setup.state['instrument'] + \
             '.' + setup.state['instrument']

    instr = importlib.import_module(module)


    result = instr.read_fits([file],setup.state['linearity_info'])

    hdr = result[2]
    
    mode = hdr[0]['MODE'][0]
    modefile = setup.state['instrument_path']+'/'+mode +'_flatinfo.fits'

    modeinfo = ps.extract.flat.read_flatcal_file(modefile)



    result = ps.extract.flat.locate_orders(result[0],
                                           modeinfo['guesspos'],
                                           modeinfo['xranges'],
                                           modeinfo['step'],
                                           modeinfo['slith_range'],
                                           modeinfo['edgedeg'],
                                           modeinfo['ybuffer'],
                                           modeinfo['flatfrac'],
                                           modeinfo['comwidth'],
                                           qa_show=False)

    assert result[0][4][0][0] > 0

