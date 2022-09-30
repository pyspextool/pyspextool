import numpy as np
from scipy import interpolate

from pyspextool.cl import config
from pyspextool.utils.math import mean_data_stack


def make_spatial_profiles():

    #
    # Check the continue variables
    #
    
    if config.state['pscontinue'] < 1 or config.state['xscontinue'] < 1:
        print('Previous reductions steps not complete.')

    profiles = []
        
    for i in range(config.state['norders']):

        # Unpack the data
        
        order = config.state['rectorders'][i]

        img = order['img']
        y = order['y']
        w = order['w']

        ncols = len(w)
        nrows = len(y)

        # Subtract the background

        medbg = np.median(img,axis=0)

        bgimg = np.tile(medbg, (nrows, 1))

        np.subtract(img, bgimg, out=img)

        if config.state['wavecalfile'] != None:

            # Do the interpolate of the atmosphere
            
            f = interpolate.interp1d(config.state['atrans_wave'],
                                     config.state['atrans_trans'],
                                     fill_value=1)
            rtrans = f(w)

            # Clip low points and create weight array
            
            rtrans = np.where(rtrans <= 0.1, rtrans, 0.1)

            weights = np.tile((1/rtrans)**2, (nrows, 1))            

            # Combine them together.  Must rotate first to work with
            # mean_data_stack

            mean, mvar, mask = mean_data_stack(np.rot90(img, 3),
                                               weights=np.rot90(weights, 3),
                                               robust=5)

        else:

            mean, mvar, mask = mean_data_stack(np.rot90(img, 3),
                                               robust=5)            

        # Normalize by the total absolute flux

        mean  = mean/np.sum(np.abs(mean))

        # Package up

        profiles.append({'order':config.state['orders'][i],'y':y,
                         'p':np.flip(mean)})

    config.state['profiles'] = profiles
