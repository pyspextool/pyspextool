import numpy as np
from scipy import interpolate

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.utils.math import mean_data_stack


def make_spatial_profiles(iplot=True, qafile=False):
    #
    # Check the continue variables
    #

    check_continue(2)

    #
    # Build the profiles
    #

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

        medbg = np.median(img, axis=0)

        bgimg = np.tile(medbg, (nrows, 1))

        np.subtract(img, bgimg, out=img)

        if config.state['wavecalfile'] is not None:

            # Do the interpolate of the atmosphere

            f = interpolate.interp1d(config.state['atrans_wave'],
                                     config.state['atrans_trans'],
                                     fill_value=1)
            rtrans = f(w)

            # Clip low points and create weight array

            rtrans = np.where(rtrans <= 0.1, rtrans, 0.1)

            weights = np.tile((1 / rtrans) ** 2, (nrows, 1))

            # Combine them together.  Must rotate first to work with
            # mean_data_stack

            mean, mvar, mask = mean_data_stack(np.rot90(img, 3),
                                               weights=np.rot90(weights, 3),
                                               robust=5)

        else:

            mean, mvar, mask = mean_data_stack(np.rot90(img, 3),
                                               robust=5)

            # Normalize by the total absolute flux

        mean = mean / np.sum(np.abs(mean))

        # Package up

        profiles.append({'order': config.state['orders'][i], 'y': y,
                         'p': np.flip(mean)})

    config.state['profiles'] = profiles

    if iplot is True:
        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      np.ones(config.state['norders'], dtype=int))

    if qafile is True:
        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
                      'filename': config.state['qafilename'] + '_profiles',
                      'extension': '.pdf'}

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      np.ones(config.state['norders'], dtype=int),
                      qafileinfo=qafileinfo)

    # Set the continue flags

    config.state['continue'] = 3
