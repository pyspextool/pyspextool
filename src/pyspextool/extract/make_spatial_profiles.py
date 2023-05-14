import numpy as np
from scipy import interpolate

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_profiles import plot_profiles
from pyspextool.utils.math import mean_data_stack


def make_spatial_profiles(qa_plot=None, qa_file=None, qa_plotsize=(6, 10),
                          verbose=None):
    """
    To create 1D "average" spatial profiles of the orders.

    Parameters 
    ----------
    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(6, 10)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_plotsize)) for the interactive plot.

    qa_file : {None, True, False}, optional
        Set to True/False to override config.state['qa_file'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    verbose : {None, True, False}, optional
        Set to True/False to override config.state['verbose'] in the
        pyspextool config file.

    Returns 
    -------
    None
    Fills the config.extract['profiles'] variable and optionally creates a 
    QA file.

    """

    #
    # Check to make sure we can proceed.
    #

    if extract.state['type_done'] is False:
        message = 'Previous steps not completed.'
        print(message)
        return

    #
    # Check the parameters
    #

    check_parameter('make_spatial_profiles', 'qa_plot', qa_plot,
                    ['NoneType', 'bool'])

    check_parameter('make_spatial_profiles', 'verbose', verbose,
                    ['NoneType', 'bool'])

    check_parameter('make_spatial_profiles', 'qa_file', qa_file,
                    ['NoneType', 'bool'])

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_file is None:
        qa_file = setup.state['qa_file']

    if qa_plot is None:
        qa_plot = setup.state['qa_plot']

    if verbose is None:
        verbose = setup.state['verbose']

        #
    # Save user inputs
    #

    extract.profiles['qaplot'] = qa_plot
    extract.profiles['qafile'] = qa_file
    extract.profiles['verbose'] = verbose

    #
    # Build the profiles
    #

    if verbose is True:
        print('Creating the spatial profiles...')

    profiles = []
    for i in range(extract.state['norders']):

        # Unpack the data

        order = extract.state['rectorders'][i]

        img = order['img']
        y = order['y']
        w = order['w']

        nrows = len(y)

        # Subtract the background

        medbg = np.median(img, axis=0)

        bgimg = np.tile(medbg, (nrows, 1))

        np.subtract(img, bgimg, out=img)

        if extract.load['wavecalfile'] is not None:

            # Do the interpolation of the atmosphere

            f = interpolate.interp1d(extract.state['atrans_wave'],
                                     extract.state['atrans_trans'],
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

        profiles.append({'order': extract.state['orders'][i], 'y': y,
                         'p': np.flip(mean)})

    extract.state['profiles'] = profiles

    if qa_plot is True:
        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      plot_size=qa_plotsize)

    if qa_file is True:
        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': setup.state['qa_path'],
                      'filename': extract.state['qafilename'] + '_profiles',
                      'extension': setup.state['qa_extension']}

        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      np.ones(extract.state['norders'], dtype=int),
                      file_info=qafileinfo)

    # Set the done variable

    extract.state['profile_done'] = True
