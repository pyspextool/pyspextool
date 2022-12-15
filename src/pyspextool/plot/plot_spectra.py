import numpy as np
import matplotlib.pyplot as pl
import os 

from pyspextool.fit.robust_savgol import robust_savgol
from pyspextool.plot.limits import get_spec_range
from pyspextool.io.check import check_parameter

def plot_spectra(spectra, orders, ap=None, xtitle=None, ytitle=None,
                 ladder=True, continuous=False, intensity=True,
                 uncertainty=False, snr=False, qafileinfo=None):

    #
    # check parameters
    #
    check_parameter('plot_spectra', 'spectra', spectra, 'list')

    check_parameter('plot_spectra', 'orders', orders, ['int', 'list',
                                                        'ndarray'])

    check_parameter('plot_spectra', 'xtitle', xtitle, ['str', 'NoneType'])

    check_parameter('plot_spectra', 'ytitle', ytitle, ['str', 'NoneType'])

    check_parameter('plot_spectra', 'ladder', ladder, 'bool')

    check_parameter('plot_spectra', 'continuous', continuous, 'bool')

    check_parameter('plot_spectra', 'qafileinfo', qafileinfo,
                     ['NoneType', 'dict'])                            

    #
    # Get some prep work done
    #
    
    # Get the number of orders and apertures

    orders = np.asarray(orders)
    norders = np.size(orders)
    naps = len(spectra)//norders

    # Must check naps versus aps if not None
    
    # Check the plot types

    sum = intensity+uncertainty+snr
    if sum > 1:

        message = 'Only one of `intensity`, `uncertainty`, `snr` can be set.'
        raise ValueError(message)

    sum = continuous+ladder
    if sum > 1:

        message = 'Only one of `ladder`, `continuous` can be set.'
        raise ValueError(message)

    # Get the figure size
        
    if qafileinfo is not None:

        figsize = qafileinfo['figsize']

    else:

        figsize = (8.5,11)
    

#    pl.figure(figsize=figsize)
#    pl.subplots_adjust(hspace=2)


    # Get aperture loop set up

    range_value = naps if ap is None else ap-1
        
    for i in range(range_value):

        






        

    
    # Get the plotting ranges

#    plot_ranges = {}
#    for i in range(nkeys):

#        all = spectra

    

    
    
    
    # Deal with the possibility that both ladder and continuous are set.

    if ladder is True and continuous is True:

        continuous = False

    if ladder is True:

        x = 1

        
    
