import numpy as np
import logging
import matplotlib.pyplot as pl
from os.path import join


from pyspextool import config as setup
from pyspextool.combine import config as combine
from pyspextool.io.check import check_parameter, check_range, check_qakeywords
from pyspextool.utils.math import scale_data_stack
from pyspextool.combine.core import plot_allorders

from pyspextool.pyspextoolerror import pySpextoolError


def scale_spectra(order:int=None,
                  wavelength_range:list=None,
                  range_fraction:int | float=0.8,
                  verbose:bool=None,
                  qa_show:bool=None,
                  qa_showscale:float | int=None,
                  qa_showblock:bool=None,
                  qa_write:bool=None):

    """
    To scale the orders to a common intensity level

    Parameters
    ----------

    order : int
        The order used to determine the scale factors.

    wavelength_range : str
        The wavelength range over which the scale factors are determined in
        order `order`.

    range_fraction : int or float
        If `range` is None, then range_fraction of the wavelength range
        centered around the center of the order are used to determine the
        scale factors.
    
    Returns
    -------
    None
    
    """

    #
    # Check the parameters and QA keywords
    #

    check_parameter('scale_spectra', 'order', order , ['int', 'NoneType'])

    check_parameter('scale_spectra', 'wavelength_range', wavelength_range ,
                    ['list', 'NoneType'])    

    check_parameter('scale_spectra', 'range_fraction', range_fraction ,
                    ['int', 'float'])

    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)    
    
    logging.info(' Scaling the spectra to a common intensity level.')
    
    #
    # First determine which order we are using to scale
    #

    if order is None:

        # Do it ourselves.  Pick the middle-most order.

        order = int(np.median(combine.state['orders']))

    else:

        # Let the user choose, but make sure it is an order we have.

        z = combine.state['orders'] == order
        if np.sum(z) == 0:

            message = 'Order '+str(order)+' not in the files.'
            raise pySpextoolError(message)

    logging.info(' Using order '+str(order)+' to determine scale factors.')

    # Store the results

    combine.state['scale_order'] = order

    #
    # Now let's figure out the wavelength range over which to do the scaling.

    z_order = combine.state['orders'] == combine.state['scale_order']    
    
    wave = np.squeeze(combine.state['wavelengths'][0, z_order, :])

    min_wave = np.nanmin(wave)
    max_wave = np.nanmax(wave)        

    
    if wavelength_range is None:

        # Determine the wavelenth range ourselves.

        median_wave = np.nanmedian(wave)

        delta_wave = (max_wave-min_wave)*range_fraction/2

        scale_range = np.array([median_wave-delta_wave,median_wave+delta_wave])

        combine.state['scale_range'] = scale_range

    else:

        # Let the user choose, but make sure it falls in range.

        check_range(wavelength_range[0],[min_wave,max_wave],'gele',
                    variable_name='wavelength_range[0]')

        check_range(wavelength_range[1],[min_wave,max_wave],'gele',
                    variable_name='wavelength_range[1]')        
               
        combine.state['wavelength_range'] = wavelength_range
    
    #
    # Loop over each aperture and order
    #

    scales = np.empty((combine.state['final_napertures'],
                       combine.state['nspectra']))

    for i in range(combine.state['final_napertures']):

        # Determine which pixels we are using for the scaling

        wave = np.squeeze(combine.state['wavelengths'][i, z_order, :])
        intensity = np.squeeze(combine.state['intensities'][i,z_order,:,:])

        z_wave = np.where((wave > combine.state['scale_range'][0]) &
                          (wave < combine.state['scale_range'][1]))

        # Get scale factors

        junk, junk, scale = scale_data_stack(intensity[:, z_wave[0]], None)

        scales[i,:] = scale

        # Now scale each order

        shape = np.shape(intensity)                
        reshape_info = (shape[0], 1)
        tile_info = (1, shape[1])        

        scale_array = np.tile(np.reshape(scale, reshape_info), tile_info)

        for j in range(combine.state['norders']):

            np.multiply(combine.state['intensities'][i,j,:,:], scale_array,
                        out=combine.state['intensities'][i,j,:,:])

            np.multiply(combine.state['uncertainties'][i,j,:,:],
                        np.sqrt(scale_array),
                        out=combine.state['uncertainties'][i,j,:,:])

    combine.state['scales'] = scales
    combine.state['spectra_scaled'] = True

    #
    # Do the QA
    #

    if qa['write'] is True:

        plot_allorders(setup.plots['combine_spectra'],
                       setup.plots['portrait_size'],
                       setup.plots['font_size'],
                       plot_scalerange=True,
                       title='Scaled Spectra')

        
        pl.savefig(join(setup.state['qa_path'],combine.state['output_name']+\
                        '_scaled'+setup.state['qa_extension']))
        pl.close()

    if qa['show'] is True:

        scaled_size = (setup.plots['portrait_size'][0]*qa['showscale'],
                       setup.plots['portrait_size'][1]*qa['showscale'])

        scaled_font = setup.plots['font_size']*qa['showscale']
        
        plot_allorders(setup.plots['combine_spectra'],
                       scaled_size,
                       scaled_font,
                       plot_scalerange=True,
                       title='Scaled Spectra')                       
        
        pl.show(block=qa['showblock'])
        if qa['showblock'] is False: pl.pause(1)
    

        


    
    
    
