import numpy as np

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.io.check import check_range
from pyspextool.plot.plot_profiles import plot_profiles


def define_aperture_parameters(aperture_radii, psf_radius=None, bg_radius=None,
                               bg_width=None, bg_regions=None, bg_fit_degree=1,
                               iplot=False, qafile=False):

    """
    To define the extraction parameters

    Parameters
    ----------
    aperture_radii : int, float, or list
        The aperture radii.
        If point source extraction, a single value giving the radius for all 
        apertures.
        If extended source extraction, (naps,) list of aperture radii. 

    psf_radius : float or None, optional
        If point source extraction, will perform an optimal extraction 
        if provided.  

    bg_radius : float or None, optional
        If point source extraction, the background radius value.

    bg_width : float or None, optional
        If point source extraction, the background width value.

    bg_regions : list, optional
        If extended source extraction, the background regions.

    bg_fit_degree : int, default=1, optional
        The polynomial degree of the fit to the background

    iplot : {False, True}, optional
        Set to plot the results interactively.

    qafile : {False, True}, optional
        Set to write the QA plot to disk.

    Returns
    -------
    None
    Updates the config.state['psfradius'], config.state['bgradius'], 
    config.state['bgwidth'] and, config.state['bgfitdeg'] variables and 
    optional plots the results.

    Notes
    -----
    None

    Examples
    --------
    define_apertures_parameters(1.5)         - point source
    define_apertures_parameters([1.5, 2])    - extended source

    """
    
    #
    # Continue status
    #

    check_continue(4)

    # Check commone parameters
    
    check_parameter('define_aperture_parameters', 'bg_fit_degree',
                    bg_fit_degree, ['int'], [1, 2])

    # Now things proceed depending on the extraction mode
    
    if config.state['exttype'] == 'ps':

        #
        # Check parameters
        #

        check_parameter('define_aperture_parameters', 'aperture_radii',
                        aperture_radii, ['int', 'float'])

        check_parameter('define_aperture_parameters', 'psf_radius',
                        psf_radius, ['int', 'float', 'NoneType'])

        check_parameter('define_aperture_parameters', 'bg_radius',
                        bg_radius, ['int', 'float', 'NoneType'])

        check_parameter('define_aperture_parameters', 'bg_width',
                        bg_width, ['int', 'float', 'NoneType'])

        #
        # Make sure the right sets of things are present
        #

        if bg_radius is not None and bg_width is None:
            message = '`bg_width` required if `bg_radius` is passed.'
            raise ValueError(message)

        if bg_width is not None and bg_radius is None:
            message = '`bg_radius` required if `bg_width` is passed.'
            raise ValueError(message)

        if psf_radius is not None and (bg_width is None or bg_radius is None):
            message = '`bg_radius` and `bg_width` required if `psf_radius` ' + \
                      'is passed.'
            raise ValueError(message)

        #
        # Now confirm values are "correct"
        #

        if bg_radius is not None:

            print('Updating this section to use check_range')
            if bg_radius <= aperture_radii:
                message = '`bg_radius` must be > `aperture_radii`.'
                raise ValueError(message)

        if psf_radius is not None:

            # Make sure it is larger than aperture_radii

            if psf_radius < aperture_radii:
                message = '`psf_radius` must be >= `aperture_radii`.'
                raise ValueError(message)

        #
        # Now store the results
        #
        
        config.state['psfradius'] = psf_radius
        config.state['bgradius'] = bg_radius
        config.state['bgwidth'] = bg_width
        config.state['bgfitdeg'] = bg_fit_degree        
            
        #
        # Get set up for plotting
        #

        aperture_radii = np.full(config.state['naps'], aperture_radii)
        doorders = config.state['psdoorders']
        if bg_radius is not None:

            psbginfo = [bg_radius, bg_width]

        else:

            psbginfo = None

        # Force xsbginfo to None 
            
        xsbginfo = None

    else:

        #
        # Check parameters
        #

        check_parameter('define_aperture_parameters', 'aperture_radii',
                        aperture_radii, ['int', 'float', 'list'])
        
        check_parameter('define_aperture_parameters', 'bg_regions',
                        bg_regions, ['list', 'str', 'NoneType'])

        #
        # Do some user checking
        #

        # Check number of apertures and radii are equal
        
        aperture_radii = np.array(aperture_radii)
        nradii = np.size(aperture_radii)

        if nradii != config.state['naps']:
            message = 'Number of aperture radii must equal number apertures.'
            raise ValueError(message)

        # Now deal with the background region
        
        if bg_regions is not None:

            # String or list?
            
            if type(bg_regions).__name__ == 'str':

                # Split on commas

                groups = bg_regions.split(',')
                xsbginfo = []
                for regions in groups:

                    ranges = regions.split('-')
                    ranges = [float(i) for i in ranges]
                    xsbginfo.append(ranges)

            else:

                # Must be a list.  Check to ensure in range.  
                
                xsbginfo = bg_regions
                
            # Check to make sure things are in range
                
            check_range(xsbginfo,
                        [0, config.state['slith_arc']], 'gele',
                        variable_name='bg_regions')

        # Store the results

        config.state['apradii'] = aperture_radii
        config.state['bgregions'] = xsbginfo
        config.state['bgfitdeg'] = bg_fit_degree
            
        # Force the psbginfo to None
            
        psbginfo = None
        
    if config.state['exttype'] == 'xs':

        doorders = config.state['xsdoorders']

    else:

        doorders = config.state['psdoorders']

    if iplot is True:

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'],
                      aperture_radii=aperture_radii, psf_radius=psf_radius,
                      psbginfo=psbginfo, xsbginfo=xsbginfo)

    if qafile is True:

        qafileinfo = {'figsize': (8.5, 11), 'filepath': config.state['qapath'],
                      'filename': config.state['qafilename'] + '_apertureparms',
                      'extension': config.state['qaextension']}

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'],
                      aperture_radii=aperture_radii, psf_radius=psf_radius,
                      psbginfo=psbginfo, xsbginfo=xsbginfo,
                      qafileinfo=qafileinfo)
