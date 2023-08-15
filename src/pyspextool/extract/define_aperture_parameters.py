import numpy as np

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.io.check import check_parameter
from pyspextool.io.check import check_range
from pyspextool.plot.plot_profiles import plot_profiles


def define_aperture_parameters(aperture_radii, psf_radius=None, bg_radius=None,
                               bg_width=None, bg_regions=None, bg_fit_degree=1,
                               qa_show=None, qa_showsize=(6, 10), qa_write=None):
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
        The polynomial degree of the fit to the background.

    qa_show : {None, True, False}, optional
        Set to True/False to override config.state['qa_show'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_showsize : tuple, default=(6,6)
        A (2,) tuple giving the plot size that is passed to matplotlib as,
        pl.figure(figsize=(qa_showsize)) for the interactive plot.

    qa_write : {None, True, False}, optional
        Set to True/False to override config.state['qa_write'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be written to disk.

    Returns
    -------
    None

    Updates various parameters in the config file in pyspextool.extract.

    """

    #
    # Check if we can proceed
    #

    if extract.state['trace_done'] is False:

        message = "extract.state['trace_done']=False.  "+\
          "Previous steps not complete."        
        raise ValueError(message)

    #
    # Check common parameters
    #

    check_parameter('define_aperture_parameters', 'bg_fit_degree',
                    bg_fit_degree, ['int'], [1, 2])

    check_parameter('define_aperture_parameters', 'qa_show', qa_show,
                    ['NoneType', 'bool'])

    check_parameter('define_aperture_parameters', 'qa_write', qa_write,
                    ['NoneType', 'bool'])

    check_parameter('define_aperture_parameters', 'qa_showsize', qa_showsize,
                    'tuple')

    #
    # Check the qa and verbose variables and set to system default if need be.
    #

    if qa_write is None:
        qa_write = setup.state['qa_write']

    if qa_show is None:
        qa_show = setup.state['qa_show']

    #
    # Store user inputs 
    #

    extract.parameters['apradii'] = aperture_radii
    extract.parameters['psfradius'] = psf_radius
    extract.parameters['bgradius'] = bg_radius
    extract.parameters['bgwidth'] = bg_width
    extract.parameters['bgregions'] = bg_regions
    extract.parameters['bgdeg'] = bg_fit_degree
    extract.parameters['qaplot'] = qa_show
    extract.parameters['qafile'] = qa_write
    extract.parameters['qaplotsize'] = qa_showsize

    #
    # Now things proceed depending on the extraction mode
    #

    if extract.state['type'] == 'ps':

        #
        # ======================= Point Source ===========================
        #

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

        extract.state['psfradius'] = psf_radius
        extract.state['bgradius'] = bg_radius
        extract.state['bgwidth'] = bg_width
        extract.state['bgfitdeg'] = bg_fit_degree

        #
        # Get the background psbginfo list together
        #

        if bg_radius is not None:

            psbginfo = [bg_radius, bg_width]

        else:

            psbginfo = None

        # Force xsbginfo to None 

        xsbginfo = None

        # Store radii 

        extract.state['apradii'] = aperture_radii

    else:

        #
        # ======================= Extended Source ===========================
        #

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

        aperture_radii = np.array(aperture_radii, ndmin=1)
        nradii = np.size(aperture_radii)

        if nradii != extract.state['naps']:
            message = 'Number of aperture radii must equal number apertures.'
            raise ValueError(message)

        extract.state['apradii'] = aperture_radii

        # Now deal with the background region

        if bg_regions is not None:

            # String or list?

            if isinstance(bg_regions, str) is True:

                # Split on commas

                groups = bg_regions.split(',')
                xsbginfo = []
                for regions in groups:
                    ranges = regions.split('-')
                    ranges = [float(i) for i in ranges]
                    xsbginfo.append(ranges)

            else:

                # Must be a list.  

                xsbginfo = bg_regions

            # Check to make sure things are in range

            check_range(xsbginfo, [0, extract.state['slith_arc']], 'gele',
                        variable_name='bg_regions')

            extract.state['bgregions'] = xsbginfo

        else:

            extract.state['bgregions'] = None

        # Store the results

        extract.state['bgfitdeg'] = bg_fit_degree

        # Force the psbginfo to None

        psbginfo = None

    #
    # Now do the plotting
    #

    if extract.state['type'] == 'xs':

        doorders = extract.state['xsdoorders']
        plot_aperture_radii = aperture_radii

    else:

        doorders = extract.state['psdoorders']
        plot_aperture_radii = np.full(extract.state['naps'], aperture_radii)

    if qa_show is True:
        
        number = plot_profiles(extract.state['profiles'],
                               extract.state['slith_arc'],
                               doorders, apertures=extract.state['apertures'],
                               aperture_radii=plot_aperture_radii,
                               psf_radius=psf_radius, ps_bginfo=psbginfo,
                               xs_bginfo=xsbginfo, plot_size=qa_showsize,
                               plot_number=extract.state['profiles_plotnum'])
        extract.state['profiles_plotnum'] = number
        
    if qa_write is True:
        qafileinfo = {'figsize': (8.5, 11),
                      'filepath': setup.state['qa_path'],
                      'filename': extract.state['qafilename'] + '_apertureparms',
                      'extension': setup.state['qa_extension']}

        plot_profiles(extract.state['profiles'], extract.state['slith_arc'],
                      doorders, apertures=extract.state['apertures'],
                      aperture_radii=plot_aperture_radii, psf_radius=psf_radius,
                      ps_bginfo=psbginfo, xs_bginfo=xsbginfo,
                      file_info=qafileinfo)

    #
    # Set the done variable
    #

    extract.state['parameters_done'] = True
