import numpy as np

from pyspextool.cl import config
from pyspextool.cl.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.plot.plot_profiles import plot_profiles

def define_aperture_parameters(aperture_radii, psf_radius=None, bg_radius=None,
                               bg_width=None, bg_regions=None, iplot=True,
                               qafile=False):

    #
    # Continue status
    #

    check_continue(4)

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

            message = '`bg_radius` and `bg_width` required if `psf_radius` '+\
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

        config.state['psfradius'] = psf_radius
        config.state['bgradius'] = bg_radius
        config.state['bgwidth'] = bg_width

        #
        # Get set up for plotting
        #

        aperture_radii = np.full(config.state['naps'], aperture_radii)
        doorders = config.state['psdoorders']
        if bg_radius is not None:

            psbginfo = [bg_radius, bg_width]

        else:

            psbginfo = None

        xsbginfo = None
                
    else:
        
        #
        # Check parameters
        #

        check_parameter('define_aperture_parameters', 'bg_regions',
                        bg_regions, ['list', 'NoneType'])




        



        # Check to make sure the aperture_radii is a single number.


        # Now we know we have both bg_width and bg_radius.
        

        
#    #
#    # Store the results
#    #

#    
#    if iplot is True or qafile is True:
#
#        # Get set up for the plotting.
#
#        if config.state['exttype'] == 'ps':
#
#            

#                
#
#        else:
#
#            doorders = config.state['xsdoorders']             
#
    if config.state['exttype'] == 'xs':
        
        doorders = config.state['xsdoorders']
        
    else:

        doorders = config.state['psdoorders'] 
        
    if iplot is True:

        plot_profiles(config.state['profiles'],config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'],
                      aperture_radii= aperture_radii, psf_radius=psf_radius,
                      psbginfo=psbginfo, xsbginfo=xsbginfo)
    

    if qafile is True:

        qafileinfo = {'figsize': (8.5,11), 'filepath':config.state['qapath'],
                      'filename':config.state['qafilename']+'_apertureparms',
                      'extension':'.pdf'}

        plot_profiles(config.state['profiles'], config.state['slith_arc'],
                      doorders, apertures=config.state['apertures'],
                      aperture_radii=aperture_radii, psf_radius=psf_radius,
                      psbginfo=psbginfo, xsbginfo=xsbginfo,
                      qafileinfo=qafileinfo)


    
