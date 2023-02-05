import numpy as np
import os

from pyspextool.extract import config
from pyspextool.extract.check_continue import check_continue
from pyspextool.io.check import check_parameter
from pyspextool.io.write_apertures_fits import write_apertures_fits
from pyspextool.spectroscopy.extract_extendedsource_1dxd import extract_extendedsource_1dxd


def extract_apertures(output_filenames=None, verbose=True):

    # Do the extraction
    
    if config.state['exttype'] == 'ps':

        #
        #========================= Point Source ============================
        #

        x = 1

    else:

        #
        #======================= Extended Source ===========================
        #

        #
        # Check parameters
        #
        check_parameter('extract_apertures', 'output_filenames',
                        output_filenames, ['str', 'list'])

        # Get the File names and make list if it isn't one.
        
        full_filenames = output_filenames
        if isinstance(full_filenames, str):
            full_filenames = [full_filenames]

        # deal with background

        if config.state['bgregions'] is not None:

            bginfo = {'regions': config.state['bgregions'],
                      'degree': config.state['bgfitdeg']}

        else:

            bginfo = None

        # Grab which orders are being extracted
            
        z = config.state['xsdoorders'] == 1

        # Do the extraction
        
        spectra = extract_extendedsource_1dxd(config.state['workimage'],
                                              config.state['varimage'],
                                              config.state['ordermask'],
                                              config.state['orders'][z],
                                              config.state['wavecal'],
                                              config.state['spatcal'],
                                              config.state['apertures'][z],
                                              config.state['apradii'],
                                              bginfo=bginfo)
        
        #
        # Write the file to disk
        #

        # Let's get set up

        if config.state['wavecalfile'] != None:

            wavecalinfo = {'file':config.state['wavecalfile'],
                           'wavecaltype':config.state['wavecaltype'],
                           'wavetype':'vacuum'}                           

        
        if config.state['reductionmode'] == 'A':

            x = 1

        
        elif config.state['reductionmode'] == 'A-B':

            # Get file information
            
            hdrinfo = config.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = config.state['hdrinfo'][1]['FILENAME'][0]
            output_fullpath = os.path.join(config.state['procpath'],
                                           full_filenames[0])

        elif config.state['reductionmode'] == 'A-Sky/Dark':

            x = 1
            
        else:

            print('Reduction Mode Unknown.')

        # Actually do it!
            
        write_apertures_fits(spectra['spectra'], config.state['xranges'][z,:],
                             aimage, sky, config.state['flatfile'],
                             config.state['naps'],
                             config.state['orders'][z], hdrinfo,
                             config.state['apertures'][z],
                             config.state['apradii'],
                             config.state['plate_scale'],
                             config.state['slith_pix'],
                             config.state['slith_arc'],
                             config.state['slitw_pix'],
                             config.state['slitw_arc'],
                             config.state['resolvingpower'],
                             'um', 'DN s-1', 'microns', 'flux',
                             config.state['version'],
                             output_fullpath,
                             wavecalinfo=wavecalinfo,
                             background_spectra=spectra['background'],
                             xsbginfo=bginfo, verbose=verbose)

