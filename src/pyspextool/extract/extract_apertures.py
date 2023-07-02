import numpy as np
import os

from pyspextool import config as setup
from pyspextool.extract import config as extract
from pyspextool.extract.extract_pointsource_1dxd import extract_pointsource_1dxd
from pyspextool.extract.extract_extendedsource_1dxd import extract_extendedsource_1dxd
from pyspextool.io.check import check_parameter
from pyspextool.io.write_apertures_fits import write_apertures_fits
from pyspextool.plot.plot_spectra import plot_spectra


def extract_apertures(qa_plot=None, qa_plotsize=(10, 6), qa_file=None,
                      verbose=None):
    """
    User function to extract spectra.

    Parameters
    ----------
    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(6,6)
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

    """

    #
    # Check if we can proceed.
    #

    if extract.state['parameters_done'] is False:
        message = 'Previous steps not completed.'
        print(message)
        return

    #
    # Check parameters.  Just do it now to not waste extracting.
    #

    check_parameter('extract_apertures', 'qa_plot', qa_plot,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'qa_plotsize', qa_plotsize,
                    'tuple')

    check_parameter('extract_apertures', 'qa_file', qa_file,
                    ['NoneType', 'bool'])

    check_parameter('extract_apertures', 'verbose',
                    verbose, ['NoneType', 'bool'])

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
    # Store user inputs
    #

    extract.extract['qafile'] = qa_file
    extract.extract['qaplot'] = qa_plot
    extract.extract['verbose'] = verbose

    #
    # The procedure depends on the extraction type, point or extended source
    #

    xsbginfo = None
    psbginfo = None

    if extract.state['type'] == 'ps':

        #
        # ========================= Point Source ============================
        #

        # Grab which orders are being extracted

        z = extract.state['psdoorders'] == 1

        # Create the background information dictionary

        if extract.state['bgradius'] is not None:
            psbginfo = {'radius': extract.state['bgradius'],
                        'width': extract.state['bgwidth'],
                        'degree': extract.state['bgfitdeg']}

        # Do the extraction

        spectra = extract_pointsource_1dxd(extract.state['workimage'],
                                           extract.state['varimage'],
                                           extract.state['ordermask'],
                                           extract.state['orders'][z],
                                           extract.state['wavecal'],
                                           extract.state['spatcal'],
                                           extract.state['tracecoeffs'],
                                           extract.state['apradii'],
                                           extract.state['apsigns'],
                                           bginfo=psbginfo, verbose=verbose)

    else:

        #
        # ======================= Extended Source ===========================
        #

        # Grab which orders are being extracted

        z = extract.state['xsdoorders'] == 1

        # Create background information dictionary

        if extract.state['bgregions'] is not None:
            xsbginfo = {'regions': extract.state['bgregions'],
                        'degree': extract.state['bgfitdeg']}

        # Do the extraction

        spectra = extract_extendedsource_1dxd(extract.state['workimage'],
                                              extract.state['varimage'],
                                              extract.state['ordermask'],
                                              extract.state['orders'][z],
                                              extract.state['wavecal'],
                                              extract.state['spatcal'],
                                              extract.state['apertures'][z],
                                              extract.state['apradii'],
                                              bginfo=xsbginfo, verbose=verbose)

    #
    # Write the results to disk
    #

    write_apertures(spectra, psbginfo=psbginfo, xsbginfo=xsbginfo,
                    qa_file=qa_file, qa_plot=qa_plot, qa_plotsize=qa_plotsize,
                    verbose=verbose)

    #
    # Set the done variable
    #

    extract.state['extract_done'] = True


def write_apertures(spectra, psbginfo=None, xsbginfo=None, qa_plot=None,
                    qa_file=None, qa_plotsize=(10, 6), verbose=None):
    """
    To write extracted spectra to disk.

    Parameters
    ----------
    spectra : dict
        The dictionary returned from any of the extract_ functions.
        
        `"spectra"` : list
            (norders*naps,) list of spectral planes of the extracted target.

        `"background"` : list
            (norders*naps,) list of spectral planes of the background of the 
            target.

    psbginfo : dict, optional
        `"radius"` : float
            The background start radius in arcseconds

        `"width"` : float
            The background annulus width in arcseconds

        `"degree"` : int
            The background polynomial fit degree.

    xsbginfo : dict, optional

        `"regions"` : str
            A comma-separated list of background regions, e.g. '1-2,14-15'.

        `"degree"` : int
            The background polynomial fit degree.

    qa_plot : {None, True, False}, optional
        Set to True/False to override config.state['qa_plot'] in the
        pyspextool config file.  If set to True, quality assurance
        plots will be interactively generated.

    qa_plotsize : tuple, default=(10,6)
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

       Writes FITS files to disk.

    """

    #
    # Get the plotting ready
    #

    qafileinfo = {'figsize': (11, 8.5), 'filepath': setup.state['qa_path'],
                  'filename': extract.state['qafilename'],
                  'extension': setup.state['qa_extension']}

    # Unpack the spectra dictionary
    #if spectra['background'] is not None:
    #    background = spectra
    #else:
    #    background = None

    spectra = spectra['spectra']

    # Get the wavelengh calibration info if used.

    if extract.load['wavecalfile'] is not None:

        wavecalinfo = {'file': extract.load['wavecalfile'],
                       'wavecaltype': extract.state['wavecaltype'],
                       'wavetype': 'vacuum'}

    else:

        wavecalinfo = None

    # Things proceed depending on the extraction mode

    if extract.state['type'] == 'ps':

        #
        # ========================= Point Source ============================
        #

        # Grab things all modes need

        xranges = extract.state['xranges'][extract.state['psdoorders']]
        orders = extract.state['orders'][extract.state['psdoorders']]
        apertures = extract.state['apertures'][extract.state['psdoorders']]
        norders = np.sum(extract.state['psdoorders'])

        # Now go mode by mode

        if extract.state['reductionmode'] == 'A':

            # Get file information

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'
            output_fullpath = extract.state['output_files'][0]

            # Write ther results

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 setup.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 psbginfo=psbginfo,
                                 verbose=verbose)

            # Plot the spectra

            filename = os.path.basename(output_fullpath)            
            if qa_plot is True:
                
                number = plot_spectra(output_fullpath + '.fits',
                                      title=filename+'.fits',                                                         plot_size=qa_plotsize,
                            plot_number=extract.state['spectra_a_plotnum'])

                extract.state['spectra_a_plotnum'] = number

            if qa_file is True:

                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-B':

            # Are there positive apertures?

            z_pos = extract.state['apsigns'] == 1
            pos_naps = int(np.sum(z_pos))
            if pos_naps != 0:

                # Get file information

                hdrinfo = extract.state['hdrinfo'][0]
                aimage = hdrinfo['FILENAME'][0]
                sky = extract.state['hdrinfo'][1]['FILENAME'][0]
                output_fullpath = extract.state['output_files'][0]

                # Now get the indices for the spectra array

                full_apsign = np.tile(extract.state['apsigns'], norders)
                z = (np.where(full_apsign == 1))[0]
                pos_spectra = [spectra[i] for i in z]

                write_apertures_fits(pos_spectra, xranges, aimage, sky,
                                     extract.load['flatfile'], pos_naps,
                                     orders, hdrinfo, apertures,
                                     extract.state['apradii'],
                                     extract.state['plate_scale'],
                                     extract.state['slith_pix'],
                                     extract.state['slith_arc'],
                                     extract.state['slitw_pix'],
                                     extract.state['slitw_arc'],
                                     extract.state['resolvingpower'],
                                     'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                     'Wavelength ($\mu$m)',
                                     'Count Rate (DN s$^{-1}$)',
                                     setup.state['version'],
                                     output_fullpath, wavecalinfo=wavecalinfo,
#                                     background_spectra=background[z],
                                     psbginfo=psbginfo,
                                     verbose=verbose)

                # Plot the spectra

                filename = os.path.basename(output_fullpath)
                if qa_plot is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                 title=filename+'.fits',
                                 plot_size=qa_plotsize,
                                 plot_number=extract.state['spectra_a_plotnum'])

                    extract.state['spectra_a_plotnum'] = number

                    
                if qa_file is True:

                    qafileinfo['filename'] = filename
                    plot_spectra(output_fullpath + '.fits',
                                 file_info=qafileinfo)

            # Are there negative apertures?

            z_neg = extract.state['apsigns'] == -1
            neg_naps = int(np.sum(z_neg))
            if neg_naps != 0:

                # Get file information

                hdrinfo = extract.state['hdrinfo'][1]
                aimage = hdrinfo['FILENAME'][1]
                sky = extract.state['hdrinfo'][0]['FILENAME'][0]
                output_fullpath = extract.state['output_files'][1]

                # Now get the indices for the spectra array

                full_apsign = np.tile(extract.state['apsigns'], norders)
                z = (np.where(full_apsign == -1))[0]
                neg_spectra = [spectra[i] for i in z]

                write_apertures_fits(neg_spectra, xranges, aimage, sky,
                                     extract.load['flatfile'],
                                     neg_naps, orders, hdrinfo, apertures,
                                     extract.state['apradii'],
                                     extract.state['plate_scale'],
                                     extract.state['slith_pix'],
                                     extract.state['slith_arc'],
                                     extract.state['slitw_pix'],
                                     extract.state['slitw_arc'],
                                     extract.state['resolvingpower'],
                                     'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                     'Wavelength ($\mu$m)',
                                     'Count Rate (DN s$^{-1}$)',
                                     setup.state['version'],
                                     output_fullpath, wavecalinfo=wavecalinfo,
#                                     background_spectra=background[z],
                                     psbginfo=psbginfo,
                                     verbose=verbose)

                # Plot the spectra

                filename = os.path.basename(output_fullpath)
                if qa_plot is True:

                    number = plot_spectra(output_fullpath + '.fits',
                                          title=filename+'.fits',
                                          plot_size=qa_plotsize,
                                plot_number=extract.state['spectra_b_plotnum'])
                    extract.state['spectra_b_plotnum'] = number
                    
                if qa_file is True:

                    qafileinfo['filename'] = filename
                    plot_spectra(output_fullpath + '.fits',
                                 file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-Sky/Dark':

            x = 1

        else:

            print('Reduction Mode Unknown.')

    else:

        #
        # ======================= Extended Source ===========================
        #

        # Grab things all modes need

        xranges = extract.state['xranges'][extract.state['psdoorders']]
        orders = extract.state['orders'][extract.state['psdoorders']]
        apertures = extract.state['apertures'][extract.state['psdoorders']]

        if extract.state['reductionmode'] == 'A':

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = 'None'
            output_fullpath = extract.state['output_files'][0]
            
            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 setup.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

            # Plot the spectra

            if qa_plot is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_plotsize)

            if qa_file is True:
                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-B':

            hdrinfo = extract.state['hdrinfo'][0]
            aimage = hdrinfo['FILENAME'][0]
            sky = extract.state['hdrinfo'][1]['FILENAME'][0]
            output_fullpath = extract.state['output_files'][0]            

            write_apertures_fits(spectra, xranges, aimage, sky,
                                 extract.load['flatfile'],
                                 extract.state['naps'],
                                 orders, hdrinfo, apertures,
                                 extract.state['apradii'],
                                 extract.state['plate_scale'],
                                 extract.state['slith_pix'],
                                 extract.state['slith_arc'],
                                 extract.state['slitw_pix'],
                                 extract.state['slitw_arc'],
                                 extract.state['resolvingpower'],
                                 'um', 'DN s-1', '$\mu$m', 'DN s$^{-1}$',
                                 'Wavelength ($\mu$m)',
                                 'Count Rate (DN s$^{-1}$)',
                                 extract.state['version'],
                                 output_fullpath, wavecalinfo=wavecalinfo,
#                                 background_spectra=background[z],
                                 xsbginfo=xsbginfo,
                                 verbose=verbose)

            # Plot the spectra

            if qa_plot is True:
                plot_spectra(output_fullpath + '.fits', title=output_fullpath,
                             plot_size=qa_plotsize)

            if qa_file is True:
                qafileinfo['filename'] = os.path.basename(output_fullpath)
                plot_spectra(output_fullpath + '.fits', file_info=qafileinfo)

        elif extract.state['reductionmode'] == 'A-Sky/Dark':

            x = 1

        else:

            print('Reduction Mode Unknown.')

    print(' ')
