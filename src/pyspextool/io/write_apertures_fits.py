import os
from astropy.io import fits
import numpy as np

from pyspextool.io.check import check_parameter
from pyspextool.utils.split_text import split_text


def write_apertures_fits(spectra, xranges, aimage, sky, flat, naps, orders,
                         header_info, aperture_positions, aperture_radii,
                         plate_scale, slith_pix, slith_arc, slitw_pix,
                         slitw_arc, resolving_power, xunits, yunits,
                         latex_xunits, latex_yunits, latex_xlabel,
                         latex_ylabel, version, output_fullpath,
                         wavecalinfo=None, psbginfo=None, xsbginfo=None,
                         optimal_info=None, badpixel_info=None,
                         lincormax=None, overwrite=True, verbose=True):

    """
    To write a spextool spectral FITS file to disk

    Parameters
    ----------

    Returns
    -------
    None.  Writes FITS files to disk.


    """

    #
    # Check parameters
    #
    check_parameter('write_apertures_fits', 'spectra', spectra, 'list')

    check_parameter('write_apertures_fits', 'xranges', xranges, 'ndarray')

    check_parameter('write_apertures_fits', 'aimage', aimage, 'str')

    check_parameter('write_apertures_fits', 'flat', flat, 'str')

    check_parameter('write_apertures_fits', 'naps', naps, 'int')

    check_parameter('write_apertures_fits', 'orders', orders,
                    ['int', 'list', 'ndarray'])

    check_parameter('write_apertures_fits', 'header_info', header_info,
                    'dict')

    check_parameter('write_apertures_fits', 'aperture_positions',
                    aperture_positions, 'ndarray')

    check_parameter('write_apertures_fits', 'aperture_radii', aperture_radii,
                    ['int', 'float', 'ndarray'])

    check_parameter('write_apertures_fits', 'plate_scale', plate_scale, 'float')

    check_parameter('write_apertures_fits', 'slith_pix', slith_pix, 'float')

    check_parameter('write_apertures_fits', 'slith_arc', slith_arc, 'float')

    check_parameter('write_apertures_fits', 'slitw_pix', slitw_pix, 'float')

    check_parameter('write_apertures_fits', 'slitw_arc', slitw_arc, 'float')

    check_parameter('write_apertures_fits', 'resolving_power',
                    resolving_power, 'float')

    check_parameter('write_apertures_fits', 'xunits', xunits, 'str')

    check_parameter('write_apertures_fits', 'yunits', yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xunits', latex_xunits, 'str')

    check_parameter('write_apertures_fits', 'latex_yunits', latex_yunits, 'str')

    check_parameter('write_apertures_fits', 'latex_xlabel', latex_xlabel, 'str')

    check_parameter('write_apertures_fits', 'latex_ylabel', latex_ylabel, 'str')

    check_parameter('write_apertures_fits', 'version', version, 'str')

    check_parameter('write_apertures_fits', 'output_fullpath',
                    output_fullpath, 'str')

    check_parameter('write_apertures_fits', 'psbginfo', psbginfo,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'xsbginfo', xsbginfo,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'optimal_info', optimal_info,
                    ['NoneType', 'dict'])

    check_parameter('write_apertures_fits', 'badpixel_info', badpixel_info,
                    ['NoneType', 'dict'])    

    check_parameter('write_apertures_fits', 'wavecalinfo', wavecalinfo,
                    ['dict', 'NoneType'])

    check_parameter('write_apertures_fits', 'lincormax', lincormax,
                    'NoneType')

    check_parameter('write_apertures_fits', 'overwrite', overwrite,
                    'bool')

    check_parameter('write_apertures_fits', 'verbose', verbose,
                    'bool')

    #    
    # Get set up
    #

    orders = np.asarray(orders)
    norders = len(orders)

    #
    # Get the list of spectra into an array format
    #

    # Determine the maximum size of each aperture array

    npixels = []
    for slice in spectra:
        npixels.append(np.shape(slice)[1])

    max_npixels = np.max(np.asarray(npixels))

    # Now create arrays into which the slices will be placed

    array = np.full((norders * naps, 4, max_npixels), np.nan)

    # Now fill in the arrays

    l = 0
    for slice in spectra:
        array[l, :, 0:npixels[l]] = slice
        l += 1

    #
    # Now write the file(s) to disk
    #

    # Create the headers

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    # Add the original file header keywords.  Cut out any history first
    # if need be.

    is_history_present = 'HISTORY' in header_info

    if is_history_present is True:

        old_history = header_info['HISTORY']

        # remove it from the header_info dictionary

        header_info.pop('HISTORY')

    keys = list(header_info.keys())
    
    for i in range(len(keys)):

        if keys[i] == 'COMMENT':

            junk = 1

        else:

            hdr[keys[i]] = (header_info[keys[i]][0], header_info[keys[i]][1])

    # Add spextool keywords

    hdr['MODULE'] = ('extract', ' Creation module')
    hdr['VERSION'] = (version, ' pySpextool version')
    hdr['AIMAGE'] = (aimage, ' A image')
    hdr['SKYORDRK'] = (sky, ' Sky or dark image')
    hdr['FLAT'] = (flat, ' Flat field image')

    if wavecalinfo is not None:
        hdr['WAVECAL'] = (wavecalinfo['file'], ' Wavecal file')
        hdr['WCTYPE'] = (wavecalinfo['wavecaltype'],
                         ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (wavecalinfo['wavetype'], ' Wavelength type')

    else:

        hdr['WAVECAL'] = (None, ' Wavecal file')
        hdr['WCTYPE'] = (None, ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (None, ' Wavelength type')
        

    hdr['NORDERS'] = (norders, ' Number of orders')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), ' Orders')
    hdr['NAPS'] = (naps, ' Number of apertures')

    hdr['PLTSCALE'] = (plate_scale, ' Plate scale (arcseconds pixel-1)')
    hdr['SLTH_ARC'] = (slith_arc, ' Nominal slit height (arcseconds)')
    hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
    hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit height (pixels)')
    hdr['SLTW_PIX'] = (slitw_pix, ' Slit width (pixels)')
    hdr['RP'] = (resolving_power, ' Average resolving power')

    # Add the aperture positions

    for i in range(norders):
        name = 'APOSO' + str(orders[i]).zfill(3)
        comment = ' Aperture positions (arcseconds) for order ' + \
                  str(orders[i]).zfill(3)
        val = ','.join([str(round(elem, 2)) for elem in aperture_positions[i]])
        hdr[name] = (val, comment)

    if optimal_info is not None:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')
        hdr['THRESH'] = (optimal_info['thresh'], ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (True, ' Optimal extraction?')
        hdr['SUMEXT'] = (False, 'Sum extraction?')
        hdr['PSFRAD'] = (optimal_info['psfradius'], ' PSF radius (arcseconds)')

                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        
        hdr['THRESH'] = (None, ' Bad pixel sigma threshold')
        hdr['OPTEXT'] = (False, ' Optimal extraction?')
        hdr['SUMEXT'] = (True, ' Sum extraction?') 
        hdr['PSFRAD'] = (None, ' PSF radius (arcseconds)')


    if badpixel_info is not None:

        hdr['BDPXFIX'] = (True, ' Bad pixels fixed?')        
        hdr['THRESH'] = (badpixel_info['thresh'], ' Bad pixel sigma threshold')
                
    else:

        hdr['BDPXFIX'] = (False, ' Bad pixels fixed?')        
                
    # Add the aperture radii

    if isinstance(aperture_radii, np.ndarray):

        val = ','.join([str(elem) for elem in aperture_radii])
        hdr['APRADII'] = (val, ' Aperture radii (arcseconds)')

    else:

        hdr['APRADII'] = (aperture_radii, ' Aperture radii (arcseconds)')

        # Add the background info

    if xsbginfo is not None:

        tmplist = []
        for val in xsbginfo['regions']:
            region = str(val[0]) + '-' + str(val[1])
            tmplist.append(region)

        hdr['BGREGS'] = (','.join(tmplist),
                         ' Background regions (arcseconds)')
        hdr['BGDEGREE'] = (xsbginfo['degree'],
                           ' Background polynomial fit degree')

    hdr['XUNITS'] = (xunits, ' Units of the x axis')
    hdr['YUNITS'] = (yunits, ' Units of the y axis')

    hdr['LXUNITS'] = (latex_xunits, ' LateX units of the x axis')
    hdr['LYUNITS'] = (latex_yunits, ' LateX units of the y axis')

    hdr['LXLABEL'] = (latex_xlabel, ' LateX x axis label')
    hdr['LYLABEL'] = (latex_ylabel, ' LateX Y axis label')

    # Now add the HISTORY

    history = 'Spextool FITS files contain an array of size ' \
              '[nwaves,4,norders*naps]. The ith image (array[*,*,i]) ' \
              'contains the data for a single extraction aperture within ' \
              'an order such that, lambda=array[*,0,i], flux=array[*,1,i], ' \
              'uncertainty=array[*,2,i],flag=array[*,3,i].  The zeroth ' \
              'image (array[*,*,0]) contains the data for the aperture in ' \
              'the order closest to the bottom of the detector that is ' \
              'closest to the bottom of the slit (i.e. also closest to the ' \
              'bottom of the detector).  Moving up the detector, the FITS ' \
              'array is filled in with subsequent extraction apertures.  ' \
              'If no orders have been deselected in the extraction process, ' \
              'the contents of the ith aperture in order j can be found as ' \
              'follows: lambda=array[*,0,{j-min(orders)}*naps + (i-1)], ' \
              'flux=array[*,1,{j-min(orders)}*naps + (i-1)], ' \
              'uncertainty=array[*,2,{j-min(orders)}*naps + (i-1)], ' \
              'flag=array[*,3,{j-min(orders)}*naps + (i-1)].'

    history = split_text(history, length=65)

    if is_history_present is True:

        # Join the two histories

        history = old_history + [' '] + history
    
    for hist in history:
        hdr['HISTORY'] = hist

    # Set the file name for the spectra file and write

    hdr['FILENAME'] = (os.path.basename(output_fullpath)+'.fits', ' File name')

    fits.writeto(output_fullpath + '.fits', array, hdr, overwrite=overwrite)

    #
    # Update the user
    #

    if verbose is True:

        if xsbginfo is None:
            print('Wrote', os.path.basename(output_fullpath)+'.fits',
                  'to disk.')

        if xsbginfo is not None:
            print('Wrote', os.path.basename(output_fullpath)+'.fits', 'and',
                  os.path.basename(output_fullpath)+'.fits', 'to disk.')

    #
    # Return the file name written to disk
    #

    return hdr['FILENAME']
