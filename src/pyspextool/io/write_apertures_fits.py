from astropy.io import fits

import numpy as np
from pyspextool.io.check import check_parameter

def write_apertures_fits(spectra, xranges, aimage, sky, flat, naps, orders, 
                         hdrinfo, aperture_positions, apradii, plate_scale,
                         slith_pix, slith_arc, slitw_pix, slitw_arc,
                         resolving_power, xunits, yunits, xtitle, ytitle,
                         version, output_fullpath, wavecalinfo=None,
                         psinfo=None, xsinfo=None, lincormax=None,
                         overwrite=True):

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

    check_parameter('write_apertures_fits', 'hdrinfo', hdrinfo,
                    'dict')

    check_parameter('write_apertures_fits', 'aperture_positions',
                    aperture_positions, 'ndarray')

    check_parameter('write_apertures_fits', 'apradii', apradii, 'ndarray')

    check_parameter('write_apertures_fits', 'plate_scale', plate_scale, 'float')
    
    check_parameter('write_apertures_fits', 'slith_pix', slith_pix, 'float')

    check_parameter('write_apertures_fits', 'slith_arc', slith_arc, 'float')    

    check_parameter('write_apertures_fits', 'slitw_pix', slitw_pix, 'float')

    check_parameter('write_apertures_fits', 'slitw_arc', slitw_arc, 'float')

    check_parameter('write_apertures_fits', 'resolving_power',
                    resolving_power, 'float')

    check_parameter('write_apertures_fits', 'xunits', xunits, 'str')

    check_parameter('write_apertures_fits', 'yunits', yunits, 'str')

    check_parameter('write_apertures_fits', 'xtitle', xtitle, 'str')

    check_parameter('write_apertures_fits', 'ytitle', ytitle, 'str')

    check_parameter('write_apertures_fits', 'version', version, 'str')

    check_parameter('write_apertures_fits', 'output_fullpath',
                    output_fullpath, 'str')

    check_parameter('write_apertures_fits', 'psinfo', psinfo, 'NoneType')

    check_parameter('write_apertures_fits', 'xsinfo', xsinfo, 'NoneType')

    check_parameter('write_apertures_fits', 'wavecalinfo', wavecalinfo,
                    ['dict', 'NoneType'])

    check_parameter('write_apertures_fits', 'lincormax', lincormax,
                    'NoneType')

    check_parameter('write_apertures_fits', 'overwrite', overwrite,
                    'bool')                                                    
    


    #    
    # Get set up
    #
    
    orders = np.asarray(orders)
    norders = len(orders)

    # Determine the maximum size of each aperture array

    npixels = []
    for slice in spectra:

        npixels.append(np.shape(slice)[1])

    max_npixels = np.max(np.asarray(npixels))

    # Now create an array into which the slices will be placed

    array = np.full((norders*naps, 4, max_npixels), np.nan)

    # Now fill this array

    
    l = 0
    for slice in spectra:

        
        shape = np.shape(slice)
        array[l, :, 0:npixels[l]] = slice
        l += 1

    # Now write the file

    phdu = fits.PrimaryHDU()
    hdr = phdu.header

    hdr['CREPROG'] = (' xspextool', ' Creation program')
    hdr['VERSION'] = (version, ' Spextool version')
    hdr['AIMAGE'] = (aimage, ' A image')
    hdr['SKYORDRK'] = (sky, ' Sky or dark image')
    hdr['FLAT'] = (flat, ' Flat field image')

    if wavecalinfo is not None:

        hdr['WAVECAL'] = (wavecalinfo['file'], ' Wavecal file')
        hdr['WCTYPE'] = (wavecalinfo['wavecaltype'],
                         ' Wavelength calibration type ')
        hdr['WAVETYPE'] = (wavecalinfo['wavetype'], ' Wavelength type')
        
    
    hdr['NORDERS'] = (norders, ' Number of orders')
    hdr['ORDERS'] = (','.join(str(o) for o in orders), ' Orders')
    hdr['NAPS'] = (naps, ' Number of apertures')


    hdr['PLTSCALE'] = (plate_scale, ' Plate scale (arcseconds pixel-1)')
    hdr['SLTH_ARC'] = (slith_arc, ' Slit height (arcseconds)')
    hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
    hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit height (pixels)')
    hdr['SLTW_PIX'] = (slitw_pix, ' Nominal slit width (pixels)')        
    
    fits.writeto(output_fullpath, array, hdr, overwrite=True)


