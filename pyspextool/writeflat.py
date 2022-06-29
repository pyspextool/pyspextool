from astropy.io import fits
import numpy as np

def writeflat(flat,var,flag,hdrlist,rotate,orders,edgecoeffs,xranges,\
              ps,slith_pix,slith_arc,slitw_pix,slitw_arc,modename,rms,rp,\
              version,history,oname,linmax=None,overwrite=True):

    '''
    To write a Spextool flat FITS file to disk.


    Input Parameters
    ----------------
    flat : numpy.ndarray of float
        The flat field image.


    var : numpy.ndarray of float
        The variance image.


    flag : numpy.ndarray of int
        A biset image.  


    rotate : {0, 1, 2, 3, 4, 5, 6, 7}, optional 
        Direction to rotate a raw image so that the dispersion direction 
        roughly aligns with the columns of the array and wavelength 
        increases to the right, and the spatial axis roughly aligns with 
        the rows of the array.  


        Direction  Transpose?  Rotation Counterclockwise
        -------------------------------------------------

        0          No          None
        1          No          90 deg
        2          No          180 deg
        3          No          270 deg
        4          Yes         None
        5          Yes         90 deg
        6          Yes         180 deg
        7          Yes         270 deg

        The directions follow the IDL rotate function convention.


    orders : list of int
        The order numbers.  orders[0] is the order number of the 
        order nearest the bottom of the image after rotation. 


    edgecoeffs : array_like of float
        (norders,ncoeffs+1,2) array giving the polynomial coefficients 
        delineating the top and bottom of each order.  edgecoeffs[0,0,:]
        gives the coefficients for the bottom of the order closest to the 
        bottom of the image and edgecoeffs[0,1,:] gives the coefficients 
        for the top of said order.  


    xranges : array_like of float
        An (norders,2) array giving the column numbers over which to 
        operate.  xranges[0,0] gives the starting column number for the 
        order nearest the bottom of the image and xranges[0,1] gives the 
        end column number for said order.


    ps : float
        The plate scale (arcseconds per pixel).


    slith_pix : float
        The slit height (pixels).


    slith_arc : float
        The nominal slit height (arcseconds).


    slitw_arc : float
        The slit width (arcseconds).


    slitw_pix : float
        The slit width (pixels).


    modename : str
        The name of the instrument mode.


    rms : list of float
        (norders,) list of RMS values for each order.


    rp: float
        The nominal resolving power.


    version : str
        The version of the pySpextool software used to create the flat.


    oname : str
        The filename of the flat field image to written to disk.  


    linmax : int, optional
        The linearity maximum used to identified pixels beyond the 
        linearity limit.


    history : list of str, optional
        The history string spliced to fit in a FITS file.


    overwrite : {True, False}, optional
        Set to True to overwrite an existing file.
        

    Returns
    -------
    None
        Writes a FITS file to PROCPATH.


    Notes
    -----
    None


    Examples
    --------
    Later


    Modification History
    --------------------
    2022-06-29 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program writeflat.pro.

    '''

# Get basic things

    norders = len(orders)
    edgedeg = edgecoeffs.shape[-1]
    
# Create the primary HDU
    
    phdu = fits.PrimaryHDU()
    hdr = phdu.header

# Add the hdrlist keywords and values
    
    keys = hdrlist.keys()
    for key in keys:

        if key == 'COMMENT':

            comments = hdrlist['COMMENT']
            
        elif key != 'HISTORY':  # probably not necessary, just in case
            
            hdr[key] = tuple(hdrlist[key])
    
# Now add new ones

        hdr['MODE'] = (modename,' Instrument Mode')
        hdr['NORDERS'] = (norders, ' Number of orders identified')
        hdr['ORDERS']=(','.join(str(o) for o in orders),'Orders identified')
        hdr['PLTSCALE'] = (ps, 'Plate scale (arcseconds per pixel)')
        hdr['SLTH_PIX'] = (slith_pix, ' Nominal slit length (pixels)')
        hdr['SLTH_ARC'] = (slith_arc, ' Slit length (arcseconds)')
        hdr['SLTW_PIX'] = (slitw_pix, ' Slit width (pixels)')
        hdr['SLTW_ARC'] = (slitw_arc, ' Slit width (arcseconds)')
        hdr['RP'] = (rp, ' Nominal resolving power')
        hdr['ROTATION'] = (rotate, ' IDL rotate value')
        hdr['VERSION'] = (version, ' Spextool version')

# Record linearity maximum if given 
        
    if linmax is not None:

        hdr['LINMAX'] = (linmax,' Linearity maximum')

# Add the RMS values.  Check to make sure not NaN

    if np.sum(np.isnan(rms)) == 0:
        
        for i in range(norders):

            name = 'OR'+str(orders[i]).zfill(3)+'_RMS'
            comment = ' RMS of normalized order '+str(orders[i]).zfill(3)
            hdr[name] = (rms[i],comment)

# Add the xranges

    for i in range(norders):

            name = 'OR'+str(orders[i]).zfill(3)+'_XR'
            comment = ' Extraction range for order '+str(orders[i]).zfill(3)
            hdr[name] = (','.join(str(x) for x in xranges[i,:]),comment)
        

# Add the edgecoeffs

    hdr['EDGEDEG']=(edgedeg,' Degree of the polynomial fit to order edges')
    for i in range(norders):

        for j in range(2):

            for k in range(edgedeg):

                if j == 0:

                    name = 'OR'+str(orders[i]).zfill(3)+'_B'+str(k+1)
                    comment = ' a'+str(k)+\
                      ' edge coefficient for bottom of order '+\
                      str(orders[i]).zfill(3)

                    hdr[name] = (edgecoeffs[i,j,k], comment)

                if j == 1:

                    name = 'OR'+str(orders[i]).zfill(3)+'_T'+str(k+1)
                    comment = ' a'+str(k)+\
                      ' edge coefficient for top of order '+\
                      str(orders[i]).zfill(3)

                    hdr[name] = (edgecoeffs[i,j,k], comment)

# Now add the comments 

    if 'comments' in locals():
                    
        for com in comments:

            hdr['COMMENT'] = com

# and then the history
            
    label = '         ============ Spextool History ============'
    hdr['HISTORY'] = label
    for hist in history:
            
        hdr['HISTORY'] = hist


# Write the results

    img_hdu = fits.ImageHDU(flat)
    var_hdu = fits.ImageHDU(var)
    flg_hdu = fits.ImageHDU(flag)
            
    hdu = fits.HDUList([phdu,img_hdu,var_hdu,flg_hdu])
    hdu.writeto(oname,overwrite=overwrite)
            
