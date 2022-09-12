import numpy as np
import os
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate

from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.plot.limits import get_spec_range
from pyspextool.utils.arrays import trim_nan
from pyspextool.utils.loop_progress import loop_progress

def find_lines_1dxd(spectra, orders, line_info, pix_thresh, qafileinfo=None,
                    clupdate=True):

    """
    To find the location of arc lines in cross-dispersed spectra.

    The "1dxd" label means the function can be used to identify lines 
    in spectra where the arrays columns are aligned with wavelength 
    (the 1d part) and can be crossed dispersed spectra (the xd part).


    Parameters
    ----------
    spectra:  dict
        An (norder,) dictionary where each key has the form ORNNN_AP01
        where NNN is the order number.  The value of each entry is a
        (4, nwave) numpy.ndarray where:
        wavelength = (0,:)
        intensity = (1,:)
        uncertainty = (2,:)
        flags = (3,:)

    orders : array_like
        An (norders,) array of int giving the order numbers.  By convention,
        orders[0] is the order closest to the bottom of the array.

    line_info : dict

    pix_thresh: int
        The threshold (in pixels) beyond which an identification is deemed 
        bad.  That is, the fit is deemed bad if 
        (abs(fit[1]-guess) > pix_thresh.


    Returns
    -------
    list
         A list of integers giving the individual file numbers


    Notes
    -----


    Examples
    --------

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program XS

    """


    # Get basic information
    
    norders = len(orders)
    nlines = len(line_info['wavelength'])

    # Setup the output arrays

    line_xpos = np.full(nlines, np.nan)
    line_fwhm = np.full(nlines, np.nan)
    line_inten = np.full(nlines, np.nan)
    line_goodbad = np.zeros(nlines, dtype=np.uint8)

    # Setup the qaplot if asked

    if qafileinfo is not None:

        pl.rcParams['font.size'] = '12'
        pl.rcParams['lines.linewidth'] = 0.75

        pdf = PdfPages(os.path.join(qafileinfo['filepath'],
                                    qafileinfo['filename'] + \
                                    '_findlines') + \
                                    qafileinfo['extension'])
    # Loop over each line

    for i in range(nlines):

        # Find the order associated with the line
        
        z  = orders == line_info['order'][i]
        if np.sum(z) == 0:
            continue

        # get spectra for this order and clip
    
        key = 'OR'+str(line_info['order'][i]).zfill(3)+'_AP01'
        x = trim_nan(spectra[key][0,:],2,trim=True)
        y = trim_nan(spectra[key][1,:],2,trim=True)        
        
        # Cut the line out

        zline = (x >=line_info['range_min_xguess'][i]) & \
                (x <=line_info['range_max_xguess'][i])
        if np.sum(zline) < line_info['num_parms'][i]:
            line_goodbad[i] = 0
            continue

        # Get fit type

        if line_info['fit_type'][i] == 'L':
            type = 'lorentzian'

        elif line_info['fit_type'][i] == 'G':
            type = 'gaussian'

        elif line_info['fit_type'][i] == 'C':
            type = 'centroid'

        else:
            print('Unknown `fittype`.')

        if type != 'centroid':


            offset = np.min(y[zline]) if line_info['num_parms'][i] == 3 else 0



            try:

                fit = fit_peak1d(x[zline],y[zline],
                                 nparms=line_info['num_parms'][i],
                                 type=type,positive=True)

            except RuntimeError:

                fit = {'parms':[np.nan, np.nan, np.nan],
                       'fit':np.full(len(x[zline]), np.nan)}


            # Store the results
    
            line_xpos[i] = fit['parms'][1]
            line_fwhm[i] = fit['parms'][2]*2.354
            line_inten[i] = fit['parms'][0]
                    
        else:

            med = np.median((y[zline][0],y[zline][-1]))

            m1 = np.sum(x[zline]*(y[zline]-med)) / np.sum(y[zline]-med)
            
            m2 = np.sum((x[zline]-m1)**2*(y[zline]-med)) / np.sum(y[zline]-med)

            f = interpolate.interp1d(x[zline],y[zline])
            inten = f(m1)

            line_xpos[i] = m1            
            line_fwhm[i] = np.sqrt(m2)*2.354
            line_inten[i] = inten

        # Now let's check to see whether it is a good find or not.


        if (line_xpos[i] <= line_info['xguess'][i]+pix_thresh) and \
           (line_xpos[i] >= line_info['xguess'][i]-pix_thresh) and \
           (line_fwhm[i] > 0) and (line_inten[i] > 0):
           line_goodbad[i] = 1

        if qafileinfo:

            fig, (axes1, axes2) = pl.subplots(2, figsize=qafileinfo['figsize'],
                                                  constrained_layout=False)

            # Plot the entire spectrum with vertical line for the line
            
            yrange = get_spec_range(y,frac=0.1)

            title = 'Order '+str(line_info['order'][i])+': $\lambda$='+\
                line_info['wavelength'][i]+' $\mu$m'
            
            axes1.step(x, y, 'black')
            axes1.set_title(title)            
            axes1.set_ylim(ymin=0, ymax=yrange[1])
            axes1.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes1.axvline(x=line_info['xguess'][i], linestyle='solid',color='r')

            # Plot the line itself with fit.
            
            yrange = get_spec_range(y[zline],frac=0.1)

            goodbad = 'Good Fit' if line_goodbad[i] == 1 else 'Bad Fit'
            
            axes2.step(x[zline], y[zline], 'black')
            axes2.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes2.axvline(x=line_info['xguess'][i], linestyle='solid',
                          color='r')
            axes2.axvline(x=line_xpos[i], linestyle='solid', color='g')            
            if type != 'centroid':
                axes2.step(x[zline], fit['fit'], 'g')

            axes2.set_title(goodbad)
            
            axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])

            pdf.savefig(fig)
            pl.close(fig)

            if clupdate is not None:
                loop_progress(i, 0, nlines)

    if qafileinfo is not None:
        pdf.close()

    # Add the results

    line_info['x'] = line_xpos
    line_info['fwhm_pix'] = line_fwhm
    line_info['intensity'] = line_inten
    line_info['goodbad'] = line_goodbad

    return line_info
        
