import numpy as np
from pyspextool.cl.setup import config
from pyspextool.utils.arrays import find_index
from pyspextool.utils.arrays import trim_nan
from pyspextool.fit.fit_peak1d import fit_peak1d
from pyspextool.utils.for_print import for_print
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
from pyspextool.plot.limits import get_spec_range
from pyspextool.utils.loop_progress import loop_progress

def find_lines_1dxd(spectra, orders, line_info, pix_thresh, qafileinfo=None):

# Get basic information
    
    norders = len(orders)
    nlines = len(line_info['strwave'])

# Setup the output arrays

    line_x = np.full(nlines, np.nan)
    line_fwhm = np.full(nlines, np.nan)
    line_inten = np.full(nlines, np.nan)
    line_goodbad = np.zeros(nlines, dtype=np.uint8)

# Setup the qaplot if asked

    if qafileinfo is not None:

        pl.rcParams['font.size'] = '12'
        pl.rcParams['lines.linewidth'] = 0.75
        pdf = PdfPages('alltogther.pdf')

    
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

        zline = (x >=line_info['xleft'][i]) & (x <=line_info['xright'][i])
        if np.sum(zline) < line_info['nterms'][i]:
            line_goodbad[i] = 0
            continue

        # Get fit type

        if line_info['fittype'][i] == 'L':
            type = 'lorentzian'

        elif line_info['fittype'][i] == 'G':
            type = 'gaussian'

        elif line_info['fittype'][i] == 'C':
            type = 'centroid'

        else:
            print('Unknown `fittype`.')

        if type != 'centroid':


            offset = np.min(y[zline]) if line_info['nterms'][i] == 3 else 0

            fit = fit_peak1d(x[zline],y[zline],nparms=line_info['nterms'][i],
                             type=type,positive=True)
        
        else:

            print('do later.')

# Store the results

        
        line_x[i] = fit['parms'][1]
        line_fwhm[i] = fit['parms'][2]*2.354
        line_inten[i] = fit['parms'][0]

# Now let's check to see whether it is a good find or not.


        if (fit['parms'][1] <= line_info['xguess'][i]+pix_thresh) and \
           (fit['parms'][1] >= line_info['xguess'][i]-pix_thresh) and \
           (fit['parms'][2] > 0) and (fit['parms'][0] > 0):
           line_goodbad[i] = 1
    
        if qafileinfo:

                        
            fig, (axes1, axes2) = pl.subplots(2, figsize=qafileinfo['figsize'],
                                                  constrained_layout=False)

            yrange = get_spec_range(y,frac=0.1)


            title = 'Order '+str(line_info['order'][i])+': $\lambda$='+\
                line_info['strwave'][i]+' $\mu$m'
            
            axes1.step(x, y, 'black')
            axes1.set_title(title)            
            axes1.set_ylim(ymin=0, ymax=yrange[1])
            axes1.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes1.axvline(x=line_info['xguess'][i], linestyle='solid',color='r')

            yrange = get_spec_range(y[zline],frac=0.1)

            goodbad = 'Good Fit' if line_goodbad[i] == 1 else 'Bad Fit'
            
            axes2.step(x[zline], y[zline], 'black')
            axes2.set(xlabel='Column Number (pixels)',ylabel='Intensity')
            axes2.axvline(x=line_info['xguess'][i], linestyle='solid',
                          color='r')
            axes2.axvline(x=fit['parms'][1], linestyle='solid', color='g')            
            axes2.step(x[zline], fit['fit'], 'g')
            axes2.set_title(goodbad)
            
            axes2.set_ylim(ymin=yrange[0], ymax=yrange[1])

            

            pdf.savefig(fig)
            pl.close(fig)

        loop_progress(i,0,nlines)
            
    if qafileinfo is not None:
        pdf.close()

        
