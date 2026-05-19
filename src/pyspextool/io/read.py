import numpy as np
import numpy.typing as npt
import re
import os
import ast

from pyspextool import config as setup
from pyspextool.io.check import check_file
from pyspextool.pyspextoolerror import pySpextoolError

def read_hlines_file(*args:str):

    """
    To load H I wavelengths and line IDs into numpy arrays.

    Parameters
    ----------
    args : str, optional
    If given, the fullname path to the pySpextool HI.dat file.  
    If not given, the default location is used.

    Returns
    -------
    ndarray
        An (nlines,) array of H I wavelengths in microns in vacuo.

    ndarray
        An (nlines,) array of H I latex line IDs.

    ndarray
        An (nlines,) array of H I ASCII line IDs.
     
    """
    
    if len(args) == 0:
        
        fullpath = os.path.join(setup.state["package_path"], "data", "HI.dat")  

    else:

        fullpath = check_file(args[0])
        

    wavelength, latex_lineid, ascii_lineid = np.loadtxt(
        fullpath, 
        comments="#", 
        unpack=True, 
        dtype="str", 
        delimiter="|")


    wavelength = np.array(wavelength).astype(float)
    
    latex_lineid = np.char.strip(latex_lineid)
    ascii_lineid = np.char.strip(ascii_lineid)
    
    return wavelength, latex_lineid, ascii_lineid


def read_telluric_file(
    fullpath:str):

    """
    To read a pySpextool MODE_telluric.dat file.

    Parameters
    ----------
    fullpath : str
        The fullname path to the pySpextool MODE_telluric.dat file

    Returns
    --------
    dict
        'VEGAMODEL' : str
            The Vega model subscript.

        'METHOD' : str
            The convolution method, deconvolution or IP.

        'DECONVOLUTION_LINE_ORDER' : int, optional
            If METHOD=deconvolution, the order with the line to deconvolve.

        'DECONVOLUTION_LINE' : float, optional
            If METHOD=deconvolution, the wavelength of the line to deconvolve.

        'DECONVOLUTION_CONTINUUM_RANGE' : tuple, optional
            If METHOD=deconvolution, the wavelength range around the line.

        'DECONVOLUTION_CONTINUUM_FIT_DEGREE' : int, optional
            If METHOD=deconvolution, the polynomial degree of the continuum fit.

        'DECONVOLUTION_LINE_FIT_TYPE' : str {'gaussian', 'lorentzian'}
            If METHOD=deconvolution, the functional form of the line fit.

        'DECONVOLUTION_RVFIT_NFWHM' : int, optional
            If METHOD=deconvolution, the number of FWHM around the line to 
            determine the radial velocity.

        'DECONVOLUTIONN_NFWHM' : int, optional 
            If METHOD=deconvolution, the number of FWHM around the line to 
            deconvolve the line.

        'PIXEL_SHIFT_RANGE=

        'EW_SCALE_INFO': list, optional
            A list of dictionaries with X keys: 

        'TELLURIC_NOISE_INFO : list, optional
            A list of dictionaries with two keys: 'Order Number' : int, 
            and 'Wavelength Range' (2,)


    """

    # Read the file into to string arrays

    labels, vals = np.loadtxt(
        fullpath,
        comments='#',
        delimiter='=',
        unpack=True,
        dtype='str')

    # Strip any edge white spaces from poor formatting of the file

    labels = np.char.strip(labels)
    vals = np.char.strip(vals)

    # Start the search for each keyword.

    output = {}

    keyword = 'TELLURIC_VEGAMODEL'
    z = _find_keyword(labels, keyword)
    output[keyword.lower()] = str(vals[z].item())

    keyword = 'TELLURIC_METHOD'
    z = _find_keyword(labels, keyword)
    output[keyword.lower()] = str(vals[z].item())

    keyword = 'DECONVOLUTION_LINE_ORDER'
    z = _find_keyword(labels, keyword, required=False)
    output[keyword.lower()] = int(vals[z].item()) if z is not None else None

    keyword = 'DECONVOLUTION_LINE'
    z = _find_keyword(labels, keyword, required=False)
    output[keyword.lower()] = float(vals[z].item()) if z is not None else None
    
    keyword = 'DECONVOLUTION_CONTINUUM_FIT_RANGE'
    z = _find_keyword(labels, keyword, required=False)

    output[keyword.lower()] = ast.literal_eval(vals[z].item()) \
        if z is not None else None
    
    keyword = 'DECONVOLUTION_CONTINUUM_FIT_DEGREE'
    z = _find_keyword(labels, keyword, required=False)
    output[keyword.lower()] = int(vals[z].item()) if z is not None else None
    
    keyword = 'DECONVOLUTION_LINE_FIT_TYPE'
    z = _find_keyword(labels, keyword, required=False)
    output[keyword.lower()] = str(vals[z].item()) if z is not None else None
    
    keyword = 'DECONVOLUTION_RVFIT_NFWHM'
    z = _find_keyword(labels, keyword, required=False)

    if z is not None:

        lis = [int(vals[z].item()),int(vals[z].item())]
        output[keyword.lower()] = lis 

    else:

        output[keyword.lower()] = [None,None]
    
    keyword = 'DECONVOLUTION_NFWHM'
    z = _find_keyword(labels, keyword, required=False)
    
    if z is not None:

        lis = [int(vals[z].item()),int(vals[z].item())]
        output[keyword.lower()] = lis 

    else:

        output[keyword.lower()] = [None, None]

    keyword = 'PIXEL_SHIFT_INFO'
    z = _find_keyword(labels, keyword, required=False)
    
    if z is not None:

        dict = {}

        tmp = np.char.strip(vals[z[0]][0].split(':'))
        
        dict['Order Number'] = int(tmp[0])
        dict['Wavelength Range'] = ast.literal_eval(tmp[1])
        
    else:
        
        dict = None

    output[keyword.lower()] = dict

    keyword = 'EW_SCALE_INFO'
    z = _find_keyword(labels, keyword, required=False)

    # Were any changes requested?
    
    if z is None:

        # No.

        fit_info = None

    else:

        # Yes.  Parse the results into a list of dictionaries equivalent 
        # to what the user can pass (eventually)

        nfits = np.size(z)

        # First, read the hydrogen line list

        waves, junk, ascii_hlineids = read_hlines_file()

        # And now start the parsing.

        all = []
        for i in range(nfits):

            dict = {}
            tmp = np.char.strip(vals[z[0]][i].split(':'))            

            dict['Order Number'] = int(tmp[0])
            dict['Fit Range'] = ast.literal_eval(tmp[1])
            dict['Line IDs'] = ast.literal_eval(tmp[2])
            idx = np.in1d(ascii_hlineids, np.array(dict['Line IDs']))
            dict['Lines'] = list(waves[idx])
            dict['Poly Degree'] = int(tmp[3])
            dict['Fit Tolerance'] = float(tmp[4])

            all.append(dict)

        # Now do the conversion to the code inputs

        orders = np.array([a['Order Number'] for a in all])
        unique_orders = np.unique(orders)
        
        by_order = []

        for order in list(unique_orders):

            z = np.where(order == orders)[0]
            tmp = []

            for idx in list(z):
                
                all[idx].pop("Order Number")
                tmp.append(all[idx])

            by_order.append(tmp)
            
        fit_info = {'Order Numbers':list(unique_orders),
                    'Fits':by_order}

    output[keyword.lower()] = fit_info

    keyword = 'TELLURIC_NOISE_INFO'
    z = _find_keyword(labels, 
        keyword, 
        required=False)
    
    if z is not None:


        result = []
        for i in range(np.size(z)):

            dict = {}
            tmp = np.char.strip(vals[z[0]][i].split(':'))

            dict['Order Number'] = int(tmp[0])
            dict['Wavelength Range'] = ast.literal_eval(tmp[1])

            result.append(dict)
    else:

        result = None
    
    output[keyword.lower()] = result

    return output


def _find_keyword(
    labels:npt.ArrayLike,
    keyword:str,
    required:bool=True):

    """



    """

    z = np.where(labels == keyword)
    if np.size(z):

        return z

    else:

        if required is True:

            message = 'Keyword '+keyword+' required and not found.'
            raise pySpextoolError(message)
        
        else:

            return None
