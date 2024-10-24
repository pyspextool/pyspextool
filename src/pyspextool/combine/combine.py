import logging

from pyspextool.io.check import check_parameter, check_qakeywords
from pyspextool.combine.load_spectra import load_spectra
from pyspextool.combine.scale_spectra import scale_spectra as scale_them
from pyspextool.combine.combine_spectra import combine_spectra

logger = logging.getLogger("pyspextool")


def combine(files:str | list,
            output_name:str,
            scale_spectra:bool=True,
            scale_order:int=None,
            scale_wavelengthrange:list=None,
            scale_wavelengthfraction:float=0.7,
            correct_spectral_shape:bool=False,
            statistic:str='robust weighted mean',
            robust_sigma:int | float=8,
            verbose:bool=None,
            qa_show:bool=None,
            qa_showscale:float | int=None,
            qa_showblock:bool=None,
            qa_write:bool=None):    

    """
    To combine raw extracted spectra

    Parameters
    ----------
    files : str or list
        If type is str, then a comma-separated string of full file names, 
        e.g. 'spc00001.a.fits, spc00002.b.fits'.

        If type is list, then a two-element list where
        files[0] is a string giving the perfix.
        files[1] is a string giving the index numbers of the files.

        e.g. ['spc', '1-2']

    output_name : str
        The output file name sans the suffix.

    scale_spectra : {True, False}
        Set to scale the spectra to a common flux level.  See below.

    scale_order : int or None
        The order number of the spectra to use to compute the scale factors.
        If None, defaults to middle-most order.

    scale_wavelenthrange : list or None
        A (2,) list giving the wavelength range, e.g. [1.1,1.15].  If None, 
        defaults to the central `scale_range_fraction` of `scale_order`.

    scale_wavelengthfraction : float, default 0.7
        A float giving the central fraction of the wavelength range in order
        `scale_order` to be used to determine the scale factors.  

    statistic : {'robust weighted mean', 'robust mean', 'weighted mean', 
                 'mean', 'median'}
        For any of the means, the uncertainty is the standard error on the 
        mean, e.g. s/np.sqrt(n), where s is the sample standard deviation and
        n is the number of data points.

        For the median, the uncertainty is given by 1.482*MAD/np.sqrt(n) where
        MAD is the median absolute deviation given by,

        1.482*median(|x_i-x_med|).

    robust_sigma : float or int, default=8
        The sigma threshold of a robust statistic.   Values are identified as
        outliers if

        |x_i - x_med|/MAD > robust_sigma,

        where x_i is the ith data point, x_med is the median of the data, and 
        MAD is the median absolute deviattion given by,

        1.482*medianm(|x_i-x_med|).
    
    verbose : {None, True, False}
        Set to True to report updates to the command line.
        Set to False to not report updates to the command line.
        Set to None to default to setup.state['verbose'].
    
    qa_show : {None, True, False}
        Set to True to show a QA plot on the screen.
        Set to False to not show a QA plot on the screen.
        Set to None to default to setup.state['qa_show'].

    qa_write : {None, True, False}
        Set to True to write a QA plot to disk
        Set to False to not write a QA plot to disk.
        Set to None to default to setup.state['qa_write'].
    
    qa_showblock : {None, True, False}
        Set to True to block the screen QA plot.
        Set to False to not block the screen QA plot.
        Set to None to default to setup.state['qa_block'].
    
    qa_showscale : float or int, default=None
        The scale factor by which to increase or decrease the default size of
        the plot window.  Set to None to default to setup.state['qa_scale'].    
    
    Returns
    -------
    None
        Writes a pyspextool FITS file to disk.

    """

    #
    # Check the parameters and QA keywords
    #

    check_parameter('combine', 'files', files, ['str', 'list'])

    check_parameter('combine', 'output_name', output_name, 'str')    


    check_parameter('combine', 'scale_spectra', scale_spectra, 'bool')

    check_parameter('combine', 'scale_order', scale_order, ['NoneType', 'int'])

    check_parameter('combine', 'scale_wavelengthrange',
                    scale_wavelengthrange, ['NoneType', 'list'])

    check_parameter('combine', 'scale_wavelengthfraction',
                    scale_wavelengthfraction, 'float')
    
    check_parameter('combine', 'correct_spectral_shape',
                    correct_spectral_shape, 'bool')

    check_parameter('combine', 'statistic', statistic, 'str')

    check_parameter('combine', 'robust_sigma', robust_sigma, ['int', 'float'])

    check_parameter('combine', 'verbose', verbose, ['NoneType', 'bool'])
    
    check_parameter('combine', 'qa_write', qa_write, ['NoneType', 'bool'])

    check_parameter('combine', 'qa_show', qa_show, ['NoneType', 'bool'])

    check_parameter('combine', 'qa_showscale', qa_showscale,
                    ['int', 'float', 'NoneType'])

    check_parameter('combine', 'qa_showblock', qa_showblock,
                    ['NoneType', 'bool'])
    
    
    qa = check_qakeywords(verbose=verbose,
                          show=qa_show,
                          showscale=qa_showscale,
                          showblock=qa_showblock,
                          write=qa_write)

    #
    # Load the spectra
    #

    load_spectra(files,
                 output_name,
                 verbose=qa['verbose'],
                 qa_show=qa['show'],
                 qa_showscale=qa['showscale'],
                 qa_showblock=qa['showblock'],
                 qa_write=qa['write'])

    #
    # Scale the spectra
    #

    if scale_spectra is True:
    
        scale_them(order=scale_order,
                   wavelength_range=scale_wavelengthrange,
                   wavelength_fraction=scale_wavelengthfraction,
                   verbose=qa['verbose'],
                   qa_show=qa['show'],
                   qa_showscale=qa['showscale'],
                   qa_showblock=qa['showblock'],
                   qa_write=qa['write'])
    
    #
    # Combine and write the results to disk
    #

    combine_spectra(statistic=statistic,
                    robust_sigma=robust_sigma,
                    verbose=qa['verbose'],
                    qa_show=qa['show'],
                    qa_showscale=qa['showscale'],
                    qa_showblock=qa['showblock'],
                    qa_write=qa['write'])
    
        
    
