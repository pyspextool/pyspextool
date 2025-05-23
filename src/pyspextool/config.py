state = {"instruments": ['uspex', 'spex'],
         "qa_extensions": ['.png', '.pdf'],
         "search_extensions": ['.fits*', '.fits','.fits.gz'],         
         "units":['W m-2 um-1', 'erg s-1 cm-2 A-1', 'W m-2 Hz-1',
                  'ergs s-1 cm-2 Hz-1', 'Jy', 'mJy', 'uJy'],
         "version": None,
         "telluric_correctiontypes":['A0 V', 'basic', 'reflectance'],
         "vega_zmag":-0.03,
         "vega_zfd":3.46e-9*10, # ergs s-1 cm-2 A-1
         "vega_zlambda":5556e-4, # microns
         "xunits":'um',
         "lxunits":'$\mu$m'}


plotwindows = {'telluric_deconvolution':None,
               'telluric_normalize':None,
               'telluric_rv':None}


plots = {"portrait_size":(7,9),
         "square_size":(7,7),
         "landscape_size":(9,7),
         "profile_size":(6,3),
         "profilestack_max":4,
         "font_size":12,
         "spectrum_linewidth":0.5,
         "zoomspectrum_linewidth":1.5,         
         "spine_linewidth":1.5,
         "subplot_size":(6,3),
         "stack_max":4,
         'flat':1,
         'locate_orders':1,
         'wavecal_image':1,
         'combine_image':1,
         'pixel_shift':2,
         'wavecal_residuals':3,
         'profiles':4,
         'abeam_spectra':5,
         'bbeam_spectra':6,
         'abeam_snr':7,
         'bbeam_shr':8,
         'combine_spectra':9,
         'normalize_order':10,
         'radial_velocity':11,
         'deconvolution':12,
         'shifts':13,
         'ewscales':14}

