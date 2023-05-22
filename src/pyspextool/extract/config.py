state = {'apertures':1,
         'apertures_done':False,
         'apradii':1,
         'apsigns':1,
         'atrans_wave':1,
         'atrans_trans':1,
         'bad_pixel_mask': None,
         'bgfitdeg':1,
         'bgradius':1,
         'bgregions':'',
         'bgwidth':1,         
         'continue':0,
         'edgecoeffs':1, 
         'edgedeg':1,
         'extract_done':False,
         'filereadmode':'',
         'flat':1,
         'hdrinfo':1,
         'keywords':[],
         'lincormax':0,
         'linearity_info':0,
         'load_done':False,
         'modename':'',
         'ncols':1,
         'norders':1,
         'nrows':1,
         'orders':1,
         'orders_done':False,
         'ordermask':1,
         'output_files':'',
         'output_prefix':'spectra',
         'parameters_done':False,
         'plate_scale':1,
         'prefix':'',
         'profiles_done':False,                             
         'psdoorders':1,
         'psfradius':1,
         'qafilename':'',
         'raw_bad_pixel_mask': None,
         'rectindices':1,
         'rectorders':1,
         'reductionmode':1,
         'resolvingpower':1,
         'rotation':1,
         'slith_arc':1,
         'slitw_arc':1,
         'slith_pix':1,
         'slitw_pix':1,
         'spatcal':1,
         'trace_done':False,
         'tracecoeffs':1,
         'type':{'type':''},
         'type_done':False,                    
         'varimage':1,           
         'wavecal':1,
         'wavecaltype':1,
         'workimage':1,
         'xranges':1,
         'xsdoorders':1}

load = {'flatfile':'', 'wavecalfile':'', 'reduction_mode':'', 'directory':'',
        'suffix':'', 'doflat':True, 'dolinearity':True, 'qafile':True,
        'qaplot':True, 'qaplotsize':(6,6), 'verbose':True}

type = {'type':'', 'verbose':True}

profiles = {'qaplot':True, 'qafile':True, 'qaplotsize':(6,10), 'verbose':True}


apertures = {'apertures':1, 'method':'', 'fwhm':0.8, 'qaplot':False,
             'qafile':False, 'qaplotsize':(6,10), 'verbose':True}

orders = {'include':'', 'exclude':'', 'include_all':False,
          'verbose':True, 'qaplot':False, 'qaplotsize':(6,10),
          'qafile':False}

    
trace = {'fitdeg':1, 'step':1, 'sumwidth':1, 'centhresh':1, 'fwhm':1,
         'verbose':True, 'qaplot':False, 'qafile':False, 'qaplotsize':(6,10)}
        
parameters = {'apradii':1, 'psfradius':1, 'bgradius':1, 'bgwidth':1,
              'bgregions':1, 'bgdeg':1, 'qaplot':False, 'qafile':False,
              'qaplotsize':(6,10)}

extract = {'verbose':True, 'qaplot':False, 'qafile':False, 'qaplotsize':(10,6)}

combine = {'beam_mode':1, 'scale_orders':False, 'background_subtraction':False, 
           'flat_field':None,'overwrite':True, 'verbose':True,
           'qaplot':False, 'qafile':False, 'qaplotsize':(10,6)}
