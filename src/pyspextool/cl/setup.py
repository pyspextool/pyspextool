import os
from pyspextool.cl import config

def setup(instrument=config.state['instruments'][0], rawpath=None,
          calpath=None, procpath=None, qapath=None, clupdate=True):
    

    """
    Set up basic information for pyspextool to run.


    Parameters
    ----------------
    instrument : str, default = config.state['instruments'][0], optional
        The name of the instrument.

    rawpath : str, optional
        The path to the raw directory.

    calpath : str, optional
        The path to the calibration directory.

    procpath : str, optional
        The path to the processed directory.

    rawpath : str, optional
        The path to the quality assurance directory.

    clupdate : bool, default = True
        Set to report the setup results.

    Returns
    -------
    None


    Notes
    -----
    Hmmm

    Examples
    --------
    later

    """

    # Get the package path

    dir = os.path.dirname(config.__file__)+'../../../../'
    packagepath = os.path.normpath(dir)
    config.state['packagepath'] = packagepath
    
    # Get the instrument directory

    test = os.path.join(packagepath,'instruments',instrument)
    if os.path.isdir(test) is True:

        config.state['instrument'] = instrument

    else:

        message = 'Unknown instrument.  Possible instruments are: '+\
          str.join(", ",config.state['instruments'])+'.'
        raise ValueError(message)

    # Now do the paths.  First we check for the .path file in the user's home
    # directory

    homedir = os.path.expanduser('~')    
    filename = os.path.join(homedir, '.pyspextool_'+\
                            config.state['instrument']+'.dat')
    filename = os.path.join(homedir,filename)

    if os.path.isfile(filename) is True:

        f = open(filename, 'r')
        paths = []
        for line in f:

            paths.append(line.strip())

        config.state['rawpath'] = paths[0]
        config.state['calpath'] = paths[1]
        config.state['procpath'] = paths[2]
        config.state['qapath'] = paths[3]


    else:
        cwd = os.path.abspath(os.getcwd())
        config.state['rawpath'] = cwd
        config.state['calpath'] = cwd
        config.state['procpath'] = cwd
        config.state['qapath'] = cwd
                                
        
# Now let's modify the paths based on the user
        
    # Is rawpath passed?
    
    if rawpath is not None:

        # Is the path real?
        
        if os.path.isdir(rawpath) is True:

            # get the absolute path 
            
            config.state['rawpath'] = os.path.abspath(rawpath)
            
        else:

            # Path is not real.
            
            message = 'rawpath not a directory.'
            raise ValueError(message)

    # Is calpath passed?
    
    if calpath is not None:

        # Is the path real?
        
        if os.path.isdir(calpath) is True:

            # get the absolute path 
            
            config.state['calpath'] = os.path.abspath(calpath)
            
        else:

            # Path is not real.
            
            message = 'calpath not a directory.'
            raise ValueError(message)

    # Is procpath passed?
    
    if procpath is not None:

        # Is the path real?
        
        if os.path.isdir(procpath) is True:

            # get the absolute path 
            
            config.state['procpath'] = os.path.abspath(procpath)
            
        else:

            # Path is not real.
            
            message = 'procpath not a directory.'
            raise ValueError(message)

    # Is qapath passed?
    
    if qapath is not None:

        # Is the path real?
        
        if os.path.isdir(qapath) is True:

            # get the absolute path 
            
            config.state['qapath'] = os.path.abspath(qapath)
            
        else:

            # Path is not real.
            
            message = 'qapath not a directory.'
            raise ValueError(message)
        
# Now write the paths to the user home directory

    f = open(os.path.join(homedir,'.pyspextool_'+config.state['instrument']+\
                          '.dat'),'w')
    f.write('%s \n' % config.state['rawpath'])   
    f.write('%s \n' % config.state['calpath'])
    f.write('%s \n' % config.state['procpath'])
    f.write('%s \n' % config.state['qapath'])
    f.close()

    if clupdate is True:

        print('Pyspextool Setup')
        print('----------------')
        print('Instrument: ',config.state['instrument'])
        print()
        print('rawpath: ',config.state['rawpath'])
        print('calpath: ',config.state['calpath'])
        print('procpath: ',config.state['procpath'])
        print('qapath: ',config.state['qapath'])
                                


            
        

        

    

