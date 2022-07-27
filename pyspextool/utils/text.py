"""Functions for text utlitiy management."""


def dictaddentry(idict,key,direction,newkey,nvalue):

    """
    Converts a sexigesmal number to a decimal


    Input Parameters
    ----------------
    idict : dict 
        a dictionary

    key : str
        ocation in `dict` to insert the new values

    direction : {'before','after'}
        insertion direction relative to `key`

    newkey : str
        name of the new key

    newvalue : any
        the value associated with `newkey`

    Returns
    --------
    dict
        the original dictionary with a new entry

    Procedure
    ---------
    #       https://stackoverflow.com/questions/44390818/how-to-insert-key-value-pair-into-dictionary-at-a-specified-position


    Examples
    --------
    > dict = {'HA':1,'PA':2,'MJD':3}
    > dictaddentry(dict,'MJD','before','new',4)

    {'HA': 1, 'PA': 2, 'new': 4, 'MJD': 3}

    > dict = {'HA':1,'PA':2,'MJD':3}
    > dictaddentry(dict,'HA','after','new',(3,4))

    {'HA': 1, 'new': (3, 4), 'PA': 2, 'MJD': 3}


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """
    
    pos   = list(idict.keys()).index(key)
    items = list(idict.items())       

    if direction == 'after':

        items.insert(pos+1,(newkey,nvalue))
        
    elif direction == 'before':   

        items.insert(pos,(newkey,nvalue))        
        
    else:

        print('Unknown direction.')
        return(-1)

    odict = dict(items)
    
    return(odict)

    
def forprint(*args):

    '''
    Prints lists of data as columns to the command line.


    Input Parameters
    ----------------
    *args : lists or 1D numpy.ndarrays
        "Arrays" of values to be printed as columns together.


    Returns
    -------
    None

    Notes
    -----
    Inputs must have the same length.

    Examples
    --------
    > x = [1,2,3]
    > y = [4,5,6]
    > forprint(x,y)
      1 4
      2 5
      3 6

    Modification History
    --------------------
    2022-06-28 - Written by M. Cushing, University of Toledo.
    Inspired by (but less than) IDL astrolib routine forprint.pro

    '''

    # Check to make sure all inputs have the same length

    nargs = len(args)

    lengths = []
    for i in range(nargs):

        lengths.append(len(args[i]))

    if lengths.count(lengths[0]) != nargs:

        raise Exception('Inputs are not the same length.')

    # Now print things
    
    for i in range(len(args[0])):

        for j in range(nargs):

            end = ' ' if j != nargs-1 else '\n'
            print(args[j][i],end=end)


def splittext(text,length=68):

    '''
    To split text into lines of a given (i.e., not longer than) length.


    Input Parameters
    ----------------
    text : str
        The text to be split.


    length : int, optional, default=68
        The requested length of each line.


    Returns
    -------
    list of str
         A list of str where element's length is <=`length`.


    Notes
    -----
    Just some parsing.  


    Examples
    --------
    > str = 'This flat was created by combining twenty images.'
    > splittext(str,length=10)
      ['This flat was ', 'created by ', 'combining twenty ', 'images.']


    Modification History
    --------------------
    2022-06-27 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_splittext.pro.

    '''
    
    leftoverlength = len(text)

    ntext = text

    i = 0

    if len(text) <= length: return(text)

    arr = []

    while len(ntext) > length:

        # Cannot cut at least than the length of a word
        
        i = 0        
        pos = -1
        while pos == -1:

            pos = ntext.rfind(' ',0,length-1+i)            
            i +=1
            
        left, ntext = ntext[:pos+1], ntext[pos+1:]        
        arr.append(left)
        
        leftoverlength = len(ntext)

    arr.append(ntext)

    return(arr)


def wherelist(lst,find):

    """
    To identify the elements of a list that match the user request.


    Input Parameters
    ----------------
    lst : list
        A list to be searched.

    find : any
        A value to search for in `lst`.

    Returns
    --------
    list
        A list of subscripts in `lst` that match `find`.  

    Procedure
    ---------
    https://www.journaldev.com/23759/python-find-string-in-list

    Example
    --------
    x = [1,'two',3]
    wherelist(x,'two')
    [1]

    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.

    """
    
    idx = []

    i = 0
    length = len(lst)

    while i < length:

        if find == lst[i]:
            idx.append(i)

        i += 1

    return(idx)