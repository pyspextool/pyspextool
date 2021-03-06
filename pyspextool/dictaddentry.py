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

    


    
