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
