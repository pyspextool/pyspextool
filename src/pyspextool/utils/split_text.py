import copy

def split_text(text, length=68):

    """
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

    """

    if isinstance(text,list):
        mtext = ''
        for t in text: mtext+=t
        text = copy.deepcopy(mtext)

    ntext = text

    if len(text) <= length:
        return text

    arr = []

    while len(ntext) > length:

        # Cannot cut at less than the length of a word
        
        i = 0        
        pos = -1
        while pos == -1:

            pos = ntext.rfind(' ', 0, length-1+i)
            i += 1
            
        left, ntext = ntext[:pos+1], ntext[pos+1:]        
        arr.append(left)

    arr.append(ntext)

    return arr
