def fsextract(string,method):

    """
    Extracts the indices or filenames from a comma-separated string

    Input Parameters
    ----------------
    string : str
        a comma separated string of either file names or file index numbers.

    method : {'index','filename'}
        'index' if the values passed are index values and 'filename' if the 
        values passed our file names

    Returns
    --------
    list
        `method`='index'
         A list of integers giving the individual file numbers

        `method`='filename'
         A list of integers giving the individual file names

    Procedure
    ---------
    index = True
    1.  separate into groups based on the comma.
    2.  loop over each group, and separate based on dash.
    3.  if no dash detected, append group (which is a number) to output list.
    3.  if dash detected, generate sequential numbers between limits and add to         output list.
    
    filename = True
    1.  separate into groups based on the comma.

    Examples
    --------
    > fsextract('1-3,5,7,10-12','index')
    [1, 2, 3, 5, 7, 10, 11, 12]

    > fsextract('spc00001.a.fits,spc00002.a.fits','filename')
    ['spc00001.a.fits', 'spc00002.a.fits']


    Modification History
    --------------------
    2022-05-24 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_fsextract.pro

    """

#  Check whether you are in index or filename mode
    
    if method == 'index':

#  Split on the comma
        
        groups = string.split(',')

#  Now loop and split on the dash, full in the missing numbers, and convert
#  to a string
        
        oarr = []
        for group in groups:
            lowupp = group.split('-')
            if len(lowupp) == 1:

# no dash, just add to output list               

                oarr.append(int(lowupp[0]))

            else:

# dash dectected, generate sequential numbers and add to output list
                
                arr = list(range(int(lowupp[0]),int(lowupp[1])+1))
                oarr+=arr

        return oarr

    elif method == 'filename':

# Just split on the comma and return
        
        return string.split(',')
        
    else:

        print('method unknown.')
        return 
    
