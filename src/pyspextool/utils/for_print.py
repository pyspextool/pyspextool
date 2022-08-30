def for_print(*args):

    """
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

    """

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
            print(args[j][i], end=end)
