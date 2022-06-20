def round_school(x):

    """
    performs "school" kind of rounding where 0.5 goes to 1

    Input Parameters
    ----------------
    x : float
        a floating point value 

    Returns
    --------
    int
        The value rounded (see procedure)

    Procedure
    ---------
    https://stackoverflow.com/questions/43851273/how-to-round-float-0-5-up-to-1-0-while-still-rounding-0-45-to-0-0-as-the-usual

    Examples
    --------
    
    > round_school(1.5)
    2

    > round_school(1.49)
    1

    > round_school(0.5)
    1

    > round_school(0.49)
    0

    > round_school(-0.49)
    0

    > round_school(-0.50)
    -1

    > round_school(-1.49)
    -1

    > round_school(-1.50)
    -2

    Modification History
    --------------------
    2022-06-06 - Written by M. Cushing, University of Toledo.


    """

    i, f = divmod(x, 1)
    return int(i + ((f >= 0.5) if (x > 0) else (f > 0.5)))
    
