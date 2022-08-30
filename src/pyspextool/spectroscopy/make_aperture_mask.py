import numpy as np
from pyspextool.utils.arrays import find_index
from pyspextool.utils.math import round


def make_aperture_mask(slit_arc, appos, apradius, psbginfo=None, xsbginfo=None):

    """
    To create an aperture mask for use in extraction.

    Input Parameters
    ----------------
    slit_arc : array_like
        An array of angular positions along the slit (typically in arcseconds).

    appos : array_like or float or int
        The aperture positions (same units as `slit_arc`).

    apradius : array_like or float or int
        The aperture radius (same units as `slit_arc`).

    psbginfo : list, optional
        (2,) list giving the angular distance from the aperture to begin the
        background region and the angular width of the background region.

    xsbginfo : list of list, optional
        (n_bg,) list where each element is a 2-element list of angular
        positions that indicate a background region.

    Returns
    -------
    numpy.ndarray
         An array of the aperture mask.  The pixel values denote their status:

         val = 0 - nothing
         val = -1 - a background pixel
         0 < val <= 1 - a pixel belonging to aperture 1
         1 < val <= 2 - a pixel belonging to aperture 2
         etc.

    Notes
    -----
    The function generates an array where each pixel value denotes its status.
    To construct the mask, the background regions are identified without
    regard to potential overlap.  If `psbginfo` is passed, then overlap is
    dealt with by clearing out from `appos`-`bgstart` to `appos`+`bgstart`.
    Finally, the apertures are labeled by their aperture number, e.g. 1, 2, etc.
    Pixels at the edge of an aperture are adjust to reflect their fractional
    pixel value, e.g. 1 -> 0.5 or 2 -> 1.6.

    Examples
    --------
    > slit_arc = np.arange(100)
    > mkapmask(slit_arc,[50,60],[4,4.1],psbginfo=[11,20])
      [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.
        0.   0.   0.   0.   0.5  1.   1.   1.   1.   1.   1.   1.   0.5  0.
        1.6  2.   2.   2.   2.   2.   2.   2.   1.6  0.   0.   0.   0.   0.
        0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.   0.   0.   0.   0.
        0.   0. ]

    > slit_arc = np.arange(100)
    > mkapmask(slit_arc,[50,60],[4,4.1],xsbginfo=[[7,15],[80,90]])
      [ 0.   0.   0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.  -1.  -1.  -1.
       -1.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.5  1.   1.   1.   1.   1.   1.   1.   0.5  0.
        1.6  2.   2.   2.   2.   2.   2.   2.   1.6  0.   0.   0.   0.   0.
        0.   0.   0.   0.   0.   0.   0.   0.   0.   0.  -1.  -1.  -1.  -1.
       -1.  -1.  -1.  -1.  -1.  -1.   0.   0.   0.   0.   0.   0.   0.   0.
        0.   0. ]

    Modification History
    --------------------
    2022-07-01 - Written by M. Cushing, University of Toledo.
    Based on Spextool IDL program mc_mkapmask.pro.

    """

    # Convert to numpy arrays

    appos = np.array(appos, dtype='float', ndmin=1)
    apradius = np.array(apradius, dtype='float', ndmin=1)

    # Now do basic things

    npix = len(slit_arc)
    naps = len(appos)

    # Now check to see if the apertures overlap

    mask = np.zeros([npix])
    for i in range(naps):
        ap = [appos[i] - apradius[i], appos[i] + apradius[i]]
        idx = find_index(slit_arc, ap)
        mask[int(idx[0]):int(idx[1])] += 1

    test = mask > 1
    if sum(test) >= 1:
        exception = 'Apertures overlap.'
        raise ValueError(exception)

    # Create the mask

    mask = np.zeros([npix])

    # Let's deal with the background regions, point source?

    if psbginfo is not None:

        bgstart = psbginfo[0]
        bgwidth = psbginfo[1]

        # Set all pixels that fall in the background region around each aperture
        # to -1 regardless of whether there is overlap.

        for i in range(naps):

            if bgstart <= apradius[i]:
                error = 'psbginfo[0] (bgstart) must be greater than ' + \
                        'the aperture radius.'
                raise ValueError(error)

            leftbg = [appos[i] - bgstart - bgwidth, appos[i] - bgstart]
            rghtbg = [appos[i] + bgstart, appos[i] + bgstart + bgwidth]

            leftidx = find_index(slit_arc, leftbg)
            rghtidx = find_index(slit_arc, rghtbg)

            mask[int(leftidx[0]):int(leftidx[1])] = -1
            mask[int(rghtidx[0]):int(rghtidx[1])] = -1

        # Now clear out from -bgstart to +bgstart around each aperture position to
        # deal with any overlap

        for i in range(naps):
            clr = [appos[i] - bgstart, appos[i] + bgstart]
            clridx = find_index(slit_arc, clr)
            mask[int(clridx[0]):int(clridx[1])] = 0

    # extended source?

    if xsbginfo is not None:

        # Let's check against a single background region

        if any(isinstance(el, list) for el in xsbginfo) is False:
            xsbginfo = [xsbginfo]

        for i in range(len(xsbginfo)):
            idx = find_index(slit_arc, xsbginfo[i])
            mask[int(idx[0]):int(idx[1])] = -1

    # Now fill in the apertures

    for i in range(naps):

        ap = [appos[i] - apradius[i], appos[i] + apradius[i]]
        idx = find_index(slit_arc, ap)
        mask[int(idx[0]):int(idx[1])] += i + 1

        # Fix end points to reflect fractional pixels

        if idx[0] - np.floor(idx[0]) >= 0.5:

            mask[int(idx[0])] = 0
            mask[int(idx[0]) + 1] = 0.5 + round(idx[0]) - idx[0] + i

        else:

            mask[int(idx[0])] = 0.5 - (idx[0] - np.floor(idx[0])) + i

        if idx[1] - np.floor(idx[1]) >= 0.5:

            mask[int(idx[0]) + 1] = 0.5 - (round(idx[1]) - idx[1]) + i

        else:

            mask[int(idx[1])] = 0.5 + (idx[1] - round(idx[1])) + i

    return mask
