"""
Pick a working SIMBAD mirror.

The main SIMBAD server (Strasbourg) is sometimes down while its mirrors are
up. These helpers check the mirrors in order and point astroquery at the
first one that responds. They follow the approach used in astrodb_utils
(astrodbtoolkit/astrodb_utils PRs #275 and #276).
"""

import logging
import time

import requests
from astroquery.simbad import Simbad
from astroquery.simbad import conf as simbad_conf

__all__ = [
    "SIMBAD_MIRRORS",
    "select_simbad_mirror",
    "set_simbad_mirror",
    "simbad_available",
]

logger = logging.getLogger(__name__)

SIMBAD_MIRRORS = ("simbad.cds.unistra.fr", "simbad.harvard.edu", "simbad.idia.ac.za")

# Result of the last mirror check: (a mirror host name or None, time.monotonic() of the check)
_simbad_status = None

# After a failed check, wait this many seconds before probing the mirrors again
SIMBAD_RECHECK_SECONDS = 60


def set_simbad_mirror(server: str):
    """
    Choose which SIMBAD mirror astroquery queries.

    The shared ``astroquery.simbad.Simbad`` object is used by SIMBAD lookups,
    so this affects every call made afterwards in the current Python session.
    New ``Simbad()`` instances do not pick this up automatically; copy the
    setting with ``sb.server = Simbad.server``.

    Parameters
    ----------
    server : str
        Host name of the mirror, without ``https://`` or a path, e.g.
        ``"simbad.cds.unistra.fr"`` (default, Strasbourg),
        ``"simbad.harvard.edu"`` (CfA) or ``"simbad.idia.ac.za"``
        (IDIA, South Africa).

    """
    server = server.replace("https://", "").replace("http://", "").strip("/")
    # astroquery only knows about some mirrors by default; register others so the setter accepts them
    if server not in simbad_conf.servers_list:
        simbad_conf.servers_list = [*simbad_conf.servers_list, server]
    Simbad.server = server
    logger.info(f"Using SIMBAD mirror: {server}")


def select_simbad_mirror(mirrors=None, timeout: float = 5):
    """
    Use the first SIMBAD mirror that responds, trying them in order.

    Parameters
    ----------
    mirrors : sequence of str, optional
        Host names to try, in order of preference. Default: the mirror
        currently in use (so a mirror chosen with `set_simbad_mirror` is kept
        while it works), then Strasbourg (CDS), Harvard (CfA) and IDIA
        (South Africa).
    timeout : float
        Seconds to wait for each mirror.

    Returns
    -------
    str or None
        The mirror now in use, or None if none responded (the current setting
        is left unchanged).

    """
    global _simbad_status
    if mirrors is None:
        # current mirror first, no duplicates
        mirrors = dict.fromkeys((Simbad.server, *SIMBAD_MIRRORS))
    logger.info("Checking which SIMBAD mirror is available...")
    for server in mirrors:
        try:
            ok = requests.get(
                f"https://{server}/simbad/sim-tap/capabilities", timeout=timeout
            ).ok
        except requests.exceptions.RequestException:
            ok = False
        if ok:
            set_simbad_mirror(server)
            _simbad_status = (server, time.monotonic())
            return server
        logger.warning(f"SIMBAD mirror {server} is not responding")
    logger.warning("No SIMBAD mirror responded")
    _simbad_status = (None, time.monotonic())
    return None


def simbad_available() -> bool:
    """
    Check (lazily) whether SIMBAD can be used, switching to a working mirror if needed.

    The first call probes the mirrors with `select_simbad_mirror`. After a
    success the answer is remembered for the rest of the session. After a
    failure it is remembered for `SIMBAD_RECHECK_SECONDS`, so a long session
    notices when SIMBAD comes back without probing on every call.

    Returns
    -------
    bool
        True if a SIMBAD mirror is reachable.

    """
    if _simbad_status is not None:
        mirror, checked = _simbad_status
        if mirror is not None or time.monotonic() - checked < SIMBAD_RECHECK_SECONDS:
            return mirror is not None
    return select_simbad_mirror() is not None
