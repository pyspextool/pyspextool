"""Tests for SIMBAD mirror selection (mocked, so no network is needed)."""

import pytest
import requests
from astroquery.simbad import conf as simbad_conf

import pyspextool.io.simbad_mirrors as mirrors
from pyspextool.io.simbad_mirrors import (
    select_simbad_mirror,
    set_simbad_mirror,
    simbad_available,
)


class FakeResponse:
    def __init__(self, ok):
        self.ok = ok


@pytest.fixture(autouse=True)
def restore_simbad_state():
    """Undo any mirror changes so other tests are unaffected."""
    old_server = mirrors.Simbad.server
    old_servers = list(simbad_conf.servers_list)
    old_status = mirrors._simbad_status
    yield
    mirrors._simbad_status = old_status
    simbad_conf.servers_list = old_servers
    mirrors.Simbad.server = old_server


def fake_get(up_hosts, not_ok_hosts=()):
    """Mock for requests.get: listed hosts respond, `not_ok_hosts` return an error status, others time out."""

    def _get(url, timeout=None):
        if any(host in url for host in up_hosts):
            return FakeResponse(True)
        if any(host in url for host in not_ok_hosts):
            return FakeResponse(False)
        raise requests.exceptions.ConnectTimeout

    return _get


def test_select_mirror_prefers_first(monkeypatch):
    set_simbad_mirror("simbad.cds.unistra.fr")  # start from the default mirror
    monkeypatch.setattr(
        mirrors.requests,
        "get",
        fake_get(["simbad.cds.unistra.fr", "simbad.harvard.edu"]),
    )
    assert select_simbad_mirror() == "simbad.cds.unistra.fr"
    assert mirrors.Simbad.server == "simbad.cds.unistra.fr"


def test_select_mirror_falls_back_after_error_status_and_timeout(monkeypatch):
    # Strasbourg answers with an error (like the 404 seen during an outage), Harvard times out
    set_simbad_mirror("simbad.cds.unistra.fr")
    monkeypatch.setattr(
        mirrors.requests,
        "get",
        fake_get(["simbad.idia.ac.za"], not_ok_hosts=["simbad.cds.unistra.fr"]),
    )
    assert select_simbad_mirror() == "simbad.idia.ac.za"
    assert mirrors.Simbad.server == "simbad.idia.ac.za"


def test_select_mirror_none_available(monkeypatch):
    monkeypatch.setattr(mirrors.requests, "get", fake_get([]))
    mirrors.Simbad.server = "simbad.harvard.edu"
    assert select_simbad_mirror() is None
    assert mirrors.Simbad.server == "simbad.harvard.edu"  # unchanged


def test_select_mirror_keeps_current_choice_while_it_works(monkeypatch):
    set_simbad_mirror("simbad.harvard.edu")
    monkeypatch.setattr(
        mirrors.requests,
        "get",
        fake_get(["simbad.cds.unistra.fr", "simbad.harvard.edu"]),
    )
    # not switched back to Strasbourg
    assert select_simbad_mirror() == "simbad.harvard.edu"


def test_set_mirror_cleans_host_and_registers_unknown_servers():
    set_simbad_mirror("https://simbad.idia.ac.za/")
    assert mirrors.Simbad.server == "simbad.idia.ac.za"
    assert "simbad.idia.ac.za" in simbad_conf.servers_list


def test_simbad_available_caches_and_rechecks(monkeypatch):
    calls = []

    def fake_select():  # behaves like select_simbad_mirror when every mirror is down
        calls.append(1)
        mirrors._simbad_status = (None, mirrors.time.monotonic())

    monkeypatch.setattr(mirrors, "select_simbad_mirror", fake_select)
    mirrors._simbad_status = None
    assert simbad_available() is False
    assert simbad_available() is False  # remembered, no second probe
    assert len(calls) == 1
    monkeypatch.setattr(mirrors, "SIMBAD_RECHECK_SECONDS", -1)  # failure memory expired
    simbad_available()
    assert len(calls) == 2


def test_simbad_available_remembers_success(monkeypatch):
    monkeypatch.setattr(mirrors.requests, "get", fake_get(["simbad.harvard.edu"]))
    mirrors._simbad_status = None
    assert simbad_available() is True
    # now everything is down, but the earlier success is remembered
    monkeypatch.setattr(mirrors.requests, "get", fake_get([]))
    assert simbad_available() is True
