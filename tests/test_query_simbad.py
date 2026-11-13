from pyspextool.io.query_simbad import query_simbad
from pyspextool.setup_utils import pyspextool_setup
from pyspextool.pyspextoolerror import pySpextoolError
import pytest



@pytest.mark.simbad
def test_query_simbad():

    result = query_simbad('HD 12345')

    assert result == {'name':'HD  12345', 'sptype':'G8III','vmag':8.770,
                      'bmag':9.520}
    

    with pytest.raises(pySpextoolError):

        result = query_simbad('HD 00000')


    result = query_simbad(['02:00:48.75','-12:52:29.9'])

    assert result == {'name':'HD  12345', 'sptype':'G8III','vmag':8.770,'bmag':9.520}


    result = query_simbad({'id':'HD  12345','sptype':'G8III','vmag':8.77, 'bmag':9.52})

    assert result == {'name':'HD  12345', 'sptype':'G8III','vmag':8.770,'bmag':9.520}



def _simbad_down():
    raise RuntimeError("No mirror responded")


def test_query_simbad_dict_does_not_need_simbad(monkeypatch):
    # A dictionary input never touches SIMBAD, so it works even when SIMBAD is down
    import pyspextool.io.query_simbad as qs
    monkeypatch.setattr(qs, "_simbad_mirror_found", False)
    monkeypatch.setattr(qs.Simbad, "find_mirror", _simbad_down)
    result = query_simbad({'id':'HD  12345','sptype':'G8III','vmag':8.77, 'bmag':9.52})
    assert result['name'] == 'HD  12345'


def test_query_simbad_raises_when_simbad_down(monkeypatch):
    import pyspextool.io.query_simbad as qs
    monkeypatch.setattr(qs, "_simbad_mirror_found", False)
    monkeypatch.setattr(qs.Simbad, "find_mirror", _simbad_down)
    with pytest.raises(pySpextoolError, match="SIMBAD is not reachable"):
        query_simbad('HD 12345')


def test_find_simbad_mirror_remembers_success(monkeypatch):
    import pyspextool.io.query_simbad as qs
    calls = []
    monkeypatch.setattr(qs, "_simbad_mirror_found", False)
    monkeypatch.setattr(qs.Simbad, "find_mirror", lambda: calls.append(1))
    assert qs.find_simbad_mirror() is True
    assert qs.find_simbad_mirror() is True
    assert len(calls) == 1  # checked only once


def test_find_simbad_mirror_retries_after_failure(monkeypatch):
    import pyspextool.io.query_simbad as qs
    calls = []

    def down():
        calls.append(1)
        raise RuntimeError("No mirror responded")

    monkeypatch.setattr(qs, "_simbad_mirror_found", False)
    monkeypatch.setattr(qs.Simbad, "find_mirror", down)
    assert qs.find_simbad_mirror() is False
    assert qs.find_simbad_mirror() is False
    assert len(calls) == 2  # a failure is not remembered
