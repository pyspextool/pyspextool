from pyspextool.io.query_simbad import query_simbad
from pyspextool.setup_utils import pyspextool_setup
from pyspextool.pyspextoolerror import pySpextoolError
import pytest



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

