##
# contains tests for dictaddentry
##

from pyspextool.dictaddentry import dictaddentry

def test_dictaddentry():
    dicta = {'HA':1,'PA':2,'MJD': 3}
    res=dictaddentry(dicta,'MJD','before','new',4)
    assert res['HA']== 1
    assert res['PA']== 2
    assert res['new']== 4
    assert res['MJD']== 3

    dicta = {'HA':1,'PA':2,'MJD':3}
    res=dictaddentry(dicta,'HA','after','new',(3,4))
    assert res['HA']== 1
    assert res['new']== (3, 4)
    assert res['PA']== 2
    assert res['MJD'] ==3
