from pyspextool.utils.add_entry import *

dict = {'HA':1, 'PA':2, 'MJD':3}
result = add_entry(dict, 'HA','before','new',4)
assert result == {'HA': 1, 'PA': 2, 'new': 4, 'MJD': 3}

dict = {'HA':1, 'PA':2, 'MJD':3}
result = add_entry(dict, 'HA','after','new',(3,4))
assert result == {'HA': 1, 'new': (3, 4), 'PA': 2, 'MJD': 3}


