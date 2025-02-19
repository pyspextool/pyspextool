from pyspextool import config as setup
from pyspextool.setup_utils import __version__
from importlib.metadata import version
from pyspextool.setup_utils import pyspextool_setup

def test_version():
    assert isinstance(__version__, str)
    assert __version__ == version("pyspextool")    
    assert setup.state["version"] == __version__
    

def test_version_setup():
    pyspextool_setup()
    assert setup.state["version"] == __version__