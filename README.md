# pySpextool
Python code to reduce data obtained with the SpeX spectrograph.

# Recommended Installation Instructions

1) [Install miniconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html#), if you don't have it already.

2) Make a `pyspextool_3.13` conda environment and install necessary dependencies, including Python 3.13.
   ```
   conda create -n pyspextool_3.13 python=3.13
   ```

3) Activate the new environment
   ```
   conda activate pyspextool_3.13
   ```

4) Install the `pyspextool` package from PyPi
   ```
   pip install pyspextool
   ```

5) Optional - Install other Python packages

   If you want to use Jupyter notebooks or other Python tools with `pyspextool`, you need to install them in the new environment. For example:
   ```
   pip install jupyterlab
   ```


## Developer Instructions

If you plan to contribute code to `pyspextool`, you should clone a fork of this repo. 
You will also need `pytest` and the test data in order to run the tests.

1) Make your own fork of this repository.

2) Clone the fork to your computer.

3) Make and activate a dedicated virtual environment, ideally with Python 3.13.

4) Make an editable install of pyspextool and install the extra packages needed for testing and developing. In the` pyspextool/` direcotry:
```
pip install -e ".[test]"

```

5) Setup and download the `test_data/` submodule
   ```
   git submodule update --init --recursive 
   ```
