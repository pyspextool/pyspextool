# pySpextool
Python code to reduce data obtained with the SpeX spectrograph.

# Recommended Installation Instructions

1) Setup your computer and install miniforge, if you haven't already.

   Follow Steps 0-2 of the [Astropy Workshop Setup instructions](https://github.com/astropy/astropy-workshop/tree/main/00-Install_and_Setup) to install WSL (for Windows) and install conda with miniforge.


2) Make a `pyspextool_3.13` conda environment with Python 3.13.
   ```
   conda create -n pyspextool_3.13 python=3.13
   ```

3) Activate the new environment
   ```
   conda activate pyspextool_3.13
   ```

4) Install the `pyspextool` package from PyPI using `pip`. This also installs the other Python packages that `pyspextool` depends upon.
   ```
   pip install pyspextool
   ```

5) Install other Python packages

   If you want to use Jupyter notebooks or other Python tools with `pyspextool`, you need to install them in the new environment. For example:
   ```
   pip install jupyterlab
   ```


## Instructions for using the example notebooks

1) The Jupyter notebook tutorials are under `notebooks/` for preupgraded and upgraded SpeX Prism, SXD, and LXD data.

2) Follow the instructions above to download the test data which will be saved under `tests/test_data/`. To get the data used by the Jupyter notebooks, clone the [test_data](https://github.com/pyspextool/test_data) repository.

   ```
   git clone git@github.com:pyspextool/test_data.git
   ```

   The raw and processed example data products are available under `tests/test_data/raw/` and `tests/test_data/processed/` respectively. The data products include the `cal` (calibration data), `proc` (processed data), and `qa` (quality assurance plots).

3) While users should be able to run the Jupyter notebook tutorials without changing the output paths, users are encouraged to move the data folder and rename the output paths. The default saving directory is `/test_pyspextool/test_data/processed/`.


## Developer Instructions

If you plan to contribute code to `pyspextool`, you should clone a fork of this repo. 
You will also need `pytest` and the test data in order to run the tests.

1) Make your own fork of this repository.

2) Clone the fork to your computer.

3) Make and activate a dedicated virtual environment, ideally with Python 3.13.

4) Install and editable version of `pyspextool` and install the extra packages needed for testing and developing. In the` pyspextool/` direcotry:
   ```
   pip install -e ".[test]"

   ```

5) Setup and download the `test_data/` submodule
   ```
   git submodule update --init --recursive 
   ```
