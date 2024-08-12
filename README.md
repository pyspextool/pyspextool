# pySpextool
Python code to reduce data obtained with the SpeX spectrograph.

# Installation Instructions
1) [Install miniconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html#), if you don't have it already.

2) Clone the pyspextool repository to your local computer
    ```
    git clone https://github.com/pyspextool/pyspextool.git
    ```

3) Navigate into the `pyspextool/` directory
   ```
   cd pyspextool/
   ```

4) Make a `pyspextool_3.11` conda environment and install necessary dependencies, including Python 3.11.
   ```
   conda env create -f environment.yml
   ```

5) Activate the new environment
   ```
   conda activate pyspextool_3.11
   ```
6) Install the `pyspextool` package
   ```
   pip install -e ./
   ```

7) Optional - Install other Python packages

   If you want to use Jupyter notebooks or other Python tools with `pyspextool`, you need to install them in the new environment. For example:
   ```
   pip install jupyterlab
   ```


## Developer Instructions

If you plan to contribute code to `pyspextool`, you'll need `pytest` and the test data in order to run the tests.

8) Install pytest
   ```
   pip install pytest
   ```

9) Setup and download the `test_data/` submodule
   ```
   git submodule update --init --recursive 
   ```
