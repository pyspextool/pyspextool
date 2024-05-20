# pySpextool
Python code to reduce data obtained with the SpeX spectrograph.

# Installation Instructions
1) [Install miniconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html#), if you don't have it already.

2) Clone the pyspextool repository to your local computer
    ```
    git clone https://github.com/pyspextool/pyspextool.git
    ```
3) navigate into the `pyspextool` directory
   ```
   cd pyspextool/
   ```

4) Make a `pyspextool` conda environment and install necessary dependencies, including Python 3.11.
   ```
   conda env create -f environment.yml
   ```
5) Activate the new environment
   ```
   conda activate pyspextool_3.11
   ```
6) Install the package
   ```
   pip install -e ./
   ```



<!---
3) Install the package

#   **For Users**
    ```
    pip install pyspextool@git+https://github.com/pyspextool/pyspextool.git
    ```

    **For Developers**
    * clone the repository to your local computer
    ```
    git clone https://github.com/pyspextool/pyspextool.git
    ```
    * navigate into the `pyspextool` directory and install with pip/setuptools
    
    ```
    cd pyspextool
    pip install -e ./
    ```
    -->
