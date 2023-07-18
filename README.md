# pySpextool
Python code to reduce data obtained with the SpeX spectrograph.

# Installation Instructions
1) [Install git-lfs](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage) 
The calibration files stored in the `data/` and `instruments/` directories
are stored using [Github Large File Storage](https://docs.github.com/en/repositories/working-with-files/managing-large-files). Git-lfs must be installed before using this package.

2) Ensure you are in a Python 3.9 environment and the packages listed in [`requirements.txt`](https://github.com/pyspextool/pyspextool/blob/main/requirements.txt) file are installed.

3) Install the package

   **For Users**
    ```
    pip install -e pyspextool@git+https://github.com/pyspextool/pyspextool.git
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

