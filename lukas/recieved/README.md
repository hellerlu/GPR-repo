# GPR-plotting

## Installation
Open your terminal and go the location where you want to store the project. Then type the following replying ```y``` if needed:

```bash
mkdir gpr
cd gpr
conda update -n base -c defaults conda
conda create -n gpr python=3.7 jupyter
conda activate gpr
conda install pandas h5py pytz 
conda install -c conda-forge obspy
pip install readgssi
```