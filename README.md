# Probabilistic Modelling of mRNA Electropherograms in Fluid Mixtures
This repository contains implementations of Gamma and Gaussian mixture model approaches to fitting the distribution of genetic marker values conditioned on present fluid pairs.

## Dependencies:
Python 3.8 (or more recent) is recommended, along with `pip` — Python’s package installer. All dependencies are listed in `requirements.txt` and can be installed with `pip install -r requirements.txt`.

**Note:** Ensure that the path of the console/environment used to execute `pip install -r requirements.txt` is that of the root of this repository.

## Programs:

- `data_preprocessing.py`: Applies all pre-processing to `./data/mixtures.csv` detailed in `report.pdf`, e.g. replacing NaNs with 0s, one-hot encoding the fluids present and dropping unneeded columns.
- `gaussianMixtures.py`: All relevant documentation is contained inside the file.
- `gammaMixtures.py`: All relevant documentation is contained inside the file.

**Note:** Ensure that the path of the console/environment used to run these programs is that of the root of this repository.

## Other:
- `./data/`: Contains the raw data as well as their pre-processed equivalents.
- `./figures/`: Contains all relevant figures corresponding to our work. Has its own README for clarity.
- `NFI_Mixtures_Results.xlsx` includes ...