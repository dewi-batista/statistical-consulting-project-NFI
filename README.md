# Probabilistic Modelling of mRNA Electropherograms in Fluid Mixtures

## Overview:
This repository contains implementations of Gausian and Gamma mixture model approaches to fitting the distribution of genetic marker values conditioned on present fluid pairs.

## Dependencies:
Most programs in this repository are written in Python. Python 3.8 or above is recommended, along with `pip` — Python’s package installer. All dependencies are listed in `requirements.txt` and can be installed with `pip install -r requirements.txt`.

**Note:** Ensure that the path of the console used to execute `pip install -r requirements.txt` is that of the root of this repository.

## Program descriptions:

Mixture model fitting:
- `gaussian_mixture.py`: Fits a bunch of Gaussian mixtures.
- `gamma_mixture.py`: 

Data exploration/pre-processing:
- `analysis_of_correlations.py`:
- `data_preprocessing.py`: Applies all pre-processing to `./data/mixtures.csv` detailed in `report.pdf`, e.g. replacing NaNs with 0s, one-hot encoding the fluids present and drops unneeded columns.
- `visuals_of_2D_histograms.py`: 

**Note:** Ensure that the path of the console/environment used to run any of these programs is that of the root of this repository.

## Miscellaneous:
- `blelele.xlsx` includes ...
- `./data/`: contains
- `./figures/`: Contains all relevant figures from our work. It has its own README detailing the figures included.