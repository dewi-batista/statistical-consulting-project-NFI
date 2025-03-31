# Probabilistic Modelling of mRNA Electropherograms in Fluid Mixtures
This repository contains implementations of Gamma and Gaussian mixture model approaches to fitting the distribution of genetic marker values conditioned on present fluid pairs.

## Dependencies:
**Python -** Python 3.8 (or more recent) is recommended, along with `pip` — Python’s package installer. All dependencies are listed in `requirements.txt` and can be installed with `pip install -r requirements.txt`.

**Note:** To execute `pip install -r requirements.txt`, ensure that the path of your console/environment is that of the root of this repository.

**R -** To install all relevant dependencies for our R program, execute `install.packages(c("readr", "dplyr", "mixtools", "evmix"))` in a console in which R is running.


## Programs:

- `data_preprocessing.py`: Applies all pre-processing to `./data/mixtures.csv` detailed in `report.pdf`, e.g. replacing NaNs with 0s, one-hot encoding the fluids present and dropping unneeded columns.
- `gaussianMixtures.py`: All relevant documentation is contained inside the file.
- `gammaMixtures.py`: All relevant documentation is contained inside the file.

**Note:** To run any of these programs, ensure that the path of your console/environment is that of the root of this repository.

## Other:
- `./data/`: Contains the raw data as well as their pre-processed equivalents.
- `./figures/`: Contains all relevant figures corresponding to our work. Has its own README for clarity.
- `NFI_Mixtures_Results.xlsx` Contains the BIC of the selected model for each marker distribution and fluid pair. Also contains the p-values corresponding to the two-sample KS test performed for said models. In the third sheet of the file, `N` denotes that the selected Gaussian mixture was the better fit (than the Gamma mixture fit) and `G` denotes that the selected Gamma mixture was the better fit.