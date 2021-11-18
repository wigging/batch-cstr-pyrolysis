# CSTR and batch reactor modeling for biomass pyrolysis

This repository contains Python code for running batch and continuous stirred tank reactor (CSTR) models to estimate thermochemical biomass conversion at fast pyrolysis conditions. Results from the reduced-order models are compared to bubbling fluidized bed reactor experiments using various biomass derived feedstocks. This work was supported by the Feedstock-Conversion Interface Consortium (FCIC).

## Data

The `data` folder provides Cantera kinetic mechanism files (`.cti`) for biomass pyrolysis. This folder also contains the `data/feedstocks.json` file which provides parameters related to each feedstock.

## Diagrams

The `diagrams` folder includes `.drawio` files. These files can be edited at https://app.diagrams.net which is an online flowchart and diagram software.

## Source code

The `src` folder contains all the Python code for the reactor models. Files prepended with `run_` can be executed from the command line. See the comments in each file for more information. Examples of running the code in a terminal are given below.

```bash
# Go to the project directory
$ cd fcic-pyrolysis

# Run a batch reactor model for a single feedstock
$ python src/run_batch_single.py

# Run a batch reactor model for all the feedstocks and compare results
$ python src/run_batch_all.py

# Determine the biomass composition for a single feedstock
$ python src/run_biocomp_single.py
```

## Documentation

Documentation for this project is generated with LaTeX. The `tex` folder contains all the LaTeX files along with the associated figure files. See the `tex/main.pdf` to read the documentation or click [here](https://github.com/wigging/batch-cstr-pyrolysis/blob/main/tex/main.pdf) to view the PDF online.

## Citation

To cite this work, use the "Cite this repository" feature available on the right side of this repository page or use the reference text given below.

> Gavin Wiggins. CSTR and batch reactor modeling for biomass pyrolsyis. GitHub repository. Available at https://github.com/wigging/batch-cstr-pyrolysis.
