# FCIC pyrolysis modeling

This repository contains Python code for modeling the pyrolysis of biomass feedstocks associated with the FCIC project. The reduced-order models represent a bubbling fluidized bed reactor operating a fast pyrolysis conditions.

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

Documentation for this project is generated with LaTeX. The `tex` folder contains all the LaTeX files along with the associated figure files. See the `tex/main.pdf` to read the documentation.
