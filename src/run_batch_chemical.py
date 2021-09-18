"""
Run a batch reactor model to compare final yield of a chemical species for
each feedstock using the Debiagi 2018 kinetics. According to Kristiina Iisa,
formaldehyde, propionaldehyde, heavy molecular weight lignin, and
acetaldehyde are detrimental to upgrading by catalytic pyrolysis. She also
suggests that in general acids are bad such as acetic acid.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import reactor as rct
from feedstock import Feedstock

# Parameters
# ----------------------------------------------------------------------------

temp = 773.15                       # reactor temperature [K]
p = 101325.0                        # reactor pressure [Pa]
time = np.linspace(0, 20.0, 100)    # reaction time steps [s]

energy = 'off'                      # reactor energy
cti = 'data/debiagi_sw_meta.cti'    # Cantera input file

# Select a chemical species used to predict final yield
# make sure `y_chemical[i]` is modified appropriately, see below
# make sure `ax.set_xlabel()` is modified appropriately, see below

# phenol is C6H5OH
# chemical = 'C6H5OH'

# furfural is FURFURAL
# chemical = 'FURFURAL'

# formaldehyde is CH2O
# chemical = 'CH2O'

# propionaldehyde is C2H5CHO
# chemical = 'C2H5CHO'

# heavy molecular weight lignin is C24H28O4
# chemical = 'C24H28O4'

# acetaldehyde is CH3CHO
chemical = 'CH3CHO'

# acetic acid is CH2OHCHO and CH3CO2H
# chemical = ('CH2OHCHO', 'CH3CO2H')

# Feedstocks
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Cantera batch reactor
# ----------------------------------------------------------------------------

n = len(feedstocks)

names = []
y_chemical = np.zeros(n)

# Run batch reactor model for each feedstock
for i in range(n):
    feedstock = feedstocks[i]
    names.append(feedstock.name)

    # Calculate optimized biomass composition (daf) and splitting parameters
    bc, splits = feedstock.calc_biocomp()
    cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

    # Get feedstock moisture content as mass fraction
    yh2o = feedstock.prox_ad[3] / 100

    # Perform batch reactor simulation
    y0 = {'CELL': cell, 'GMSW': hemi, 'LIGC': ligc, 'LIGH': ligh, 'LIGO': ligo,
          'TANN': tann, 'TGL': tgl, 'ACQUA': yh2o}

    states = rct.run_batch_simulation(cti, p, temp, time, y0, energy)

    y_chemical[i] = states(chemical).Y[-1]
    # y_chemical[i] = states(*chemical).Y.sum(axis=1)[-1]

# Plot
# ----------------------------------------------------------------------------

y = np.arange(n)

_, ax = plt.subplots(tight_layout=True)
ax.barh(y, y_chemical)
ax.set_xlabel('Acetaldehyde, mass fraction [-]')
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_frame_on(False)
ax.set_axisbelow(True)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

plt.show()
