"""
Run a batch reactor model for a single FCIC feedstock using the Debiagi 2018
biomass pyrolysis kinetics. All reactions and chemical species are defined in
Cantera `cti` files where `debiagi_sw.cti` is for softwood reactions,
`debiagi_hw.cti` is for hardwood reactions, and `debiagi_gr.cti` is for grass
reactions.

Note
----
The `debiagi_sw.cti` uses the original reaction rates from Debiagi 2018 while
`debiagi_sw_meta.cti` uses modified rates for the metaplastic reactions.

Reference
---------
P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, and T. Faravelli.
A predictive model of biochar formation and characterization. Journal of
Analytical and Applied Pyrolysis, vol. 134, pp. 326-335, 2018.
"""

import cantera as ct
import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()

# Parameters
# ----------------------------------------------------------------------------

temp = 773.15                       # reactor temperature [K]
p = 101325.0                        # reactor pressure [Pa]
time = np.linspace(0, 20.0, 100)    # reaction time steps [s]

energy = 'off'                      # reactor energy
cti = 'data/debiagi_sw_meta.cti'    # Cantera input file

# Feedstock
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstock = Feedstock(fdata[0])     # change index to choose feedstock

# Biomass composition
# ----------------------------------------------------------------------------

# Calculate optimized biomass composition (daf) and splitting parameters
bc, splits = feedstock.calc_biocomp()
cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

# Get feedstock moisture content as mass fraction
yh2o = feedstock.prox_ad[3] / 100

# Cantera batch reactor
# ----------------------------------------------------------------------------

y0 = {'CELL': cell, 'GMSW': hemi, 'LIGC': ligc, 'LIGH': ligh, 'LIGO': ligo,
      'TANN': tann, 'TGL': tgl, 'ACQUA': yh2o}

gas = ct.Solution(cti)
gas.TPY = temp, p, y0

r = ct.IdealGasReactor(gas, energy=energy)
sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

for t in time:
    sim.advance(t)
    states.append(r.thermo.state, t=t)

# Chemical species representing gas, liquid, solid, metaplastic phases
sp_gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2')

sp_liquids = (
    'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
    'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H',
    'CH3OH', 'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12',
    'VANILLIN'
)

sp_solids = (
    'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
    'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR', 'ACQUA'
)

sp_metaplastics = (
    'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
    'GH2', 'GC2H6'
)

# Sum of chemical species mass fractions at each time step
y_gases = states(*sp_gases).Y.sum(axis=1)
y_liquids = states(*sp_liquids).Y.sum(axis=1)
y_solids = states(*sp_solids).Y.sum(axis=1)
y_metaplastics = states(*sp_metaplastics).Y.sum(axis=1)

# Print
# ----------------------------------------------------------------------------

print(
    f'\n{" Batch reactor model ":*^70}\n'
    '\nParameters for reactor model\n'
    f'feedstock     {feedstock.name + ", Cycle " + str(feedstock.cycle)}\n'
    f'final time    {time[-1]} s\n'
    f'temperature   {temp} K\n'
    f'pressure      {p:,} Pa\n'
    f'energy        {energy}\n'
    f'cti file      {cti}'
)

print(
    '\nOptimized splitting parameters\n'
    f'α         {splits[0]:.4f}\n'
    f'β         {splits[1]:.4f}\n'
    f'γ         {splits[2]:.4f}\n'
    f'δ         {splits[3]:.4f}\n'
    f'ε         {splits[4]:.4f}\n'
    '\nChemical analysis, wt. ％ daf\n'
    '           exp      opt\n'
    f'cell      {feedstock.chem_bc[0]:.2f}    {cell * 100:.2f}\n'
    f'hemi      {feedstock.chem_bc[1]:.2f}    {hemi * 100:.2f}\n'
    f'lignin    {feedstock.chem_bc[2]:.2f}    {(ligc + ligh + ligo) * 100:.2f}\n'
    '\nBiomass composition, wt. ％ daf\n'
    f'cell      {cell * 100:.2f}\n'
    f'hemi      {hemi * 100:.2f}\n'
    f'lig-c     {ligc * 100:.2f}\n'
    f'lig-h     {ligh * 100:.2f}\n'
    f'lig-o     {ligo * 100:.2f}\n'
    f'tann      {tann * 100:.2f}\n'
    f'tgl       {tgl * 100:.2f}'
)

print(
    '\nLumped species final yields as wt. %\n'
    f'gases         {y_gases[-1] * 100:.2f}\n'
    f'liquids       {y_liquids[-1] * 100:.2f}\n'
    f'solids        {y_solids[-1] * 100:.2f}\n'
    f'metaplastics  {y_metaplastics[-1] * 100:.2f}\n'
    f'total solids  {(y_solids[-1] + y_metaplastics[-1]) * 100:.2f}'
)


# Plot
# ----------------------------------------------------------------------------

def style(ax, xlabel, ylabel, loc=None):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if loc:
        ax.legend(loc=loc)


# Initial biomass composition species
fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELL').Y[:, 0], label='CELL')
ax.plot(states.t, states('GMSW').Y[:, 0], label='GMSW')
# ax.plot(states.t, states('XYHW').Y[:, 0], label='XYHW')
# ax.plot(states.t, states('XYGR').Y[:, 0], label='XYGR')
ax.plot(states.t, states('LIGC').Y[:, 0], label='LIGC')
ax.plot(states.t, states('LIGH').Y[:, 0], label='LIGH')
ax.plot(states.t, states('LIGO').Y[:, 0], label='LIGO')
ax.plot(states.t, states('TANN').Y[:, 0], label='TANN')
ax.plot(states.t, states('TGL').Y[:, 0], label='TGL')
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]', loc='best')

# Intermediate species
fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELLA').Y[:, 0], label='CELLA')
ax.plot(states.t, states('HCE1').Y[:, 0], label='HCE1')
ax.plot(states.t, states('HCE2').Y[:, 0], label='HCE2')
ax.plot(states.t, states('LIGCC').Y[:, 0], label='LIGCC')
ax.plot(states.t, states('LIGOH').Y[:, 0], label='LIGOH')
ax.plot(states.t, states('LIG').Y[:, 0], label='LIG')
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]', loc='best')

# Lumped species
_, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, y_gases, label='gases')
ax.plot(states.t, y_liquids, label='liquids')
ax.plot(states.t, y_solids, label='solids')
ax.plot(states.t, y_metaplastics, label='metaplastics')
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]', loc='best')

# Reactor temperature
_, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states.T, color='m')
style(ax, xlabel='Time [s]', ylabel='Temperature [K]')

# Final yield for all chemical species
species = states.species_names
ys = [states(sp).Y[:, 0][-1] for sp in species]
ypos = np.arange(len(species))

fig, ax = plt.subplots(figsize=(6.4, 8), tight_layout=True)
ax.barh(ypos, ys, align='center')
ax.set_xlabel('Mass fraction [-]')
ax.set_ylim(min(ypos) - 1, max(ypos) + 1)
ax.set_yticks(ypos)
ax.set_yticklabels(species)
ax.invert_yaxis()
ax.set_axisbelow(True)
ax.set_frame_on(False)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

# Final yield for lumped species as gases, liquids, solids, metaplastics
y_gas = np.arange(len(sp_gases))
x_gas = states(*sp_gases).Y[-1]

y_liquid = np.arange(len(sp_liquids))
x_liquid = states(*sp_liquids).Y[-1]

y_solid = np.arange(len(sp_solids))
x_solid = states(*sp_solids).Y[-1]

y_meta = np.arange(len(sp_metaplastics))
x_meta = states(*sp_metaplastics).Y[-1]

_, axes = plt.subplots(2, 2, figsize=(10, 9.6), sharex='col', tight_layout=True)

axes[0, 0].barh(y_gas, x_gas, color='C6')
axes[0, 0].set_title('Gases')
axes[0, 0].set_ylim(min(y_gas) - 1, max(y_gas) + 1)
axes[0, 0].set_yticks(y_gas)
axes[0, 0].set_yticklabels(list(sp_gases))
axes[0, 0].invert_yaxis()
axes[0, 0].set_axisbelow(True)
axes[0, 0].set_frame_on(False)
axes[0, 0].tick_params(color='0.8')
axes[0, 0].xaxis.grid(True, color='0.8')

axes[0, 1].barh(y_liquid, x_liquid, color='C4')
axes[0, 1].set_title('Liquids')
axes[0, 1].set_ylim(min(y_liquid) - 1, max(y_liquid) + 1)
axes[0, 1].set_yticks(y_liquid)
axes[0, 1].set_yticklabels(list(sp_liquids))
axes[0, 1].invert_yaxis()
axes[0, 1].set_axisbelow(True)
axes[0, 1].set_frame_on(False)
axes[0, 1].tick_params(color='0.8')
axes[0, 1].xaxis.grid(True, color='0.8')

axes[1, 0].barh(y_solid, x_solid, color='C2')
axes[1, 0].set_xlabel('Mass fraction [-]')
axes[1, 0].set_title('Solids')
axes[1, 0].set_ylim(min(y_solid) - 1, max(y_solid) + 1)
axes[1, 0].set_yticks(y_solid)
axes[1, 0].set_yticklabels(list(sp_solids))
axes[1, 0].invert_yaxis()
axes[1, 0].set_axisbelow(True)
axes[1, 0].set_frame_on(False)
axes[1, 0].tick_params(color='0.8')
axes[1, 0].xaxis.grid(True, color='0.8')

axes[1, 1].barh(y_meta, x_meta, color='C1')
axes[1, 1].set_xlabel('Mass fraction [-]')
axes[1, 1].set_title('Metaplastics')
axes[1, 1].set_ylim(min(y_meta) - 1, max(y_meta) + 1)
axes[1, 1].set_yticks(y_meta)
axes[1, 1].set_yticklabels(list(sp_metaplastics))
axes[1, 1].invert_yaxis()
axes[1, 1].set_axisbelow(True)
axes[1, 1].set_frame_on(False)
axes[1, 1].tick_params(color='0.8')
axes[1, 1].xaxis.grid(True, color='0.8')

plt.show()
