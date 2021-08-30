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

# Feedstock
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstock = Feedstock(fdata[1])     # change index to choose feedstock

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

states = rct.run_batch_simulation(cti, p, temp, time, y0, energy='off')

# Chemical species representing each phase
sp_gases = rct.sp_gases
sp_liquids = rct.sp_liquids
sp_solids = rct.sp_solids
sp_metaplastics = rct.sp_metaplastics

# Mass fractions of the phases at each time step
y_gas = states(*sp_gases).Y.sum(axis=1)
y_liquid = states(*sp_liquids).Y.sum(axis=1)
y_solid = states(*sp_solids).Y.sum(axis=1)
y_metaplastic = states(*sp_metaplastics).Y.sum(axis=1)

# Carbon, hydrogen, and oxygen fractions
yc_gases, yh_gases, yo_gases = rct.get_ycho_gases(states)
yc_liquids, yh_liquids, yo_liquids = rct.get_ycho_liquids(states)

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
    '\nFinal mixture yields, wt. %\n'
    f'gas           {y_gas[-1] * 100:.2f}\n'
    f'liquid        {y_liquid[-1] * 100:.2f}\n'
    f'solid         {y_solid[-1] * 100:.2f}\n'
    f'metaplastic   {y_metaplastic[-1] * 100:.2f}\n'
    f'total solids  {(y_solid[-1] + y_metaplastic[-1]) * 100:.2f}'
)

print(
    '\nFinal mixture CHO fractions\n'
    '       gas   liquid\n'
    f'yc    {yc_gases.sum(axis=1)[-1]:.2f}     {yc_liquids.sum(axis=1)[-1]:.2f}\n'
    f'yh    {yh_gases.sum(axis=1)[-1]:.2f}     {yh_liquids.sum(axis=1)[-1]:.2f}\n'
    f'yo    {yo_gases.sum(axis=1)[-1]:.2f}     {yo_liquids.sum(axis=1)[-1]:.2f}'
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


# ---

_, ax = plt.subplots()
ax.plot(states.t, yc_gases.sum(axis=1), label='carbon')
ax.plot(states.t, yh_gases.sum(axis=1), label='hydrogen')
ax.plot(states.t, yo_gases.sum(axis=1), label='oxygen')
style(ax, xlabel='Time [s]', ylabel='Gas mixture CHO fractions [-]', loc='best')

# ---

_, ax = plt.subplots()
ax.plot(states.t, yc_liquids.sum(axis=1), label='carbon')
ax.plot(states.t, yh_liquids.sum(axis=1), label='hydrogen')
ax.plot(states.t, yo_liquids.sum(axis=1), label='oxygen')
style(ax, xlabel='Time [s]', ylabel='Liquid mixture CHO fractions [-]', loc='best')

# ---

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

# ---

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELLA').Y[:, 0], label='CELLA')
ax.plot(states.t, states('HCE1').Y[:, 0], label='HCE1')
ax.plot(states.t, states('HCE2').Y[:, 0], label='HCE2')
ax.plot(states.t, states('LIGCC').Y[:, 0], label='LIGCC')
ax.plot(states.t, states('LIGOH').Y[:, 0], label='LIGOH')
ax.plot(states.t, states('LIG').Y[:, 0], label='LIG')
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]', loc='best')

# ---

_, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, y_gas, label='gas')
ax.plot(states.t, y_liquid, label='liquid')
ax.plot(states.t, y_solid, label='solid')
ax.plot(states.t, y_metaplastic, label='metaplastic')
style(ax, xlabel='Time [s]', ylabel='Mixture mass fraction [-]', loc='best')

# ---

_, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states.T, color='m')
style(ax, xlabel='Time [s]', ylabel='Temperature [K]')

# ---

y_gas = np.arange(len(sp_gases))
x_gas = states(*sp_gases).Y[-1]

y_liquid = np.arange(len(sp_liquids))
x_liquid = states(*sp_liquids).Y[-1]

y_solid = np.arange(len(sp_solids))
x_solid = states(*sp_solids).Y[-1]

y_meta = np.arange(len(sp_metaplastics))
x_meta = states(*sp_metaplastics).Y[-1]

_, ax = plt.subplots(tight_layout=True)
ax.barh(y_gas, x_gas, color='C6')
ax.set_xlabel('Gas species mass fraction [-]')
ax.set_yticks(y_gas)
ax.set_yticklabels(list(sp_gases))
ax.set_frame_on(False)
ax.set_axisbelow(True)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

_, ax = plt.subplots(tight_layout=True)
ax.barh(y_liquid, x_liquid, color='C4')
ax.set_xlabel('Liquid species mass fraction [-]')
ax.set_yticks(y_liquid)
ax.set_yticklabels(list(sp_liquids))
ax.set_frame_on(False)
ax.set_axisbelow(True)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

_, ax = plt.subplots(tight_layout=True)
ax.barh(y_solid, x_solid, color='C2')
ax.set_xlabel('Solid species mass fraction [-]')
ax.set_yticks(y_solid)
ax.set_yticklabels(list(sp_solids))
ax.set_frame_on(False)
ax.set_axisbelow(True)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

_, ax = plt.subplots(tight_layout=True)
ax.barh(y_meta, x_meta, color='C1')
ax.set_xlabel('Metaplastic species mass fraction [-]')
ax.set_yticks(y_meta)
ax.set_yticklabels(list(sp_metaplastics))
ax.set_frame_on(False)
ax.set_axisbelow(True)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

plt.show()
