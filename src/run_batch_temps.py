"""
Run a batch reactor model for a single FCIC feedstock at different
temperatures using the Debiagi 2018 biomass pyrolysis kinetics. All reactions
and chemical species are defined in a Cantera `cti` file.

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

temps = [723.15, 773.15, 823.15, 873.15]    # reactor temperatures [K]
p = 101325.0                                # reactor pressure [Pa]
time = np.linspace(0, 20.0, 100)            # reaction time steps [s]

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

# Cantera batch reactor
# ----------------------------------------------------------------------------

# Store results for each tempature
ntemps = len(temps)
ntime = len(time)
yt_gas = np.zeros((ntemps, ntime))
yt_liquid = np.zeros((ntemps, ntime))
yt_solid = np.zeros((ntemps, ntime))
yt_metaplastic = np.zeros((ntemps, ntime))

# Chemical species representing each phase
sp_gases = rct.sp_gases
sp_liquids = rct.sp_liquids
sp_solids = rct.sp_solids
sp_metaplastics = rct.sp_metaplastics

# Get feedstock moisture content as mass fraction
y0_h2o = feedstock.prox_ad[3] / 100

y0 = {'CELL': cell, 'GMSW': hemi, 'LIGC': ligc, 'LIGH': ligh, 'LIGO': ligo,
      'TANN': tann, 'TGL': tgl, 'ACQUA': y0_h2o}

for i in range(ntemps):
    temp = temps[i]
    states = rct.run_batch_simulation(cti, p, temp, time, y0, energy)

    # Mass fractions of the phases at each time step for a defined temperature
    yt_gas[i] = states(*sp_gases).Y.sum(axis=1)
    yt_liquid[i] = states(*sp_liquids).Y.sum(axis=1)
    yt_solid[i] = states(*sp_solids).Y.sum(axis=1)
    yt_metaplastic[i] = states(*sp_metaplastics).Y.sum(axis=1)

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

# Plot
# ----------------------------------------------------------------------------


def style(ax, xlabel, ylabel, loc=None, title=None):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if loc:
        ax.legend(loc=loc)
    if title:
        ax.set_title(title)


_, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(8, 4.8), sharey=True, tight_layout=True)

for i in range(ntemps):
    if temps[i] == 773.15:
        ax1.plot(time, yt_gas[i], 'k--')
        ax2.plot(time, yt_liquid[i], 'k--')
        ax3.plot(time, yt_solid[i] + yt_metaplastic[i], 'k--', label=f'{temps[i]} K')
    else:
        ax1.plot(time, yt_gas[i])
        ax2.plot(time, yt_liquid[i])
        ax3.plot(time, yt_solid[i] + yt_metaplastic[i], label=f'{temps[i]} K')

style(ax1, 'Time [s]', 'Mass fraction [-]', title='Gases')
style(ax2, 'Time [s]', '', title='Liquids')
style(ax3, 'Time [s]', '', loc='best', title='Solids')

plt.show()
