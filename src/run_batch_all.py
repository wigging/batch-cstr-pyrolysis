"""
Run a batch reactor model for all FCIC feedstocks using the Debiagi 2018
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

# Feedstocks
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Cantera batch reactor
# ----------------------------------------------------------------------------

# Store results for each feedstock
names = []
n = len(feedstocks)
yf_gas = np.zeros(n)
yf_liquid = np.zeros(n)
yf_solid = np.zeros(n)
yf_metaplastic = np.zeros(n)
exp_gases = np.zeros(n)
exp_liquids = np.zeros(n)
exp_solids = np.zeros(n)
exp_ash = np.zeros(n)

# Run batch reactor model for each feedstock
for i, feedstock in enumerate(feedstocks):
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

    # Final mixture mass fractions
    yf_gas[i] = y_gas[-1]
    yf_liquid[i] = y_liquid[-1]
    yf_solid[i] = y_solid[-1]
    yf_metaplastic[i] = y_metaplastic[-1]

    # Experiment lumped yields
    # exp_gases[i] = feedstock.lump_yield[0]
    # exp_liquids[i] = feedstock.lump_yield[1]
    # exp_solids[i] = feedstock.lump_yield[2]
    exp_gases[i] = feedstock.lump2_yield[0]
    exp_liquids[i] = feedstock.lump2_yield[1]
    exp_solids[i] = feedstock.lump2_yield[2]
    exp_ash[i] = feedstock.prox_ad[2]

# Print
# ----------------------------------------------------------------------------

print(
    f'\n{" Batch reactor model ":*^70}\n'
    '\nParameters for reactor model\n'
    f'final time    {time[-1]} s\n'
    f'temperature   {temp} K\n'
    f'pressure      {p:,} Pa\n'
    f'energy        {energy}\n'
    f'cti file      {cti}'
)


# Plot
# ----------------------------------------------------------------------------

def style_barh(ax):
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.set_frame_on(False)
    ax.tick_params(color='0.8')
    ax.xaxis.grid(True, color='0.8')


y = np.arange(n)
yf_gas = yf_gas * 100
yf_liquid = yf_liquid * 100
yf_solid = yf_solid * 100
yf_metaplastic = yf_metaplastic * 100
yf_solidmeta = yf_solid + yf_metaplastic

_, ax = plt.subplots(tight_layout=True, figsize=(9, 4.8))
b1 = ax.barh(y, yf_gas, color='C6', label='Gases')
b2 = ax.barh(y, yf_liquid, left=yf_gas, color='C4', label='Liquids')
b3 = ax.barh(y, yf_solid, left=yf_gas + yf_liquid, color='C2', label='Solids')
b4 = ax.barh(y, yf_metaplastic, left=yf_gas + yf_liquid + yf_solid, color='C1', label='Metaplastics')
ax.bar_label(b1, color='w', label_type='center', fmt='%.1f')
ax.bar_label(b2, color='w', label_type='center', fmt='%.1f')
ax.bar_label(b3, color='w', label_type='center', fmt='%.1f')
ax.bar_label(b4, color='w', label_type='center', fmt='%.1f')
ax.legend(bbox_to_anchor=[0.5, 1.02], loc='center', ncol=4, frameon=False)
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Final yield [wt. %]')
style_barh(ax)

# ---

h = 0.4

_, ax = plt.subplots(tight_layout=True, figsize=(9, 8))

b1 = ax.barh(y, yf_gas, edgecolor='k', height=h, color='C6', label='Gases')
b2 = ax.barh(y, yf_liquid, edgecolor='k', left=yf_gas, height=h, color='C4', label='Liquids')
b3 = ax.barh(y, yf_solidmeta, edgecolor='k', left=yf_gas + yf_liquid, height=h, color='C2', label='Solids')

e1 = ax.barh(y + h, exp_gases, edgecolor='k', height=h, color='C6')
e2 = ax.barh(y + h, exp_liquids, edgecolor='k', height=h, left=exp_gases, color='C4')
e3 = ax.barh(y + h, exp_solids, edgecolor='k', height=h, left=exp_gases + exp_liquids, color='C2')

ax.bar_label(b1, color='w', label_type='center', fmt='%.1f')
ax.bar_label(e1, color='w', label_type='center', fmt='%.1f')
ax.bar_label(b2, color='w', label_type='center', fmt='%.1f')
ax.bar_label(e2, color='w', label_type='center', fmt='%.1f')
ax.bar_label(b3, color='w', label_type='center', fmt='%.1f')
ax.bar_label(e3, color='w', label_type='center', fmt='%.1f')

ax.legend(bbox_to_anchor=[0.5, 1.02], loc='center', ncol=3, frameon=False)
ax.set_yticks(y + h / 2)
ax.set_yticklabels(names)
ax.set_xlabel('Final yield [wt. %]')
style_barh(ax)


# ---

def style(ax):
    ax.grid(color='0.9')
    ax.tick_params(color='0.9')


pfit_gas_model = np.poly1d(np.polyfit(exp_ash, yf_gas, 1))
pfit_gas_exp = np.poly1d(np.polyfit(exp_ash, exp_gases, 1))

pfit_liq_model = np.poly1d(np.polyfit(exp_ash, yf_liquid, 1))
pfit_liq_exp = np.poly1d(np.polyfit(exp_ash, exp_liquids, 1))

pfit_sld_model = np.poly1d(np.polyfit(exp_ash, yf_solid, 1))
pfit_sld_exp = np.poly1d(np.polyfit(exp_ash, exp_solids, 1))

_, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(4.8, 8), sharex=True, tight_layout=True)

ax1.plot(exp_ash, yf_gas, 'ro', label='model')
ax1.plot(exp_ash, pfit_gas_model(exp_ash), 'r')
ax1.plot(exp_ash, exp_gases, 'ko', label='exp')
ax1.plot(exp_ash, pfit_gas_exp(exp_ash), 'k')
ax1.set_ylabel('Gases [wt. %]')
style(ax1)

ax2.plot(exp_ash, yf_liquid, 'ro', label='model')
ax2.plot(exp_ash, pfit_liq_model(exp_ash), 'r')
ax2.plot(exp_ash, exp_liquids, 'ko', label='exp')
ax2.plot(exp_ash, pfit_liq_exp(exp_ash), 'k')
ax2.set_ylabel('Liquids [wt. %]')
ax2.legend(loc='best')
style(ax2)

ax3.plot(exp_ash, yf_solid, 'ro', label='model')
ax3.plot(exp_ash, pfit_sld_model(exp_ash), 'r')
ax3.plot(exp_ash, exp_solids, 'ko', label='exp')
ax3.plot(exp_ash, pfit_sld_exp(exp_ash), 'k')
ax3.set_xlabel('Ash [wt. %]')
ax3.set_ylabel('Solids [wt. %]')
style(ax3)

plt.show()
