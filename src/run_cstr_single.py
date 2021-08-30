"""
Run a series of CSTR models for a single FCIC feedstock using the Debiagi 2018
biomass pyrolysis kinetics. All reactions and chemical species are defined in
Cantera `cti` files where `debiagi_sw.cti` is for softwood reactions,
`debiagi_hw.cti` is for hardwood reactions, and `debiagi_gr.cti` is for grass
reactions.

Note
----
The `debiagi_sw.cti` uses the original reaction rates from Debiagi 2018 while
`debiagi_sw_metan2.cti` uses modified rates for the metaplastic reactions
such that b=1 instead of b=0 where b is the temperature exponent Tᵇ.

Reference
---------
P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, and T. Faravelli.
A predictive model of biochar formation and characterization. Journal of
Analytical and Applied Pyrolysis, vol. 134, pp. 326-335, 2018.
"""

import chemics as cm
import json
import matplotlib.pyplot as plt
import reactor as rct
from feedstock import Feedstock

# Parameters
# ----------------------------------------------------------------------------

diam = 0.0525    # inner reactor diameter [m], 5.25 cm
length = 0.4318  # total reactor length or height [m], 43.18 cm
temp = 773.15    # reactor temperature [K]
p = 101_325      # reactor absolute pressure [Pa]

ghr_bio = 420    # biomass inlet feedrate [g/hr]
slm_n2 = 14      # inlet nitrogen gas flowrate [SLM]

n_cstrs = 1000                      # number of CSTRs in series
tau = 20 / n_cstrs                  # residence time [s]
energy = 'off'                      # reactor energy
cti = 'data/debiagi_sw_metan2.cti'  # Cantera input file

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
y0_h2o = feedstock.prox_ad[3] / 100

# Reactor inputs
# ----------------------------------------------------------------------------

# Convert biomass feedrate from g/hr to kg/s
mf_bio = ghr_bio / 1000 / 3600

# Convert nitrogen gas flowrate from SLM to kg/s
# 1 m³/s = 60,000 liter/minute
# N₂ gas molecular weight = 28 g/mol
# N₂ gas density at STP = 1.2506 kg/m³
lpm_n2 = cm.slm_to_lpm(slm_n2, p / 1000, temp)
rhog_n2 = cm.rhog(28, 101325, 773.15)
mf_n2 = lpm_n2 / 60_000 * rhog_n2

# Mass fraction of nitrogen gas into reactor
y0_n2 = mf_n2 / (mf_n2 + mf_bio)

# All mass fractions for reactor input
y0 = {'N2': y0_n2, 'CELL': cell, 'GMSW': hemi, 'LIGC': ligc, 'LIGH': ligh,
      'LIGO': ligo, 'TANN': tann, 'TGL': tgl, 'ACQUA': y0_h2o}

# Cantera CSTR model
# ----------------------------------------------------------------------------

states = rct.run_cstr_simulation(cti, diam, length, n_cstrs, p, tau, temp, y0)

# Mass fractions of phases and nitrogen gas in each CSTR (N₂ basis)
y_n2 = states('N2').Y[:, 0]
y_gas_n2 = states(*rct.sp_gases_n2).Y.sum(axis=1)
y_liquid_n2 = states(*rct.sp_liquids).Y.sum(axis=1)
y_solid_n2 = states(*rct.sp_solids).Y.sum(axis=1)
y_metaplastic_n2 = states(*rct.sp_metaplastics).Y.sum(axis=1)

# Mass fractions of the phases excluding nitrogen gas in each CSTR (N₂ free basis)
sum_no_n2 = y_gas_n2 + y_liquid_n2 + y_solid_n2 + y_metaplastic_n2 - y_n2
y_gas = (y_gas_n2 - y_n2) / sum_no_n2
y_liquid = y_liquid_n2 / sum_no_n2
y_solid = y_solid_n2 / sum_no_n2
y_metaplastic = y_metaplastic_n2 / sum_no_n2

# Carbon, hydrogen, and oxygen fractions
yc_gases, yh_gases, yo_gases = rct.get_ycho_gases(states)
yc_liquids, yh_liquids, yo_liquids = rct.get_ycho_liquids(states)

# Print
# ----------------------------------------------------------------------------

print(f'\n{" CSTR model ":*^70}\n')

print(
    'Parameters for reactor model\n'
    f'feed       {feedstock.name + ", Cycle " + str(feedstock.cycle)}\n'
    f'diam       {diam} m\n'
    f'length     {length} m\n'
    f'temp       {temp} K\n'
    f'p          {p:,} Pa\n'
    f'tau        {tau} s\n'
    f'n_cstrs    {n_cstrs}\n'
    f'energy     {energy}\n'
    f'cti file   {cti}\n'
)

print(
    'Calculated parameters\n'
    f'mf_bio   {mf_bio:.4g} kg/s\n'
    f'mf_n2    {mf_n2:.4g} kg/s\n'
)

y_total_n2 = y_gas_n2[-1] + y_liquid_n2[-1] + y_solid_n2[-1] + y_metaplastic_n2[-1]
y_total = y_gas[-1] + y_liquid[-1] + y_solid[-1] + y_metaplastic[-1]

print(
    'Exit mass fractions\n'
    f'              N₂ basis    N₂ free\n'
    f'gas          {y_gas_n2[-1]:>8.4f} {y_gas[-1]:>10.4f}\n'
    f'liquid       {y_liquid_n2[-1]:>8.4f} {y_liquid[-1]:>10.4f}\n'
    f'solid        {y_solid_n2[-1]:>8.4f} {y_solid[-1]:>10.4f}\n'
    f'metaplastic  {y_metaplastic_n2[-1]:>8.4f} {y_metaplastic[-1]:>10.4f}\n'
    f'total        {y_total_n2:>8.2f} {y_total:>10.2f}\n'
)

print(
    'Exit yields as wt. %\n'
    f'                N₂ basis   N₂ free\n'
    f'gas           {y_gas_n2[-1] * 100:>8.2f} {y_gas[-1] * 100:>10.2f}\n'
    f'liquid        {y_liquid_n2[-1] * 100:>8.2f} {y_liquid[-1] * 100:>10.2f}\n'
    f'solid         {y_solid_n2[-1] * 100:>8.2f} {y_solid[-1] * 100:>10.2f}\n'
    f'metaplastic   {y_metaplastic_n2[-1] * 100:>8.2f} {y_metaplastic[-1] * 100:>10.2f}\n'
    f'total         {y_total_n2 * 100:>8.2f} {y_total * 100:>10.2f}\n\n'
    f'total solids  {(y_solid_n2[-1] + y_metaplastic_n2[-1]) * 100:>8.2f} {(y_solid[-1] + y_metaplastic[-1]) * 100:>10.2f}'
)

print(
    '\nFinal mixture CHO fractions (N₂ basis)\n'
    '       gas   liquid\n'
    f'yc    {yc_gases.sum(axis=1)[-1]:.2f}     {yc_liquids.sum(axis=1)[-1]:.2f}\n'
    f'yh    {yh_gases.sum(axis=1)[-1]:.2f}     {yh_liquids.sum(axis=1)[-1]:.2f}\n'
    f'yo    {yo_gases.sum(axis=1)[-1]:.2f}     {yo_liquids.sum(axis=1)[-1]:.2f}'
)


# Plot
# ----------------------------------------------------------------------------

def _config(ax, xlabel, ylabel, title=None, legend=None):
    """
    Configure and style the plot figure.
    """
    ax.grid(True, color='0.9')
    ax.set_frame_on(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(color='0.9')

    if legend == 'side':
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    elif legend == 'best':
        ax.legend(loc='best', frameon=False)


# Yields for N₂ basis
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(10, 4.8), tight_layout=True)

ax1.plot(y_gas_n2, label='gas')
ax1.plot(y_liquid_n2, label='liquid')
ax1.plot(y_solid_n2, label='solid')
ax1.plot(y_metaplastic_n2, label='metaplastic')
_config(ax1, xlabel='CSTR [-]', ylabel='Mass fraction [-]', legend='best')

ax2.plot(states.T)
_config(ax2, xlabel='CSTR [-]', ylabel='Temperature [K]')

ax3.plot(states.P / 1000)
_config(ax3, xlabel='CSTR [-]', ylabel='Pressure [kPa]')

# Yields for N₂ free basis
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(10, 4.8), tight_layout=True)

ax1.plot(y_gas, label='gas')
ax1.plot(y_liquid, label='liquid')
ax1.plot(y_solid, label='solid')
ax1.plot(y_metaplastic, label='metaplastic')
_config(ax1, xlabel='CSTR [-]', ylabel='Mass fraction [-]', legend='best')

ax2.plot(states.T)
_config(ax2, xlabel='CSTR [-]', ylabel='Temperature [K]')

ax3.plot(states.P / 1000)
_config(ax3, xlabel='CSTR [-]', ylabel='Pressure [kPa]')

plt.show()
