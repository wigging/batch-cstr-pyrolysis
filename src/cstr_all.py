"""
Run a series of CSTR models for all FCIC feedstocks using the Debiagi 2018
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

import cantera as ct
import chemics as cm
import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()

# Parameters
# ----------------------------------------------------------------------------

diam = 0.0525    # inner reactor diameter [m], 5.25 cm
length = 0.4318  # total reactor length or height [m], 43.18 cm
temp = 773.15    # reactor temperature [K]
pabs = 101_325   # reactor absolute pressure [Pa]
tau = 20         # residence time [s]

ghr_bio = 420    # biomass inlet feedrate [g/hr]
slm_n2 = 14      # inlet nitrogen gas flowrate [SLM]

n_cstrs = 1000                      # number of CSTRs in series
energy = 'off'                      # reactor energy
cti = 'data/debiagi_sw_metan2.cti'  # Cantera input file

# Feedstocks
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Reactor inputs
# ----------------------------------------------------------------------------

# Convert biomass feedrate from g/hr to kg/s
mf_bio = ghr_bio / 1000 / 3600

# Convert nitrogen gas flowrate from SLM to kg/s
# 1 m³/s = 60,000 liter/minute
# N₂ gas molecular weight = 28 g/mol
# N₂ gas density at STP = 1.2506 kg/m³
lpm_n2 = cm.slm_to_lpm(slm_n2, pabs / 1000, temp)
rhog_n2 = cm.rhog(28, 101325, 773.15)
mf_n2 = lpm_n2 / 60_000 * rhog_n2

# Mass fraction of nitrogen gas into reactor
y0_n2 = mf_n2 / (mf_n2 + mf_bio)

# Cantera CSTR model
# ----------------------------------------------------------------------------

# Store results for each feedstock
names = []
nf = len(feedstocks)
yf_gases = np.zeros(nf)
yf_liquids = np.zeros(nf)
yf_solids = np.zeros(nf)
yf_metas = np.zeros(nf)
exp_gases = np.zeros(nf)
exp_liquids = np.zeros(nf)
exp_solids = np.zeros(nf)
exp_ash = np.zeros(nf)

# Run CSTR reactor model for each feedstock
for i, feedstock in enumerate(feedstocks):
    print('Run', i, feedstock.name)

    # Calculate optimized biomass composition (daf) and splitting parameters
    bc, splits = feedstock.calc_biocomp()
    cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

    # Get feedstock moisture content as mass fraction
    yh2o = feedstock.prox_ad[3] / 100

    # All mass fractions for reactor input
    y0_all = {'N2': y0_n2, 'CELL': cell, 'GMSW': hemi, 'LIGC': ligc, 'LIGH': ligh,
              'LIGO': ligo, 'TANN': tann, 'TGL': tgl, 'ACQUA': yh2o}

    # Setup gas phase
    gas = ct.Solution(cti)
    gas.TPY = temp, pabs, y0_all

    # Create a stirred tank reactor (CSTR)
    dz = length / n_cstrs
    area = np.pi * (diam**2)
    vol = area * dz
    cstr = ct.IdealGasConstPressureReactor(gas, energy=energy, volume=vol)

    # Reservoirs for the inlet and outlet of each CSTR
    inlet = ct.Reservoir(gas)
    outlet = ct.Reservoir(gas)

    # Mass flow rate into the CSTR
    mf = ct.MassFlowController(upstream=inlet, downstream=cstr, mdot=cstr.mass / tau)

    # Determine pressure in the CSTR
    ct.PressureController(upstream=cstr, downstream=outlet, master=mf, K=1e-5)

    # Create a reactor network for performing the simulation
    sim = ct.ReactorNet([cstr])

    # Store the results for each CSTR
    states = ct.SolutionArray(cstr.thermo)
    states.append(cstr.thermo.state)

    for n in range(n_cstrs):
        gas.TPY = cstr.thermo.TPY
        inlet.syncState()
        sim.reinitialize()
        sim.advance_to_steady_state()
        states.append(cstr.thermo.state)

    # Species representing each phase
    gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2', 'N2')

    liquids = (
        'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
        'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H',
        'CH3OH', 'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12',
        'VANILLIN'
    )

    solids = (
        'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
        'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR', 'ACQUA'
    )

    metaplastics = (
        'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
        'GH2', 'GC2H6'
    )

    # Mass fractions including nitrogen and biomass (N₂ basis)
    y_n2 = states('N2').Y[:, 0]
    y_gas = states(*gases).Y.sum(axis=1)
    y_liquid = states(*liquids).Y.sum(axis=1)
    y_solid = states(*solids).Y.sum(axis=1)
    y_meta = states(*metaplastics).Y.sum(axis=1)

    # Mass fractions from only biomass (no nitrogen gas, N₂ free basis)
    sum_bio = y_gas + y_liquid + y_solid + y_meta - y_n2
    y_gas_bio = (y_gas - y_n2) / sum_bio
    y_liquid_bio = y_liquid / sum_bio
    y_solid_bio = y_solid / sum_bio
    y_meta_bio = y_meta / sum_bio

    # Store exit yields for feedstock
    names.append(feedstock.name)
    yf_gases[i] = y_gas_bio[-1]
    yf_liquids[i] = y_liquid_bio[-1]
    yf_solids[i] = y_solid_bio[-1]
    yf_metas[i] = y_meta_bio[-1]

    # Experiment lumped yields
    exp_gases[i] = feedstock.lump_yield[0]
    exp_liquids[i] = feedstock.lump_yield[1]
    exp_solids[i] = feedstock.lump_yield[2]
    exp_ash[i] = feedstock.prox_ad[2]

# Print
# ----------------------------------------------------------------------------

print(f'\n{" CSTR model ":*^70}\n')

print(
    'Parameters for reactor model\n'
    f'diam       {diam} m\n'
    f'length     {length} m\n'
    f'temp       {temp} K\n'
    f'pabs       {pabs:,} Pa\n'
    f'tau        {tau} s\n'
    f'n_cstrs    {n_cstrs}\n'
    f'energy     {energy}\n'
    f'cti file   {cti}'
)


# Plot
# ----------------------------------------------------------------------------

def style_barh(ax):
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.set_frame_on(False)
    ax.tick_params(color='0.8')
    ax.xaxis.grid(True, color='0.8')


y = np.arange(nf)

_, ax = plt.subplots(tight_layout=True, figsize=(9, 4.8))
b1 = ax.barh(y, yf_gases * 100, color='C6', label='Gases')
b2 = ax.barh(y, yf_liquids * 100, left=yf_gases * 100, color='C4', label='Liquids')
b3 = ax.barh(y, yf_solids * 100, left=(yf_gases + yf_liquids) * 100, color='C2', label='Solids')
b4 = ax.barh(y, yf_metas * 100, left=(yf_gases + yf_liquids + yf_solids) * 100, color='C1', label='Metaplastics')
ax.bar_label(b1, label_type='center', fmt='%.1f')
ax.bar_label(b2, label_type='center', fmt='%.1f')
ax.bar_label(b3, label_type='center', fmt='%.1f')
ax.bar_label(b4, label_type='center', fmt='%.1f')
ax.legend(bbox_to_anchor=[0.5, 1.02], loc='center', ncol=4, frameon=False)
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Final yield [wt. %]')
style_barh(ax)

# ---

h = 0.4
yf_gases = yf_gases * 100
yf_liquids = yf_liquids * 100
yf_solids = yf_solids * 100
yf_metas = yf_metas * 100
yf_solidmeta = yf_solids + yf_metas

_, ax = plt.subplots(tight_layout=True, figsize=(9, 8))

b1 = ax.barh(y, yf_gases, edgecolor='k', height=h, color='C6', label='Gases')
b2 = ax.barh(y, yf_liquids, edgecolor='k', left=yf_gases, height=h, color='C4', label='Liquids')
b3 = ax.barh(y, yf_solidmeta, edgecolor='k', left=yf_gases + yf_liquids, height=h, color='C2', label='Solids')

e1 = ax.barh(y + h, exp_gases, edgecolor='k', height=h, color='C6')
e2 = ax.barh(y + h, exp_liquids, edgecolor='k', height=h, left=exp_gases, color='C4')
e3 = ax.barh(y + h, exp_solids, edgecolor='k', height=h, left=exp_gases + exp_liquids, color='C2')

ax.bar_label(b1, label_type='center', fmt='%.1f')
ax.bar_label(e1, label_type='center', fmt='%.1f')
ax.bar_label(b2, label_type='center', fmt='%.1f')
ax.bar_label(e2, label_type='center', fmt='%.1f')
ax.bar_label(b3, label_type='center', fmt='%.1f')
ax.bar_label(e3, label_type='center', fmt='%.1f')

ax.legend(bbox_to_anchor=[0.5, 1.02], loc='center', ncol=3, frameon=False)
ax.set_yticks(y + h / 2)
ax.set_yticklabels(names)
ax.set_xlabel('Final yield [wt. %]')
style_barh(ax)


# ---

def style(ax):
    ax.grid(color='0.9')
    ax.tick_params(color='0.9')


pfit_gas_model = np.poly1d(np.polyfit(exp_ash, yf_gases, 1))
pfit_gas_exp = np.poly1d(np.polyfit(exp_ash, exp_gases, 1))

pfit_liq_model = np.poly1d(np.polyfit(exp_ash, yf_liquids, 1))
pfit_liq_exp = np.poly1d(np.polyfit(exp_ash, exp_liquids, 1))

pfit_sld_model = np.poly1d(np.polyfit(exp_ash, yf_solids, 1))
pfit_sld_exp = np.poly1d(np.polyfit(exp_ash, exp_solids, 1))

_, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(4.8, 8), sharex=True, tight_layout=True)

ax1.plot(exp_ash, yf_gases, 'ro', label='model')
ax1.plot(exp_ash, pfit_gas_model(exp_ash), 'r')
ax1.plot(exp_ash, exp_gases, 'ko', label='exp')
ax1.plot(exp_ash, pfit_gas_exp(exp_ash), 'k')
ax1.set_ylabel('Gases [wt. %]')
style(ax1)

ax2.plot(exp_ash, yf_liquids, 'ro', label='model')
ax2.plot(exp_ash, pfit_liq_model(exp_ash), 'r')
ax2.plot(exp_ash, exp_liquids, 'ko', label='exp')
ax2.plot(exp_ash, pfit_liq_exp(exp_ash), 'k')
ax2.set_ylabel('Liquids [wt. %]')
ax2.legend(loc='best')
style(ax2)

ax3.plot(exp_ash, yf_solids, 'ro', label='model')
ax3.plot(exp_ash, pfit_sld_model(exp_ash), 'r')
ax3.plot(exp_ash, exp_solids, 'ko', label='exp')
ax3.plot(exp_ash, pfit_sld_exp(exp_ash), 'k')
ax3.set_xlabel('Ash [wt. %]')
ax3.set_ylabel('Solids [wt. %]')
style(ax3)

plt.show()
