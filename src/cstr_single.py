"""
Run a series of CSTR models for a single FCIC feedstock using the Debiagi 2018
biomass pyrolysis kinetics. All reactions and chemical species are defined in
Cantera `cti` files where `debiagi_sw.cti` is for softwood reactions,
`debiagi_hw.cti` is for hardwood reactions, and `debiagi_gr.cti` is for grass
reactions.

Note
----
The `debiagi_sw.cti` uses the original reaction rates from Debiagi 2018 while
`debiagi_sw2.cti` uses modified rates for the metaplastic reactions such that
b=1 instead of b=0 where b is the temperature exponent Tᵇ.

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
from objfunc import objfunc
from scipy.optimize import minimize

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

# Feedstock
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstock = Feedstock(fdata[0])     # change index to choose feedstock

# Biomass composition
# ----------------------------------------------------------------------------

# Get feedstock name, ultimate analysis data, and chemical analysis data
ult_daf = feedstock.calc_ult_daf()
ult_cho = feedstock.calc_ult_cho(ult_daf)
chem_daf = feedstock.calc_chem_daf()
chem_bc = feedstock.calc_chem_bc(chem_daf) / 100

# C and H mass fractions from ultimate analysis data
yc = ult_cho[0] / 100
yh = ult_cho[1] / 100

# Determine optimized splitting parameters using default values for `x0`
# where each parameter is bound within 0 to 1
x0 = [0.6, 0.8, 0.8, 1, 1]
bnds = ((0, 1), (0, 1), (0, 1), (0, 1), (0, 1))
res = minimize(objfunc, x0, args=(yc, yh, chem_bc), method='L-BFGS-B', bounds=bnds)

# Calculate biomass composition as dry ash-free basis (daf)
bc = cm.biocomp(yc, yh, alpha=res.x[0], beta=res.x[1], gamma=res.x[2], delta=res.x[3], epsilon=res.x[4])
cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

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
y_n2 = mf_n2 / (mf_n2 + mf_bio)

# All mass fractions for reactor input
y_all = {
    'N2': y_n2,
    'CELL': cell,
    'GMSW': hemi,
    'LIGC': ligc,
    'LIGH': ligh,
    'LIGO': ligo,
    'TANN': tann,
    'TGL': tgl
}

# Cantera CSTR model
# ----------------------------------------------------------------------------

# Setup gas phase
gas = ct.Solution(cti)
gas.TPY = temp, pabs, y_all

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
    f'feed       {feedstock.name + ", Cycle " + str(feedstock.cycle)}\n'
    f'cti file   {cti}\n'
)

print(
    'Calculated parameters\n'
    f'mf_bio   {mf_bio:.4g} kg/s\n'
    f'mf_n2    {mf_n2:.4g} kg/s\n'
)

y_total = y_gas[-1] + y_liquid[-1] + y_solid[-1] + y_meta[-1]
y_total_bio = y_gas_bio[-1] + y_liquid_bio[-1] + y_solid_bio[-1] + y_meta_bio[-1]

print(
    'Exit mass fractions\n'
    f'                N₂ basis   N₂ free\n'
    f'gases         {y_gas[-1]:>8.4f} {y_gas_bio[-1]:>10.4f}\n'
    f'liquids       {y_liquid[-1]:>8.4f} {y_liquid_bio[-1]:>10.4f}\n'
    f'solids        {y_solid[-1]:>8.4f} {y_solid_bio[-1]:>10.4f}\n'
    f'metaplastics  {y_meta[-1]:>8.4f} {y_meta_bio[-1]:>10.4f}\n'
    f'total         {y_total:>8.2f} {y_total_bio:>10.2f}\n'
)

print(
    'Exit yields wt. %\n'
    f'                N₂ basis   N₂ free\n'
    f'gases         {y_gas[-1] * 100:>8.2f} {y_gas_bio[-1] * 100:>10.2f}\n'
    f'liquids       {y_liquid[-1] * 100:>8.2f} {y_liquid_bio[-1] * 100:>10.2f}\n'
    f'solids        {y_solid[-1] * 100:>8.2f} {y_solid_bio[-1] * 100:>10.2f}\n'
    f'metaplastics  {y_meta[-1] * 100:>8.2f} {y_meta_bio[-1] * 100:>10.2f}\n'
    f'total         {y_total * 100:>8.2f} {y_total_bio * 100:>10.2f}\n\n'
    f'total solids  {(y_solid[-1] + y_meta[-1]) * 100:>8.2f} {(y_solid_bio[-1] + y_meta_bio[-1]) * 100:>10.2f}'
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

ax1.plot(y_gas, label='gas')
ax1.plot(y_liquid, label='liquid')
ax1.plot(y_solid, label='solid')
ax1.plot(y_meta, label='meta')
_config(ax1, xlabel='CSTR [-]', ylabel='Mass fraction [-]', legend='best')

ax2.plot(states.T)
_config(ax2, xlabel='CSTR [-]', ylabel='Temperature [K]')

ax3.plot(states.P / 1000)
_config(ax3, xlabel='CSTR [-]', ylabel='Pressure [kPa]')

# Yields for N₂ free basis
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(10, 4.8), tight_layout=True)

ax1.plot(y_gas_bio, label='gas')
ax1.plot(y_liquid_bio, label='liquid')
ax1.plot(y_solid_bio, label='solid')
ax1.plot(y_meta_bio, label='meta')
_config(ax1, xlabel='CSTR [-]', ylabel='Mass fraction [-]', legend='best')

ax2.plot(states.T)
_config(ax2, xlabel='CSTR [-]', ylabel='Temperature [K]')

ax3.plot(states.P / 1000)
_config(ax3, xlabel='CSTR [-]', ylabel='Pressure [kPa]')

plt.show()
