"""
here
"""

import cantera as ct
import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()

# Parameters
# ----------------------------------------------------------------------------

diam = 0.0525    # inner reactor diameter [m], 5.25 cm
length = 0.4318  # total reactor length or height [m], 43.18 cm
temp = 773.15    # reactor temperature [K]
pabs = 101_325   # reactor absolute pressure [Pa]
ghr_bio = 420    # biomass inlet feedrate [g/hr]
slm_n2 = 14      # inlet nitrogen gas flowrate [SLM]

n_cstrs = 20     # number of CSTRs in series
energy = 'off'   # reactor energy

# biomass composition [cell, hemi, ligc, ligh, ligo, tann, tgl]
bc = [0.2898, 0.2202, 0.0058, 0.0879, 0.2716, 0.0160, 0.1088]

# Calculated parameters
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
    'CELL': bc[0],
    'GMSW': bc[1],
    'LIGC': bc[2],
    'LIGH': bc[3],
    'LIGO': bc[4],
    'TANN': bc[5],
    'TGL': bc[6]
}

# CSTR model
# ----------------------------------------------------------------------------

# Setup gas phase
gas = ct.Solution('data/debiagi_sw_n2.cti')
gas.TPY = temp, pabs, y_all

# Create a stirred tank reactor (CSTR)
cstr = ct.IdealGasConstPressureReactor(gas, energy=energy)

dz = length / n_cstrs
area = np.pi * (diam**2)
cstr.volume = area * dz

# Reservoirs for the inlet and outlet of each CSTR
inlet = ct.Reservoir(gas)
outlet = ct.Reservoir(gas)

# Mass flow rate into the CSTR
mf = ct.MassFlowController(upstream=inlet, downstream=cstr, mdot=mf_n2 + mf_bio)

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

# Mass fractions including nitrogen and biomass
y_n2 = states('N2').Y[:, 0]
y_gas = states(*gases).Y.sum(axis=1)
y_liquid = states(*liquids).Y.sum(axis=1)
y_solid = states(*solids).Y.sum(axis=1)
y_meta = states(*metaplastics).Y.sum(axis=1)

# Mass fractions from only biomass (no nitrogen gas)
sum_bio = y_gas + y_liquid + y_solid + y_meta - y_n2
y_gas_bio = (y_gas - y_n2) / sum_bio
y_liquid_bio = y_liquid / sum_bio
y_solid_bio = y_solid / sum_bio
y_meta_bio = y_meta / sum_bio

# Print parameters and results
# ----------------------------------------------------------------------------

print(
    f'\n{" Parameters ":*^70}\n\n'
    f'diam     {diam} m\n'
    f'length   {length} m\n'
    f'temp     {temp} K\n'
    f'pabs     {pabs:,} Pa\n'
    f'n_cstrs  {n_cstrs}'
)

print(
    f'\n{" Calculated parameters ":*^70}\n\n'
    f'mf_bio   {mf_bio:.4g} kg/s\n'
    f'mf_n2    {mf_n2:.4g} kg/s'
)

print(
    f'\n{" Final mass fractions (N₂ basis) ":*^70}\n\n'
    f'gas      {y_gas[-1]:.3f}\n'
    f'liquid   {y_liquid[-1]:.3f}\n'
    f'solid    {y_solid[-1]:.3f}\n'
    f'meta     {y_meta[-1]:.3f}\n'
    f'sum      {y_gas[-1] + y_liquid[-1] + y_solid[-1] + y_meta[-1]:.3f}'
)

print(
    f'\n{" Final mass fractions (N₂ free basis) ":*^70}\n\n'
    f'gas      {y_gas_bio[-1]:.3f}\n'
    f'liquid   {y_liquid_bio[-1]:.3f}\n'
    f'solid    {y_solid_bio[-1]:.3f}\n'
    f'meta     {y_meta_bio[-1]:.3f}\n'
    f'sum      {y_gas_bio[-1] + y_liquid_bio[-1] + y_solid_bio[-1] + y_meta_bio[-1]:.3f}'
)


# Plot results
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
