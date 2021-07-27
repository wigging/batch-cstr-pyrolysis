"""
Run a batch reactor model for a single FCIC feedstock using the Debiagi 2018
biomass pyrolysis kinetics. All reactions and chemical species are defined in
Cantera `cti` files where `debiagi_sw.cti` is for softwood reactions,
`debiagi_hw.cti` is for hardwood reactions, and `debiagi_gr.cti` is for grass
reactions.

Note
----
The `debiagi_sw.cti` uses the original reaction rates from Debiagi 2018 while
`debiagi_sw2.cti` uses modified rates for the metaplastic reactions such that
b=1 instead of b=0 where b is the temperature exponent Tᵇ.

References
----------
P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, and T. Faravelli.
A predictive model of biochar formation and characterization. Journal of
Analytical and Applied Pyrolysis, vol. 134, pp. 326-335, 2018.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

ct.suppress_thermo_warnings()

# Parameters
# ----------------------------------------------------------------------------

tk = 773.15                         # reactor temperature [K]
p = 101325.0                        # reactor pressure [Pa]
time = np.linspace(0, 10.0, 100)    # reaction time steps [s]

# Residues feedstock

label = 'Residues, Cycle 1'
y0 = 'CELL:0.2846 GMSW:0.2228 LIGC:0.0014 LIGH:0.1160 LIGO:0.2509 TANN:0.0201 TGL:0.1043'
gas = ct.Solution('data/debiagi_sw.cti')

# Stem wood feedstock

# label = 'Stem wood, Cycle 2'
# y0 = 'CELL:0.3948 GMSW:0.2528 LIGC:0.0 LIGH:0.2589 LIGO:0.0437 TANN:0.0035 TGL:0.0463'
# gas = ct.Solution('data/debiagi_sw2.cti')

# Bark feedstock

# label = 'Bark, Cycle 3'
# y0 = 'CELL:0.3383 GMSW:0.1838 LIGC:0.3190 LIGH:0.0023 LIGO:0.0252 TANN:0.0909 TGL:0.0405'
# gas = ct.Solution('data/debiagi_sw2.cti')

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas.TPY = tk, p, y0
r = ct.IdealGasReactor(gas, energy='off')

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
    f'\n{" Batch reactor model ":*^60}\n\n'
    f'Feedstock     {label}\n'
    f'Final time    {time[-1]} s\n'
    f'Temperature   {tk} K\n'
    f'Pressure      {p:,} Pa'
)

print('\nChemical species final yields as mass fraction, dry ash-free basis')
print(f'Number of species = {len(states.species_names)}')
for sp in states.species_names:
    print(f"{sp:11} {states(sp).Y[:, 0][-1]:.4f}")

print(
    '\nLumped species final yields as mass fraction, dry ash-free basis\n'
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
