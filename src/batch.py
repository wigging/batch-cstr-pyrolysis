"""
Batch reactor example with Debiagi 2018 biomass pyrolysis kinetics.

This example uses Cantera to calculate conversion in a batch reactor
environment. All reactions and species are defined in the `cti` file where
`debiagi_sw.cti` is for softwood, `debiagi_hw.cti` is for hardwood, and
finally `debiagi_gr.cti` is for grass. Examples of initial mass fractions for
types of biomass such as softwood, hardwood, and grass are given below. Don't
forget to uncomment the relevant lines depending on the type of biomass.

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

tk = 773.15     # reactor temperature [K]
p = 101325.0    # reactor pressure [Pa]

# Residues feedstock
# label = 'Residues'
# y0 = 'CELL:0.2799 GMSW:0.2280 LIGC:0.0014 LIGH:0.1160 LIGO:0.2504 TANN:0.02 TGL:0.1043'
# gas = ct.Solution('data/debiagi_sw.cti')

# Stem wood feedstock
label = 'Stem wood'
y0 = 'CELL:0.3754 GMSW:0.2713 LIGC:0.0 LIGH:0.3167 LIGO:0.0023 TANN:0.0 TGL:0.0343'
gas = ct.Solution('data/debiagi_sw2.cti')

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 10.0, 100)

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas.TPY = tk, p, y0
r = ct.IdealGasReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

for t in time:
    sim.advance(t)
    states.append(r.thermo.state, t=t)

# species representing gases
sp_gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2')

# species representing liquids (tars)
sp_liquids = (
    'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
    'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H',
    'CH3OH', 'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12',
    'VANILLIN', 'ACQUA'
)

# species representing solids
sp_solids = (
    'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
    'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR'
)

# species representing metaplastics
sp_metaplastics = (
    'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
    'GH2', 'GC2H6'
)

# sum of species mass fractions for gases, liquids, solids, metaplastics
y_gases = states(*sp_gases).Y.sum(axis=1)
y_liquids = states(*sp_liquids).Y.sum(axis=1)
y_solids = states(*sp_solids).Y.sum(axis=1)
y_metaplastics = states(*sp_metaplastics).Y.sum(axis=1)

# Print
# ----------------------------------------------------------------------------

print(f'--- {label} Final Yields ---')
print(f'Number of species = {len(states.species_names)}')
for sp in states.species_names:
    print(f"{sp:11} {states(sp).Y[:, 0][-1]:.4f}")

print(
    f'\n{" " + label + " lumped yields ":*^70}\n\n'
    f'gases         {y_gases[-1] * 100:.2f}\n'
    f'liquids       {y_liquids[-1] * 100:.2f}\n'
    f'solids        {y_solids[-1] * 100:.2f}\n'
    f'metaplastics  {y_metaplastics[-1] * 100:.2f}\n'
)


# Plot
# ----------------------------------------------------------------------------

def style(ax, xlabel, ylabel):
    ax.grid(True, color='0.9')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.set_frame_on(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(color='0.9')


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
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELLA').Y[:, 0], label='CELLA')
ax.plot(states.t, states('HCE1').Y[:, 0], label='HCE1')
ax.plot(states.t, states('HCE2').Y[:, 0], label='HCE2')
ax.plot(states.t, states('LIGCC').Y[:, 0], label='LIGCC')
ax.plot(states.t, states('LIGOH').Y[:, 0], label='LIGOH')
ax.plot(states.t, states('LIG').Y[:, 0], label='LIG')
style(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.8), tight_layout=True)
ax1.plot(states.t, y_gases, label='gases')
ax1.plot(states.t, y_liquids, label='liquids')
ax1.plot(states.t, y_solids, label='solids')
ax1.plot(states.t, y_metaplastics, label='metaplastics')
style(ax1, xlabel='Time [s]', ylabel='Mass fraction [-]')

ax2.plot(states.t, states.T, color='m')
style(ax2, xlabel='Time [s]', ylabel='Temperature [K]')

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

plt.show()
