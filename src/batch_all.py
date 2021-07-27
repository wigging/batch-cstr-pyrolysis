"""
Run a batch reactor model for every FCIC feedstock using the Debiagi 2018
biomass pyrolysis kinetics. All reactions and chemical species are defined in
Cantera `cti` files where `debiagi_sw.cti` is for softwood reactions,
`debiagi_hw.cti` is for hardwood reactions, and `debiagi_gr.cti` is for grass
reactions. Biomass compositions for each feedstock were determined from the
`biocomp.py` file.

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

# Feedstock biomass compositions

feedstocks = [
    {'name': 'residues', 'biocomp': [28.46, 22.28, 0.14, 11.6, 25.09, 2.01, 10.43]},
    {'name': 'stem wood', 'biocomp': [39.48, 25.28, 0.0, 25.89, 4.37, 0.35, 4.63]},
    {'name': 'bark', 'biocomp': [33.83, 18.38, 31.90, 0.23, 2.52, 9.09, 4.05]},
]

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

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas = ct.Solution('data/debiagi_sw.cti')

n = len(feedstocks)
names = []
yf_gases = np.zeros(n)
yf_liquids = np.zeros(n)
yf_solids = np.zeros(n)
yf_metas = np.zeros(n)

print('\nLumped species final yields as mass fraction, dry ash-free basis\n')
print('             Gases    Liquids  Solids   Metaplastics')

for i, feed in enumerate(feedstocks):

    bc = feed['biocomp']
    y0 = np.array(bc) / 100

    y00 = f'CELL:{y0[0]} GMSW:{y0[1]} LIGC:{y0[2]} LIGH:{y0[3]} LIGO:{y0[4]} TANN:{y0[5]} TGL:{y0[6]}'
    gas.TPY = tk, p, y00
    r = ct.IdealGasReactor(gas, energy='off')

    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])

    for t in time:
        sim.advance(t)
        states.append(r.thermo.state, t=t)

    y_gases = states(*sp_gases).Y.sum(axis=1)
    y_liquids = states(*sp_liquids).Y.sum(axis=1)
    y_solids = states(*sp_solids).Y.sum(axis=1)
    y_metaplastics = states(*sp_metaplastics).Y.sum(axis=1)

    names.append(feed['name'])
    yf_gases[i] = y_gases[-1]
    yf_liquids[i] = y_liquids[-1]
    yf_solids[i] = y_solids[-1]
    yf_metas[i] = y_metaplastics[-1]

    name = feed['name']
    print(f'{name:10} {y_gases[-1]:8.4f} {y_liquids[-1]:8.4f} {y_solids[-1]:8.4f} {y_metaplastics[-1]:8.4f}')


# Plot
# ----------------------------------------------------------------------------

def style_barh(ax):
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.set_frame_on(False)
    ax.tick_params(color='0.8')
    ax.xaxis.grid(True, color='0.8')


y = np.arange(n)

_, ax = plt.subplots(tight_layout=True)
ax.barh(y, yf_gases, color='C6')
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Mass fraction [-]')
ax.set_title('Gases')
style_barh(ax)

_, ax = plt.subplots(tight_layout=True)
ax.barh(y, yf_liquids, color='C4')
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Mass fraction [-]')
ax.set_title('Liquids')
style_barh(ax)

_, ax = plt.subplots(tight_layout=True)
ax.barh(y, yf_solids, color='C2')
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Mass fraction [-]')
ax.set_title('Solids')
style_barh(ax)

_, ax = plt.subplots(tight_layout=True)
ax.barh(y, yf_metas, color='C1')
ax.set_yticks(y)
ax.set_yticklabels(names)
ax.set_xlabel('Mass fraction [-]')
ax.set_title('Metaplastics')
style_barh(ax)

plt.show()
