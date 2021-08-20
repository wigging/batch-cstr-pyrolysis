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

import cantera as ct
import chemics as cm
import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock
from objfunc import objfunc
from scipy.optimize import minimize

ct.suppress_thermo_warnings()

# Parameters
# ----------------------------------------------------------------------------

temp = 773.15                       # reactor temperature [K]
p = 101325.0                        # reactor pressure [Pa]
time = np.linspace(0, 20.0, 100)    # reaction time steps [s]
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
yf_gases = np.zeros(n)
yf_liquids = np.zeros(n)
yf_solids = np.zeros(n)
yf_metas = np.zeros(n)
exp_gases = np.zeros(n)
exp_liquids = np.zeros(n)
exp_solids = np.zeros(n)
exp_ash = np.zeros(n)

# Run batch reactor model for each feedstock
for i, f in enumerate(feedstocks):

    # Get feedstock name, ultimate analysis data, and chemical analysis data
    ult_daf = f.calc_ult_daf()
    ult_cho = f.calc_ult_cho(ult_daf)
    chem_daf = f.calc_chem_daf()
    chem_bc = f.calc_chem_bc(chem_daf) / 100

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

    # Perform batch reactor simulation
    y0 = f'CELL:{cell} GMSW:{hemi} LIGC:{ligc} LIGH:{ligh} LIGO:{ligo} TANN:{tann} TGL:{tgl}'

    gas = ct.Solution(cti)
    gas.TPY = temp, p, y0

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

    # Store results for feedstock
    names.append(f.name)
    yf_gases[i] = y_gases[-1]
    yf_liquids[i] = y_liquids[-1]
    yf_solids[i] = y_solids[-1]
    yf_metas[i] = y_metaplastics[-1]

    # Experiment lumped yields
    exp_yields, norm_yields = f.calc_yields()
    exp_lumps, _ = f.calc_lump_yields(exp_yields, norm_yields)
    exp_gases[i] = exp_lumps[0]
    exp_liquids[i] = exp_lumps[1]
    exp_solids[i] = exp_lumps[2]
    exp_ash[i] = f.prox[2]

# Print
# ----------------------------------------------------------------------------

print(
    f'\n{" Batch reactor model ":*^70}\n'
    '\nParameters for reactor model\n'
    f'final time    {time[-1]} s\n'
    f'temperature   {temp} K\n'
    f'pressure      {p:,} Pa\n'
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

z = np.polyfit(exp_ash, yf_liquids, 1)
p = np.poly1d(z)

_, ax = plt.subplots()
ax.plot(exp_ash, yf_liquids, 'o')
ax.plot(exp_ash, p(exp_ash))
ax.set_xlabel('Ash [wt. %]')
ax.set_ylabel('Liquids [wt. %]')

plt.show()
