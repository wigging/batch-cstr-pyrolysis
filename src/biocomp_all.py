"""
Determine the biomass composition for all FCIC feedstocks using ultimate
analysis data and chemical analysis data. The splitting parameter values are
obtained using an optimization function and comparing the results to the
measured chemical analysis data.
"""

import chemics as cm
import json
from feedstock import Feedstock
from objfunc import objfunc
from scipy.optimize import minimize

# Feedstock parameters
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Biomass composition
# ----------------------------------------------------------------------------

for f in feedstocks:

    # Get feedstock name, ultimate analysis data, and chemical analysis data
    name = f.name
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

    # Print results
    print(
        f'\n{" " + name + ", Cycle " + str(f.cycle) + " ":*^70}\n'
        '\nUltimate analysis, wt. ％ CHO basis\n'
        f'C         {yc * 100:.2f}\n'
        f'H         {yh * 100:.2f}\n'
        '\nOptimized splitting parameters\n'
        f'α         {res.x[0]:.4f}\n'
        f'β         {res.x[1]:.4f}\n'
        f'γ         {res.x[2]:.4f}\n'
        f'δ         {res.x[3]:.4f}\n'
        f'ε         {res.x[4]:.4f}\n'
        '\nChemical analysis, wt. ％ daf\n'
        '           exp      opt\n'
        f'cell      {chem_bc[0] * 100:.2f}    {cell * 100:.2f}\n'
        f'hemi      {chem_bc[1] * 100:.2f}    {hemi * 100:.2f}\n'
        f'lignin    {chem_bc[2] * 100:.2f}    {(ligc + ligh + ligo) * 100:.2f}\n'
        '\nBiomass composition, wt. ％ daf\n'
        f'cell      {cell * 100:.2f}\n'
        f'hemi      {hemi * 100:.2f}\n'
        f'lig-c     {ligc * 100:.2f}\n'
        f'lig-h     {ligh * 100:.2f}\n'
        f'lig-o     {ligo * 100:.2f}\n'
        f'tann      {tann * 100:.2f}\n'
        f'tgl       {tgl * 100:.2f}'
    )
