"""
Determine the biomass composition for all FCIC feedstocks using ultimate
analysis data and chemical analysis data. The splitting parameter values are
obtained using an optimization function which compares the results to the
measured chemical analysis data.
"""

import json
from feedstock import Feedstock

# Feedstock
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Biomass composition
# ----------------------------------------------------------------------------

for feedstock in feedstocks:

    # Calculate optimized biomass composition (daf) and splitting parameters
    bc, splits = feedstock.calc_biocomp()
    cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

    # Print results
    print(
        f'\n{" " + feedstock.name + ", Cycle " + str(feedstock.cycle) + " ":*^70}\n'
        '\nOptimized splitting parameters\n'
        f'α         {splits[0]:.4f}\n'
        f'β         {splits[1]:.4f}\n'
        f'γ         {splits[2]:.4f}\n'
        f'δ         {splits[3]:.4f}\n'
        f'ε         {splits[4]:.4f}\n'
        '\nChemical analysis, wt. ％ daf\n'
        '           exp      opt\n'
        f'cell      {feedstock.chem_bc[0]:.2f}    {cell * 100:.2f}\n'
        f'hemi      {feedstock.chem_bc[1]:.2f}    {hemi * 100:.2f}\n'
        f'lignin    {feedstock.chem_bc[2]:.2f}    {(ligc + ligh + ligo) * 100:.2f}\n'
        '\nBiomass composition, wt. ％ daf\n'
        f'cell      {cell * 100:.2f}\n'
        f'hemi      {hemi * 100:.2f}\n'
        f'lig-c     {ligc * 100:.2f}\n'
        f'lig-h     {ligh * 100:.2f}\n'
        f'lig-o     {ligo * 100:.2f}\n'
        f'tann      {tann * 100:.2f}\n'
        f'tgl       {tgl * 100:.2f}'
    )
