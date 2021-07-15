"""
Analyze the chemical analysis values from the FCIC feedstock data spreadsheet
and convert the data to other bases. Compare the max wt. % difference of
as-determined values for each feedstock
"""

# FIX: need to get chemical analysis data for cycle 15 and cycle 16, last two feedstocks

import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock

np.set_printoptions(precision=4, suppress=True)

# Get feedstock data from JSON file and create Feedstock objects

with open("feedstocks.json") as json_file:
    fdata = json.load(json_file)

feeds = [Feedstock(f) for f in fdata]

# Chemical analysis values for each feedstock

n = len(feeds)
chems = np.zeros((n, 12))

for i, fd in enumerate(feeds):
    chems[i] = fd.chem
    chem_daf = fd.calc_chem_daf()
    biocomp = fd.calc_bio_comp(chem_daf)

    tot_d = sum(fd.chem)
    tot_daf = sum(chem_daf)
    tot_bio = sum(biocomp)

    print(f'\n{" " + fd.name + ", Cycle " + str(fd.cycle) + " ":*^70}\n')
    print(
        'Chemical analysis wt. %        d      daf \n'
        f'structural inorganics     {fd.chem[0]:>8} \n'
        f'non-structural inorganics {fd.chem[1]:>8} \n'
        f'water extractives         {fd.chem[2]:>8}{chem_daf[2]:>8.2f} \n'
        f'ethanol extractives       {fd.chem[3]:>8}{chem_daf[3]:>8.2f} \n'
        f'acetone extractives       {fd.chem[4]:>8}{chem_daf[4]:>8.2f} \n'
        f'lignin                    {fd.chem[5]:>8}{chem_daf[5]:>8.2f} \n'
        f'glucan                    {fd.chem[6]:>8}{chem_daf[6]:>8.2f} \n'
        f'xylan                     {fd.chem[7]:>8}{chem_daf[7]:>8.2f} \n'
        f'galactan                  {fd.chem[8]:>8}{chem_daf[8]:>8.2f} \n'
        f'arabinan                  {fd.chem[9]:>8}{chem_daf[9]:>8.2f} \n'
        f'mannan                    {fd.chem[10]:>8}{chem_daf[10]:>8.2f} \n'
        f'acetyl                    {fd.chem[11]:>8}{chem_daf[11]:>8.2f} \n'
        f'total                     {tot_d:>8.2f}{tot_daf:>8.2f}'
    )

    print(
        '\nBiomass composition wt. %     daf \n'
        f'cellulose                 {biocomp[0]:>8.2f} \n'
        f'hemicellulose             {biocomp[1]:>8.2f} \n'
        f'lignin                    {biocomp[2]:>8.2f} \n'
        f'total                     {tot_bio:>8.2f}'
    )

# Max weight percent difference for dry basis

# FIX: need to get chemical analysis data for cycle 15 and cycle 16, last two feedstocks

wt_max = [max(col) - min(col) for col in chems[:-2].T]

print(f'\n{" Max wt. ％ difference (d) ":*^70}\n')
print(
    f'structural inorganics      {wt_max[0]:.2f} \n'
    f'non-structural inorganics  {wt_max[1]:.2f} \n'
    f'water extractives          {wt_max[2]:.2f} \n'
    f'ethanol extractives        {wt_max[3]:.2f} \n'
    f'acetone extractives        {wt_max[4]:.2f} \n'
    f'lignin                     {wt_max[5]:.2f} \n'
    f'glucan                     {wt_max[6]:.2f} \n'
    f'xylan                      {wt_max[7]:.2f} \n'
    f'galactan                   {wt_max[8]:.2f} \n'
    f'arabinan                   {wt_max[9]:.2f} \n'
    f'mannan                     {wt_max[10]:.2f} \n'
    f'acetyl                     {wt_max[11]:.2f} \n'
)


# Plot comparison of chemical analysis values

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


labels = [
    'struct. inorg.', 'non-struct. inorg.', 'water ext.', 'ethanol ext.', 'acetone ext.',
    'lignin', 'glucan', 'xylan', 'galactan', 'arabinan', 'mannan', 'acetyl'
]

_, ax = plt.subplots(tight_layout=True)

for c in chems:
    ax.plot(c, 'o')
ax.set_ylabel('Weight % (as-determined)')
ax.set_xlabel('Chemical analysis')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)
ax.tick_params(axis='x', rotation=90)
style(ax)

plt.show()
