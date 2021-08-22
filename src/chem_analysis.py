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

# Feedstock
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Chemical analysis
# ----------------------------------------------------------------------------

n = len(feedstocks)
chems_d = np.zeros((n, 12))

for i, f in enumerate(feedstocks):
    chems_d[i] = f.chem_d
    chem_daf = f.chem_daf
    biocomp = f.chem_bc

    total_d = sum(f.chem_d)
    total_daf = sum(chem_daf)
    total_bio = sum(biocomp)

    print(f'\n{" " + f.name + ", Cycle " + str(f.cycle) + " ":*^70}\n')
    print(
        'Chemical analysis wt. %        d      daf \n'
        f'structural inorganics     {f.chem_d[0]:>8} \n'
        f'non-structural inorganics {f.chem_d[1]:>8} \n'
        f'water extractives         {f.chem_d[2]:>8}{chem_daf[2]:>8.2f} \n'
        f'ethanol extractives       {f.chem_d[3]:>8}{chem_daf[3]:>8.2f} \n'
        f'acetone extractives       {f.chem_d[4]:>8}{chem_daf[4]:>8.2f} \n'
        f'lignin                    {f.chem_d[5]:>8}{chem_daf[5]:>8.2f} \n'
        f'glucan                    {f.chem_d[6]:>8}{chem_daf[6]:>8.2f} \n'
        f'xylan                     {f.chem_d[7]:>8}{chem_daf[7]:>8.2f} \n'
        f'galactan                  {f.chem_d[8]:>8}{chem_daf[8]:>8.2f} \n'
        f'arabinan                  {f.chem_d[9]:>8}{chem_daf[9]:>8.2f} \n'
        f'mannan                    {f.chem_d[10]:>8}{chem_daf[10]:>8.2f} \n'
        f'acetyl                    {f.chem_d[11]:>8}{chem_daf[11]:>8.2f} \n'
        f'total                     {total_d:>8.2f}{total_daf:>8.2f}'
    )

    print(
        '\nBiomass composition wt. %     daf \n'
        f'cellulose                 {biocomp[0]:>8.2f} \n'
        f'hemicellulose             {biocomp[1]:>8.2f} \n'
        f'lignin                    {biocomp[2]:>8.2f} \n'
        f'total                     {total_bio:>8.2f}'
    )

# Max weight percent difference for dry basis
wt_max = [max(col) - min(col) for col in chems_d[:-2].T]

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


# Plot
# ----------------------------------------------------------------------------

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


labels = [
    'struct. inorg.', 'non-struct. inorg.', 'water ext.', 'ethanol ext.', 'acetone ext.',
    'lignin', 'glucan', 'xylan', 'galactan', 'arabinan', 'mannan', 'acetyl'
]

_, ax = plt.subplots(tight_layout=True)

for c in chems_d:
    ax.plot(c, 'o')
ax.set_ylabel('Weight % (as-determined)')
ax.set_xlabel('Chemical analysis')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)
ax.tick_params(axis='x', rotation=90)
style(ax)

plt.show()
