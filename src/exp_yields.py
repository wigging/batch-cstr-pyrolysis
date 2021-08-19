"""
Experimental and lumped yields for all FCIC feedstocks.
"""

import json
from feedstock import Feedstock

# Feedstock parameters
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Experimental yields
# ----------------------------------------------------------------------------

for f in feedstocks:
    exp_yields, norm_yields = f.calc_yields()
    exp_lumps, norm_lumps = f.calc_lump_yields(exp_yields, norm_yields)

    print(
        f'\n{" " + f.name + ", Cycle " + str(f.cycle) + " ":*^70}\n\n'
        f'Yields as wt. % wet     exp     norm\n'
        f'oil                {exp_yields[0]:>8} {norm_yields[0]:>8.2f}\n'
        f'condensables       {exp_yields[1]:>8} {norm_yields[1]:>8.2f}\n'
        f'light gas          {exp_yields[2]:>8} {norm_yields[2]:>8.2f}\n'
        f'water vapor        {exp_yields[3]:>8} {norm_yields[3]:>8.2f}\n'
        f'char               {exp_yields[4]:>8} {norm_yields[4]:>8.2f}\n'
        f'total              {sum(exp_yields):>8.2f} {sum(norm_yields):>8.2f}'
    )

    print(
        f'\nLump yields as wt. %    exp     norm\n'
        f'gases              {exp_lumps[0]:>8.2f} {norm_lumps[0]:>8.2f}\n'
        f'liquids            {exp_lumps[1]:>8.2f} {norm_lumps[1]:>8.2f}\n'
        f'char               {exp_lumps[2]:>8.2f} {norm_lumps[2]:>8.2f}\n'
        f'total              {sum(exp_lumps):>8.2f} {sum(norm_lumps):>8.2f}'
    )
