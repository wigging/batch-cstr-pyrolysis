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
    yields, norm_yields = f.calc_yields()
    lumped = f.calc_lump_yields(norm_yields)

    print(
        f'\n{" " + f.name + ", Cycle " + str(f.cycle) + " ":*^70}\n\n'
        f'Yields as wt. % wet     exp     norm\n'
        f'oil yield          {yields[0]:>8} {norm_yields[0]:>8.2f}\n'
        f'condensables yield {yields[1]:>8} {norm_yields[1]:>8.2f}\n'
        f'light gas yield    {yields[2]:>8} {norm_yields[2]:>8.2f}\n'
        f'water vapor yield  {yields[3]:>8} {norm_yields[3]:>8.2f}\n'
        f'char yield         {yields[4]:>8} {norm_yields[4]:>8.2f}\n'
        f'ash                {yields[5]:>8} {norm_yields[5]:>8.2f}\n'
        f'total              {sum(yields):>8.2f} {sum(norm_yields):>8.2f}'
    )

    print(
        f'\nLumped yields as wt. %\n'
        f'gases     {lumped[0]:.2f}\n'
        f'liquids   {lumped[1]:.2f}\n'
        f'solids    {lumped[2]:.2f}\n'
        f'total     {sum(lumped):.2f}'
    )
