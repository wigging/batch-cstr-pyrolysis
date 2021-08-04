"""
Experimental yields for all FCIC feedstocks.
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
    exp_norm = f.calc_norm_yield()

    print(
        f'\n{" " + f.name + ", Cycle " + str(f.cycle) + " ":*^70}\n'
        f'Given as wt. % wet      exp     norm\n'
        f'oil yield          {f.exp_yield[0]:>8} {exp_norm[0]:>8.2f}\n'
        f'condensables yield {f.exp_yield[1]:>8} {exp_norm[1]:>8.2f}\n'
        f'light gas yield    {f.exp_yield[2]:>8} {exp_norm[2]:>8.2f}\n'
        f'water vapor yield  {f.exp_yield[3]:>8} {exp_norm[3]:>8.2f}\n'
        f'char yield         {f.exp_yield[4]:>8} {exp_norm[4]:>8.2f}\n'
        f'total              {sum(f.exp_yield):>8.2f} {sum(exp_norm):>8.2f}'
    )
