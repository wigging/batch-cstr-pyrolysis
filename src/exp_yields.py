"""
Experimental and lumped yields for all FCIC feedstocks.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock

# Feedstock parameters
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(fd) for fd in fdata]

# Experimental yields
# ----------------------------------------------------------------------------

ash = []
char = []
oil = []
liquids = []

for f in feedstocks:
    exp_yields, norm_yields = f.calc_yields()
    exp_lumps, norm_lumps = f.calc_lump_yields(exp_yields, norm_yields)

    ash.append(f.prox[2])
    char.append(exp_yields[4])
    oil.append(exp_yields[0])
    liquids.append(exp_lumps[1])

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


# Plot
# ----------------------------------------------------------------------------

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


z1 = np.polyfit(ash, oil, 1)
p1 = np.poly1d(z1)

z2 = np.polyfit(ash, char, 1)
p2 = np.poly1d(z2)

z3 = np.polyfit(ash, liquids, 1)
p3 = np.poly1d(z3)

_, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(9, 4.8), sharey=True, tight_layout=True)

ax1.plot(ash, oil, 'o')
ax1.plot(ash, p1(ash))
ax1.set_ylabel('Experiment yield [wt. %]')
ax1.set_title('Oil')
style(ax1)

ax2.plot(ash, char, 'o')
ax2.plot(ash, p2(ash))
ax2.set_xlabel('Experiment ash [wt. %]')
ax2.set_title('Char')
style(ax2)

ax3.plot(ash, liquids, 'o')
ax3.plot(ash, p3(ash))
ax3.set_title('Liquids')
style(ax3)

plt.show()
