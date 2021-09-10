"""
Experimental yield data and calculated lumped yields for all FCIC feedstocks.
Compare experiment yields to feedstock ash content.
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

for feedstock in feedstocks:
    exp_yield = feedstock.exp_yield
    normexp_yield = feedstock.normexp_yield
    lump_yield = feedstock.lump_yield
    normlump_yield = feedstock.normlump_yield

    ash.append(feedstock.prox_ad[2])
    char.append(exp_yield[4])
    oil.append(exp_yield[0])
    liquids.append(lump_yield[1])

    print(
        f'\n{" " + feedstock.name + ", Cycle " + str(feedstock.cycle) + " ":*^70}\n\n'
        f'Yields as wt. % wet     exp     norm\n'
        f'oil                {exp_yield[0]:>8} {normexp_yield[0]:>8.2f}\n'
        f'condensables       {exp_yield[1]:>8} {normexp_yield[1]:>8.2f}\n'
        f'light gas          {exp_yield[2]:>8} {normexp_yield[2]:>8.2f}\n'
        f'water vapor        {exp_yield[3]:>8} {normexp_yield[3]:>8.2f}\n'
        f'char               {exp_yield[4]:>8} {normexp_yield[4]:>8.2f}\n'
        f'total              {sum(exp_yield):>8.2f} {sum(normexp_yield):>8.2f}'
    )

    print(
        f'\nLump yields as wt. %    exp     norm\n'
        f'gases              {lump_yield[0]:>8.2f} {normlump_yield[0]:>8.2f}\n'
        f'liquids            {lump_yield[1]:>8.2f} {normlump_yield[1]:>8.2f}\n'
        f'char               {lump_yield[2]:>8.2f} {normlump_yield[2]:>8.2f}\n'
        f'total              {sum(lump_yield):>8.2f} {sum(normlump_yield):>8.2f}'
    )


# Plot
# ----------------------------------------------------------------------------

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


# Trendlines
z1 = np.polyfit(ash, oil, 1)
p1 = np.poly1d(z1)

z2 = np.polyfit(ash, char, 1)
p2 = np.poly1d(z2)

z3 = np.polyfit(ash, liquids, 1)
p3 = np.poly1d(z3)

# Yields vs ash content and show trendline
_, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(9, 4.8), sharey=True, tight_layout=True)

ax1.plot(ash, oil, 'o')
ax1.plot(ash, p1(ash))
ax1.set_ylabel('Experiment yield [wt. %]')
ax1.set_title('Oil')
style(ax1)

ax2.plot(ash, char, 'o')
ax2.plot(ash, p2(ash))
ax2.set_xlabel('Feedstock ash [wt. %]')
ax2.set_title('Char')
style(ax2)

ax3.plot(ash, liquids, 'o')
ax3.plot(ash, p3(ash))
ax3.set_title('Liquids')
style(ax3)

# Yields for each feedstock
labels = ['oil', 'condensables', 'light gas', 'water vapor', 'char']

_, ax = plt.subplots(tight_layout=True)
for f in feedstocks:
    ax.plot(f.exp_yield, 'o')
ax.set_ylabel('Weight % (wet basis)')
ax.set_xlabel('Experiment yield')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)
style(ax)

plt.show()
