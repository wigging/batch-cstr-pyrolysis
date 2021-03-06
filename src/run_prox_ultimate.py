"""
Analyze proximate and ultimate analysis values from the FCIC feedstock data
spreadsheet and convert the as-determined basis (ad) to other bases. Compare
the max wt. % difference of as-determined values for each feedstock.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from feedstock import Feedstock

np.set_printoptions(precision=4, suppress=True)

# Feedstocks
# ----------------------------------------------------------------------------

with open("data/feedstocks.json") as json_file:
    fdata = json.load(json_file)

feedstocks = [Feedstock(f) for f in fdata]

# Proximate and ultimate analysis data
# ----------------------------------------------------------------------------

n = len(feedstocks)
proxs_ad = np.zeros((n, 4))
ults_ad = np.zeros((n, 7))

for i, feedstock in enumerate(feedstocks):
    proxs_ad[i] = feedstock.prox_ad
    ults_ad[i] = feedstock.ult_ad

    prox_ad = feedstock.prox_ad
    prox_ar = feedstock.prox_ar
    prox_d = feedstock.prox_d
    prox_daf = feedstock.prox_daf

    ult_ad = feedstock.ult_ad
    ult_ar = feedstock.ult_ar
    ult_d = feedstock.ult_d
    ult_daf = feedstock.ult_daf
    ult_cho = feedstock.ult_cho

    sum_ult_ad = sum(ult_ad) - ult_ad[6]  # exclude moisture content

    print(f'\n{" " + feedstock.name + ", Cycle " + str(feedstock.cycle) + " ":*^70}\n')

    print('Proximate analysis wt. %')
    print(f'{"ad":>13} {"ar":>10} {"d":>10}{"daf":>11}')
    print(f'FC {prox_ad[0]:>10.2f} {prox_ar[0]:>10.2f} {prox_d[0]:>10.2f} {prox_daf[0]:>10.2f}')
    print(f'VM {prox_ad[1]:>10.2f} {prox_ar[1]:>10.2f} {prox_d[1]:>10.2f} {prox_daf[1]:>10.2f}')
    print(f'ash {prox_ad[2]:>9.2f} {prox_ar[2]:>10.2f} {prox_d[2]:>10.2f}')
    print(f'M {prox_ad[3]:>11.2f} {prox_ar[3]:>10.2f}')
    print(f'Σ {sum(prox_ad):>11.2f} {sum(prox_ar):>10.2f} {sum(prox_d):>10.2f} {sum(prox_daf):>10.2f}')

    print('\nUltimate analysis wt. %')
    print('★ reported H and O for ad-basis excludes H and O in moisture')
    print(f'{"ad":>10} {"ar":>10} {"d":>10} {"daf":>10} {"cho":>10}')
    print(f'C {ult_ad[0]:>9.2f} {ult_ar[0]:>10.2f} {ult_d[0]:>10.2f} {ult_daf[0]:>10.2f} {ult_cho[0]:>10.2f}')
    print(f'H {ult_ad[1]:>9.2f} {ult_ar[1]:>10.2f} {ult_d[1]:>10.2f} {ult_daf[1]:>10.2f} {ult_cho[1]:>10.2f}')
    print(f'O {ult_ad[2]:>9.2f} {ult_ar[2]:>10.2f} {ult_d[2]:>10.2f} {ult_daf[2]:>10.2f} {ult_cho[2]:>10.2f}')
    print(f'N {ult_ad[3]:>9.2f} {ult_ar[3]:>10.2f} {ult_d[3]:>10.2f} {ult_daf[3]:>10.2f}')
    print(f'S {ult_ad[4]:>9.2f} {ult_ar[4]:>10.2f} {ult_d[4]:>10.2f} {ult_daf[4]:>10.2f}')
    print(f'ash {ult_ad[5]:>7.2f} {ult_ar[5]:>10.2f} {ult_d[5]:>10.2f}')
    print(f'M {ult_ad[6]:>9.2f}★ {ult_ar[6]:>9.2f}')
    print(f'Σ {sum_ult_ad:>9.2f} {sum(ult_ar):>10.2f} {sum(ult_d):>10.2f} {sum(ult_daf):>10.2f} {sum(ult_cho):>10.2f}')

# Max weight percent difference for FC, VM, ash, moisture for as-determined basis
wt_max = [max(col) - min(col) for col in proxs_ad.T]

print(f'\n{" Max wt. ％ difference for all FC, VM, ash, moisture (ad) ":*^70}\n')
print(f'FC       {wt_max[0]:.2f}')
print(f'VM       {wt_max[1]:.2f}')
print(f'ash      {wt_max[2]:.2f}')
print(f'moisture {wt_max[3]:.2f}')

# Max weight percent difference for C, H, O, N, S for as-determined basis
wt_max = [max(col) - min(col) for col in ults_ad.T]

print(f'\n{" Max wt. ％ difference for all C, H, O, N, S (ad) ":*^70}\n')
print(f'C  {wt_max[0]:.2f}')
print(f'H  {wt_max[1]:.2f}')
print(f'O  {wt_max[2]:.2f}')
print(f'N  {wt_max[3]:.2f}')
print(f'S  {wt_max[4]:.2f}')


# Plot
# ----------------------------------------------------------------------------

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


_, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.8), tight_layout=True)

labels = ['FC', 'VM', 'ash', 'moisture']

for p in proxs_ad:
    ax1.plot(p, 'o')
ax1.set_ylabel('Weight % (as-determined)')
ax1.set_xlabel('Proximate analysis')
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels)
style(ax1)

labels = ['C', 'H', 'O', 'N', 'S']

for u in ults_ad:
    ax2.plot(u[:5], 'o')
ax2.set_xlabel('Ultimate analysis')
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels)
style(ax2)

plt.show()
