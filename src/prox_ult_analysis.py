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

# Get feedstock data from JSON file and create Feedstock objects

with open("feedstocks.json") as json_file:
    fdata = json.load(json_file)

feeds = [Feedstock(f) for f in fdata]

# Compare proximate and ultimate analysis values for each feedstock

n = len(feeds)
proxs = np.zeros((n, 4))
ults = np.zeros((n, 7))

for i, fd in enumerate(feeds):
    proxs[i] = fd.prox
    ults[i] = fd.ult

    prox = fd.prox
    prox_ar = fd.calc_prox_ar()
    prox_d = fd.calc_prox_dry()
    prox_daf = fd.calc_prox_daf()

    ult = fd.ult
    ult_ar = fd.calc_ult_ar()
    ult_d = fd.calc_ult_dry()
    ult_daf = fd.calc_ult_daf()
    ult_cho = fd.calc_ult_cho(ult_daf)

    tot_ad = sum(ult[:-1])
    tot_ar = sum(ult_ar)
    tot_d = sum(ult_d)
    tot_daf = sum(ult_daf)
    tot_cho = sum(ult_cho)

    print(f'\n{" " + fd.name + ", Cycle " + str(fd.cycle) + " ":*^70}\n')

    print('Proximate analysis wt. %')
    print(f'{"ad":>13} {"ar":>10} {"d":>10}{"daf":>11}')
    print(f'FC {prox[0]:>10.2f} {prox_ar[0]:>10.2f} {prox_d[0]:>10.2f} {prox_daf[0]:>10.2f}')
    print(f'VM {prox[1]:>10.2f} {prox_ar[1]:>10.2f} {prox_d[1]:>10.2f} {prox_daf[1]:>10.2f}')
    print(f'ash {prox[2]:>9.2f} {prox_ar[2]:>10.2f} {prox_d[2]:>10.2f}')
    print(f'M {prox[3]:>11.2f} {prox_ar[3]:>10.2f}')
    print(f'Σ {sum(prox):>11.2f} {sum(prox_ar):>10.2f} {sum(prox_d):>10.2f} {sum(prox_daf):>10.2f}')

    print('\nUltimate analysis wt. %')
    print('reported H and O for ad-basis excludes H and O in moisture')
    print(f'{"ad":>11} {"ar":>10} {"d":>10} {"daf":>10} {"cho":>10}')
    print(f'C {ult[0]:>9.2f} {ult_ar[0]:>10.2f} {ult_d[0]:>10.2f} {ult_daf[0]:>10.2f} {ult_cho[0]:>10.2f}')
    print(f'H {ult[1]:>9.2f} {ult_ar[1]:>10.2f} {ult_d[1]:>10.2f} {ult_daf[1]:>10.2f} {ult_cho[1]:>10.2f}')
    print(f'O {ult[2]:>9.2f} {ult_ar[2]:>10.2f} {ult_d[2]:>10.2f} {ult_daf[2]:>10.2f} {ult_cho[2]:>10.2f}')
    print(f'N {ult[3]:>9.2f} {ult_ar[3]:>10.2f} {ult_d[3]:>10.2f} {ult_daf[3]:>10.2f}')
    print(f'S {ult[4]:>9.2f} {ult_ar[4]:>10.2f} {ult_d[4]:>10.2f} {ult_daf[4]:>10.2f}')
    print(f'ash {ult[5]:>7.2f} {ult_ar[5]:>10.2f} {ult_d[5]:>10.2f}')
    print(f'M {ult[6]:>9.2f} {ult_ar[6]:>10.2f}')
    print(f'Σ {tot_ad:>9.2f} {tot_ar:>10.2f} {tot_d:>10.2f} {tot_daf:>10.2f} {tot_cho:>10.2f}')

# Max weight percent difference for FC, VM, ash, moisture for as-determined basis

wt_max = [max(col) - min(col) for col in proxs.T]

print(f'\n{" Max wt. ％ difference for all FC, VM, ash, moisture (ad) ":*^70}\n')
print(f'FC       {wt_max[0]:.2f}')
print(f'VM       {wt_max[1]:.2f}')
print(f'ash      {wt_max[2]:.2f}')
print(f'moisture {wt_max[3]:.2f}')

# Max weight percent difference for C, H, O, N, S for as-determined basis

wt_max = [max(col) - min(col) for col in ults.T]

print(f'\n{" Max wt. ％ difference for all C, H, O, N, S (ad) ":*^70}\n')
print(f'C  {wt_max[0]:.2f}')
print(f'H  {wt_max[1]:.2f}')
print(f'O  {wt_max[2]:.2f}')
print(f'N  {wt_max[3]:.2f}')
print(f'S  {wt_max[4]:.2f}')


# Plot comparison of proximate and ultiamte analysis values

def style(ax):
    ax.grid(color='0.9')
    ax.set_frame_on(False)
    ax.tick_params(color='0.9')


_, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.8), tight_layout=True)

labels = ['FC', 'VM', 'ash', 'moisture']

for p in proxs:
    ax1.plot(p, 'o')
ax1.set_ylabel('Weight % (as-determined)')
ax1.set_xlabel('Proximate analysis')
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels)
style(ax1)

labels = ['C', 'H', 'O', 'N', 'S']

for u in ults:
    ax2.plot(u[:5], 'o')
ax2.set_xlabel('Ultimate analysis')
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels)
style(ax2)

plt.show()
