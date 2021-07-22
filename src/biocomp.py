"""
Determine the biomass composition for a given FCIC feedstock using ultimate
analysis data. Splitting parameters should be adjusted such that cellulose,
hemicellulose, and total lignin are similar to the chemical analysis data.
"""

import matplotlib.pyplot as plt
import chemics as cm

# Feedstock parameters
# Default splitting parameters are α = 0.6, β = 0.8, γ = 0.8, δ = 1, ε = 1

# name = 'Residues'
# yc = 0.5331
# yh = 0.0641
# alpha = 0.51
# beta = 0.98
# gamma = 1
# delta = 0.7
# epsilon = 0.9

# name = 'Stem wood'
# yc = 0.5094
# yh = 0.0639
# alpha = 0.56
# beta = 1
# gamma = 1
# delta = 0.92
# epsilon = 0.9

name = 'Bark'
yc = 0.5569
yh = 0.0589
alpha = 0.6
beta = 0.05
gamma = 0.05
delta = 0.7
epsilon = 0.8

# Calculate biomass composition

bc = cm.biocomp(yc, yh, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon)

cell = bc['y_daf'][0]
hemi = bc['y_daf'][1]
ligc = bc['y_daf'][2]
ligh = bc['y_daf'][3]
ligo = bc['y_daf'][4]
tann = bc['y_daf'][5]
tgl = bc['y_daf'][6]
total = cell + hemi + ligc + ligh + ligo + tann + tgl
total_lig = sum(bc['y_daf'][2:5])

# Print parameters and biomass composition as mass fraction, dry ash-free basis (daf)

print(
    f'\nParameters for {name} \n'
    f'C = {yc * 100} \n'
    f'H = {yh * 100} \n'
    f'α = {alpha} \n'
    f'β = {beta} \n'
    f'γ = {gamma} \n'
    f'δ = {delta} \n'
    f'ε = {epsilon} \n'
)

print(
    'Biomass composition    daf \n'
    f'cellulose             {cell * 100:8>.2f} \n'
    f'hemicellulose         {hemi * 100:8>.2f} \n'
    f'lignin-c              {ligc * 100:8>.2f} \n'
    f'lignin-h              {ligh * 100:8>.2f} \n'
    f'lignin-o              {ligo * 100:8>.2f} \n'
    f'tannins               {tann * 100:8>.2f} \n'
    f'triglycerides         {tgl * 100:8>.2f} \n'
    f'total                 {total * 100:8>.2f} \n'
    f'total lignin          {total_lig * 100:8>.2f}'
)

# Plot the reference mixtures

fig, ax = plt.subplots(tight_layout=True)
cm.plot_biocomp(ax, yc, yh, bc['y_rm1'], bc['y_rm2'], bc['y_rm3'])
ax.set_title(name, y=1, pad=-16)

plt.show()
