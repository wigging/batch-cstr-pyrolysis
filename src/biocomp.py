import matplotlib.pyplot as plt
import chemics as cm

# Carbon and hydrogen mass fractions from ultimate analysis
yc = 0.5331
yh = 0.0641

# Calculate biomass composition and print results
bc = cm.biocomp(yc, yh, alpha=0.45, delta=0.9, epsilon=0.8, printcomp=True)

# Print composition as mass fraction, dry ash-free basis (daf)
cell = bc['y_daf'][0]
hemi = bc['y_daf'][1]
lignin = sum(bc['y_daf'][2:5])
tann = bc['y_daf'][5]
tgl = bc['y_daf'][6]
total = cell + hemi + lignin + tann + tgl

print(
    '\nComponent daf \n'
    f'cell     {cell:8>.4f} \n'
    f'hemi     {hemi:8>.4f} \n'
    f'lignin   {lignin:8>.4f} \n'
    f'tann     {tann:8>.4f} \n'
    f'tgl      {tgl:8>.4f} \n'
    f'total    {total:8>.2f}'
)

# Plot the reference mixtures
fig, ax = plt.subplots(tight_layout=True)
cm.plot_biocomp(ax, yc, yh, bc['y_rm1'], bc['y_rm2'], bc['y_rm3'])

plt.show()
