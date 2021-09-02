"""
Compare feedstock conversion time and tar yields using the Papadikis 2010
kinetics which includes primary and secondary reactions.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

temp = 773.15                       # reactor temperature [K]
p = 101325.0                        # reactor pressure [Pa]
time = np.linspace(0, 20.0, 100)    # reaction time steps [s]

energy = 'off'                      # reactor energy
cti = 'data/papadikis.cti'          # Cantera input file

y0 = {'wood': 1, 'gas': 0, 'tar': 0, 'char': 0}

# Feedstocks and associated density [kg/m³]
feedstocks = ['Residues', 'Stem wood', 'Bark', 'Needles', 'Air classified (10 Hz)', 'Air classified (28 Hz)']
rho = [593, 598, 457, 621, 462, 535]

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas = ct.Solution(cti)
gas.TPY = temp, p, y0

r = ct.IdealGasReactor(gas, energy=energy)
sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

for t in time:
    sim.advance(t)
    states.append(r.thermo.state, t=t)

# Plot
# ----------------------------------------------------------------------------

_, ax = plt.subplots(tight_layout=True)
for i in range(len(feedstocks)):
    ax.plot(states.t, states('wood').Y * rho[i], label=feedstocks[i])
ax.set_xlabel('Time [s]')
ax.set_ylabel('Wood concentration [kg/m³]')
ax.grid(color='0.9')
ax.set_frame_on(False)
ax.tick_params(color='0.9')
ax.legend(loc='best')

_, ax = plt.subplots(tight_layout=True)
for i in range(len(feedstocks)):
    ax.plot(states.t, states('tar').Y * rho[i], label=feedstocks[i])
ax.set_xlabel('Time [s]')
ax.set_ylabel('Tar concentration [kg/m³]')
ax.grid(color='0.9')
ax.set_frame_on(False)
ax.tick_params(color='0.9')
ax.legend(loc='best')

_, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('tar').Y)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Tar mass fraction [-]')
ax.grid(color='0.9')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

plt.show()
