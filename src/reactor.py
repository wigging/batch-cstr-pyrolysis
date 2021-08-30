"""
Reactor module for batch and CSTR simulations.
"""

import cantera as ct
import numpy as np

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()

# Chemical species representing gas phase
sp_gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2')

# Chemical species representing gas phase including nitrogen for fluidization
sp_gases_n2 = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2', 'N2')

# Chemical species representing liquid phase
sp_liquids = (
    'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
    'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H', 'CH3OH',
    'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12', 'VANILLIN')

# Chemical species representing solid phase
sp_solids = (
    'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
    'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR', 'ACQUA')

# Chemical species representing metaplastic phase
sp_metaplastics = (
    'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
    'GH2', 'GC2H6')


def run_batch_simulation(cti, pressure, temp, time, y0, energy='off'):
    """
    Run the batch reactor simulation.
    """

    gas = ct.Solution(cti)
    gas.TPY = temp, pressure, y0

    r = ct.IdealGasReactor(gas, energy=energy)
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])

    for t in time:
        sim.advance(t)
        states.append(r.thermo.state, t=t)

    return states


def run_cstr_simulation(cti, diam, length, n_cstrs, pressure, tau, temp, y0, energy='off'):
    """
    Run the CSTR simulation.
    """

    # Setup gas phase
    gas = ct.Solution(cti)
    gas.TPY = temp, pressure, y0

    # Create a continuously stirred tank reactor (CSTR)
    dz = length / n_cstrs
    area = np.pi * (diam**2)
    vol = area * dz
    cstr = ct.IdealGasConstPressureReactor(gas, energy=energy, volume=vol)

    # Reservoirs for the inlet and outlet of each CSTR
    inlet = ct.Reservoir(gas)
    outlet = ct.Reservoir(gas)

    # Mass flow rate into the CSTR
    mf = ct.MassFlowController(upstream=inlet, downstream=cstr, mdot=cstr.mass / tau)

    # Determine pressure in the CSTR
    ct.PressureController(upstream=cstr, downstream=outlet, master=mf, K=1e-5)

    # Create a reactor network for performing the simulation
    sim = ct.ReactorNet([cstr])

    # Store the results for each CSTR
    states = ct.SolutionArray(cstr.thermo)
    states.append(cstr.thermo.state)

    for n in range(n_cstrs):
        gas.TPY = cstr.thermo.TPY
        inlet.syncState()
        sim.reinitialize()
        sim.advance_to_steady_state()
        states.append(cstr.thermo.state)

    return states


def get_ycho_gases(states):
    """
    Carbon, hydrogen, and oxygen fractions at each time step for each gas
    species. Fraction values obtained from output of `ycho_fractions.py`.
    """
    c_fracs = [0.86, 0.80, 0.40, 0.75, 0.43, 0.27, 0]
    h_fracs = [0.14, 0.20, 0.07, 0.25, 0, 0, 1]
    o_fracs = [0, 0, 0.53, 0, 0.57, 0.73, 0]

    yc = states(*sp_gases).Y * c_fracs
    yh = states(*sp_gases).Y * h_fracs
    yo = states(*sp_gases).Y * o_fracs

    return yc, yh, yo


def get_ycho_liquids(states):
    """
    Carbon, hydrogen, and oxygen fractions at each time step for each liquid
    species. Fraction values obtained from output of `ycho_fractions.py`.
    """
    c_fracs = [
        0.64, 0.62, 0.52, 0.45, 0.44, 0.78, 0.77, 0.57, 0.76, 0.49, 0.40,
        0.55, 0.40, 0.37, 0.41, 0.78, 0.63, 0.00, 0.26, 0.77, 0.74, 0.63
    ]

    h_fracs = [
        0.07, 0.10, 0.13, 0.06, 0.06, 0.07, 0.06, 0.05, 0.07, 0.08, 0.07,
        0.09, 0.07, 0.13, 0.03, 0.07, 0.04, 0.11, 0.04, 0.12, 0.11, 0.05
    ]

    o_fracs = [
        0.29, 0.28, 0.35, 0.48, 0.49, 0.15, 0.17, 0.38, 0.17, 0.43, 0.53,
        0.36, 0.53, 0.50, 0.55, 0.15, 0.33, 0.89, 0.70, 0.11, 0.15, 0.32
    ]

    yc = states(*sp_liquids).Y * c_fracs
    yh = states(*sp_liquids).Y * h_fracs
    yo = states(*sp_liquids).Y * o_fracs

    return yc, yh, yo
