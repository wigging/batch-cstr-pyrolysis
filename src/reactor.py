"""
Reactor module for batch and CSTR simulations.
"""

import cantera as ct
import numpy as np

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()

# Chemical species representing gas, liquid, solid, metaplastic phases
# Notice `N2` is added to the gases to include inlet fluidization gas flow
sp_gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2', 'N2')

sp_liquids = (
    'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
    'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H', 'CH3OH',
    'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12', 'VANILLIN')

sp_solids = (
    'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
    'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR', 'ACQUA')

sp_metaplastics = (
    'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
    'GH2', 'GC2H6')


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
