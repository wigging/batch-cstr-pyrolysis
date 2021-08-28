"""
Batch reactor model class.
"""

import cantera as ct

# Disable warnings about discontinuity at polynomial mid-point in thermo data.
# Remove this line to show the warnings.
ct.suppress_thermo_warnings()


class BatchReactor:

    def __init__(self, y0, temp, pressure, time, cti, energy='off'):
        self.y0 = y0
        self.temp = temp
        self.pressure = pressure
        self.time = time
        self.cti = cti
        self.energy = energy

        # Chemical species representing gas, liquid, solid, metaplastic phases
        self.sp_gases = ('C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2')

        self.sp_liquids = (
            'C2H3CHO', 'C2H5CHO', 'C2H5OH', 'C5H8O4', 'C6H10O5', 'C6H5OCH3', 'C6H5OH',
            'C6H6O3', 'C24H28O4', 'CH2OHCH2CHO', 'CH2OHCHO', 'CH3CHO', 'CH3CO2H',
            'CH3OH', 'CHOCHO', 'CRESOL', 'FURFURAL', 'H2O', 'HCOOH', 'MLINO', 'U2ME12',
            'VANILLIN')

        self.sp_solids = (
            'CELL', 'CELLA', 'GMSW', 'HCE1', 'HCE2', 'ITANN', 'LIG', 'LIGC', 'LIGCC',
            'LIGH', 'LIGO', 'LIGOH', 'TANN', 'TGL', 'CHAR', 'ACQUA')

        self.sp_metaplastics = (
            'GCH2O', 'GCO2', 'GCO', 'GCH3OH', 'GCH4', 'GC2H4', 'GC6H5OH', 'GCOH2',
            'GH2', 'GC2H6')

        self.states = None

    def run_simulation(self):
        """
        Run the batch reactor simulation.
        """
        cti = self.cti
        temp = self.temp
        pressure = self.pressure
        y0 = self.y0
        time = self.time

        gas = ct.Solution(cti)
        gas.TPY = temp, pressure, y0

        r = ct.IdealGasReactor(gas, energy=self.energy)
        sim = ct.ReactorNet([r])
        states = ct.SolutionArray(gas, extra=['t'])

        for t in time:
            sim.advance(t)
            states.append(r.thermo.state, t=t)

        self.states = states

    def get_ycho_gases(self):
        """
        Carbon, hydrogen, and oxygen fractions at each time step for each gas
        species. Fraction values obtained from output of `ycho_fractions.py`.
        """
        states = self.states
        sp_gases = self.sp_gases

        c_fracs = [0.86, 0.80, 0.40, 0.75, 0.43, 0.27, 0]
        h_fracs = [0.14, 0.20, 0.07, 0.25, 0, 0, 1]
        o_fracs = [0, 0, 0.53, 0, 0.57, 0.73, 0]

        yc = states(*sp_gases).Y * c_fracs
        yh = states(*sp_gases).Y * h_fracs
        yo = states(*sp_gases).Y * o_fracs

        return yc, yh, yo

    def get_ycho_liquids(self):
        """
        Carbon, hydrogen, and oxygen fractions at each time step for each liquid
        species. Fraction values obtained from output of `ycho_fractions.py`.
        """
        states = self.states
        sp_liquids = self.sp_liquids

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

    def get_y_gasmix(self):
        """
        Mass fraction of the gas mixture at each time step.
        """
        states = self.states
        sp_gases = self.sp_gases
        y_gases = states(*sp_gases).Y.sum(axis=1)
        return y_gases

    def get_y_liqmix(self):
        """
        Mass fraction of the liquid mixture at each time step.
        """
        states = self.states
        sp_liquids = self.sp_liquids
        y_liquids = states(*sp_liquids).Y.sum(axis=1)
        return y_liquids

    def get_y_sldmix(self):
        """
        Mass fraction of the solid mixture at each time step.
        """
        states = self.states
        sp_solids = self.sp_solids
        y_solids = states(*sp_solids).Y.sum(axis=1)
        return y_solids

    def get_y_metamix(self):
        """
        Mass fraction of the metaplastic mixture at each time step.
        """
        states = self.states
        sp_metaplastics = self.sp_metaplastics
        y_metaplastics = states(*sp_metaplastics).Y.sum(axis=1)
        return y_metaplastics
