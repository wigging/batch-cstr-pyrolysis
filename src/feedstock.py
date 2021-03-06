import chemics as cm
import numpy as np
from objfunc import objfunc
from scipy.optimize import minimize


class Feedstock:
    """
    Feedstock properties.

    Parameters
    ----------
    data : dict
        Dictionary containing feedstock data. Required keys are name, cycle,
        proximate, ultimate, chemical, and yield.

    Attributes
    ----------
    name : str
        Name of the feedstock.
    cycle : int
        Cycle (experiment) number corresponding to feedstock.
    prox_ad : ndarray
        Proximate analysis data reported on an as-determined basis (ad). Values
        are [FC, VM, ash, moisture] in units of weight percent (wt. %).
    prox_ar : ndarray
        Proximate analysis as-received basis (ar). Values are [FC, VM, ash,
        moisture] in units of weight percent (wt. %).
    prox_d : ndarray
        Proximate analysis dry basis (d). Values are [FC, VM, ash] in units of
        weight percent (wt. %).
    prox_daf : ndarray
        Proximate analysis dry ash-free basis (daf). Values are [FC, VM] in
        units of weight percent (wt. %).
    ult_ad : ndarray
        Ultimate analysis data reported on an as-determined basis (ad). Values
        are [C, H, O, N, S, ash, moisture] in units of weight percent
        (wt. %). Reported values for H and O exclude the H and O from
        moisture content.
    ult_ar : ndarray
        Ultimate analysis as-received basis (ar). Values are [C, H, O, N, S,
        ash, moisture] in units of weight percent (wt. %).
    ult_d : ndarray
        Ultimate analysis dry basis (d). Values are [C, H, O, N, S, ash] in
        units of weight percent (wt. %).
    ult_daf : ndarray
        Ultimate analysis dry ash-free basis (daf). Values are [C, H, O, N, S]
        in units of weight percent (wt. %).
    ult_cho : ndarray
        Ultimate analysis CHO basis. Values are [C, H, O] in units of weight
        percent (wt. %).
    chem_d : ndarray
        Chemical analysis data reported on a dry basis (d). Values are
        [structural organics, non-structural organics, water extractable,
        ethanol extractives, acetone extractives, lignin, glucan, xylan,
        galactan, arabinan, mannan, acetyl] in units of weight percent
        (wt. %).
    chem_daf : ndarray
        Chemical analysis dry ash-free basis (daf). Values are [0, 0, water
        extractable, ethanol extractives, acetone extractives, lignin,
        glucan, xylan, galactan, arabinan, mannan, acetyl] in units of weight
        percent(wt. %).
    chem_bc : ndarray
        Biomass composition calculated from the reported chemical analysis
        data on a dry ash-free basis (daf). Values are given as
        [cellulose, hemicellulose, lignin] in units of weight percent
        (wt. %).
    ADL : int
        Assume 22 wt. % for air-dry loss (ADL) when calculating as-received
        basis (ar) from the as-determined basis (ad).
    exp_yield : ndarray
        Experimental yield data reported on a wet basis. Values are
        [oil, condensables, light gas, water vapor, char] in units of weight
        percent (wt. %).
    normexp_yield : ndarray
        Normalized experimental yield data on a wet basis. Values are
        [oil, condensables, light gas, water vapor, char] in units of weight
        percent (wt. %).
    lump_yield : ndarray
        Lumped yields from measured experiment yield data. Values are
        [gases, liquids, solids] in units of weight percent (wt. %).
    lump2_yield : ndarray
        Lumped yields from measured experiment yield data. Values are
        [gases, liquids, solids] in units of weight percent (wt. %).
    normlump_yield : ndarray
        Normalized lumped yields from normalized experiment yield data. Values
        are [gases, liquids, solids] in units of weight percent (wt. %).
    residence_time : optional float
        Residence of the feedstock in the reactor [s].
    """

    def __init__(self, data):
        self.name = data['name']
        self.cycle = data['cycle']
        self.prox_ad = np.array(data['proximate'])
        self.prox_ar = None
        self.prox_d = None
        self.prox_daf = None
        self.ult_ad = np.array(data['ultimate'])
        self.ult_ar = None
        self.ult_d = None
        self.ult_daf = None
        self.ult_cho = None
        self.chem_d = np.array(data['chemical'])
        self.chem_daf = None
        self.chem_bc = None
        self.ADL = 22
        self.exp_yield = np.array(data['yield'])
        self.normexp_yield = np.array(data['yield']) / sum(np.array(data['yield'])) * 100
        self.lump_yield = None
        self.lump2_yield = None
        self.normlump_yield = None

        # Assign residence time if it is available
        if 'residenceTime' in data:
            self.residence_time = data['residenceTime']

        # Calculate proximate and ultimate analysis bases
        self._prox_bases()
        self._ult_bases()

        # Calculate lumped yields from experiment yield data
        self._lump_yields()

        # Calculate chemical analysis bases
        self._chem_bases()

    def _prox_bases(self):
        # Get as-determined (ad) values
        FC_ad = self.prox_ad[0]
        VM_ad = self.prox_ad[1]
        ash_ad = self.prox_ad[2]
        M_ad = self.prox_ad[3]
        ADL = self.ADL

        # Calculate proximate analysis as-received basis (ar) from
        # as-determined basis (ad) in units of weight percent (wt. %)
        M_ar = (M_ad * (100 - ADL) / 100) + ADL
        FC_ar = FC_ad * (100 - M_ar) / (100 - M_ad)
        VM_ar = VM_ad * (100 - M_ar) / (100 - M_ad)
        ash_ar = ash_ad * (100 - M_ar) / (100 - M_ad)
        self.prox_ar = np.array([FC_ar, VM_ar, ash_ar, M_ar])

        # Calculate proximate analysis dry basis (d) from as-determined basis
        # (ad) in units of weight percent (wt. %)
        FC_d = FC_ad * 100 / (100 - M_ad)
        VM_d = VM_ad * 100 / (100 - M_ad)
        ash_d = ash_ad * 100 / (100 - M_ad)
        self.prox_d = np.array([FC_d, VM_d, ash_d])

        # Calculate proximate analysis dry ash-free basis (daf) from
        # as-determined basis (ad) in units of weight percent (wt. %)
        FC_daf = FC_ad * 100 / (100 - M_ad - ash_ad)
        VM_daf = VM_ad * 100 / (100 - M_ad - ash_ad)
        self.prox_daf = np.array([FC_daf, VM_daf])

    def _ult_bases(self):
        # Get as-determined (ad) values
        C_ad = self.ult_ad[0]
        H_ad = self.ult_ad[1]
        O_ad = self.ult_ad[2]
        N_ad = self.ult_ad[3]
        S_ad = self.ult_ad[4]
        ash_ad = self.ult_ad[5]
        M_ad = self.ult_ad[6]
        ADL = self.ADL

        # Calculate ultimate analysis as-received basis (ar) from
        # as-determined basis (ad) in units of weight percent (wt. %)
        M_ar = (M_ad * ((100 - ADL) / 100)) + ADL
        C_ar = C_ad * (100 - M_ar) / (100 - M_ad)
        H_ar = (H_ad - 0.1119 * M_ad) * (100 - M_ar) / (100 - M_ad)
        O_ar = (O_ad - 0.8881 * M_ad) * (100 - M_ar) / (100 - M_ad)
        N_ar = N_ad * (100 - M_ar) / (100 - M_ad)
        S_ar = S_ad * (100 - M_ar) / (100 - M_ad)
        ash_ar = ash_ad * (100 - M_ar) / (100 - M_ad)
        self.ult_ar = np.array([C_ar, H_ar, O_ar, N_ar, S_ar, ash_ar, M_ar])

        # Calculate ultimate analysis dry basis (d) from as-determined basis
        # (ad) in units of weight percent (wt. %)
        C_d = C_ad * 100 / (100 - M_ad)
        H_d = (H_ad - 0.1119 * M_ad) * 100 / (100 - M_ad)
        O_d = (O_ad - 0.8881 * M_ad) * 100 / (100 - M_ad)
        N_d = N_ad * 100 / (100 - M_ad)
        S_d = S_ad * 100 / (100 - M_ad)
        ash_d = ash_ad * 100 / (100 - M_ad)
        self.ult_d = np.array([C_d, H_d, O_d, N_d, S_d, ash_d])

        # Calculate ultimate analysis dry ash-free basis (daf) from
        # as-determined basis (ad) in units of weight percent (wt. %)
        C_daf = C_ad * 100 / (100 - M_ad - ash_ad)
        H_daf = (H_ad - 0.1119 * M_ad) * 100 / (100 - M_ad - ash_ad)
        O_daf = (O_ad - 0.8881 * M_ad) * 100 / (100 - M_ad - ash_ad)
        N_daf = N_ad * 100 / (100 - M_ad - ash_ad)
        S_daf = S_ad * 100 / (100 - M_ad - ash_ad)
        self.ult_daf = np.array([C_daf, H_daf, O_daf, N_daf, S_daf])

        # Calculate ultimate analysis CHO basis from daf basis in units of
        # weight percent (wt. %)
        C_cho = C_daf * 100 / (100 - N_daf - S_daf)
        H_cho = H_daf * 100 / (100 - N_daf - S_daf)
        O_cho = O_daf * 100 / (100 - N_daf - S_daf)
        self.ult_cho = np.array([C_cho, H_cho, O_cho])

    def _lump_yields(self):
        # Calculate lumped yields from measured experiment yield data where
        # experiment yields = [oil, condensables, light gas, water vapor, char]
        # gases = light gas
        # liquids = oil + condensables + water vapor
        # solids = char
        gases = self.exp_yield[2]
        liquids = self.exp_yield[0] + self.exp_yield[1] + self.exp_yield[3]
        solids = self.exp_yield[4]
        self.lump_yield = np.array([gases, liquids, solids])

        # Another approach to lump the experiment yields as the following
        # gases = light gas + condensables + water vapor
        # liquids = oil
        # solids = char
        gases = self.exp_yield[2] + self.exp_yield[1] + self.exp_yield[3]
        liquids = self.exp_yield[0]
        solids = self.exp_yield[4]
        self.lump2_yield = np.array([gases, liquids, solids])

        # Calculate normalized lumped yields from normalized experiment yield data.
        norm_gases = self.normexp_yield[2]
        norm_liquids = self.normexp_yield[0] + self.normexp_yield[1] + self.normexp_yield[3]
        norm_solids = self.normexp_yield[4]
        self.normlump_yield = np.array([norm_gases, norm_liquids, norm_solids])

    def _chem_bases(self):
        # Calculate the chemical analysis dry ash-free basis (daf) from the
        # dry basis in units of weight percent (wt. %).
        chem_d = self.chem_d
        total_d = sum(chem_d)
        chem_daf = chem_d * 100 / (total_d - chem_d[0] - chem_d[1])
        chem_daf[0:2] = 0
        self.chem_daf = chem_daf

        # Calculate the biomass composition from chemical analysis dry
        # ash-free basis (daf) values where
        # cellulose = glucan
        # hemicellulose = xylan + galactan + arabinan + mannan + acetyl
        cell = chem_daf[6]
        hemi = chem_daf[7] + chem_daf[8] + chem_daf[9] + chem_daf[10] + chem_daf[11]
        lig = chem_daf[5]
        self.chem_bc = np.array([cell, hemi, lig])

    def calc_biocomp(self):
        """
        Calculate the optimized splitting parameters and associated biomass
        composition (daf) of the feedstock.
        """

        # Get mass fractions for ultimate analysis data
        # Get mass fractions of biomass composition from chemical analysis data
        yc = self.ult_cho[0] / 100
        yh = self.ult_cho[1] / 100
        ybc = self.chem_bc / 100

        # Determine optimized splitting parameters using default values for `x0`
        # where each parameter is bound within 0 to 1
        x0 = [0.6, 0.8, 0.8, 1, 1]
        bnds = ((0, 1), (0, 1), (0, 1), (0, 1), (0, 1))
        res = minimize(objfunc, x0, args=(yc, yh, ybc), method='L-BFGS-B', bounds=bnds)

        # Calculate biomass composition as dry ash-free basis (daf)
        bc = cm.biocomp(yc, yh, alpha=res.x[0], beta=res.x[1], gamma=res.x[2], delta=res.x[3], epsilon=res.x[4])
        cell, hemi, ligc, ligh, ligo, tann, tgl = bc['y_daf']

        # Optimized splitting parameters in order of [alpha, beta, gamma, delta, epsilon]
        splits = np.array([res.x[0], res.x[1], res.x[2], res.x[3], res.x[4]])

        return bc, splits
