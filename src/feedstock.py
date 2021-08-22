import numpy as np


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
        Chemical analysis data given as weight percent (wt. %) reported on a
        dry basis (d). Values are [structural organics, non-structural
        organics, water extractable, ethanol extractives, acetone
        extractives, lignin, glucan, xylan, galactan, arabinan, mannan,
        acetyl].
    ADL : int
        Assume 22 wt. % for air-dry loss (ADL) when calculating as-received
        basis (ar) from the as-determined basis (ad).
    exp_yield: ndarray
        Experimental yield data given as weight percent (wt. %) wet basis.
        Yeild values are [oil, condensables, light gas, water vapor, char].
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
        self.ADL = 22
        self.exp_yield = np.array(data['yield'])

        # Calculate proximate and ultimate analysis bases
        self._prox_bases()
        self._ult_bases()

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

    def calc_chem_daf(self):
        """
        Calculate the chemical analysis dry ash-free basis (daf) from the dry
        basis in units of weight percent (wt. %).
        """
        chem_d = self.chem
        tot_d = sum(chem_d)
        chem_daf = chem_d * 100 / (tot_d - chem_d[0] - chem_d[1])
        chem_daf[0:2] = 0

        return chem_daf

    @staticmethod
    def calc_chem_bc(chem_daf):
        """
        Calculate the biomass composition from the chemical analysis data.
        cellulose = glucan
        hemicellulose = xylan + galactan + arabinan + mannan + acetyl
        """
        cell = chem_daf[6]
        hemi = chem_daf[7] + chem_daf[8] + chem_daf[9] + chem_daf[10] + chem_daf[11]
        lig = chem_daf[5]

        return np.array([cell, hemi, lig])

    def calc_yields(self):
        """
        Return experimental yields array and calculate normalized yield values
        from the experimental yield data. Yield values returned in order of
        [oil, condensables, light gas, water vapor, char].
        """
        oil = self.exp_yield[0]
        condensable = self.exp_yield[1]
        lightgas = self.exp_yield[2]
        watervap = self.exp_yield[3]
        char = self.exp_yield[4]

        exp_yields = np.array([oil, condensable, lightgas, watervap, char])
        norm_yields = exp_yields / sum(exp_yields) * 100

        return exp_yields, norm_yields

    @staticmethod
    def calc_lump_yields(exp_yields, norm_yields):
        """
        Calculate the lumped yields for comparing to the reactor model yields.
        Values returned in order of [gases, liquids, char] where
            gases = light gases
            liquids = oil + condensables + water vapor
            char = char
        """

        # Lumped groups using measured experiment yields
        exp_gases = exp_yields[2]
        exp_liquids = exp_yields[0] + exp_yields[1] + exp_yields[3]
        exp_char = exp_yields[4]
        exp_lumps = np.array([exp_gases, exp_liquids, exp_char])

        # Lumped groups using normalized experiment yields
        norm_gases = norm_yields[2]
        norm_liquids = norm_yields[0] + norm_yields[1] + norm_yields[3]
        norm_char = norm_yields[4]
        norm_lumps = np.array([norm_gases, norm_liquids, norm_char])

        return exp_lumps, norm_lumps
