import numpy as np


class Feedstock:
    """
    Attributes
    ----------
    name : str
        Name of the feedstock.
    cycle : int
        Cycle (experiment) number corresponding to feedstock.
    prox : ndarray
        Proximate analysis FC, VM, ash, moisture as weight percent
        (wt. %) reported on an as-determined basis (ad).
    ult : ndarray
        Ultimate analysis C, H, O, N, S, ash, moisture as weight percent
        (wt. %) reported on an as-determined basis (ad). Reported values for
        H and O exclude H and O from moisture content.
    """

    def __init__(self, data):
        self.name = data['name']
        self.cycle = data['cycle']
        self.prox = np.array(data['proximate'])
        self.ult = np.array(data['ultimate'])
        self.chem = np.array(data['chemical'])

        # Assume 22 wt. % for air-dry loss (ADL) when calculating as-received
        # basis (ar) from the as-determined basis (ad)
        self.ADL = 22

    def calc_prox_ar(self):
        """
        Calculate proximate analysis as-received basis (ar) from as-determined
        basis (ad) in units of weight percent (wt. %).
        """
        FC_ad = self.prox[0]
        VM_ad = self.prox[1]
        ash_ad = self.prox[2]
        M_ad = self.prox[3]
        ADL = self.ADL

        M_ar = (M_ad * (100 - ADL) / 100) + ADL

        FC_ar = FC_ad * (100 - M_ar) / (100 - M_ad)
        VM_ar = VM_ad * (100 - M_ar) / (100 - M_ad)
        ash_ar = ash_ad * (100 - M_ar) / (100 - M_ad)

        return np.array([FC_ar, VM_ar, ash_ar, M_ar])

    def calc_prox_dry(self):
        """
        Calculate proximate analysis dry basis (d) from as-determined basis
        (ad) in units of weight percent (wt. %).
        """
        FC_ad = self.prox[0]
        VM_ad = self.prox[1]
        ash_ad = self.prox[2]
        M_ad = self.prox[3]

        FC_d = FC_ad * 100 / (100 - M_ad)
        VM_d = VM_ad * 100 / (100 - M_ad)
        ash_d = ash_ad * 100 / (100 - M_ad)

        return np.array([FC_d, VM_d, ash_d])

    def calc_prox_daf(self):
        """
        Calculate proximate analysis dry ash-free basis (daf) from
        as-determined basis (ad) in units of weight percent (wt. %).
        """
        FC_ad = self.prox[0]
        VM_ad = self.prox[1]
        ash_ad = self.prox[2]
        M_ad = self.prox[3]

        FC_daf = FC_ad * 100 / (100 - M_ad - ash_ad)
        VM_daf = VM_ad * 100 / (100 - M_ad - ash_ad)

        return np.array([FC_daf, VM_daf])

    def calc_ult_ar(self):
        """
        Calculate ultimate analysis as-received basis (ar) from as-determined
        basis (ad) in units of weight percent (wt. %).
        """
        C_ad = self.ult[0]
        H_ad = self.ult[1]
        O_ad = self.ult[2]
        N_ad = self.ult[3]
        S_ad = self.ult[4]
        ash_ad = self.ult[5]
        M_ad = self.ult[6]
        ADL = self.ADL

        M_ar = (M_ad * ((100 - ADL) / 100)) + ADL

        C_ar = C_ad * (100 - M_ar) / (100 - M_ad)
        H_ar = (H_ad - 0.1119 * M_ad) * (100 - M_ar) / (100 - M_ad)
        O_ar = (O_ad - 0.8881 * M_ad) * (100 - M_ar) / (100 - M_ad)
        N_ar = N_ad * (100 - M_ar) / (100 - M_ad)
        S_ar = S_ad * (100 - M_ar) / (100 - M_ad)
        ash_ar = ash_ad * (100 - M_ar) / (100 - M_ad)

        return np.array([C_ar, H_ar, O_ar, N_ar, S_ar, ash_ar, M_ar])

    def calc_ult_dry(self):
        """
        Calculate ultimate analysis dry basis (d) from as-determined basis
        (ad) in units of weight percent (wt. %).
        """
        C_ad = self.ult[0]
        H_ad = self.ult[1]
        O_ad = self.ult[2]
        N_ad = self.ult[3]
        S_ad = self.ult[4]
        ash_ad = self.ult[5]
        M_ad = self.ult[6]

        C_d = C_ad * 100 / (100 - M_ad)
        H_d = (H_ad - 0.1119 * M_ad) * 100 / (100 - M_ad)
        O_d = (O_ad - 0.8881 * M_ad) * 100 / (100 - M_ad)
        N_d = N_ad * 100 / (100 - M_ad)
        S_d = S_ad * 100 / (100 - M_ad)
        ash_d = ash_ad * 100 / (100 - M_ad)

        return np.array([C_d, H_d, O_d, N_d, S_d, ash_d])

    def calc_ult_daf(self):
        """
        Calculate ultimate analysis dry ash-free basis (daf) from
        as-determined basis (ad) in units of weight percent (wt. %).
        """
        C_ad = self.ult[0]
        H_ad = self.ult[1]
        O_ad = self.ult[2]
        N_ad = self.ult[3]
        S_ad = self.ult[4]
        ash_ad = self.ult[5]
        M_ad = self.ult[6]

        C_daf = C_ad * 100 / (100 - M_ad - ash_ad)
        H_daf = (H_ad - 0.1119 * M_ad) * 100 / (100 - M_ad - ash_ad)
        O_daf = (O_ad - 0.8881 * M_ad) * 100 / (100 - M_ad - ash_ad)
        N_daf = N_ad * 100 / (100 - M_ad - ash_ad)
        S_daf = S_ad * 100 / (100 - M_ad - ash_ad)

        return np.array([C_daf, H_daf, O_daf, N_daf, S_daf])

    @staticmethod
    def calc_ult_cho(ult_daf):
        """
        Calculate ultimate analysis CHO basis from daf basis in units of
        weight percent (wt. %).
        """
        C_daf = ult_daf[0]
        H_daf = ult_daf[1]
        O_daf = ult_daf[2]
        N_daf = ult_daf[3]
        S_daf = ult_daf[4]

        C_cho = C_daf * 100 / (100 - N_daf - S_daf)
        H_cho = H_daf * 100 / (100 - N_daf - S_daf)
        O_cho = O_daf * 100 / (100 - N_daf - S_daf)

        return np.array([C_cho, H_cho, O_cho])

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
    def calc_bio_comp(chem_daf):
        """
        Calculate the biomass composition such that
        cellulose = glucan
        hemicellulose = xylan + galactan + arabinan + mannan + acetyl
        """
        cell = chem_daf[6]
        hemi = chem_daf[7] + chem_daf[8] + chem_daf[9] + chem_daf[10] + chem_daf[11]
        lig = chem_daf[5]

        return np.array([cell, hemi, lig])
