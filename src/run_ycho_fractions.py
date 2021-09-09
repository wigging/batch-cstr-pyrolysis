"""
Calculate the C, H, and O mass fractions for the gas and liquid species. The
calculated fractions are used in the reactor models.
"""

import chemics as cm


def get_mw_total(species):
    mw_total = 0
    for i in range(len(species)):
        sp = species[i]
        mw_sp = cm.mw(sp)
        mw_total += mw_sp
    return mw_total


mw_c = cm.mw('C')
mw_h = cm.mw('H')
mw_o = cm.mw('O')

# Gas species
# ----------------------------------------------------------------------------

sp_gases = ['C2H4', 'C2H6', 'CH2O', 'CH4', 'CO', 'CO2', 'H2']
mol_c = [2, 2, 1, 1, 1, 1, 0]
mol_h = [4, 6, 2, 4, 0, 0, 2]
mol_o = [0, 0, 1, 0, 1, 2, 0]

print('Gases       yc       yh       yo')
for i in range(len(sp_gases)):
    sp = sp_gases[i]
    mw_sp = cm.mw(sp)
    yc = mol_c[i] * mw_c / mw_sp
    yh = mol_h[i] * mw_h / mw_sp
    yo = mol_o[i] * mw_o / mw_sp
    total = yc + yh + yo

    print(f'{sp:6} {yc:>8.2f} {yh:>8.2f} {yo:>8.2f} {total:>8.2f}')

# Liquid species
# ----------------------------------------------------------------------------

sp_liquids = [
    'C3H4O', 'C3H6O', 'C2H6O', 'C5H8O4', 'C6H10O5', 'C7H8O', 'C6H6O', 'C6H6O3',
    'C24H28O4', 'C3H6O2', 'C2H4O2', 'C2H4O', 'C2H4O2', 'CH4O', 'C2H2O2', 'C7H8O',
    'C5H4O2', 'H2O', 'CH2O2', 'C19H34O2', 'C13H22O2', 'C8H8O3'
]

mol_c = [3, 3, 2, 5, 6, 7, 6, 6, 24, 3, 2, 2, 2, 1, 2, 7, 5, 0, 1, 19, 13, 8]
mol_h = [4, 6, 6, 8, 10, 8, 6, 6, 28, 6, 4, 4, 4, 4, 2, 8, 4, 2, 2, 34, 22, 8]
mol_o = [1, 1, 1, 4, 5, 1, 1, 3, 4, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 3]

print('\nLiquids       yc       yh       yo')
for i in range(len(sp_liquids)):
    sp = sp_liquids[i]
    mw_sp = cm.mw(sp)
    yc = mol_c[i] * mw_c / mw_sp
    yh = mol_h[i] * mw_h / mw_sp
    yo = mol_o[i] * mw_o / mw_sp
    total = yc + yh + yo

    print(f'{sp:8} {yc:>8.2f} {yh:>8.2f} {yo:>8.2f} {total:>8.2f}')
