# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
Dictionary of common elements present for a given nominal mass

Those data are deprecated as they are stored in the sqlite database
"""

def getOrganicAt(m0):
    import numpy as np
    from . import get_mass
    res = []
    for C in range(1, 1+m0//12):
       for N in range(1+(m0-12*C)//15):
        for S in range(1+(m0-12*C-15*N)//32):
            for O in range(1+(m0-12*C-15*N-32*S)//16):
                H = m0-12*C-15*N-32*S-16*O
                elt = ''
                for e in 'CNSOH':
                    if locals()[e]>0:
                        elt += e
                        if locals()[e]>1:
                            elt += str(locals()[e])
                if int(np.round(get_mass(elt), 0))==m0:
                    res.append(elt)
    return res
    
elts = {
    1: ['H'],
    2: ['^2H','H2'],
    3: ['H3'],
    6: ['^6Li'],
    7: ['Li'],
    10:['^10B'],
    11:['B'],
    12:['C'],
    13:['CH'],
    14:['N','CH2'],
    15:['CH3'],
    16:['O','NH2','^13CH3'],
    17:['OH','NH3'],
    18:['NH4','H2O'],
    19:['H3O','CLi'],
    23:['Na'],
    24:['Mg','C2'],
    25:['^25Mg','C2H'],
    26:['^26Mg','^25MgH','MgH2','C2H2','CN'],
    27:['C2H3','Al','BO','CHN','BCH4','B2H5','BNH2'],
    28:['Si','AlH','CO','N2','CH2N','C2H4','^13CCH3','BNH3'],
    29:['^29Si','SiH','CHO','N2H','C2H5','CH3N'],
    30:['CH4N','C2H6'],
    31:['CF','CH3O','N2H3','CH5N','^13CH4N'],
    32:['S','O2','PH','CHF','CH4O','N2H4','CH6N'],
    33:['HS','O2H','NOH3','CH5O','N2H5','^13CH6N'],
    34:['^34S','H2S','NOH4','N2H6'],
    35:['H3S','N2H7'],
    36:['C3'],
    37:['C3H','NaN','C^25Mg'],
    38:['C3H2'],
    39:['C2HN','HF2','C3H3','K'],
    40:['C2H2N','C3H4','C2O','CN2','Ca'],
    41:['C2H3N','C3H5','CaH','^41K','C2HO'],
    42:['SiN','C2H4N','C3H6','CH2N2','CNO'],
    43:['C2H5N','C3H7','CH3N2','C2H3O'],
    44:['C3H8','C2H6N','C2H4O','CH4N2','CH2NO','SiO','AlOH','CO2'],
    45:['CHS','CHO2','C2H5O','CH5N2','C2H7N','C3H9'],
    46:['Na2','CH2S','NO2','CH4NO','C2H8N','^29SiHO'],
    47:['CH^34S','CH3S','CH3O2','CH7N2','C^13CH8N','SiF'],
    48:['C4','CH4S','SNH2','CH6NO','CH2^34S','SO','Ti'],
    49:['C4H','SOH','O3H','CH5S','SNH3','N2OH5','CH2Cl','CH5O2','^13CH4O2','NO2H3'],
    50:['O3H2','SOH2','^34SNH2','C4H2','C3N','CH6O2'],
    51:['C3HN','C4H3','CHF2'],
    52:['C3H2N','C4H4','C2N2','C3O'],
    53:['C4H5','C3H3N','C3HO','H3^34SO','C2HN2','CH2K'],
    54:['C3H4N','C3H2O','C4H6'],
    55:['C3H3O','C4H7','C3H5N','C2H3N2'],
    56:['Fe','C3H6N','C3H4O','C4H8','C2H2NO','C2H4N2'],
    57:['C2H3NO','C3H5O','C3H7N','C4H9','C2HS','C2H5N2','CaOH'],
    58:['CSN','C2H2S','N3O','C2H2O2','CH2N2O','N4H2','C3H8N','C2H4NO','C2H6N2','C3H6O'],
    59:['C3H7O','C2H5NO','C2H7N2','C3H9N','C2H3S'],
    60:['CSO','SN2','CH2SN','C2H4S','CH2NO2','C5','CH4N2O','C2H6NO','CH6N3','C3H8O','C2H8N2','C3H10N','C2H4O2','SiO2'],
    61:['CHSO','SN2H','CHO3','CH3SN', 'N2O2H','C5H','C2H5S','CH3NO2','N3OH3','C2H5O2','CH5N2O','N4H5','C2H7NO','C2H9N2'],
    62:['SNO','CH2SO','NO3','C5H2','C2H6S','CH4NO2','C2H8NO','C3H10O','^46TiO'],
    63:['^47TiO','C5H3','Cu','C4HN','CH3O3','Na2OH','CH5NO2','C2H7O2','CH7N2O'],
    64:['TiO'],
    65:['^49TiO','^65Cu'],
    66:['^50TiO'],
    69:['^69Ga'],
    70:['C5H10','C4H8N','C3H6N2','C2H4N3','C4H6O','C3H2O2','C2NO2'],
    71:['^71Ga'],
    74:['C3H8NO','CH2N2O2','C2H4NO2','C3H6O2','C3H10N2','C4H12N','C4H10O','C2H8N3','C2H6N2O','C3H6S','C6H2','C2H4SN','C5N','C2H2O3'],
    78:['^46TiO_2'],
    79:['^47TiO_2'],
    80:['TiO_2'],
    81:['^49TiO_2','C4HS','C4H3NO','C3H3N3','C5H5O','C4H5N2','C5H7N','C6H9'],
    82:['^50TiO_2','C4H6N2','C5H6O','C3H4N3','C5H8N','C6H10'],
    83:['C6H11','C5H9N','C4H7N2','C5H7O','C3H5N3','C4H5NO','C3H3N2O','C4H3O2','C2HN3O','C3HNO2'],
    84:['Si3','C3O3','C7','C2H2N3O','C4H4O2','C4H6NO','C5H8O','C3H6N3','C4H8N2','C5H10N','C6H12'],
    85:['C3H5N2O','C4H7NO','C3H7N3','C5H9O','C4H9N2','C5H11N','C6H13','C4H5O2','C2H3N3O','C3H3NO2','C7H','C2HN2O2','C3HO3'],
    95:['C4HNO2'],
    98:['C7H14','C6H12N','C5H10N2','C4H8N3','C3H6N4','C6H10O','C5H8NO','C5H6O2','C4H4NO2','C4H6N2O'],
    107:['^107Ag'],
    109:['^109Ag'],
    112:['C8H16','C7H14N','C6H12N2','C5H10N3','C4H8N4','C7H12O','C6H10NO','C5H8N2O','C6H8O2','C5H6NO2','C4H4N2O2'],
    113:['^113In'],
    115:['In'],
    140:['C10H20','C9H18N','C8H16N2','C9H16O','C8H14NO','C8H12O2','C7H10NO2','C7H12N2O','C7H14N3'],
    182: ['C15H2','C13H10O','C9H12NO3','C10H14O3','C13H12N','C14H14','C13H26','C12H24N','C11H22N2',
        'C10H20N3','C12H22O','C11H20NO','C11H18O2','C10H16NO2'],
    197:['Au'],
   }
   
neg_elts = {
    15: ['NH'],
    19: ['F'],
    31: ['P','^30SiH','^29SiH2'],
    32: ['S','O2'],
    35: ['Cl'],
    37: ['^37Cl'],
    45: ['PN','SiHO','CHS','CH2P','CHO2','CNF','N2OH','C2H2F','CH3NO','C2H5O'],
    47: ['PO'],
    52: ['CHK'],
    64: ['PO2H','SO2','CHSF','C4O','C2H5OF','CH6NO2'],
    80: ['PO3H','SO3','S2NH2','SiC3O','C4S','PSNH3','SNO2H2','SiC2N2','C4O2','CHO3F','C3N2O'],
    81: ['SO3H','C4HS','C2F3','C4HO2'],
    82: ['H2S2O','CHNOK','C3SN','SNOHF','C4H2S','CH3SOF','C3NO2','C2HF3','C2N3O','C4H2O2','CH3O3F'],
    83: ['C3H9F2','C5H4F','C2H2F3','C4H3S','CHO2F2','C4OF','SNOH2F','C3HSN','H3SO3','O4F',
        'CH2NOK','C4Cl','H3S2O','SiOK','CHCl2','SO2F'],
    84: ['SNOH3F','C3H2SN','CH3NOK','C2SN2','C4HCl','C3SO','CH2Cl2','CHSK','SClOH'],
    85: ['C5H6F','C4H5S','C4H2OF','COF3','C3H3SN','C2HSN2','C3HSO'],
    95: ['PO4','CH5SNO2','C5H3S','SNO3H'],
    96: ['SO4','C4SO','PO4H'],
    97: ['CH6PO3','CH2SOCl','SO4H','PH2O4','PO4H2'],
    98: ['H2SO4','C3SNO','C2H4Cl2','C3SNO','C2SN3','CO3F2','C4H2SO','C3NO3','C3H2SN2','C5F2','C2HOF3','C2HOF3','Al3OH','CHSNK','Si3N'],
    197: ['Au'],
}