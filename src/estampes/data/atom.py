"""Module on atoms-related data.

This module provides basic data related to atoms.

Attributes
----------
ELEMENTS : list
    atomic number -> atomic Label list

"""

from estampes.base import TypeAtData, TypeAtLab


# ================
# Module Constants
# ================

ELEMENTS = ('',  # NOQA
    'H' ,                                                                                'He',  # NOQA
    'Li','Be'                                                  ,'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',  # NOQA
    'Na','Mg'                                                  ,'Al','Si','P' ,'S' ,'Cl','Ar',  # NOQA
    'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',  # NOQA
    'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',  # NOQA
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',  # NOQA
        'Hf','Ta',' W','Re','Os','Ir','Pt', 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',  # NOQA
    'Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',  # NOQA
        'Rf','Db','Sg','Bh','Hs','Mt','Ds')  # NOQA


# ================
# Module Functions
# ================

def atomic_data(*atoms: TypeAtLab) -> TypeAtData:
    """Generates atomic data.

    Generates a dictionary containing atomic data for each atom given in
    argument.  Each item of the returned dictionary contains the
    following elements:

    `name`
        Full name of the atom.
    `num`
        Atomic number.
    `mass`
        Atomic mass (in amu).
    `rcov`
        Covalent radius (in Ang).
    `rvdw`
        Van der Waals radiuus (in Ang).
    `rvis`
        Visualization-related radius (in Ang).
    `rgb`
        Color of the atom (as tuple of integer 0-255).

    Parameters
    ----------
    *atoms
        List of atomic symbols (or numbers).

    Returns
    -------
    dict
        Dictionary of atomic data, containing the elements listed above
        grouped by atomic symbols.

    Raises
    ------
    KeyError
        Unrecognized atomic symbol.

    Notes
    -----
    * The function has a very basic support of symbols and atomic
      numbers.  A more robust procedure should rely on a first
      conversion by estampes.tools.convert_labsymb
    """
    at_data = {}
    for atom in set(atoms):
        try:
            at_symb = atom.title()
            at_idx = at_symb
        except AttributeError:
            at_symb = ELEMENTS[atom]
            at_idx = atom
        if at_symb == 'H':
            at_data[at_idx] = {
                'symb': 'H',
                'name': 'Hydrogen',
                'num': 1,
                'mass': 1.00790,
                'rcov': (32, None, None),
                'rvdw': 0.430,
                'rvis': 0.200,
                # 'rgb': (230, 230, 230),
                'rgb': (255, 255, 255),
            }
        elif at_symb == 'He':
            at_data[at_idx] = {
                'symb': 'He',
                'name': 'Helium',
                'num': 2,
                'mass': 4.00260,
                'rcov': (46, None, None),
                'rvdw': 0.741,
                'rvis': 0.286,
                'rgb': (217, 255, 255),
            }
        elif at_symb == 'Li':
            at_data[at_idx] = {
                'symb': 'Li',
                'name': 'Lithium',
                'num': 3,
                'mass': 6.94000,
                'rcov': (133, 124, None),
                'rvdw': 0.880,
                'rvis': 0.340,
                'rgb': (204, 128, 255),
            }
        elif at_symb == 'Be':
            at_data[at_idx] = {
                'symb': 'Be',
                'name': 'Beryllium',
                'num': 4,
                'mass': 9.01218,
                'rcov': (102, 90, 84),
                'rvdw': 0.550,
                'rvis': 0.589,
                'rgb': (194, 255, 0),
            }
        elif at_symb == 'B':
            at_data[at_idx] = {
                'symb': 'B',
                'name': 'Boron',
                'num': 5,
                'mass': 10.81000,
                'rcov': (85, 78, 73),
                'rvdw': 1.030,
                'rvis': 0.415,
                'rgb': (255, 181, 181),
            }
        elif at_symb == 'C':
            at_data[at_idx] = {
                'symb': 'C',
                'name': 'Carbon',
                'num': 6,
                'mass': 12.01100,
                'rcov': (75, 67, 60),
                'rvdw': 0.900,
                'rvis': 0.400,
                # 'rgb': (102, 102, 102),
                'rgb': (144, 144, 144),
            }
        elif at_symb == 'N':
            at_data[at_idx] = {
                'symb': 'N',
                'name': 'Nitrogen',
                'num': 7,
                'mass': 14.00670,
                'rcov': (71, 60, 54),
                'rvdw': 0.880,
                'rvis': 0.400,
                # 'rgb': (30, 76, 175),
                'rgb': (48, 80, 248),
            }
        elif at_symb == 'O':
            at_data[at_idx] = {
                'symb': 'O',
                'name': 'Oxygen',
                'num': 8,
                'mass': 15.99940,
                'rcov': (63, 57, 53),
                'rvdw': 0.880,
                'rvis': 0.400,
                # 'rgb': (255, 0, 0),
                'rgb': (255, 13, 13),
            }
        elif at_symb == 'F':
            at_data[at_idx] = {
                'symb': 'F',
                'name': 'Fluorine',
                'num': 9,
                'mass': 18.99840,
                'rcov': (64, 59, 53),
                'rvdw': 0.840,
                'rvis': 0.320,
                'rgb': (144, 224,  80),  # 50, 205,  50
            }
        elif at_symb == 'Ne':
            at_data[at_idx] = {
                'symb': 'Ne',
                'name': 'Neon',
                'num': 10,
                'mass': 20.17900,
                'rcov': (67, 96, None),
                'rvdw': 0.815,
                'rvis': 0.423,
                'rgb': (179, 227, 245),
            }
        elif at_symb == 'Na':
            at_data[at_idx] = {
                'symb': 'Na',
                'name': 'Sodium',
                'num': 11,
                'mass': 22.98977,
                'rcov': (155, 160, None),
                'rvdw': 1.170,
                'rvis': 0.485,
                'rgb': (171, 92, 242),
            }
        elif at_symb == 'Mg':
            at_data[at_idx] = {
                'symb': 'Mg',
                'name': 'Magnesium',
                'num': 12,
                'mass': 24.30500,
                'rcov': (139, 132, 127),
                'rvdw': 1.300,
                'rvis': 0.550,
                'rgb': (138, 255, 0),
            }
        elif at_symb == 'Al':
            at_data[at_idx] = {
                'symb': 'Al',
                'name': 'Aluminium',
                'num': 13,
                'mass': 26.98154,
                'rcov': (126, 113, 111),
                'rvdw': 1.550,
                'rvis': 0.675,
                'rgb': (191, 166, 166),
            }
        elif at_symb == 'Si':
            at_data[at_idx] = {
                'symb': 'Si',
                'name': 'Silicon',
                'num': 14,
                'mass': 28.08550,
                'rcov': (116, 107, 102),
                'rvdw': 1.400,
                'rvis': 0.600,
                'rgb': (240, 200, 160),
            }
        elif at_symb == 'P':
            at_data[at_idx] = {
                'symb': 'P',
                'name': 'Phosphorus',
                'num': 15,
                'mass': 30.97376,
                'rcov': (111, 102, 94),
                'rvdw': 1.250,
                'rvis': 0.525,
                # 'rgb': (128, 32, 255),
                'rgb': (255, 128, 0),
            }
        elif at_symb == 'S':
            at_data[at_idx] = {
                'symb': 'S',
                'name': 'Sulfur',
                'num': 16,
                'mass': 32.06000,
                'rcov': (103, 94, 95),
                'rvdw': 1.220,
                'rvis': 0.510,
                # 'rgb': (255, 191, 0),
                'rgb': (255, 255, 48),
            }
        elif at_symb == 'Cl':
            at_data[at_idx] = {
                'symb': 'Cl',
                'name': 'Chlorine',
                'num': 17,
                'mass': 35.45300,
                'rcov': (99, 95, 93),
                'rvdw': 1.190,
                'rvis': 0.495,
                'rgb': (31, 240, 31),
            }
        elif at_symb == 'Ar':
            at_data[at_idx] = {
                'symb': 'Ar',
                'name': 'Argon',
                'num': 18,
                'mass': 39.94800,
                'rcov': (96, 107, 96),
                'rvdw': 0.995,
                'rvis': 0.508,
                'rgb': (128, 209, 227),
            }
        elif at_symb == 'K':
            at_data[at_idx] = {
                'symb': 'K',
                'name': 'Potassium',
                'num': 19,
                'mass': 39.09830,
                'rcov': (196, 193, None),
                'rvdw': 1.530,
                'rvis': 0.665,
                'rgb': (143, 64, 212),
            }
        elif at_symb == 'Ca':
            at_data[at_idx] = {
                'symb': 'Ca',
                'name': 'Calcium',
                'num': 20,
                'mass': 40.08000,
                'rcov': (171, 147, 133),
                'rvdw': 1.190,
                'rvis': 0.495,
                'rgb': (61, 255, 0),
            }
        elif at_symb == 'Sc':
            at_data[at_idx] = {
                'symb': 'Sc',
                'name': 'Scandium',
                'num': 21,
                'mass': 44.95590,
                'rcov': (148, 116, 114),
                'rvdw': 1.670,
                'rvis': 0.735,
                'rgb': (230, 230, 230),
            }
        elif at_symb == 'Ti':
            at_data[at_idx] = {
                'symb': 'Ti',
                'name': 'Titanium',
                'num': 22,
                'mass': 47.90000,
                'rcov': (136, 117, 108),
                'rvdw': 1.640,
                'rvis': 0.720,
                'rgb': (191, 194, 199),
            }
        elif at_symb == 'V':
            at_data[at_idx] = {
                'symb': 'V',
                'name': 'Vanadium',
                'num': 23,
                'mass': 50.94150,
                'rcov': (134, 112, 106),
                'rvdw': 1.530,
                'rvis': 0.665,
                'rgb': (166, 166, 171),
            }
        elif at_symb == 'Cr':
            at_data[at_idx] = {
                'symb': 'Cr',
                'name': 'Chromium',
                'num': 24,
                'mass': 51.99600,
                'rcov': (122, 111, 103),
                'rvdw': 1.550,
                'rvis': 0.675,
                'rgb': (138, 153, 199),
            }
        elif at_symb == 'Mn':
            at_data[at_idx] = {
                'symb': 'Mn',
                'name': 'Manganese',
                'num': 25,
                'mass': 54.93800,
                'rcov': (119, 105, 103),
                'rvdw': 1.550,
                'rvis': 0.675,
                'rgb': (156, 122, 199),
            }
        elif at_symb == 'Fe':
            at_data[at_idx] = {
                'symb': 'Fe',
                'name': 'Iron',
                'num': 26,
                'mass': 55.84700,
                'rcov': (116, 109, 102),
                'rvdw': 1.540,
                'rvis': 0.670,
                'rgb': (224, 102, 51),
            }
        elif at_symb == 'Co':
            at_data[at_idx] = {
                'symb': 'Co',
                'name': 'Cobalt',
                'num': 27,
                'mass': 58.93320,
                'rcov': (111, 103, 96),
                'rvdw': 1.530,
                'rvis': 0.615,
                'rgb': (240, 144, 160),
            }
        elif at_symb == 'Ni':
            at_data[at_idx] = {
                'symb': 'Ni',
                'name': 'Nickel',
                'num': 28,
                'mass': 58.71000,
                'rcov': (110, 101, 101),
                'rvdw': 1.700,
                'rvis': 0.750,
                'rgb': (80, 208, 80),
            }
        elif at_symb == 'Cu':
            at_data[at_idx] = {
                'symb': 'Cu',
                'name': 'Copper',
                'num': 29,
                'mass': 63.54600,
                'rcov': (112, 115, 120),
                'rvdw': 1.720,
                'rvis': 0.760,
                'rgb': (200, 128, 51),
            }
        elif at_symb == 'Zn':
            at_data[at_idx] = {
                'symb': 'Zn',
                'name': 'Zinc',
                'num': 30,
                'mass': 65.38000,
                'rcov': (118, 120, None),
                'rvdw': 1.650,
                'rvis': 0.725,
                'rgb': (125, 128, 176),
            }
        elif at_symb == 'Ga':
            at_data[at_idx] = {
                'symb': 'Ga',
                'name': 'Gallium',
                'num': 31,
                'mass': 69.73500,
                'rcov': (124, 117, 121),
                'rvdw': 1.420,
                'rvis': 0.610,
                'rgb': (194, 143, 143),
            }
        elif at_symb == 'Ge':
            at_data[at_idx] = {
                'symb': 'Ge',
                'name': 'Germanium',
                'num': 32,
                'mass': 72.59000,
                'rcov': (121, 111, 114),
                'rvdw': 1.370,
                'rvis': 0.585,
                'rgb': (102, 143, 143),
            }
        elif at_symb == 'As':
            at_data[at_idx] = {
                'symb': 'As',
                'name': 'Arsenic',
                'num': 33,
                'mass': 74.92160,
                'rcov': (121, 114, 106),
                'rvdw': 1.410,
                'rvis': 0.605,
                'rgb': (189, 128, 227),
            }
        elif at_symb == 'Se':
            at_data[at_idx] = {
                'symb': 'Se',
                'name': 'Selenium',
                'num': 34,
                'mass': 78.96000,
                'rcov': (116, 107, 107),
                'rvdw': 1.420,
                'rvis': 0.610,
                'rgb': (255, 161, 0),
            }
        elif at_symb == 'Br':
            at_data[at_idx] = {
                'symb': 'Br',
                'name': 'Bromine',
                'num': 35,
                'mass': 79.90400,
                'rcov': (114, 109, 110),
                'rvdw': 1.410,
                'rvis': 0.605,
                'rgb': (166,  41,  41),
            }
        elif at_symb == 'Kr':
            at_data[at_idx] = {
                'symb': 'Kr',
                'name': 'Krypton',
                'num': 36,
                'mass': 83.80000,
                'rcov': (117, 121, 108),
                'rvdw': 1.069,
                'rvis': 0.524,
                'rgb': (92, 184, 209),
            }
        elif at_symb == 'Rb':
            at_data[at_idx] = {
                'symb': 'Rb',
                'name': 'Rubidium',
                'num': 37,
                'mass': 85.46780,
                'rcov': (210, 202, None),
                'rvdw': 1.670,
                'rvis': 0.735,
                'rgb': (112, 46, 176),
            }
        elif at_symb == 'Sr':
            at_data[at_idx] = {
                'symb': 'Sr',
                'name': 'Strontium',
                'num': 38,
                'mass': 87.62000,
                'rcov': (185, 157, 139),
                'rvdw': 1.320,
                'rvis': 0.560,
                'rgb': (0, 255, 0),
            }
        elif at_symb == 'Y':
            at_data[at_idx] = {
                'symb': 'Y',
                'name': 'Yttrium',
                'num': 39,
                'mass': 88.90590,
                'rcov': (163, 130, 124),
                'rvdw': 1.980,
                'rvis': 0.890,
                'rgb': (148, 255, 255),
            }
        elif at_symb == 'Zr':
            at_data[at_idx] = {
                'symb': 'Zr',
                'name': 'Zirconium',
                'num': 40,
                'mass': 91.22000,
                'rcov': (154, 127, 121),
                'rvdw': 1.760,
                'rvis': 0.780,
                'rgb': (148, 224, 224),
            }
        elif at_symb == 'Nb':
            at_data[at_idx] = {
                'symb': 'Nb',
                'name': 'Niobium',
                'num': 41,
                'mass': 92.90640,
                'rcov': (147, 125, 116),
                'rvdw': 1.680,
                'rvis': 0.740,
                'rgb': (115, 194, 201),
            }
        elif at_symb == 'Mo':
            at_data[at_idx] = {
                'symb': 'Mo',
                'name': 'Molybdenum',
                'num': 42,
                'mass': 95.94000,
                'rcov': (138, 121, 113),
                'rvdw': 1.670,
                'rvis': 0.735,
                'rgb': (84, 181, 181),
            }
        elif at_symb == 'Tc':
            at_data[at_idx] = {
                'symb': 'Tc',
                'name': 'Technetium',
                'num': 43,
                'mass': 98.90620,
                'rcov': (128, 120, 110),
                'rvdw': 1.550,
                'rvis': 0.675,
                'rgb': (59, 158, 158),
            }
        elif at_symb == 'Ru':
            at_data[at_idx] = {
                'symb': 'Ru',
                'name': 'Ruthenium',
                'num': 44,
                'mass': 101.0700,
                'rcov': (125, 114, 103),
                'rvdw': 1.600,
                'rvis': 0.700,
                'rgb': (36, 143, 143),
            }
        elif at_symb == 'Rh':
            at_data[at_idx] = {
                'symb': 'Rh',
                'name': 'Rhodium',
                'num': 45,
                'mass': 102.9055,
                'rcov': (125, 110, 106),
                'rvdw': 1.650,
                'rvis': 0.725,
                'rgb': (10, 125, 140),
            }
        elif at_symb == 'Pd':
            at_data[at_idx] = {
                'symb': 'Pd',
                'name': 'Palladium',
                'num': 46,
                'mass': 106.4000,
                'rcov': (120, 117, 112),
                'rvdw': 1.700,
                'rvis': 0.750,
                'rgb': (0, 105, 133),
            }
        elif at_symb == 'Ag':
            at_data[at_idx] = {
                'symb': 'Ag',
                'name': 'Silver',
                'num': 47,
                'mass': 107.8680,
                'rcov': (128, 139, 137),
                'rvdw': 1.790,
                'rvis': 0.795,
                'rgb': (192, 192, 192),
            }
        elif at_symb == 'Cd':
            at_data[at_idx] = {
                'symb': 'Cd',
                'name': 'Cadmium',
                'num': 48,
                'mass': 112.4100,
                'rcov': (136, 144, None),
                'rvdw': 1.890,
                'rvis': 0.845,
                'rgb': (255, 217, 143),
            }
        elif at_symb == 'In':
            at_data[at_idx] = {
                'symb': 'In',
                'name': 'Indium',
                'num': 49,
                'mass': 114.8200,
                'rcov': (142, 136, 146),
                'rvdw': 1.830,
                'rvis': 0.815,
                'rgb': (166, 117, 115),
            }
        elif at_symb == 'Sn':
            at_data[at_idx] = {
                'symb': 'Sn',
                'name': 'Tin',
                'num': 50,
                'mass': 118.6900,
                'rcov': (140, 130, 132),
                'rvdw': 1.660,
                'rvis': 0.730,
                'rgb': (102, 128, 128),
            }
        elif at_symb == 'Sb':
            at_data[at_idx] = {
                'symb': 'Sb',
                'name': 'Antimony',
                'num': 51,
                'mass': 121.7500,
                'rcov': (140, 133, 127),
                'rvdw': 1.660,
                'rvis': 0.730,
                'rgb': (158, 99, 181),
            }
        elif at_symb == 'Te':
            at_data[at_idx] = {
                'symb': 'Te',
                'name': 'Tellurium',
                'num': 52,
                'mass': 127.6000,
                'rcov': (136, 128, 121),
                'rvdw': 1.670,
                'rvis': 0.735,
                'rgb': (212, 122, 0),
            },
        elif at_symb == 'I':
            at_data[at_idx] = {
                'symb': 'I',
                'name': 'Iodine',
                'num': 53,
                'mass': 126.9045,
                'rcov': (133, 129, 125),
                'rvdw': 1.600,
                'rvis': 0.700,
                'rgb': (148, 0, 148),
            }
        elif at_symb == 'Xe':
            at_data[at_idx] = {
                'symb': 'Xe',
                'name': 'Xenon',
                'num': 54,
                'mass': 131.3000,
                'rcov': (131, 135, 122),
                'rvdw': 1.750,
                'rvis': 0.577,
                'rgb': (66, 158, 176),
            }
        elif at_symb == 'Cs':
            at_data[at_idx] = {
                'symb': 'Cs',
                'name': 'Caesium',
                'num': 55,
                'mass': 132.9054,
                'rcov': (232, 209, None),
                'rvdw': 1.870,
                'rvis': 0.835,
                'rgb': (87, 23, 143),
            }
        elif at_symb == 'Ba':
            at_data[at_idx] = {
                'symb': 'Ba',
                'name': 'Barium',
                'num': 56,
                'mass': 137.3300,
                'rcov': (196, 161, 149),
                'rvdw': 1.540,
                'rvis': 0.670,
                'rgb': (0, 201, 0),
            }
        elif at_symb == 'La':
            at_data[at_idx] = {
                'symb': 'La',
                'name': 'Lanthanum',
                'num': 57,
                'mass': 138.9055,
                'rcov': (180, 139, 139),
                'rvdw': 2.070,
                'rvis': 0.935,
                'rgb': (112, 212, 255),
            }
        elif at_symb == 'Ce':
            at_data[at_idx] = {
                'symb': 'Ce',
                'name': 'Cerium',
                'num': 58,
                'mass': 140.1200,
                'rcov': (163, 137, 131),
                'rvdw': 2.030,
                'rvis': 0.915,
                'rgb': (255, 255, 199),
            }
        elif at_symb == 'Pr':
            at_data[at_idx] = {
                'symb': 'Pr',
                'name': 'Praseodymium',
                'num': 59,
                'mass': 140.9077,
                'rcov': (176, 138, 128),
                'rvdw': 2.020,
                'rvis': 0.910,
                'rgb': (217, 255, 199),
            }
        elif at_symb == 'Nd':
            at_data[at_idx] = {
                'symb': 'Nd',
                'name': 'Neodymium',
                'num': 60,
                'mass': 144.2400,
                'rcov': (174, 132, None),
                'rvdw': 2.010,
                'rvis': 0.905,
                'rgb': (199, 255, 199),
            }
        elif at_symb == 'Pm':
            at_data[at_idx] = {
                'symb': 'Pm',
                'name': 'Promethium',
                'num': 61,
                'mass': None,
                'rcov': (173, 135, None),
                'rvdw': 2.000,
                'rvis': 0.900,
                'rgb': (163, 255, 199),
            }
        elif at_symb == 'Sm':
            at_data[at_idx] = {
                'symb': 'Sm',
                'name': 'Samarium',
                'num': 62,
                'mass': 150.3600,
                'rcov': (172, 134, None),
                'rvdw': 2.000,
                'rvis': 0.900,
                'rgb': (143, 255, 199),
            }
        elif at_symb == 'Eu':
            at_data[at_idx] = {
                'symb': 'Eu',
                'name': 'Europium',
                'num': 63,
                'mass': 151.9600,
                'rcov': (168, 134, None),
                'rvdw': 2.190,
                'rvis': 0.995,
                'rgb': (97, 255, 199),
            }
        elif at_symb == 'Gd':
            at_data[at_idx] = {
                'symb': 'Gd',
                'name': 'Gadolinium',
                'num': 64,
                'mass': 157.2500,
                'rcov': (169, 135, 132),
                'rvdw': 1.990,
                'rvis': 0.895,
                'rgb': (69, 255, 199),
            }
        elif at_symb == 'Tb':
            at_data[at_idx] = {
                'symb': 'Tb',
                'name': 'Terbium',
                'num': 65,
                'mass': 158.9253,
                'rcov': (168, 135, None),
                'rvdw': 1.960,
                'rvis': 0.880,
                'rgb': (48, 255, 199),
            }
        elif at_symb == 'Dy':
            at_data[at_idx] = {
                'symb': 'Dy',
                'name': 'Dysprosium',
                'num': 66,
                'mass': 162.5000,
                'rcov': (167, 133, None),
                'rvdw': 1.950,
                'rvis': 0.875,
                'rgb': (31, 255, 199),
            }
        elif at_symb == 'Ho':
            at_data[at_idx] = {
                'symb': 'Ho',
                'name': 'Holmium',
                'num': 67,
                'mass': 164.9303,
                'rcov': (166, 133, None),
                'rvdw': 1.940,
                'rvis': 0.870,
                'rgb': (0, 255, 156),
            }
        elif at_symb == 'Er':
            at_data[at_idx] = {
                'symb': 'Er',
                'name': 'Erbium',
                'num': 68,
                'mass': 167.2600,
                'rcov': (165, 133, None),
                'rvdw': 1.930,
                'rvis': 0.865,
                'rgb': (0, 230, 117),
            }
        elif at_symb == 'Tm':
            at_data[at_idx] = {
                'symb': 'Tm',
                'name': 'Thulium',
                'num': 69,
                'mass': 168.9342,
                'rcov': (164, 131, None),
                'rvdw': 1.920,
                'rvis': 0.860,
                'rgb': (0, 212, 82),
            }
        elif at_symb == 'Yb':
            at_data[at_idx] = {
                'symb': 'Yb',
                'name': 'Ytterbium',
                'num': 70,
                'mass': 173.0500,
                'rcov': (170, 129, None),
                'rvdw': 2.140,
                'rvis': 0.970,
                'rgb': (0, 191, 56),
            }
        elif at_symb == 'Lu':
            at_data[at_idx] = {
                'symb': 'Lu',
                'name': 'Lutetium',
                'num': 71,
                'mass': 174.9670,
                'rcov': (162, 131, 131),
                'rvdw': 1.920,
                'rvis': 0.860,
                'rgb': (0, 171, 36),
            }
        elif at_symb == 'Hf':
            at_data[at_idx] = {
                'symb': 'Hf',
                'name': 'Hafnium',
                'num': 72,
                'mass': 178.4900,
                'rcov': (152, 128, 122),
                'rvdw': 1.770,
                'rvis': 0.785,
                'rgb': (77, 194, 255),
            }
        elif at_symb == 'Ta':
            at_data[at_idx] = {
                'symb': 'Ta',
                'name': 'Tantalum',
                'num': 73,
                'mass': 180.9479,
                'rcov': (146, 126, 119),
                'rvdw': 1.630,
                'rvis': 0.715,
                'rgb': (77, 166, 255),
            }
        elif at_symb == 'W':
            at_data[at_idx] = {
                'symb': 'W',
                'name': 'Tungsten',
                'num': 74,
                'mass': 183.8400,
                'rcov': (137, 120, 115),
                'rvdw': 1.570,
                'rvis': 0.685,
                'rgb': (33, 148, 214),
            }
        elif at_symb == 'Re':
            at_data[at_idx] = {
                'symb': 'Re',
                'name': 'Rhenium',
                'num': 75,
                'mass': 186.2070,
                'rcov': (131, 119, 110),
                'rvdw': 1.550,
                'rvis': 0.675,
                'rgb': (38, 125, 171),
            }
        elif at_symb == 'Os':
            at_data[at_idx] = {
                'symb': 'Os',
                'name': 'Osmium',
                'num': 76,
                'mass': 190.2000,
                'rcov': (129, 116, 109),
                'rvdw': 1.570,
                'rvis': 0.685,
                'rgb': (38, 102, 150),
            }
        elif at_symb == 'Ir':
            at_data[at_idx] = {
                'symb': 'Ir',
                'name': 'Iridium',
                'num': 77,
                'mass': 192.2200,
                'rcov': (122, 115, 107),
                'rvdw': 1.520,
                'rvis': 0.660,
                'rgb': (23, 84, 135),
            }
        elif at_symb == 'Pt':
            at_data[at_idx] = {
                'symb': 'Pt',
                'name': 'Platinum',
                'num': 78,
                'mass': 195.0900,
                'rcov': (123, 112, 110),
                'rvdw': 1.700,
                'rvis': 0.750,
                'rgb': (208, 208, 224),
            }
        elif at_symb == 'Au':
            at_data[at_idx] = {
                'symb': 'Au',
                'name': 'Gold',
                'num': 79,
                'mass': 196.9665,
                'rcov': (124, 121, 123),
                'rvdw': 1.700,
                'rvis': 0.750,
                'rgb': (255, 209, 35),
            }
        elif at_symb == 'Hg':
            at_data[at_idx] = {
                'symb': 'Hg',
                'name': 'Mercury',
                'num': 80,
                'mass': 200.5900,
                'rcov': (133, 142, None),
                'rvdw': 1.900,
                'rvis': 0.850,
                'rgb': (184, 184, 208),
            }
        elif at_symb == 'Tl':
            at_data[at_idx] = {
                'symb': 'Tl',
                'name': 'Thallium',
                'num': 81,
                'mass': 204.3700,
                'rcov': (144, 142, 150),
                'rvdw': 1.750,
                'rvis': 0.775,
                'rgb': (166, 84, 77),
            }
        elif at_symb == 'Pb':
            at_data[at_idx] = {
                'symb': 'Pb',
                'name': 'Lead',
                'num': 82,
                'mass': 207.2000,
                'rcov': (144, 135, 137),
                'rvdw': 1.740,
                'rvis': 0.770,
                'rgb': (87, 89, 97),
            }
        elif at_symb == 'Bi':
            at_data[at_idx] = {
                'symb': 'Bi',
                'name': 'Bismuth',
                'num': 83,
                'mass': 208.9804,
                'rcov': (151, 141, 135),
                'rvdw': 1.740,
                'rvis': 0.770,
                'rgb': (158, 79, 181),
            }
        elif at_symb == 'Po':
            at_data[at_idx] = {
                'symb': 'Po',
                'name': 'Polonium',
                'num': 84,
                'mass': None,
                'rcov': (145, 135, 129),
                'rvdw': 1.880,
                'rvis': 0.840,
                'rgb': (171, 92, 0),
            }
        elif at_symb == 'At':
            at_data[at_idx] = {
                'symb': 'At',
                'name': 'Astatine',
                'num': 85,
                'mass': None,
                'rcov': (147, 138, 138),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (117, 79, 69),
            }
        elif at_symb == 'Rn':
            at_data[at_idx] = {
                'symb': 'Rn',
                'name': 'Radon',
                'num': 86,
                'mass': None,
                'rcov': (142, 145, 133),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (66, 130, 150),
            }
        elif at_symb == 'Fr':
            at_data[at_idx] = {
                'symb': 'Fr',
                'name': 'Francium',
                'num': 87,
                'mass': None,
                'rcov': (223, 218, None),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (66, 0, 102),
            }
        elif at_symb == 'Ra':
            at_data[at_idx] = {
                'symb': 'Ra',
                'name': 'Radium',
                'num': 88,
                'mass': None,
                'rcov': (201, 173, 159),
                'rvdw': 2.100,
                'rvis': 0.950,
                'rgb': (0, 125, 0),
            }
        elif at_symb == 'Ac':
            at_data[at_idx] = {
                'symb': 'Ac',
                'name': 'Actinium',
                'num': 89,
                'mass': None,
                'rcov': (186, 153, 140),
                'rvdw': 2.080,
                'rvis': 0.940,
                'rgb': (112, 171, 250),
            }
        elif at_symb == 'Th':
            at_data[at_idx] = {
                'symb': 'Th',
                'name': 'Thorium',
                'num': 90,
                'mass': 232.038,
                'rcov': (175, 143, 136),
                'rvdw': 1.990,
                'rvis': 0.895,
                'rgb': (0, 186, 255),
            }
        elif at_symb == 'Pa':
            at_data[at_idx] = {
                'symb': 'Pa',
                'name': 'Protactinium',
                'num': 91,
                'mass': 231.0359,
                'rcov': (169, 138, 129),
                'rvdw': 1.810,
                'rvis': 0.805,
                'rgb': (0, 161, 255),
            }
        elif at_symb == 'U':
            at_data[at_idx] = {
                'symb': 'U',
                'name': 'Uranium',
                'num': 92,
                'mass': 238.0289,
                'rcov': (170, 134, 118),
                'rvdw': 1.780,
                'rvis': 0.790,
                'rgb': (0, 143, 255),
            }
        elif at_symb == 'Np':
            at_data[at_idx] = {
                'symb': 'Np',
                'name': 'Neptunium',
                'num': 93,
                'mass': None,
                'rcov': (171, 133, 116),
                'rvdw': 1.750,
                'rvis': 0.775,
                'rgb': (0, 128, 255),
            }
        elif at_symb == 'Pu':
            at_data[at_idx] = {
                'symb': 'Pu',
                'name': 'Plutonium',
                'num': 94,
                'mass': None,
                'rcov': (172, 135, None),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (0, 107, 255),
            }
        elif at_symb == 'Am':
            at_data[at_idx] = {
                'symb': 'Am',
                'name': 'Americium',
                'num': 95,
                'mass': None,
                'rcov': (166, 135, None),
                'rvdw': 1.710,
                'rvis': 0.755,
                'rgb': (84, 92, 242),
            }
        elif at_symb == 'Cm':
            at_data[at_idx] = {
                'symb': 'Cm',
                'name': 'Curium',
                'num': 96,
                'mass': None,
                'rcov': (166, 136, None),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (120, 92, 227),
            }
        elif at_symb == 'Bk':
            at_data[at_idx] = {
                'symb': 'Bk',
                'name': 'Berkelium',
                'num': 97,
                'mass': None,
                'rcov': (168, 139, None),
                'rvdw': 0.200,
                'rvis': 1.000,
                'rgb': (138, 79, 227),
            }
        elif at_symb == 'Cf':
            at_data[at_idx] = {
                'symb': 'Cf',
                'name': 'Californium',
                'num': 98,
                'mass': None,
                'rcov': (168, 140, None),
                'rvdw': 1.730,
                'rvis': 0.765,
                'rgb': (161, 54, 212),
            }
        elif at_symb == 'Es':
            at_data[at_idx] = {
                'symb': 'Es',
                'name': 'Einsteinium',
                'num': 99,
                'mass': None,
                'rcov': (165, 140, None),
                'rvdw': 0.100,
                'rvis': 0.100,
                'rgb': (179, 31, 212),
            }
        elif at_symb == 'Fm':
            at_data[at_idx] = {
                'symb': 'Fm',
                'name': 'Fermium',
                'num': 100,
                'mass': None,
                'rcov': (167, None, None),
                'rvdw': 0.200,
                'rvis': 0.900,
                'rgb': (179, 31, 186),
            }
        elif at_symb == 'Md':
            at_data[at_idx] = {
                'symb': 'Md',
                'name': 'Mendelevium',
                'num': 101,
                'mass': None,
                'rcov': (173, 139, None),
                'rvdw': None,
                'rvis': None,
                'rgb': (179, 13, 166),
            }
        elif at_symb == 'No':
            at_data[at_idx] = {
                'symb': 'No',
                'name': 'Nobelium',
                'num': 102,
                'mass': None,
                'rcov': (176, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': (189, 13, 135),
            }
        elif at_symb == 'Lr':
            at_data[at_idx] = {
                'symb': 'Lr',
                'name': 'Lawrencium',
                'num': 103,
                'mass': None,
                'rcov': (161, 141, None),
                'rvdw': None,
                'rvis': None,
                'rgb': (199, 0, 102),
            }
        elif at_symb == 'Rf':
            at_data[at_idx] = {
                'symb': 'Rf',
                'name': 'Rutherfordium',
                'num': 104,
                'mass': None,
                'rcov': (157, 140, 131),
                'rvdw': None,
                'rvis': None,
                'rgb': (204, 0, 89),
            }
        elif at_symb == 'Db':
            at_data[at_idx] = {
                'symb': 'Db',
                'name': 'Dubnium',
                'num': 105,
                'mass': None,
                'rcov': (149, 136, 126),
                'rvdw': None,
                'rvis': None,
                'rgb': (209, 0, 79),
            }
        elif at_symb == 'Sg':
            at_data[at_idx] = {
                'symb': 'Sg',
                'name': 'Seaborgium',
                'num': 106,
                'mass': None,
                'rcov': (143, 128, 121),
                'rvdw': None,
                'rvis': None,
                'rgb': (217, 0, 69),
            }
        elif at_symb == 'Bh':
            at_data[at_idx] = {
                'symb': 'Bh',
                'name': 'Bohrium',
                'num': 107,
                'mass': None,
                'rcov': (141, 128, 119),
                'rvdw': None,
                'rvis': None,
                'rgb': (224, 0, 56),
            }
        elif at_symb == 'Hs':
            at_data[at_idx] = {
                'symb': 'Hs',
                'name': 'Hassium',
                'num': 108,
                'mass': None,
                'rcov': (134, 125, 118),
                'rvdw': None,
                'rvis': None,
                'rgb': (230, 0, 46),
            }
        elif at_symb == 'Mt':
            at_data[at_idx] = {
                'symb': 'Mt',
                'name': 'Meitnerium',
                'num': 109,
                'mass': None,
                'rcov': (129, 125, 118),
                'rvdw': None,
                'rvis': None,
                'rgb': (235, 0, 38),
            }
        elif at_symb == 'Ds':
            at_data[at_idx] = {
                'symb': 'Ds',
                'name': 'Darmstadtium',
                'num': 110,
                'mass': None,
                'rcov': (128, 116, 112),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Rg':
            at_data[at_idx] = {
                'symb': 'Rg',
                'name': 'Roentgenium',
                'num': 111,
                'mass': None,
                'rcov': (121, 116, 118),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Cn':
            at_data[at_idx] = {
                'symb': 'Cn',
                'name': 'Copernicium',
                'num': 112,
                'mass': None,
                'rcov': (122, 116, 118),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Nh':
            at_data[at_idx] = {
                'symb': 'Nh',
                'name': 'Nihonium',
                'num': 113,
                'mass': None,
                'rcov': (136, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Fl':
            at_data[at_idx] = {
                'symb': 'Fl',
                'name': 'Flerovium',
                'num': 114,
                'mass': None,
                'rcov': (143, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Mc':
            at_data[at_idx] = {
                'symb': 'Mc',
                'name': 'Moscovium',
                'num': 115,
                'mass': None,
                'rcov': (162, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Lv':
            at_data[at_idx] = {
                'symb': 'Lv',
                'name': 'Livermorium',
                'num': 116,
                'mass': None,
                'rcov': (175, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Ts':
            at_data[at_idx] = {
                'symb': 'Ts',
                'name': 'Tennessine',
                'num': 117,
                'mass': None,
                'rcov': (165, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Og':
            at_data[at_idx] = {
                'symb': 'Og',
                'name': 'Oganesson',
                'num': 118,
                'mass': None,
                'rcov': (157, None, None),
                'rvdw': None,
                'rvis': None,
                'rgb': None,
            }
        else:
            raise KeyError('Atomic symbol not recognized')
    return at_data
