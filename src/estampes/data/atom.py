"""Module on atoms-related data.

This module provides basic data related to atoms.

Attributes
----------
ELEMENTS : list
    atomic number -> atomic Label list

"""

from estampes.base import AtDatType, AtLabType


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

def atomic_data(*atoms: AtLabType) -> AtDatType:
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
        Covalent radius (in Ang), as a tuple (single, double, triple).
    `rvdw`
        Van der Waals radius (in Ang), as dictionary, with keys:
        "bondi64"
          A. Bondi, J. Phys. Chem. A 1964 (68) 441
          https://doi.org/10.1021/j100785a001
        "truhlar09:"
          M. Mantina, A.C. Chamberlin, R. Valero, C.J. Cramer,
          D.G. Truhlar, J. Phys. Chem. A 2009 (113) 5809.
          https://doi.org/10.1021/jp8111556
          Extension of Bondi's set with some additional atoms from
          main group, but some values from Bondi were ignored.
        "alvarez13"
          S. Alvarez, Dalt. Trans. 2013 (42) 8617
          https://dx.doi.org/10.1039/c3dt50599e
          Statistical analysis from Cambridge Structure Database
        "rahm16"
          M. Rahm, R. Hoffmann, N.W. Ashcroft,
          Chem. Eur. J. 2016 (22) 14625.
          https://doi.org/10.1002/chem.201602949
          Radii built by considering a threshold in density of
          0.001 e.bohr^-3, computed at the PBE0/ANO-RCC level
        "truhlar_ext"
          Extended basis set considering values of bondi64 if not
          provided in original work of Trular and coworkers.
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
                'rcov': (0.3200, None, None),
                'rvdw': {
                    'bondi64': 1.2000,
                    'truhlar09': 1.1000,
                    'alvarez13': 1.2000,
                    'truhlar_ext': 1.1000,
                },
                'rvis': 0.2000,
                'rgb': (255, 255, 255),
            }
        elif at_symb == 'He':
            at_data[at_idx] = {
                'symb': 'He',
                'name': 'Helium',
                'num': 2,
                'mass': 4.00260,
                'rcov': (0.4600, None, None),
                'rvdw': {
                    'bondi64': 1.4000,
                    'truhlar09': 1.4000,
                    'alvarez13': 1.4300,
                    'truhlar_ext': 1.4000,
                },
                'rvis': 0.2860,
                'rgb': (217, 255, 255),
            }
        elif at_symb == 'Li':
            at_data[at_idx] = {
                'symb': 'Li',
                'name': 'Lithium',
                'num': 3,
                'mass': 6.94000,
                'rcov': (1.3300, 1.2400, None),
                'rvdw': {
                    'bondi64': 1.8100,
                    'truhlar09': 1.8100,
                    'alvarez13': 2.1200,
                    'truhlar_ext': 1.8100,
                },
                'rvis': 0.3400,
                'rgb': (204, 128, 255),
            }
        elif at_symb == 'Be':
            at_data[at_idx] = {
                'symb': 'Be',
                'name': 'Beryllium',
                'num': 4,
                'mass': 9.01218,
                'rcov': (1.0200, 0.9000, 0.8400),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 1.5300,
                    'alvarez13': 1.9800,
                    'truhlar_ext': 1.5300,
                },
                'rvis': 0.5890,
                'rgb': (194, 255, 0),
            }
        elif at_symb == 'B':
            at_data[at_idx] = {
                'symb': 'B',
                'name': 'Boron',
                'num': 5,
                'mass': 10.81000,
                'rcov': (0.8500, 0.7800, 0.7300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 1.9200,
                    'alvarez13': 1.9100,
                    'truhlar_ext': 1.9200,
                },
                'rvis': 0.4150,
                'rgb': (255, 181, 181),
            }
        elif at_symb == 'C':
            at_data[at_idx] = {
                'symb': 'C',
                'name': 'Carbon',
                'num': 6,
                'mass': 12.01100,
                'rcov': (0.7500, 0.6700, 0.6000),
                'rvdw': {
                    'bondi64': 1.7000,
                    'truhlar09': 1.7000,
                    'alvarez13': 1.7700,
                    'truhlar_ext': 1.7000,
                },
                'rvis': 0.4000,
                'rgb': (144, 144, 144),
            }
        elif at_symb == 'N':
            at_data[at_idx] = {
                'symb': 'N',
                'name': 'Nitrogen',
                'num': 7,
                'mass': 14.00670,
                'rcov': (0.7100, 0.6000, 0.5400),
                'rvdw': {
                    'bondi64': 1.5500,
                    'truhlar09': 1.5500,
                    'alvarez13': 1.6600,
                    'truhlar_ext': 1.5500,
                },
                'rvis': 0.4000,
                'rgb': (48, 80, 248),
            }
        elif at_symb == 'O':
            at_data[at_idx] = {
                'symb': 'O',
                'name': 'Oxygen',
                'num': 8,
                'mass': 15.99940,
                'rcov': (0.6300, 0.5700, 0.5300),
                'rvdw': {
                    'bondi64': 1.5200,
                    'truhlar09': 1.5200,
                    'alvarez13': 1.5000,
                    'truhlar_ext': 1.5200,
                },
                'rvis': 0.4000,
                'rgb': (255, 13, 13),
            }
        elif at_symb == 'F':
            at_data[at_idx] = {
                'symb': 'F',
                'name': 'Fluorine',
                'num': 9,
                'mass': 18.99840,
                'rcov': (0.6400, 0.5900, 0.5300),
                'rvdw': {
                    'bondi64': 1.4700,
                    'truhlar09': 1.4700,
                    'alvarez13': 1.4600,
                    'truhlar_ext': 1.4700,
                },
                'rvis': 0.3200,
                'rgb': (144, 224, 80),
            }
        elif at_symb == 'Ne':
            at_data[at_idx] = {
                'symb': 'Ne',
                'name': 'Neon',
                'num': 10,
                'mass': 20.17900,
                'rcov': (0.6700, 0.9600, None),
                'rvdw': {
                    'bondi64': 1.5400,
                    'truhlar09': 1.5400,
                    'alvarez13': 1.5800,
                    'truhlar_ext': 1.5400,
                },
                'rvis': 0.4230,
                'rgb': (179, 227, 245),
            }
        elif at_symb == 'Na':
            at_data[at_idx] = {
                'symb': 'Na',
                'name': 'Sodium',
                'num': 11,
                'mass': 22.98977,
                'rcov': (1.5500, 1.6000, None),
                'rvdw': {
                    'bondi64': 2.2700,
                    'truhlar09': 2.2700,
                    'alvarez13': 2.5000,
                    'truhlar_ext': 2.2700,
                },
                'rvis': 0.4850,
                'rgb': (171, 92, 242),
            }
        elif at_symb == 'Mg':
            at_data[at_idx] = {
                'symb': 'Mg',
                'name': 'Magnesium',
                'num': 12,
                'mass': 24.30500,
                'rcov': (1.3900, 1.3200, 1.2700),
                'rvdw': {
                    'bondi64': 1.7300,
                    'truhlar09': 1.7300,
                    'alvarez13': 2.5100,
                    'truhlar_ext': 1.7300,
                },
                'rvis': 0.5500,
                'rgb': (138, 255, 0),
            }
        elif at_symb == 'Al':
            at_data[at_idx] = {
                'symb': 'Al',
                'name': 'Aluminium',
                'num': 13,
                'mass': 26.98154,
                'rcov': (1.2600, 1.1300, 1.1100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 1.8400,
                    'alvarez13': 2.2500,
                    'truhlar_ext': 1.8400,
                },
                'rvis': 0.6750,
                'rgb': (191, 166, 166),
            }
        elif at_symb == 'Si':
            at_data[at_idx] = {
                'symb': 'Si',
                'name': 'Silicon',
                'num': 14,
                'mass': 28.08550,
                'rcov': (1.1600, 1.0700, 1.0200),
                'rvdw': {
                    'bondi64': 2.1000,
                    'truhlar09': 2.1000,
                    'alvarez13': 2.1900,
                    'truhlar_ext': 2.1000,
                },
                'rvis': 0.6000,
                'rgb': (240, 200, 160),
            }
        elif at_symb == 'P':
            at_data[at_idx] = {
                'symb': 'P',
                'name': 'Phosphorus',
                'num': 15,
                'mass': 30.97376,
                'rcov': (1.1100, 1.0200, 0.9400),
                'rvdw': {
                    'bondi64': 1.8000,
                    'truhlar09': 1.8000,
                    'alvarez13': 1.9000,
                    'truhlar_ext': 1.8000,
                },
                'rvis': 0.5250,
                'rgb': (255, 128, 0),
            }
        elif at_symb == 'S':
            at_data[at_idx] = {
                'symb': 'S',
                'name': 'Sulfur',
                'num': 16,
                'mass': 32.06000,
                'rcov': (1.0300, 0.9400, 0.9500),
                'rvdw': {
                    'bondi64': 1.8000,
                    'truhlar09': 1.8000,
                    'alvarez13': 1.8900,
                    'truhlar_ext': 1.8000,
                },
                'rvis': 0.5100,
                'rgb': (255, 255, 48),
            }
        elif at_symb == 'Cl':
            at_data[at_idx] = {
                'symb': 'Cl',
                'name': 'Chlorine',
                'num': 17,
                'mass': 35.45300,
                'rcov': (0.9900, 0.9500, 0.9300),
                'rvdw': {
                    'bondi64': 1.7500,
                    'truhlar09': 1.7500,
                    'alvarez13': 1.8200,
                    'truhlar_ext': 1.7500,
                },
                'rvis': 0.4950,
                'rgb': (31, 240, 31),
            }
        elif at_symb == 'Ar':
            at_data[at_idx] = {
                'symb': 'Ar',
                'name': 'Argon',
                'num': 18,
                'mass': 39.94800,
                'rcov': (0.9600, 1.0700, 0.9600),
                'rvdw': {
                    'bondi64': 1.8800,
                    'truhlar09': 1.8800,
                    'alvarez13': 1.8300,
                    'truhlar_ext': 1.8800,
                },
                'rvis': 0.5080,
                'rgb': (128, 209, 227),
            }
        elif at_symb == 'K':
            at_data[at_idx] = {
                'symb': 'K',
                'name': 'Potassium',
                'num': 19,
                'mass': 39.09830,
                'rcov': (1.9600, 1.9300, None),
                'rvdw': {
                    'bondi64': 2.7500,
                    'truhlar09': 2.7500,
                    'alvarez13': 2.7300,
                    'truhlar_ext': 2.7500,
                },
                'rvis': 0.6650,
                'rgb': (143, 64, 212),
            }
        elif at_symb == 'Ca':
            at_data[at_idx] = {
                'symb': 'Ca',
                'name': 'Calcium',
                'num': 20,
                'mass': 40.08000,
                'rcov': (1.7100, 1.4700, 1.3300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.3100,
                    'alvarez13': 2.6200,
                    'truhlar_ext': 2.3100,
                },
                'rvis': 0.4950,
                'rgb': (61, 255, 0),
            }
        elif at_symb == 'Sc':
            at_data[at_idx] = {
                'symb': 'Sc',
                'name': 'Scandium',
                'num': 21,
                'mass': 44.95590,
                'rcov': (1.4800, 1.1600, 1.1400),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.5800,
                    'truhlar_ext': None,
                },
                'rvis': 0.7350,
                'rgb': (230, 230, 230),
            }
        elif at_symb == 'Ti':
            at_data[at_idx] = {
                'symb': 'Ti',
                'name': 'Titanium',
                'num': 22,
                'mass': 47.90000,
                'rcov': (1.3600, 1.1700, 1.0800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4600,
                    'truhlar_ext': None,
                },
                'rvis': 0.7200,
                'rgb': (191, 194, 199),
            }
        elif at_symb == 'V':
            at_data[at_idx] = {
                'symb': 'V',
                'name': 'Vanadium',
                'num': 23,
                'mass': 50.94150,
                'rcov': (1.3400, 1.1200, 1.0600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4200,
                    'truhlar_ext': None,
                },
                'rvis': 0.6650,
                'rgb': (166, 166, 171),
            }
        elif at_symb == 'Cr':
            at_data[at_idx] = {
                'symb': 'Cr',
                'name': 'Chromium',
                'num': 24,
                'mass': 51.99600,
                'rcov': (1.2200, 1.1100, 1.0300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4500,
                    'truhlar_ext': None,
                },
                'rvis': 0.6750,
                'rgb': (138, 153, 199),
            }
        elif at_symb == 'Mn':
            at_data[at_idx] = {
                'symb': 'Mn',
                'name': 'Manganese',
                'num': 25,
                'mass': 54.93800,
                'rcov': (1.1900, 1.0500, 1.0300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4500,
                    'truhlar_ext': None,
                },
                'rvis': 0.6750,
                'rgb': (156, 122, 199),
            }
        elif at_symb == 'Fe':
            at_data[at_idx] = {
                'symb': 'Fe',
                'name': 'Iron',
                'num': 26,
                'mass': 55.84700,
                'rcov': (1.1600, 1.0900, 1.0200),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4400,
                    'truhlar_ext': None,
                },
                'rvis': 0.6700,
                'rgb': (224, 102, 51),
            }
        elif at_symb == 'Co':
            at_data[at_idx] = {
                'symb': 'Co',
                'name': 'Cobalt',
                'num': 27,
                'mass': 58.93320,
                'rcov': (1.1100, 1.0300, 0.9600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4000,
                    'truhlar_ext': None,
                },
                'rvis': 0.6150,
                'rgb': (240, 144, 160),
            }
        elif at_symb == 'Ni':
            at_data[at_idx] = {
                'symb': 'Ni',
                'name': 'Nickel',
                'num': 28,
                'mass': 58.71000,
                'rcov': (1.1000, 1.0100, 1.0100),
                'rvdw': {
                    'bondi64': 1.6300,
                    'truhlar09': None,
                    'alvarez13': 2.4000,
                    'truhlar_ext': 1.6300,
                },
                'rvis': 0.7500,
                'rgb': (80, 208, 80),
            }
        elif at_symb == 'Cu':
            at_data[at_idx] = {
                'symb': 'Cu',
                'name': 'Copper',
                'num': 29,
                'mass': 63.54600,
                'rcov': (1.1200, 1.1500, 1.2000),
                'rvdw': {
                    'bondi64': 1.4000,
                    'truhlar09': None,
                    'alvarez13': 2.3800,
                    'truhlar_ext': 1.4000,
                },
                'rvis': 0.7600,
                'rgb': (200, 128, 51),
            }
        elif at_symb == 'Zn':
            at_data[at_idx] = {
                'symb': 'Zn',
                'name': 'Zinc',
                'num': 30,
                'mass': 65.38000,
                'rcov': (1.1800, 1.2000, None),
                'rvdw': {
                    'bondi64': 1.3900,
                    'truhlar09': None,
                    'alvarez13': 2.3900,
                    'truhlar_ext': 1.3900,
                },
                'rvis': 0.7250,
                'rgb': (125, 128, 176),
            }
        elif at_symb == 'Ga':
            at_data[at_idx] = {
                'symb': 'Ga',
                'name': 'Gallium',
                'num': 31,
                'mass': 69.73500,
                'rcov': (1.2400, 1.1700, 1.2100),
                'rvdw': {
                    'bondi64': 1.8700,
                    'truhlar09': 1.8700,
                    'alvarez13': 2.3200,
                    'truhlar_ext': 1.8700,
                },
                'rvis': 0.6100,
                'rgb': (194, 143, 143),
            }
        elif at_symb == 'Ge':
            at_data[at_idx] = {
                'symb': 'Ge',
                'name': 'Germanium',
                'num': 32,
                'mass': 72.59000,
                'rcov': (1.2100, 1.1100, 1.1400),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.1100,
                    'alvarez13': 2.2900,
                    'truhlar_ext': 2.1100,
                },
                'rvis': 0.5850,
                'rgb': (102, 143, 143),
            }
        elif at_symb == 'As':
            at_data[at_idx] = {
                'symb': 'As',
                'name': 'Arsenic',
                'num': 33,
                'mass': 74.92160,
                'rcov': (1.2100, 1.1400, 1.0600),
                'rvdw': {
                    'bondi64': 1.8500,
                    'truhlar09': 1.8500,
                    'alvarez13': 1.8800,
                    'truhlar_ext': 1.8500,
                },
                'rvis': 0.6050,
                'rgb': (189, 128, 227),
            }
        elif at_symb == 'Se':
            at_data[at_idx] = {
                'symb': 'Se',
                'name': 'Selenium',
                'num': 34,
                'mass': 78.96000,
                'rcov': (1.1600, 1.0700, 1.0700),
                'rvdw': {
                    'bondi64': 1.9000,
                    'truhlar09': 1.9000,
                    'alvarez13': 1.8200,
                    'truhlar_ext': 1.9000,
                },
                'rvis': 0.6100,
                'rgb': (255, 161, 0),
            }
        elif at_symb == 'Br':
            at_data[at_idx] = {
                'symb': 'Br',
                'name': 'Bromine',
                'num': 35,
                'mass': 79.90400,
                'rcov': (1.1400, 1.0900, 1.1000),
                'rvdw': {
                    'bondi64': 1.8300,
                    'truhlar09': 1.8300,
                    'alvarez13': 1.8600,
                    'truhlar_ext': 1.8300,
                },
                'rvis': 0.6050,
                'rgb': (166, 41, 41),
            }
        elif at_symb == 'Kr':
            at_data[at_idx] = {
                'symb': 'Kr',
                'name': 'Krypton',
                'num': 36,
                'mass': 83.80000,
                'rcov': (1.1700, 1.2100, 1.0800),
                'rvdw': {
                    'bondi64': 2.0200,
                    'truhlar09': 2.0200,
                    'alvarez13': 2.2500,
                    'truhlar_ext': 2.0200,
                },
                'rvis': 0.5240,
                'rgb': (92, 184, 209),
            }
        elif at_symb == 'Rb':
            at_data[at_idx] = {
                'symb': 'Rb',
                'name': 'Rubidium',
                'num': 37,
                'mass': 85.46780,
                'rcov': (2.1000, 2.0200, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 3.0300,
                    'alvarez13': 3.2100,
                    'truhlar_ext': 3.0300,
                },
                'rvis': 0.7350,
                'rgb': (112, 46, 176),
            }
        elif at_symb == 'Sr':
            at_data[at_idx] = {
                'symb': 'Sr',
                'name': 'Strontium',
                'num': 38,
                'mass': 87.62000,
                'rcov': (1.8500, 1.5700, 1.3900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.4900,
                    'alvarez13': 2.8400,
                    'truhlar_ext': 2.4900,
                },
                'rvis': 0.5600,
                'rgb': (0, 255, 0),
            }
        elif at_symb == 'Y':
            at_data[at_idx] = {
                'symb': 'Y',
                'name': 'Yttrium',
                'num': 39,
                'mass': 88.90590,
                'rcov': (1.6300, 1.3000, 1.2400),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.7500,
                    'truhlar_ext': None,
                },
                'rvis': 0.8900,
                'rgb': (148, 255, 255),
            }
        elif at_symb == 'Zr':
            at_data[at_idx] = {
                'symb': 'Zr',
                'name': 'Zirconium',
                'num': 40,
                'mass': 91.22000,
                'rcov': (1.5400, 1.2700, 1.2100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.5200,
                    'truhlar_ext': None,
                },
                'rvis': 0.7800,
                'rgb': (148, 224, 224),
            }
        elif at_symb == 'Nb':
            at_data[at_idx] = {
                'symb': 'Nb',
                'name': 'Niobium',
                'num': 41,
                'mass': 92.90640,
                'rcov': (1.4700, 1.2500, 1.1600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.5600,
                    'truhlar_ext': None,
                },
                'rvis': 0.7400,
                'rgb': (115, 194, 201),
            }
        elif at_symb == 'Mo':
            at_data[at_idx] = {
                'symb': 'Mo',
                'name': 'Molybdenum',
                'num': 42,
                'mass': 95.94000,
                'rcov': (1.3800, 1.2100, 1.1300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4500,
                    'truhlar_ext': None,
                },
                'rvis': 0.7350,
                'rgb': (84, 181, 181),
            }
        elif at_symb == 'Tc':
            at_data[at_idx] = {
                'symb': 'Tc',
                'name': 'Technetium',
                'num': 43,
                'mass': 98.90620,
                'rcov': (1.2800, 1.2000, 1.1000),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4400,
                    'truhlar_ext': None,
                },
                'rvis': 0.6750,
                'rgb': (59, 158, 158),
            }
        elif at_symb == 'Ru':
            at_data[at_idx] = {
                'symb': 'Ru',
                'name': 'Ruthenium',
                'num': 44,
                'mass': 101.07000,
                'rcov': (1.2500, 1.1400, 1.0300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4600,
                    'truhlar_ext': None,
                },
                'rvis': 0.7000,
                'rgb': (36, 143, 143),
            }
        elif at_symb == 'Rh':
            at_data[at_idx] = {
                'symb': 'Rh',
                'name': 'Rhodium',
                'num': 45,
                'mass': 102.90550,
                'rcov': (1.2500, 1.1000, 1.0600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4400,
                    'truhlar_ext': None,
                },
                'rvis': 0.7250,
                'rgb': (10, 125, 140),
            }
        elif at_symb == 'Pd':
            at_data[at_idx] = {
                'symb': 'Pd',
                'name': 'Palladium',
                'num': 46,
                'mass': 106.40000,
                'rcov': (1.2000, 1.1700, 1.1200),
                'rvdw': {
                    'bondi64': 1.6300,
                    'truhlar09': None,
                    'alvarez13': 2.1500,
                    'truhlar_ext': 1.6300,
                },
                'rvis': 0.7500,
                'rgb': (0, 105, 133),
            }
        elif at_symb == 'Ag':
            at_data[at_idx] = {
                'symb': 'Ag',
                'name': 'Silver',
                'num': 47,
                'mass': 107.86800,
                'rcov': (1.2800, 1.3900, 1.3700),
                'rvdw': {
                    'bondi64': 1.7200,
                    'truhlar09': None,
                    'alvarez13': 2.5300,
                    'truhlar_ext': 1.7200,
                },
                'rvis': 0.7950,
                'rgb': (192, 192, 192),
            }
        elif at_symb == 'Cd':
            at_data[at_idx] = {
                'symb': 'Cd',
                'name': 'Cadmium',
                'num': 48,
                'mass': 112.41000,
                'rcov': (1.3600, 1.4400, None),
                'rvdw': {
                    'bondi64': 1.5800,
                    'truhlar09': None,
                    'alvarez13': 2.4900,
                    'truhlar_ext': 1.5800,
                },
                'rvis': 0.8450,
                'rgb': (255, 217, 143),
            }
        elif at_symb == 'In':
            at_data[at_idx] = {
                'symb': 'In',
                'name': 'Indium',
                'num': 49,
                'mass': 114.82000,
                'rcov': (1.4200, 1.3600, 1.4600),
                'rvdw': {
                    'bondi64': 1.9300,
                    'truhlar09': 1.9300,
                    'alvarez13': 2.4300,
                    'truhlar_ext': 1.9300,
                },
                'rvis': 0.8150,
                'rgb': (166, 117, 115),
            }
        elif at_symb == 'Sn':
            at_data[at_idx] = {
                'symb': 'Sn',
                'name': 'Tin',
                'num': 50,
                'mass': 118.69000,
                'rcov': (1.4000, 1.3000, 1.3200),
                'rvdw': {
                    'bondi64': 2.1700,
                    'truhlar09': 2.1700,
                    'alvarez13': 2.4200,
                    'truhlar_ext': 2.1700,
                },
                'rvis': 0.7300,
                'rgb': (102, 128, 128),
            }
        elif at_symb == 'Sb':
            at_data[at_idx] = {
                'symb': 'Sb',
                'name': 'Antimony',
                'num': 51,
                'mass': 121.75000,
                'rcov': (1.4000, 1.3300, 1.2700),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.0600,
                    'alvarez13': 2.4700,
                    'truhlar_ext': 2.0600,
                },
                'rvis': 0.7300,
                'rgb': (158, 99, 181),
            }
        elif at_symb == 'Te':
            at_data[at_idx] = {
                'symb': 'Te',
                'name': 'Tellurium',
                'num': 52,
                'mass': 127.60000,
                'rcov': (1.3600, 1.2800, 1.2100),
                'rvdw': {
                    'bondi64': 2.0600,
                    'truhlar09': 2.0600,
                    'alvarez13': 1.9900,
                    'truhlar_ext': 2.0600,
                },
                'rvis': 0.7350,
                'rgb': (212, 122, 0),
            }
        elif at_symb == 'I':
            at_data[at_idx] = {
                'symb': 'I',
                'name': 'Iodine',
                'num': 53,
                'mass': 126.90450,
                'rcov': (1.3300, 1.2900, 1.2500),
                'rvdw': {
                    'bondi64': 1.9800,
                    'truhlar09': 1.9800,
                    'alvarez13': 2.0400,
                    'truhlar_ext': 1.9800,
                },
                'rvis': 0.7000,
                'rgb': (148, 0, 148),
            }
        elif at_symb == 'Xe':
            at_data[at_idx] = {
                'symb': 'Xe',
                'name': 'Xenon',
                'num': 54,
                'mass': 131.30000,
                'rcov': (1.3100, 1.3500, 1.2200),
                'rvdw': {
                    'bondi64': 2.1600,
                    'truhlar09': 2.1600,
                    'alvarez13': 2.0600,
                    'truhlar_ext': 2.1600,
                },
                'rvis': 0.5770,
                'rgb': (66, 158, 176),
            }
        elif at_symb == 'Cs':
            at_data[at_idx] = {
                'symb': 'Cs',
                'name': 'Caesium',
                'num': 55,
                'mass': 132.90540,
                'rcov': (2.3200, 2.0900, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 3.4300,
                    'alvarez13': 3.4800,
                    'truhlar_ext': 3.4300,
                },
                'rvis': 0.8350,
                'rgb': (87, 23, 143),
            }
        elif at_symb == 'Ba':
            at_data[at_idx] = {
                'symb': 'Ba',
                'name': 'Barium',
                'num': 56,
                'mass': 137.33000,
                'rcov': (1.9600, 1.6100, 1.4900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.6800,
                    'alvarez13': 3.0300,
                    'truhlar_ext': 2.6800,
                },
                'rvis': 0.6700,
                'rgb': (0, 201, 0),
            }
        elif at_symb == 'La':
            at_data[at_idx] = {
                'symb': 'La',
                'name': 'Lanthanum',
                'num': 57,
                'mass': 138.90550,
                'rcov': (1.8000, 1.3900, 1.3900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9800,
                    'truhlar_ext': None,
                },
                'rvis': 0.9350,
                'rgb': (112, 212, 255),
            }
        elif at_symb == 'Ce':
            at_data[at_idx] = {
                'symb': 'Ce',
                'name': 'Cerium',
                'num': 58,
                'mass': 140.12000,
                'rcov': (1.6300, 1.3700, 1.3100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8800,
                    'truhlar_ext': None,
                },
                'rvis': 0.9150,
                'rgb': (255, 255, 199),
            }
        elif at_symb == 'Pr':
            at_data[at_idx] = {
                'symb': 'Pr',
                'name': 'Praseodymium',
                'num': 59,
                'mass': 140.90770,
                'rcov': (1.7600, 1.3800, 1.2800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9200,
                    'truhlar_ext': None,
                },
                'rvis': 0.9100,
                'rgb': (217, 255, 199),
            }
        elif at_symb == 'Nd':
            at_data[at_idx] = {
                'symb': 'Nd',
                'name': 'Neodymium',
                'num': 60,
                'mass': 144.24000,
                'rcov': (1.7400, 1.3200, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9500,
                    'truhlar_ext': None,
                },
                'rvis': 0.9050,
                'rgb': (199, 255, 199),
            }
        elif at_symb == 'Pm':
            at_data[at_idx] = {
                'symb': 'Pm',
                'name': 'Promethium',
                'num': 61,
                'mass': None,
                'rcov': (1.7300, 1.3500, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9000,
                    'truhlar_ext': None,
                },
                'rvis': 0.9000,
                'rgb': (163, 255, 199),
            }
        elif at_symb == 'Sm':
            at_data[at_idx] = {
                'symb': 'Sm',
                'name': 'Samarium',
                'num': 62,
                'mass': 150.36000,
                'rcov': (1.7200, 1.3400, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9000,
                    'truhlar_ext': None,
                },
                'rvis': 0.9000,
                'rgb': (143, 255, 199),
            }
        elif at_symb == 'Eu':
            at_data[at_idx] = {
                'symb': 'Eu',
                'name': 'Europium',
                'num': 63,
                'mass': 151.96000,
                'rcov': (1.6800, 1.3400, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8700,
                    'truhlar_ext': None,
                },
                'rvis': 0.9950,
                'rgb': (97, 255, 199),
            }
        elif at_symb == 'Gd':
            at_data[at_idx] = {
                'symb': 'Gd',
                'name': 'Gadolinium',
                'num': 64,
                'mass': 157.25000,
                'rcov': (1.6900, 1.3500, 1.3200),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8300,
                    'truhlar_ext': None,
                },
                'rvis': 0.8950,
                'rgb': (69, 255, 199),
            }
        elif at_symb == 'Tb':
            at_data[at_idx] = {
                'symb': 'Tb',
                'name': 'Terbium',
                'num': 65,
                'mass': 158.92530,
                'rcov': (1.6800, 1.3500, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.7900,
                    'truhlar_ext': None,
                },
                'rvis': 0.8800,
                'rgb': (48, 255, 199),
            }
        elif at_symb == 'Dy':
            at_data[at_idx] = {
                'symb': 'Dy',
                'name': 'Dysprosium',
                'num': 66,
                'mass': 162.50000,
                'rcov': (1.6700, 1.3300, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8700,
                    'truhlar_ext': None,
                },
                'rvis': 0.8750,
                'rgb': (31, 255, 199),
            }
        elif at_symb == 'Ho':
            at_data[at_idx] = {
                'symb': 'Ho',
                'name': 'Holmium',
                'num': 67,
                'mass': 164.93030,
                'rcov': (1.6600, 1.3300, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8100,
                    'truhlar_ext': None,
                },
                'rvis': 0.8700,
                'rgb': (0, 255, 156),
            }
        elif at_symb == 'Er':
            at_data[at_idx] = {
                'symb': 'Er',
                'name': 'Erbium',
                'num': 68,
                'mass': 167.26000,
                'rcov': (1.6500, 1.3300, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8300,
                    'truhlar_ext': None,
                },
                'rvis': 0.8650,
                'rgb': (0, 230, 117),
            }
        elif at_symb == 'Tm':
            at_data[at_idx] = {
                'symb': 'Tm',
                'name': 'Thulium',
                'num': 69,
                'mass': 168.93420,
                'rcov': (1.6400, 1.3100, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.7900,
                    'truhlar_ext': None,
                },
                'rvis': 0.8600,
                'rgb': (0, 212, 82),
            }
        elif at_symb == 'Yb':
            at_data[at_idx] = {
                'symb': 'Yb',
                'name': 'Ytterbium',
                'num': 70,
                'mass': 173.05000,
                'rcov': (1.7000, 1.2900, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8000,
                    'truhlar_ext': None,
                },
                'rvis': 0.9700,
                'rgb': (0, 191, 56),
            }
        elif at_symb == 'Lu':
            at_data[at_idx] = {
                'symb': 'Lu',
                'name': 'Lutetium',
                'num': 71,
                'mass': 174.96700,
                'rcov': (1.6200, 1.3100, 1.3100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.7400,
                    'truhlar_ext': None,
                },
                'rvis': 0.8600,
                'rgb': (0, 171, 36),
            }
        elif at_symb == 'Hf':
            at_data[at_idx] = {
                'symb': 'Hf',
                'name': 'Hafnium',
                'num': 72,
                'mass': 178.49000,
                'rcov': (1.5200, 1.2800, 1.2200),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.6300,
                    'truhlar_ext': None,
                },
                'rvis': 0.7850,
                'rgb': (77, 194, 255),
            }
        elif at_symb == 'Ta':
            at_data[at_idx] = {
                'symb': 'Ta',
                'name': 'Tantalum',
                'num': 73,
                'mass': 180.94790,
                'rcov': (1.4600, 1.2600, 1.1900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.5300,
                    'truhlar_ext': None,
                },
                'rvis': 0.7150,
                'rgb': (77, 166, 255),
            }
        elif at_symb == 'W':
            at_data[at_idx] = {
                'symb': 'W',
                'name': 'Tungsten',
                'num': 74,
                'mass': 183.84000,
                'rcov': (1.3700, 1.2000, 1.1500),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.5700,
                    'truhlar_ext': None,
                },
                'rvis': 0.6850,
                'rgb': (33, 148, 214),
            }
        elif at_symb == 'Re':
            at_data[at_idx] = {
                'symb': 'Re',
                'name': 'Rhenium',
                'num': 75,
                'mass': 186.20700,
                'rcov': (1.3100, 1.1900, 1.1000),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4900,
                    'truhlar_ext': None,
                },
                'rvis': 0.6750,
                'rgb': (38, 125, 171),
            }
        elif at_symb == 'Os':
            at_data[at_idx] = {
                'symb': 'Os',
                'name': 'Osmium',
                'num': 76,
                'mass': 190.20000,
                'rcov': (1.2900, 1.1600, 1.0900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4800,
                    'truhlar_ext': None,
                },
                'rvis': 0.6850,
                'rgb': (38, 102, 150),
            }
        elif at_symb == 'Ir':
            at_data[at_idx] = {
                'symb': 'Ir',
                'name': 'Iridium',
                'num': 77,
                'mass': 192.22000,
                'rcov': (1.2200, 1.1500, 1.0700),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.4100,
                    'truhlar_ext': None,
                },
                'rvis': 0.6600,
                'rgb': (23, 84, 135),
            }
        elif at_symb == 'Pt':
            at_data[at_idx] = {
                'symb': 'Pt',
                'name': 'Platinum',
                'num': 78,
                'mass': 195.09000,
                'rcov': (1.2300, 1.1200, 1.1000),
                'rvdw': {
                    'bondi64': 1.7200,
                    'truhlar09': None,
                    'alvarez13': 2.2900,
                    'truhlar_ext': 1.7200,
                },
                'rvis': 0.7500,
                'rgb': (208, 208, 224),
            }
        elif at_symb == 'Au':
            at_data[at_idx] = {
                'symb': 'Au',
                'name': 'Gold',
                'num': 79,
                'mass': 196.96650,
                'rcov': (1.2400, 1.2100, 1.2300),
                'rvdw': {
                    'bondi64': 1.6600,
                    'truhlar09': None,
                    'alvarez13': 2.3200,
                    'truhlar_ext': 1.6600,
                },
                'rvis': 0.7500,
                'rgb': (255, 209, 35),
            }
        elif at_symb == 'Hg':
            at_data[at_idx] = {
                'symb': 'Hg',
                'name': 'Mercury',
                'num': 80,
                'mass': 200.59000,
                'rcov': (1.3300, 1.4200, None),
                'rvdw': {
                    'bondi64': 1.7000,
                    'truhlar09': None,
                    'alvarez13': 2.4500,
                    'truhlar_ext': 1.7000,
                },
                'rvis': 0.8500,
                'rgb': (184, 184, 208),
            }
        elif at_symb == 'Tl':
            at_data[at_idx] = {
                'symb': 'Tl',
                'name': 'Thallium',
                'num': 81,
                'mass': 204.37000,
                'rcov': (1.4400, 1.4200, 1.5000),
                'rvdw': {
                    'bondi64': 1.9600,
                    'truhlar09': 1.9600,
                    'alvarez13': 2.4700,
                    'truhlar_ext': 1.9600,
                },
                'rvis': 0.7750,
                'rgb': (166, 84, 77),
            }
        elif at_symb == 'Pb':
            at_data[at_idx] = {
                'symb': 'Pb',
                'name': 'Lead',
                'num': 82,
                'mass': 207.20000,
                'rcov': (1.4400, 1.3500, 1.3700),
                'rvdw': {
                    'bondi64': 2.0200,
                    'truhlar09': 2.0200,
                    'alvarez13': 2.6000,
                    'truhlar_ext': 2.0200,
                },
                'rvis': 0.7700,
                'rgb': (87, 89, 97),
            }
        elif at_symb == 'Bi':
            at_data[at_idx] = {
                'symb': 'Bi',
                'name': 'Bismuth',
                'num': 83,
                'mass': 208.98040,
                'rcov': (1.5100, 1.4100, 1.3500),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.0700,
                    'alvarez13': 2.5400,
                    'truhlar_ext': 2.0700,
                },
                'rvis': 0.7700,
                'rgb': (158, 79, 181),
            }
        elif at_symb == 'Po':
            at_data[at_idx] = {
                'symb': 'Po',
                'name': 'Polonium',
                'num': 84,
                'mass': None,
                'rcov': (1.4500, 1.3500, 1.2900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 1.9700,
                    'alvarez13': None,
                    'truhlar_ext': 1.9700,
                },
                'rvis': 0.8400,
                'rgb': (171, 92, 0),
            }
        elif at_symb == 'At':
            at_data[at_idx] = {
                'symb': 'At',
                'name': 'Astatine',
                'num': 85,
                'mass': None,
                'rcov': (1.4700, 1.3800, 1.3800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.0200,
                    'alvarez13': None,
                    'truhlar_ext': 2.0200,
                },
                'rvis': 1.0000,
                'rgb': (117, 79, 69),
            }
        elif at_symb == 'Rn':
            at_data[at_idx] = {
                'symb': 'Rn',
                'name': 'Radon',
                'num': 86,
                'mass': None,
                'rcov': (1.4200, 1.4500, 1.3300),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.2000,
                    'alvarez13': None,
                    'truhlar_ext': 2.2000,
                },
                'rvis': 1.0000,
                'rgb': (66, 130, 150),
            }
        elif at_symb == 'Fr':
            at_data[at_idx] = {
                'symb': 'Fr',
                'name': 'Francium',
                'num': 87,
                'mass': None,
                'rcov': (2.2300, 2.1800, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 3.4800,
                    'alvarez13': None,
                    'truhlar_ext': 3.4800,
                },
                'rvis': 1.0000,
                'rgb': (66, 0, 102),
            }
        elif at_symb == 'Ra':
            at_data[at_idx] = {
                'symb': 'Ra',
                'name': 'Radium',
                'num': 88,
                'mass': None,
                'rcov': (2.0100, 1.7300, 1.5900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': 2.8300,
                    'alvarez13': None,
                    'truhlar_ext': 2.8300,
                },
                'rvis': 0.9500,
                'rgb': (0, 125, 0),
            }
        elif at_symb == 'Ac':
            at_data[at_idx] = {
                'symb': 'Ac',
                'name': 'Actinium',
                'num': 89,
                'mass': None,
                'rcov': (1.8600, 1.5300, 1.4000),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8000,
                    'truhlar_ext': None,
                },
                'rvis': 0.9400,
                'rgb': (112, 171, 250),
            }
        elif at_symb == 'Th':
            at_data[at_idx] = {
                'symb': 'Th',
                'name': 'Thorium',
                'num': 90,
                'mass': 232.03800,
                'rcov': (1.7500, 1.4300, 1.3600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.9300,
                    'truhlar_ext': None,
                },
                'rvis': 0.8950,
                'rgb': (0, 186, 255),
            }
        elif at_symb == 'Pa':
            at_data[at_idx] = {
                'symb': 'Pa',
                'name': 'Protactinium',
                'num': 91,
                'mass': 231.03590,
                'rcov': (1.6900, 1.3800, 1.2900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8800,
                    'truhlar_ext': None,
                },
                'rvis': 0.8050,
                'rgb': (0, 161, 255),
            }
        elif at_symb == 'U':
            at_data[at_idx] = {
                'symb': 'U',
                'name': 'Uranium',
                'num': 92,
                'mass': 238.02890,
                'rcov': (1.7000, 1.3400, 1.1800),
                'rvdw': {
                    'bondi64': 1.8600,
                    'truhlar09': None,
                    'alvarez13': 2.7100,
                    'truhlar_ext': 1.8600,
                },
                'rvis': 0.7900,
                'rgb': (0, 143, 255),
            }
        elif at_symb == 'Np':
            at_data[at_idx] = {
                'symb': 'Np',
                'name': 'Neptunium',
                'num': 93,
                'mass': None,
                'rcov': (1.7100, 1.3300, 1.1600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8200,
                    'truhlar_ext': None,
                },
                'rvis': 0.7750,
                'rgb': (0, 128, 255),
            }
        elif at_symb == 'Pu':
            at_data[at_idx] = {
                'symb': 'Pu',
                'name': 'Plutonium',
                'num': 94,
                'mass': None,
                'rcov': (1.7200, 1.3500, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8100,
                    'truhlar_ext': None,
                },
                'rvis': 1.0000,
                'rgb': (0, 107, 255),
            }
        elif at_symb == 'Am':
            at_data[at_idx] = {
                'symb': 'Am',
                'name': 'Americium',
                'num': 95,
                'mass': None,
                'rcov': (1.6600, 1.3500, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.8300,
                    'truhlar_ext': None,
                },
                'rvis': 0.7550,
                'rgb': (84, 92, 242),
            }
        elif at_symb == 'Cm':
            at_data[at_idx] = {
                'symb': 'Cm',
                'name': 'Curium',
                'num': 96,
                'mass': None,
                'rcov': (1.6600, 1.3600, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 3.0500,
                    'truhlar_ext': None,
                },
                'rvis': 1.0000,
                'rgb': (120, 92, 227),
            }
        elif at_symb == 'Bk':
            at_data[at_idx] = {
                'symb': 'Bk',
                'name': 'Berkelium',
                'num': 97,
                'mass': None,
                'rcov': (1.6800, 1.3900, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 3.4000,
                    'truhlar_ext': None,
                },
                'rvis': 1.0000,
                'rgb': (138, 79, 227),
            }
        elif at_symb == 'Cf':
            at_data[at_idx] = {
                'symb': 'Cf',
                'name': 'Californium',
                'num': 98,
                'mass': None,
                'rcov': (1.6800, 1.4000, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 3.0500,
                    'truhlar_ext': None,
                },
                'rvis': 0.7650,
                'rgb': (161, 54, 212),
            }
        elif at_symb == 'Es':
            at_data[at_idx] = {
                'symb': 'Es',
                'name': 'Einsteinium',
                'num': 99,
                'mass': None,
                'rcov': (1.6500, 1.4000, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': 2.7000,
                    'truhlar_ext': None,
                },
                'rvis': 0.1000,
                'rgb': (179, 31, 212),
            }
        elif at_symb == 'Fm':
            at_data[at_idx] = {
                'symb': 'Fm',
                'name': 'Fermium',
                'num': 100,
                'mass': None,
                'rcov': (1.6700, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': 0.9000,
                'rgb': (179, 31, 186),
            }
        elif at_symb == 'Md':
            at_data[at_idx] = {
                'symb': 'Md',
                'name': 'Mendelevium',
                'num': 101,
                'mass': None,
                'rcov': (1.7300, 1.3900, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (179, 13, 166),
            }
        elif at_symb == 'No':
            at_data[at_idx] = {
                'symb': 'No',
                'name': 'Nobelium',
                'num': 102,
                'mass': None,
                'rcov': (1.7600, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (189, 13, 135),
            }
        elif at_symb == 'Lr':
            at_data[at_idx] = {
                'symb': 'Lr',
                'name': 'Lawrencium',
                'num': 103,
                'mass': None,
                'rcov': (1.6100, 1.4100, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (199, 0, 102),
            }
        elif at_symb == 'Rf':
            at_data[at_idx] = {
                'symb': 'Rf',
                'name': 'Rutherfordium',
                'num': 104,
                'mass': None,
                'rcov': (1.5700, 1.4000, 1.3100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (204, 0, 89),
            }
        elif at_symb == 'Db':
            at_data[at_idx] = {
                'symb': 'Db',
                'name': 'Dubnium',
                'num': 105,
                'mass': None,
                'rcov': (1.4900, 1.3600, 1.2600),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (209, 0, 79),
            }
        elif at_symb == 'Sg':
            at_data[at_idx] = {
                'symb': 'Sg',
                'name': 'Seaborgium',
                'num': 106,
                'mass': None,
                'rcov': (1.4300, 1.2800, 1.2100),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (217, 0, 69),
            }
        elif at_symb == 'Bh':
            at_data[at_idx] = {
                'symb': 'Bh',
                'name': 'Bohrium',
                'num': 107,
                'mass': None,
                'rcov': (1.4100, 1.2800, 1.1900),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (224, 0, 56),
            }
        elif at_symb == 'Hs':
            at_data[at_idx] = {
                'symb': 'Hs',
                'name': 'Hassium',
                'num': 108,
                'mass': None,
                'rcov': (1.3400, 1.2500, 1.1800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (230, 0, 46),
            }
        elif at_symb == 'Mt':
            at_data[at_idx] = {
                'symb': 'Mt',
                'name': 'Meitnerium',
                'num': 109,
                'mass': None,
                'rcov': (1.2900, 1.2500, 1.1800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': (235, 0, 38),
            }
        elif at_symb == 'Ds':
            at_data[at_idx] = {
                'symb': 'Ds',
                'name': 'Darmstadtium',
                'num': 110,
                'mass': None,
                'rcov': (1.2800, 1.1600, 1.1200),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Rg':
            at_data[at_idx] = {
                'symb': 'Rg',
                'name': 'Roentgenium',
                'num': 111,
                'mass': None,
                'rcov': (1.2100, 1.1600, 1.1800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Cn':
            at_data[at_idx] = {
                'symb': 'Cn',
                'name': 'Copernicium',
                'num': 112,
                'mass': None,
                'rcov': (1.2200, 1.1600, 1.1800),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Nh':
            at_data[at_idx] = {
                'symb': 'Nh',
                'name': 'Nihonium',
                'num': 113,
                'mass': None,
                'rcov': (1.3600, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Fl':
            at_data[at_idx] = {
                'symb': 'Fl',
                'name': 'Flerovium',
                'num': 114,
                'mass': None,
                'rcov': (1.4300, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Mc':
            at_data[at_idx] = {
                'symb': 'Mc',
                'name': 'Moscovium',
                'num': 115,
                'mass': None,
                'rcov': (1.6200, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Lv':
            at_data[at_idx] = {
                'symb': 'Lv',
                'name': 'Livermorium',
                'num': 116,
                'mass': None,
                'rcov': (1.7500, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Ts':
            at_data[at_idx] = {
                'symb': 'Ts',
                'name': 'Tennessine',
                'num': 117,
                'mass': None,
                'rcov': (1.6500, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        elif at_symb == 'Og':
            at_data[at_idx] = {
                'symb': 'Og',
                'name': 'Oganesson',
                'num': 118,
                'mass': None,
                'rcov': (1.5700, None, None),
                'rvdw': {
                    'bondi64': None,
                    'truhlar09': None,
                    'alvarez13': None,
                    'truhlar_ext': None,
                },
                'rvis': None,
                'rgb': None,
            }
        else:
            raise KeyError('Atomic symbol not recognized')
    return at_data
