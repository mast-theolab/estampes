"""Toolbox for atoms

Module providing tools to manipulate data relative to atoms.

Methods
-------
convert_labsymb
    Converts between atomic symbols and labels
convert_atoms_id
    Converts between atoms identifiers.
"""

import typing as tp

from estampes.data.atom import ELEMENTS


# ==============
# MODULE METHODS
# ==============

def convert_labsymb(to_symb: bool,
                    *atoms: tp.List[tp.Union[str, int]]
                    ) -> tp.Union[tp.List[str], tp.List[int]]:
    """Converts between atomic symbols and labels

    Converts to atomic symbols (`to_symb`=True) or labels (False).
    The function is voluntarily very permissive accepting mixed data sets.

    Parameters
    ----------
    to_symb
        Converts to atomic symbols (True) or labels (False).
    *atoms
        List of atomic symbols and/or labels.

    Returns
    -------
    list
        List of atomic symbols or labels.

    Raises
    ------
    ValueError
        Missing atom(s).
    IndexError
        Incorrect atom label/symbol.
    """
    if not atoms:
        raise ValueError('No atoms')
    new_list = []
    for atom in atoms:
        if isinstance(atom, int):
            if 1 <= atom <= len(ELEMENTS):
                if to_symb:
                    new_list.append(ELEMENTS[atom])
                else:
                    new_list.append(atom)
            else:
                raise IndexError('Incorrect label')
        elif isinstance(atom, str):
            if atom.isdecimal():
                if 1 <= int(atom) <= len(ELEMENTS):
                    if to_symb:
                        new_list.append(ELEMENTS[int(atom)])
                    else:
                        new_list.append(int(atom))
                else:
                    raise IndexError('Incorrect label')
            else:
                if atom.capitalize() in ELEMENTS[1:]:
                    if to_symb:
                        new_list.append(atom)
                    else:
                        new_list.append(ELEMENTS.index(atom))
                else:
                    raise IndexError('Incorrect symbol')
        else:
            raise IndexError('Unrecognized atomic label/symbol')

    if len(new_list) == 1:
        return new_list[0]
    else:
        return new_list


def convert_atoms_id(list_atoms: tp.List[tp.Union[int, str]],
                     to: tp.Optional[str] = 'atnum'
                     ) -> tp.List[tp.Union[int, str]]:
    """Converts between atoms identifiers.

    Converts between atoms numbers (`atnum`) and labels (`atlab`).

    Attributes
    ----------
    list_atoms
        List of atoms numbers or labels.
    to
        Identifier to which the list of atoms should be converted.

    Returns
    -------
    list
        List of atoms with the chosen type of identifers.
    """
    atoms = []
    if to.lower() not in ('atnum', 'atlab'):
        raise ValueError('Unsupported type of identifier')
    if to == 'atnum':
        for atom in list_atoms:
            if isinstance(atom, int):
                atoms.append(atom)
            else:
                atoms.append(ELEMENTS.index(atom.title()))
    else:
        for atom in list_atoms:
            if isinstance(atom, str):
                atoms.append(atom.title())
            elif isinstance(atom, int):
                atoms.append(ELEMENTS[atom])
            else:
                fmt = 'Unrecognized atomic identifier {}'
                raise KeyError(fmt.format(atom))
