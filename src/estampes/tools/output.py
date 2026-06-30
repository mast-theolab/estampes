"""Provide tools to facilitate some output operations.

The module provides some simple tools for common output operations.
"""
import re
import typing as tp
from collections.abc import Iterable, Sequence

from estampes.base import ArgumentError
from estampes.data.property import QBaseInfo


def fortran_fmt_D(number_s: float | Sequence[float], fmt: str) -> str:
    """Print simulating the double precision format.

    Simulate the printing a value `num` with a double precision format
    specification.
    Fortran double precision format has the form: Dw.d[Ee] with:

    - w: total width of the field.
    - d: number of decimal digits.
    - e: number of digits for the exponentiation.

    with the integer part always 0.

    If multiple values are provided, the resulting string is a
    concatenation of each one.

    Parameters
    ----------
    number_s
        Number to represent.
    fmt
        Fortran format.

    Returns
    -------
    str
        Representation of the content `number_s` with the proper format.
    """
    def num_to_repr(num: float) -> str:
        n_sign = 1 if num < 0 else 0
        if comps['e'] is not None:
            # sign + 0 + dot + decimal_part + E + sign_E + length_e
            w_min = n_sign + 2 + len_d + 2 + len_e
        else:
            # sign + 0 + dot + decimal_part + 4 (E+dd or +ddd)
            w_min = n_sign + 2 + len_d + 4
        if len_w < w_min:
            return len_w*'*'

        # Python represents number with integer part > 0 so we must correct d
        # by reducing it by -1
        py_fmt = f'{{:.{len_d-1}e}}'
        py_res = py_fmt.format(abs(num))
        p_num, p_exp = py_res.split('e')
        if num > 0:
            p_num = '0.' + p_num.replace('.', '')
        else:
            p_num = '-0.' + p_num.replace('.', '')
        if comps['e'] is not None:
            if len_e + 1 < len(p_exp):
                return len_w*'*'
            else:
                py_fmt = f'{{:+0{len_e}d}}'
                p_exp = 'D' + py_fmt.format(int(p_exp)+1)
        else:
            if len(p_exp) == 3:
                p_exp = f'D{int(p_exp)+1:+03d}'
            elif len(p_exp) == 4:
                p_exp = f'{int(p_exp)+1:+04d}'
            else:
                return len_w*'*'

        py_fmt = f'{{:>{len_w}s}}'
        return py_fmt.format(p_num+p_exp)

    parse = re.match(r'D(?P<w>\d+)\.(?P<d>\d+)(?:E(?P<e>\d+))?', fmt, re.I)
    if parse is None:
        raise ArgumentError('fmt', 'Incorrect Fortran format')

    comps = parse.groupdict()
    len_w = int(comps['w'])
    len_d = int(comps['d'])
    len_e = int(comps['e'] if comps['e'] is not None else 0)
    if 0 < len_e < 2:
        raise NotImplementedError('Small specs for D not implemented.')

    if isinstance(number_s, float):
        return num_to_repr(number_s)
    elif isinstance(number_s, (Sequence, Iterable)):
        repr_s = []
        for item in number_s:
            repr_s.append(num_to_repr(item))
        return ''.join(repr_s)
    else:
        print(type)
        raise TypeError('Unsupported representation of number_s')


def pstruct_to_labels(qinfo: QBaseInfo,
                      natoms: int | None = None,
                      linearize: bool = False,
                      sep_atom: str = ' ',
                      sep_xyz: str = '') -> str | tuple[tp.Any, ...]:
    """Construct sequence of labels from property structure information.

    Constructs a sequence of labels following the dimensions information
    or linearized, based on the property structure description.

    Parameters
    ----------
    qinfo
        Quantity information
    natoms
        Number of atoms, only relevant for properties dependent on atoms.
    linearize
        Linearize structure instead of keeping the assumed shape.
    sep_atom
        Separator character(s) to use between atoms or atom-coordinates.
        If not separator is provided, at least one space is added
        between atoms labels to ensure that they are readable.
    sep_xyz
        Separator to use between coordinates of different dimensions.

    Returns
    -------
    str or tuple
        Sequence of labels.  An empty string is returned for scalars.
    """
    def proc_atom_spec(struct: str, nat: int | None) -> tuple[tp.Any, ...]:
        if nat is None:
            raise ArgumentError('natoms', 'Missing number of atoms')
        if nat <= 0:
            raise ValueError(
                'Number of atoms cannot be negative or null')
        fmt_lab = f'{{:>{len(str(nat))}c}}'
        if struct == 'atom':
            return tuple((fmt_lab.format(i+1) for i in range(nat)))
        elif struct == 'atomLT':
            if not isinstance(sep_atom, str):
                raise TypeError('sep_atom must be a string')
            return tuple((fmt_lab.format(i+1) + sep_atom + fmt_lab.format(j+1)
                          for i in range(nat) for j in range(i+1)))
        else:
            raise ValueError('Unrecognized atom-based structure')

    xyz = ('x', 'y', 'z')
    comp_xyz = {
        'xyz': xyz,
        'xyz2D': (lab1+lab2 for lab1 in xyz for lab2 in xyz),
        'xyz3D': (lab1+lab2+lab3
                  for lab1 in xyz for lab2 in xyz for lab3 in xyz),
        'xyzLT': (lab1+lab2
                  for i, lab1 in enumerate(xyz, start=1)
                  for lab2 in xyz[:i]),
        'xyzLC': (lab1+lab2+lab3
                  for i, lab1 in enumerate(xyz, start=1)
                  for j, lab2 in enumerate(xyz[:i], start=1)
                  for lab3 in xyz[:j])
    }

    match (qinfo.form):
        case ('scalar'):
            return ''
        case str():
            if qinfo.form.startswith('atom'):
                return proc_atom_spec(qinfo.form, natoms)
            else:
                return tuple(comp_xyz[qinfo.form])
        case tuple() | list():
            comps = []
            is_xyz = []
            for item in qinfo.form:
                if isinstance(item, str):
                    if item.startswith('atom'):
                        comps.append(proc_atom_spec(item, natoms))
                        is_xyz.append(False)
                    else:
                        comps.append(comp_xyz[item])
                        is_xyz.append(True)
                elif isinstance(item, (tuple, list)):
                    # We assume there is only one max level of nesting
                    # So at most a structure tuple(tuple(scalar))
                    # which corresponds to most standard cases.
                    for subitem in item:
                        if subitem.startswith('atom'):
                            comps.append(proc_atom_spec(subitem, natoms))
                            is_xyz.append(False)
                        else:
                            comps.append(comp_xyz[subitem])
                            is_xyz.append(True)
            # We build list of components, now construct extended structure.
            # We could do a generic but in practice we will have 2 or 3
            # dimensions, so we make this explicit for simplicity.
            if any(is_xyz):
                if not isinstance(sep_xyz, str):
                    raise TypeError('String expected as separator between '
                                    + 'Cartesian components')
            if not all(is_xyz):
                if not isinstance(sep_atom, str):
                    raise TypeError('String expected as separator between '
                                    + 'atomic components')
            if linearize:
                match (len(comps)):
                    case 2:
                        sep = sep_xyz if is_xyz[0] else sep_atom
                        return tuple((lab1+sep+lab2
                                      for lab1 in comps[0]
                                      for lab2 in comps[1]))
                    case 3:
                        sep1 = sep_xyz if is_xyz[0] else sep_atom
                        sep2 = sep_xyz if is_xyz[1] else sep_atom
                        return tuple((lab1+sep1+lab2+sep2+lab3
                                      for lab1 in comps[0]
                                      for lab2 in comps[1]
                                      for lab3 in comps[2]))
                    case _:
                        raise ValueError('Unsupported component structure')
            else:
                match (len(comps)):
                    case 2:
                        sep = sep_xyz if is_xyz[0] else sep_atom
                        return tuple(((lab1+sep+lab2
                                       for lab1 in comps[0])
                                      for lab2 in comps[1]))
                    case 3:
                        sep1 = sep_xyz if is_xyz[0] else sep_atom
                        sep2 = sep_xyz if is_xyz[1] else sep_atom
                        return tuple((((lab1+sep1+lab2+sep2+lab3
                                        for lab1 in comps[0])
                                       for lab2 in comps[1])
                                      for lab3 in comps[2]))
                    case _:
                        raise ValueError('Unsupported component structure')


def sec_header(output: tp.IO, level: int, title: str, lead_spaces: int = 1):
    """Write section header in output.
    
    Writes a section header in the file object defined with output.
    
    Parameters
    ----------
    output
        File object where the section header is printed.
    level
        Hierarchical level of the section:
        
        - `-1`: Main title
        - `0`: Chapter
        - `1`: Section / Header1
        - `2`: Subsection / Header2
        - `3`: Subsubsection / Header3
        - `4`: Paragraph / Header4

    title
        Header title.
    lead_spaces
        Leading shift for the section.
    """
    ltitle = len(title.strip())
    if ltitle == 0:
        return
    
    if level == -1:
        length = max(76, ltitle+8)
        fmt = f'''\
{lead_spaces*' '}/{length*'-'}\\
{lead_spaces*' '}|{length*' '}|
{lead_spaces*' '}|{{:^{length}s}}|
{lead_spaces*' '}|{length*' '}|
{lead_spaces*' '}\\{length*'-'}/
'''
    elif level == 0:
        length = max(78, ltitle+4)
        fmt = f'''

{lead_spaces*' '}{length*'*'}

{lead_spaces*' '}{{:^{length}s}}

{lead_spaces*' '}{length*'*'}'''
    elif level == 1:
        fmt = f'''

{lead_spaces*' '}{{:{ltitle}}}
{lead_spaces*' '}{ltitle*'='}'''
    elif level == 2:
        fmt = f'''
{lead_spaces*' '}{{:{ltitle}s}}
{lead_spaces*' '}{ltitle*'-'}'''
    elif level == 3:
        fmt = f'''
{lead_spaces*' '}{{:{ltitle}s}}
{lead_spaces*' '}{ltitle*'^'}'''
    elif level >= 4:
        fmt = f'''
{lead_spaces*' '}### {{:{ltitle}^s}} ###'''
    else:
        raise ArgumentError('level', 'Unsupported header level')
    output.write(fmt.format(title.strip()) + '\n')
