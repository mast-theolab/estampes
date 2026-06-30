"""Program DERIVEUR - Shift geometry and build input module.

The module contains the routines related to the construction of
displaced geometries and the associated input files.
"""

import os
import re
import sys
import argparse
import typing as tp

import numpy as np
import numpy.typing as npt

from estampes.base import (
    AtsLabType, AtsMasType, InternalError, QuantityError)
from estampes.parser import DataFile
from estampes.data.atom import atomic_data
from estampes.tools.atom import convert_labsymb
from estampes.tools.char import convert_range_spec
from estampes.tools.math import nint
from estampes.tools.vib import orient_modes

from estampes.extras.derive_core import (
    INFO_PM_STEPS, KWORD_CHK, KWORD_DESC, KWORD_DIR, KWORD_DISP, KWORD_GEOM,
    KWORD_NUM, KWORD_STEP, KWORD_XYZ, QLABS_ATGEOM,
    check_fname_pattern, gen_fname_from_tmpl)


#
# TEMPLATE FILES
# ==============
def parse_template_file(fname: str | None = None) -> str:
    """Get template to generate content.

    Parses the content of a template file or generate a generic one.
    The parsing is very basic and only checks that the required
    placeholders are present.

    Parameters
    ----------
    fname
        Filename with the template information.

    Returns
    -------
    str
        Template for the content.

    Raises
    ------
    FileNotFoundError
        Provided file does not exist.
    KeyError
        Missing mandatory placeholders or redundant.
    ValueError
        File is too large.
    """
    # Maximum file size set to 32 MB for the template file
    max_size = 32 * 1024**2

    if fname is not None:
        if not os.path.exists(fname):
            raise FileNotFoundError('Template file does not exist')
        if os.path.getsize(fname) > max_size:
            raise ValueError('Template file is too large.')
        with open(fname, encoding='utf-8') as fobj:
            txt = fobj.read()
        n_places = txt.count(KWORD_GEOM)
        if n_places == 0:
            if KWORD_GEOM.lower() in txt.lower():
                print(f'{KWORD_GEOM} found with an incorrect case, fixing.')
                pattern = f'\\{KWORD_GEOM[:-1]}\\]'
                txt = re.sub(pattern, KWORD_GEOM, txt, flags=re.I)
            else:
                raise KeyError(KWORD_GEOM + ' not found in template file')
        elif n_places > 1:
            msg = 'Too many occurrences of {} found.  This is confusing.'
            raise KeyError(msg.format(KWORD_GEOM))
        n_places = txt.count(KWORD_DESC)
        if n_places == 0:
            if KWORD_DESC.lower() in txt.lower():
                print(f'{KWORD_DESC} found with an incorrect case, fixing.')
                pattern = f'\\{KWORD_DESC[:-1]}\\]'
                txt = re.sub(pattern, KWORD_DESC, txt, flags=re.I)
        txt_lcase = txt.lower()
        for kword in (KWORD_DISP, KWORD_STEP, KWORD_NUM, KWORD_XYZ, KWORD_DIR,
                      KWORD_CHK):
            num = txt_lcase.count(kword.lower())
            if num != txt.count(kword):
                print(f'Occurrences of {kword} found with incorrect case.',
                      'fixing.')
                pattern = rf'\{kword[:-1]}\]'
                txt = re.sub(pattern, kword, txt, flags=re.I)

    else:
        txt = f"""
Displaced geometry: {KWORD_DESC}

{KWORD_GEOM}
"""

    return txt


#
# SHIFT ROUTINES
# ==============
def build_geoms(at_crd_ref: npt.NDArray,
                at_symb: AtsLabType,
                at_mass: AtsMasType,
                displ_idx: tp.Iterable[int],
                displ_num: int,
                displ_step: float,
                displ_cart: bool,
                au2unit: float,
                tmpl_text: str,
                tmpl_files: str | None = None,
                ref_fname: str | None = None,
                displ_vecs: npt.NDArray | None = None):
    """Build displaced geometries.

    Builds displaced geometries from the reference coordinates.

    Parameters
    ----------
    at_crd_ref
        Coordinates of the reference geometry.
    at_symb
        Atomic symbols.
    at_mass
        Atomic masses.
    displ_idx
        Displacement indexes (atoms or vectors).
    displ_num
        Number of displacements along each coordinates.
    displ_step
        Displacement step, in the unit of the displacement.
    displ_cart
        If true, displacement is done along Cart. coordinates.
    au2unit
        Conversion factor from bohrs to output unit.
    tmpl_text
        Template for the content, as a single character string.
    tmpl_files
        Template of the file names.  If `None`, print in output.
    ref_fname
        Filename for the reference geometry.
    displ_vecs
        List of displacement vectors if not along Cart. coords.

    Notes
    -----
    * Checks should be done beforehand
    """
    fmt_geom = '{at:<10s}{xyz[0]:21.14f}{xyz[1]:21.14f}{xyz[2]:21.14f}'
    natoms = at_crd_ref.shape[0]

    # Let us check if all masses are consistent:
    at_data = atomic_data(*at_symb)
    at_lab = []
    for ia, mass in enumerate(at_mass):
        # We multiply by 2 to manage cases of fractional masses (Cl)
        if nint(mass*2) != nint(2*at_data[at_symb[ia]]['mass']):
            at_lab.append(f'{at_symb[ia]}(iso={nint(mass)})')
        else:
            at_lab.append(f'{at_symb[ia]}')

    if displ_cart:
        for i in displ_idx:
            for ix, axis in enumerate(('X', 'Y', 'Z')):
                displ = np.zeros(at_crd_ref.shape)
                for step, step_data in INFO_PM_STEPS.items():
                    for istep in range(1, displ_num+1):
                        displ[i, ix] = istep*step_data[1]*displ_step
                        crd = (at_crd_ref + displ)*au2unit
                        txt_keys = gen_fname_from_tmpl(
                            natoms, displ_num, i+1, istep, axis, step)
                        if tmpl_files is not None:
                            fname_inp = tmpl_files\
                                .replace(KWORD_DISP, txt_keys[KWORD_DISP]) \
                                .replace(KWORD_STEP, txt_keys[KWORD_STEP]) \
                                .replace(KWORD_NUM, txt_keys[KWORD_NUM]) \
                                .replace(KWORD_XYZ, txt_keys[KWORD_XYZ]) \
                                .replace(KWORD_DIR, txt_keys[KWORD_DIR])
                            fobj = open(fname_inp, 'w', encoding='utf-8')
                        else:
                            fobj = sys.stdout
                            fobj.write(72*'-'+'\n')
                            fname_inp = ''
                        txt_desc = f'Step {step_data[0]} along {axis} ' \
                            + f'axis from atom {i+1}'
                        txt_geom = '\n'.join(
                            [fmt_geom.format(at=at_lab[ia], xyz=crd[ia, :])
                             for ia in range(natoms)]
                        )
                        if fname_inp:
                            txt_chk = os.path.splitext(fname_inp)[0] + '.chk'
                        elif ref_fname is not None:
                            txt_chk = f'{os.path.splitext(ref_fname)[0]}' \
                                + f'_{txt_keys["KWORD_DISP"]}.chk'
                        else:
                            txt_chk = f'generic_{txt_keys["KWORD_DISP"]}.chk'
                        txt = tmpl_text.replace(KWORD_GEOM, txt_geom)\
                            .replace(KWORD_DESC, txt_desc)\
                            .replace(KWORD_CHK, txt_chk)\
                            .replace(KWORD_DISP, txt_keys[KWORD_DISP]) \
                            .replace(KWORD_STEP, txt_keys[KWORD_STEP]) \
                            .replace(KWORD_NUM, txt_keys[KWORD_NUM]) \
                            .replace(KWORD_XYZ, txt_keys[KWORD_XYZ]) \
                            .replace(KWORD_DIR, txt_keys[KWORD_DIR])
                        print(txt, file=fobj)
                        if tmpl_files is not None:
                            fobj.close()
    else:
        if displ_vecs is None:
            raise InternalError(
                'Missing displacement vectors to build geometries.')
        nvib = displ_vecs.shape[0]
        for i in displ_idx:
            displ = np.reshape(displ_vecs[i, :]*displ_step, (natoms, 3))
            for step, step_data in INFO_PM_STEPS.items():
                for istep in range(1, displ_num+1):
                    crd = (at_crd_ref + step_data[1]*displ)*au2unit
                    txt_keys = gen_fname_from_tmpl(nvib, displ_num, i+1, istep,
                                                   label_dir=step)
                    if tmpl_files is not None:
                        fname_inp = tmpl_files\
                            .replace(KWORD_DISP, txt_keys[KWORD_DISP]) \
                            .replace(KWORD_STEP, txt_keys[KWORD_STEP]) \
                            .replace(KWORD_NUM, txt_keys[KWORD_NUM]) \
                            .replace(KWORD_DIR, txt_keys[KWORD_DIR])
                        fobj = open(fname_inp, 'w', encoding='utf-8')
                    else:
                        fobj = sys.stdout
                        fobj.write(72*'-'+'\n')
                        fname_inp = ''
                    txt_desc = f'{step_data[0]} step along mode {i+1}'
                    txt_geom = '\n'.join(
                        [fmt_geom.format(at=at_lab[ia], xyz=crd[ia, :])
                            for ia in range(natoms)]
                    )
                    if fname_inp:
                        txt_chk = os.path.splitext(fname_inp)[0] + '.chk'
                    elif ref_fname is not None:
                        txt_chk = f'{os.path.splitext(ref_fname)[0]}' \
                            + f'_{txt_keys[KWORD_DISP]}.chk'
                    else:
                        txt_chk = f'generic_{txt_keys[KWORD_DISP]}.chk'
                    txt = tmpl_text.replace(KWORD_GEOM, txt_geom)\
                        .replace(KWORD_DESC, txt_desc)\
                        .replace(KWORD_CHK, txt_chk)\
                        .replace(KWORD_DISP, txt_keys[KWORD_DISP]) \
                        .replace(KWORD_STEP, txt_keys[KWORD_STEP]) \
                        .replace(KWORD_NUM, txt_keys[KWORD_NUM]) \
                        .replace(KWORD_DIR, txt_keys[KWORD_DIR])
                    print(txt, file=fobj)
                    if tmpl_files is not None:
                        fobj.close()


#
# MAIN
# ====
def main_shift(fname_ref: str,
               fname_tmpl: str | None,
               der_coord: str,
               n_steps: int,
               der_step: float,
               au2unit: float,
               opts: argparse.Namespace):
    """Run the geometry builder.

    The routine handles operation to build displace geometries and
    create files.

    Parameters
    ----------
    fname_ref
        Filename for the reference geometry.
    fname_tmpl
        Template file for the filename.
    der_coord
        Derivation coordinate.
    n_steps
        Number of steps in each direction.
    der_step
        Size of the displacement step for the differentation.
    au2unit
        Conversion factor from bohrs to output unit.
    opts
        User options.
    """
    fref = DataFile(fname_ref)
    # Extract atom numbers
    qdata = fref.get_data(error_noqty=True, **QLABS_ATGEOM)
    if qdata['atnum'] is None:
        print('ERROR: Could not parse atomic numbers')
        sys.exit(1)
    else:
        atnum = qdata['atnum'].data
    atsmb = convert_labsymb(True, *atnum)
    natoms = len(atsmb)
    if qdata['atcrd'] is None:
        print('ERROR: Could not parse atomic coordinates')
        sys.exit(1)
    else:
        atcrd = np.array(qdata['atcrd'].data).reshape((natoms, 3))
    if qdata['atmas'] is None:
        print('ERROR: Could not parse atomic masses')
        sys.exit(1)
    else:
        atmas = np.array(qdata['atmas'].data)

    # Set displacement coordinates
    match der_coord:
        case 'Q' | 'q':
            displ_cart = False
            try:
                hessdat = fref.get_hess_data(get_evec=False,
                                             get_lweigh='lwmat',
                                             norm_lweigh=True)
            except QuantityError as err:
                print('Failed to extract normal-modes coordinates')
                print(err)
                sys.exit(2)
            lwmat = orient_modes(hessdat['lwmat'])
            nvib = hessdat['lwmat'].shape[0]
            if opts.indexes is None:
                indexes = range(nvib)
            else:
                indexes = convert_range_spec(opts.indexes, py_index=True)
            if max(indexes) >= nvib:
                print('Incorrect index specification, outside valid range.')
                print(f'The system has {nvib} modes.')
                sys.exit(1)
        case 'X':
            displ_cart = True
            if opts.indexes is None:
                indexes = range(natoms)
            else:
                indexes = convert_range_spec(opts.indexes, py_index=True)
            if max(indexes) >= natoms:
                print('Incorrect index specification, outside valid range.')
                print(f'The system has {natoms} atoms.')
                sys.exit(1)
            lwmat = None
        case _:
            print(f'ERROR: Unsupported steps along {der_coord} for now')
            sys.exit(2)

    tmpl_text = parse_template_file(fname_tmpl)
    if opts.pattern is not None:
        tmpl_file = check_fname_pattern(opts.pattern, displ_cart)
        if n_steps > 1 and not (KWORD_STEP in tmpl_file
                                or KWORD_DISP in tmpl_file):
            raise ValueError('Missing step specification in file pattern, '
                             + 'mandatory for multi-steps displacements.')
    else:
        tmpl_file = None

    build_geoms(atcrd, atsmb, atmas, indexes, n_steps, der_step, displ_cart,
                au2unit, tmpl_text, tmpl_file, fname_ref, lwmat)
