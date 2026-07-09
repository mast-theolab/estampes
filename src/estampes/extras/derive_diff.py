"""Program DERIVEUR - Numerical differentiation module.

The module contains the routines specific to construct numerical
derivatives.
"""

import os
import re
import sys
import argparse
import typing as tp
from math import sqrt
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from estampes.base import (
    AtsLabType, AtsMasType, QLabel,
    InternalError, ParsingError, QuantityError)
from estampes.parser import DataFile
from estampes.data.physics import phys_fact
from estampes.data.property import QBaseInfo, property_data as pinfo
from estampes.tools.atom import convert_labsymb
from estampes.tools.char import convert_range_spec
from estampes.tools.math import square_ltmat, superpose
from estampes.tools.output import fortran_fmt_D, pstruct_to_labels, sec_header
from estampes.tools.vib import orient_modes

from estampes.extras.derive_core import (
    USE_FORTRAN_FMT, LABEL_MSTEP, LABEL_PSTEP,
    QLABS_ATGEOM, THRESH_NUMERR, THRESH_SMALL,
    get_tmpl_fmt_from_file)


#
# USER OPTIONS
# ============
def parse_qty(qlabel: str) -> dict[int, list[int]]:
    """Parse the quantity label given by the user.

    Parses a string containing the quantities and derivative orders of
    interest, and returns a dictionary containing the data to process.
    An error is raised if the label is incorrectly built.

    Parameters
    ----------
    qlabel
        string containing the quantity labels (see HELP_QTY for details)

    Returns
    -------
    dict
        dictionary containing:
        - quantity labels as keys (compatible with parse_qlabel)
        - derivative orders as a tuple (0: reference structure)

    Raises
    ------
    KeyError
        Wrong structure for the label
    """
    key = re.compile(r'(?P<dinfo>d(?P<dord>\d+-\d+|\d+))?(?P<qty>.+)',
                     flags=re.I)
    qdata = {}
    for item in qlabel.split(','):
        res = key.match(item)
        if res is None:
            fmt = 'ERROR: Unrecognized quantity specification "{}"'
            raise KeyError(fmt.format(item))
        qtag = res.groupdict()
        qkey = qtag['qty'].strip().lower()
        if qkey in ('energy', 'v'):
            qty = 1
        elif qkey in ('elecdip', 'ed', 'u'):
            qty = 101
        elif qkey in ('magdip', 'md', 'm'):
            qty = 102
        elif qkey in ('poltens', 'pt', 'a'):
            qty = 103
        elif qkey in ('fdalpha', 'fda', 'a(w)'):
            qty = 301
        elif qkey in ('fdoptrot', 'fdor', 'o(w)'):
            qty = 302
        elif qkey in ('fddipquad', 'fddq', 'q(w)'):
            qty = 304
        else:
            raise KeyError('Unrecognized quantity label')
        if qty not in qdata:
            qdata[qty] = []
        dkey = qtag['dord']
        if '-' in dkey:
            dmin, dmax = [int(item) for item in dkey.split('-')]
            qdata[qty] = sorted(set(qdata[qty] + list(range(dmin, dmax+1))))
        else:
            dord = int(dkey)
            if dkey not in qdata[qty]:
                qdata[qty] = sorted([dord] + qdata[qty])
    return qdata


#
# PRINTING
# ========
def print_header(output: tp.IO,
                 atcrd_ref: npt.NDArray,
                 atlab: AtsLabType,
                 evec_ref: npt.NDArray | None):
    """Print reference header.

    Prints reference header data in the output.

    Parameters
    ----------
    atcrd_ref
        Reference coordinates.
    atlab
        Atomic label.
    evec_ref
        Reference eigenvectors matrix.
    output
        File object where data must be printed
    """
    print(' REF_GEOM', file=output)
    fmt = '  {at:<3s}{xyz[0]:19.13f}{xyz[1]:19.13f}{xyz[2]:19.13f}'
    atcrd = np.asarray(atcrd_ref)
    for i, at in enumerate(atlab):
        print(fmt.format(at=at, xyz=atcrd[i, :]), file=output)
    if evec_ref is not None:
        print(' REF_NM', file=output)
        maxcols = 5
        tmpvec = evec_ref.reshape((-1, ))
        nitems = len(tmpvec)
        for i in range(0, nitems, maxcols):
            ncols = min(maxcols, nitems-i)
            fmt = ncols*'{:17.8E}'
            print(fmt.format(*[tmpvec[i+j] for j in range(ncols)]),
                  file=output)


def print_header_indata(output: tp.IO,
                        atcrd_ref: npt.NDArray,
                        atlab: AtsLabType,
                        atmas: AtsMasType,
                        evec_ref: npt.NDArray,
                        freq_ref: npt.NDArray,
                        nm_list: Sequence[int]):
    """Print reference header.

    Prints reference header data in the output.

    Parameters
    ----------
    atcrd_ref
        Reference coordinates.
    atlab
        Atomic label.
    atmas
        Atomic masses.
    evec_ref
        Reference eigenvectors matrix.
    freq_ref
        Reference frequencies (cm-1).
    nm_list
        List of active normal modes for the differentiation.
    output
        File object where data must be printed
    """
    MAX_COLS_R = 5
    MAX_COLS_I = 10
    FMT_D = 'D17.8'
    FMT_E = '{:17.8E}'

    print(' REF_GEOM', file=output)
    fmt = '  {at:<3s}{xyz[0]:19.13f}{xyz[1]:19.13f}{xyz[2]:19.13f}'
    atcrd = np.asarray(atcrd_ref)
    for i, at in enumerate(atlab):
        print(fmt.format(at=at, xyz=atcrd[i, :]), file=output)
    print(' REF_ATMASS', file=output)
    natoms = len(atmas)
    for i in range(0, natoms, MAX_COLS_R):
        nmax = min(i+MAX_COLS_R, natoms)
        if USE_FORTRAN_FMT:
            print(fortran_fmt_D(list(atmas[i:nmax]), FMT_D), file=output)
        else:
            fmt = (nmax-i)*FMT_E
            print(fmt.format(*atmas[i:nmax]), file=output)
    print(' REF_NM', file=output)
    tmpvec = evec_ref.reshape((-1, ))
    nitems = len(tmpvec)
    for i in range(0, nitems, MAX_COLS_R):
        nmax = min(i+MAX_COLS_R, nitems)
        if USE_FORTRAN_FMT:
            print(fortran_fmt_D(list(tmpvec[i:nmax]), FMT_D), file=output)
        else:
            fmt = (nmax-i)*FMT_E
            print(fmt.format(*tmpvec[i:nmax]), file=output)
    print(' REF_FREQ', file=output)
    nvib = len(freq_ref)
    for i in range(0, nvib, MAX_COLS_R):
        nmax = min(i+MAX_COLS_R, nvib)
        fmt = (nmax-i)*'{:12.5f}'
        print(fmt.format(*freq_ref[i:nmax]), file=output)
    print(' REF_NMLIST', file=output)
    nm_indexes = [1 if i in nm_list else 0 for i in range(nvib)]
    for i in range(0, nvib, MAX_COLS_I):
        nmax = min(i+MAX_COLS_I, nvib)
        fmt = (nmax-i)*'{:7d}'
        print(fmt.format(*nm_indexes[i:nmax]), file=output)


def print_indatanm(qtag: int,
                   qinfo: QBaseInfo,
                   ddata: dict[int, tp.Any],
                   output: tp.IO):
    """Print data in InDataNM format.

    Prints quantity related data in the InDataNM format.

    Parameters
    ----------
    qtag
        Quantity label.
    qinfo
        Quantity information.
    ddata
        Data for each derivative order.
    output
        Output file object.
    """
    if qtag == 1:
        label = 'ENERGY'
    elif qtag == 101:
        label = 'ELDIP'
    elif qtag == 102:
        label = 'MAGDIP'
    elif qtag == 103:
        label = 'POLAR'
    else:
        label = qinfo.name.replace(' ', '_').upper()

    if qinfo.size < 0:
        raise NotImplementedError('Atom-dependent properties NYI for output.')

    FMT_D = 'D17.8'
    FMT_E = qinfo.size*'{:17.8E}'
    if qinfo.size == 1:
        if USE_FORTRAN_FMT:
            def dat_to_str(dat):
                return fortran_fmt_D(dat, FMT_D)
        else:
            @tp.override
            def dat_to_str(dat):
                return FMT_E.format(dat)

        def chk_val(x: float) -> bool:  # pyright: ignore[reportRedeclaration]
            return abs(x) > THRESH_SMALL
    else:
        if USE_FORTRAN_FMT:
            def dat_to_str(dat):
                return fortran_fmt_D(dat, FMT_D)
        else:
            @tp.override
            def dat_to_str(dat):
                return FMT_E.format(*dat)

        def chk_val(x: npt.NDArray) -> bool:
            return max(np.abs(x)) > THRESH_SMALL

    for derord in sorted(ddata.keys()):
        if derord == 0:
            output.write(f' {label}_0\n')
            output.write(f'{dat_to_str(ddata[derord])}\n')
        elif derord == 1:
            output.write(f' {label}_D1N\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                if chk_val(row):
                    output.write(f'{modes+1:7d} {dat_to_str(row)}\n')
        elif derord == 2:
            output.write(f' {label}_D2N\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                if chk_val(row):
                    ii, ij = modes
                    output.write(f'{ii+1:7d} {ij+1:7d} {dat_to_str(row)}\n')
        elif derord == 3:
            output.write(f' {label}_D3N\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                if chk_val(row):
                    ii, ij, ik = modes
                    output.write(f'{ii+1:7d} {ij+1:7d} {ik+1:7d} '
                                 + f'{dat_to_str(row)}\n')
        elif derord == 4:
            output.write(f' {label}_D4N\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                if chk_val(row):
                    ii, ij, ik, il = modes
                    output.write(f'{ii+1:7d} {ij+1:7d} {ik+1:7d} {il+1:7d}'
                                 + f'{dat_to_str(row)}\n')


def print_indatax(qtag: int,
                  qinfo: QBaseInfo,
                  natoms: int,
                  ddata: dict[int, tp.Any],
                  dstat: dict[int, int],
                  output: tp.IO):
    """Print data in InDataX format.

    Prints quantity related data in the InDataX format.

    Parameters
    ----------
    qtag
        Quantity label.
    qinfo
        Quantity information.
    natoms
        Number of atoms.
    ddata
        Data for each derivative order.
    dstat
        Derivative calculation status.

        1: available analytically
        2: done numerically (at least partially).

    output
        Output file object.
    """
    if qtag == 1:
        label = 'ENERGY'
    elif qtag == 101:
        label = 'ELDIP'
    elif qtag == 102:
        label = 'MAGDIP'
    elif qtag == 103:
        label = 'POLAR'
    else:
        label = qinfo.name.replace(' ', '_').upper()

    if qinfo.size < 0:
        raise NotImplementedError('Atom-dependent properties NYI for output.')

    FMT_D = 'D17.8'
    FMT_E = qinfo.size*'{:17.8E}'
    xyz = ('x', 'y', 'z')

    if USE_FORTRAN_FMT:
        def dat_to_str(dat):
            return fortran_fmt_D(dat, FMT_D)
    else:
        @tp.override
        def dat_to_str(dat):
            if isinstance(dat, float):
                return FMT_E.format(dat)
            else:
                return FMT_E.format(*dat)

    for derord in sorted(ddata.keys()):
        if derord == 0:
            output.write(f' {label}_0\n')
            output.write(f'{dat_to_str(ddata[derord])}\n')
        elif derord == 1:
            if dstat == 2:
                output.write(f' {label}_D1N\n')
                for i, row in ddata[derord].items():
                    output.write(f'{i+1:6d} {dat_to_str(row)}\n')
            else:
                output.write(f' {label}_D1X\n')
                for ia in range(natoms):
                    for ix, labxyz in enumerate(xyz):
                        output.write(f'{ia+1:6d} {labxyz} '
                                     + f'{dat_to_str(ddata[derord][ix])}\n')
        elif derord == 2:
            if dstat == 2:
                do_prt = True
                for modes in sorted(ddata[derord]):
                    row = ddata[derord][modes]
                    if isinstance(modes, tuple):
                        if do_prt:
                            output.write(f' {label}_D2N\n')
                            do_prt = False
                        i, j = modes
                        output.write(f'{i+1:6d} {j+1:6d} {dat_to_str(row)}\n')
                    else:
                        if do_prt:
                            output.write(f' {label}_D2X\n')
                            do_prt = False
                        ind = 0
                        for ia in range(natoms):
                            for ix in xyz:
                                output.write(f'{modes+1:6d} {ia+1:6d} {ix} '
                                             + f'{dat_to_str(row[ind])}\n')
                                ind += 1
            else:
                ind = 0
                for ia in range(natoms):
                    for ix in xyz:
                        for ja in range(ia+1):
                            for jx in xyz:
                                output.write(
                                    f'{ia+1:6d} {ix} {ja+1:6d} {jx} '
                                    + f'{dat_to_str(ddata[derord][ind])}\n')
                                ind += 1
        elif derord == 3:
            output.write(f' {label}_D3X\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                if isinstance(modes, tuple):
                    i, j = modes
                    ind = 0
                    for ia in range(natoms):
                        for ix in xyz:
                            output.write(f'{i+1:6d} {j+1:6d} {ia+1:6d} {ix} '
                                         + f'{dat_to_str(row[ind])}\n')
                            ind += 1
                else:
                    ind = 0
                    for ia in range(natoms):
                        for ix in xyz:
                            for ja in range(ia+1):
                                for jx in xyz:
                                    output.write(
                                        f'{modes+1:6d} {ia+1:6d} {ix} '
                                        + f'{ja+1:6d} '
                                        + f'{jx} {dat_to_str(row[ind])}\n')
                                    ind += 1
        elif derord == 4:
            output.write(f' {label}_D4X\n')
            for modes in sorted(ddata[derord]):
                row = ddata[derord][modes]
                i, j = modes
                ind = 0
                for ia in range(natoms):
                    for ix in xyz:
                        for ja in range(ia+1):
                            for jx in xyz:
                                output.write(
                                    f'{i+1:6d} {j+1:6d} {ia+1:6d} {ix} '
                                    + f'{ja+1:6d} {jx}'
                                    + f' {dat_to_str(row[ind])}\n')
                                ind += 1


#
# DIFFERENTIATION ROUTINES
# ========================
def build_2dq_num_diff(qlabel: QLabel,
                       derord: int,
                       qinfo: QBaseInfo,
                       file_pattern: str,
                       nvib: int,
                       indexes: Sequence[int],
                       der_step: float,
                       do_d1: bool,
                       do_d2: bool,
                       fref: DataFile,
                       lmat: npt.NDArray | None,
                       freq: npt.NDArray | None,
                       error_on_missing: bool = True,
                       print_warn: bool = True):
    """Build derivatives from 2-steps numerical differentiation.

    Builds derivatives from 2-steps numerical differentiation.

    Parameter
    ---------
    qlabel
        Quantity label.
    derord
        Analytical derivative order, for Cartesian/normal conversion.
    qinfo
        Quantity information.
    file_pattern
        Pattern for the filename for the displaced geometries, as format
        specification.
    nvib
        Total number of normal modes.
    der_step
        Size of the displacement step for the differentation.
    indexes
        Differentiation indexes to take into account.
    do_d1
        Do first numerical differentiation.
    do_d2
        Do second numerical differentiation.
    fref
        Datafile for the reference geometry, used for `do_d2 = True`.
    lmat
        Conversion matrix Cartesian -> normal coordinates.
        If not present, no conversion of analytic derivatives is done.
    freq
        Harmonic wavenumbers, used for analysis of derivatives.
    error_on_missing
        Raise an error if a differentiation point is missing.
    print_warn
        Print warning messages if missing files.
    """
    data = {}
    if do_d1:
        data['d1'] = {}
    if do_d2:
        data['d2'] = {}

    qdata0 = None
    if do_d2:
        qdata0 = fref.get_data(qlabel)[str(qlabel)]
        if qdata0 is None:
            raise ParsingError(
                f'Failed to extract {qlabel} from central point')
        qdata0 = np.array(qdata0.data)

    dq = 2.0 * der_step
    dq_sq = der_step**2

    for imode in indexes:
        diffdata: dict[str, npt.NDArray] = {}
        fplus = file_pattern.format(coord=imode+1, dir=LABEL_PSTEP)
        fminus = file_pattern.format(coord=imode+1, dir=LABEL_MSTEP)
        if not os.path.exists(fplus) or not os.path.exists(fminus):
            msg = f'Missing differentiation files along mode {imode+1}'
            if print_warn:
                print(msg)
            if error_on_missing:
                print(msg + f'for {qinfo.name}')
                raise QuantityError(msg)
            continue
        dfile = DataFile(fplus)
        dplus = dfile.get_data(qlabel)[str(qlabel)]
        if dplus is None:
            msg = f'Failed to extract {qlabel} from +dQ_{imode}.'
            if print_warn:
                print(msg)
            if error_on_missing:
                raise ParsingError(msg)
            continue
        dfile = DataFile(fminus)
        dminus = dfile.get_data(qlabel)[str(qlabel)]
        if dminus is None:
            msg = f'Failed to extract {qlabel} from -dQ_{imode}.'
            if print_warn:
                print(msg)
            if error_on_missing:
                raise ParsingError(msg)
            continue
        if do_d1:
            diffdata['d1'] = \
                (np.array(dplus.data) - np.array(dminus.data))/dq
        if do_d2:
            if qdata0 is None:
                raise InternalError('reference data should have been loaded.')
            diffdata['d2'] = \
                (np.array(dplus.data) + np.array(dminus.data)
                 - 2*qdata0)/dq_sq

        if lmat is None or derord == 0:
            if do_d1:
                data['d1'][imode] = diffdata['d1'].copy()
            if do_d2:
                data['d2'][imode] = diffdata['d2'].copy()
        else:
            # We convert to normal modes, we need to check order of analytic
            # derivatives
            if derord == 2:
                if qlabel.label == 1:
                    if do_d1:
                        diffdata['d1'] = np.einsum(
                            'ij,jk,lk', lmat, square_ltmat(diffdata['d1']),
                            lmat)
                    if do_d2:
                        diffdata['d2'] = np.einsum(
                            'ij,jk,lk', lmat, square_ltmat(diffdata['d2']),
                            lmat)
                else:
                    raise NotImplementedError(
                        'Unknown analytic second derivs.')
            elif derord == 1:
                if qinfo.size == 1:
                    shape = (-1, )
                    formula = 'ij,j->i'
                else:
                    if qinfo.size < 0:
                        raise NotImplementedError(
                            'Atom-dependent properties NYI')
                    shape = (-1, qinfo.dim)
                    formula = 'ij,jk->ik'
                # Type ignore below because Pyright is just too dumb to
                # understand the structure of diffdata even when explained...
                if do_d1:
                    diffdata['d1'] = np.einsum(
                        formula, lmat,
                        diffdata['d1'].reshape(shape))  # type: ignore
                if do_d2:
                    diffdata['d2'] = np.einsum(
                        formula, lmat,
                        diffdata['d2'].reshape(shape))  # type: ignore
            elif derord != 0:
                raise NotImplementedError('Unsupported analytic der. order')

            if qinfo.d1q == 'q':
                # All indexes can permute
                if derord == 2:
                    for i in range(nvib):
                        for j in range(i+1):
                            ijk = tuple(sorted((imode, i, j),
                                        reverse=True))
                            iijk = tuple(sorted((imode, imode, i, j),
                                                reverse=True))
                            if do_d1:
                                if ijk not in data['d1']:
                                    data['d1'][ijk] = []
                                data['d1'][ijk].append(diffdata['d1'][i][j])

                            if do_d2:
                                if iijk not in data['d2']:
                                    data['d2'][iijk] = []
                                data['d2'][iijk].append(diffdata['d2'][i][j])

                elif derord == 1:
                    for i in range(nvib):
                        ij = tuple(sorted((imode, i), reverse=True))
                        iij = tuple(sorted((imode, imode, i),
                                           reverse=True))
                        if do_d1:
                            if ij not in data['d1']:
                                data['d1'][ij] = []
                            data['d1'][ij].append(diffdata['d1'][i])

                        if do_d2:
                            data['d2'][iij] = diffdata['d2'][i]
            elif qinfo.d1q == 'p':
                # Indexes of differentiation cannot permute with others
                if derord == 2:
                    for i in range(nvib):
                        for j in range(i+1):
                            ijk = tuple([imode] + sorted((i, j),
                                                         reverse=True))
                            iijk = tuple([imode, imode]
                                         + sorted((i, j), reverse=True))
                            if do_d1:
                                data['d1'][ijk] = diffdata['d1'][i][j]

                            if do_d2:
                                if iijk not in data['d2']:
                                    data['d2'][iijk] = []
                                data['d2'][iijk] = diffdata['d2'][i][j]
                elif derord == 1:
                    for i in range(nvib):
                        ij = (imode, i)
                        iij = (imode, imode, i)
                        if do_d1:
                            data['d1'][ij] = diffdata['d1'][i]

                        if do_d2:
                            data['d2'][iij] = diffdata['d2'][i]

            else:
                raise ValueError(
                    'Unknown type of coordinates for analytic derivatives.')

    if lmat is not None and qinfo.d1q == 'q':
        fmt_mode = f'{{:{len(str(nvib))}d}}'
        fmt_err_ij = '{dqty} | δQ({qi}): {fi:15.8e} | ' \
            + 'δQ({qj}): {fj:15.8e} | Error: {err:.2%}'
        fmt_err_ij_x = '{dqty} | {comp} | δQ({qi}): {fi:15.8e} | ' \
            + 'δQ({qj}): {fj:15.8e} | Error: {err:.2%}'
        fmt_err_ijk = '{dqty} | δQ({qi}): {fi:15.8e} | ' \
            + 'δQ({qj}): {fj:15.8e} | δQ({qk}): {fk:15.8e} | ' \
            + 'Error: {err:.2%}'

        sec_header(sys.stdout, 3, 'Inconsistency Check')
        if freq is None:
            raise ValueError('Missing harmonic frequencies')
        if derord == 2:
            if do_d2:
                dord = 4
                fmt_qty = 'D4(' + (dord-1)*(fmt_mode+',') + fmt_mode + ')'
                for ijkl, row in data['d2'].items():
                    i, _, j, _ = ijkl
                    wi = freq[i]
                    lab_i = fmt_mode.format(i+1)
                    wj = freq[j]
                    lab_j = fmt_mode.format(j+1)
                    lab_qty = fmt_qty.format(i+1, i+1, j+1, j+1)
                    nvals = len(row)
                    fact = phys_fact('mwq2q')**4/sqrt(wi**2*wj**2)
                    if nvals == 2:
                        vals = sorted(row, key=abs)
                        err = abs((vals[0]-vals[-1])/vals[-1])
                        if (err > THRESH_NUMERR
                                and abs(vals[-1]) > THRESH_SMALL):
                            print(fmt_err_ij.format(
                                dqty=lab_qty, qi=lab_j, qj=lab_i,
                                fi=fact*row[0], fj=fact*row[1],
                                err=err))
                        data['d2'][ijkl] = sum(row)/2.
                    elif nvals == 1:
                        data['d2'][ijkl] = row[0]
            if do_d1:
                dord = 3
                fmt_qty = 'D3(' + (dord-1)*(fmt_mode+',') + fmt_mode + ')'
                for ijk, row in data['d1'].items():
                    i, j, k = ijk
                    wi = freq[i]
                    lab_i = fmt_mode.format(i+1)
                    wj = freq[j]
                    lab_j = fmt_mode.format(j+1)
                    wk = freq[k]
                    lab_k = fmt_mode.format(k+1)
                    lab_qty = fmt_qty.format(i+1, j+1, k+1)
                    fact = phys_fact('mwq2q')**3/sqrt(wi*wj*wk)
                    nvals = len(row)
                    if nvals == 3:
                        vals = sorted(row, key=abs)
                        err = max(abs((vals[0]-vals[-1])/vals[-1]),
                                  abs((vals[1]-vals[-1])/vals[-1]))
                        if (err > THRESH_NUMERR and
                                abs(vals[-1]) > THRESH_SMALL):
                            print(fmt_err_ijk.format(
                                dqty=lab_qty, qi=lab_k, qj=lab_j, qk=lab_i,
                                fi=fact*row[0], fj=fact*row[1],
                                fk=fact*row[2], err=err))
                        data['d1'][ijk] = sum(row)/3.
                    elif nvals == 2:
                        vals = sorted(row, key=abs)
                        err = abs((vals[0]-vals[-1])/vals[-1])
                        if (err > THRESH_NUMERR and
                                abs(vals[-1]) > THRESH_SMALL):
                            if j == i:
                                lab2 = lab_i
                                lab1 = lab_k
                            elif k == j:
                                lab2 = lab_i
                                lab1 = lab_j
                            elif i not in indexes:
                                lab2 = lab_j
                                lab1 = lab_k
                            elif j not in indexes:
                                lab2 = lab_i
                                lab1 = lab_k
                            elif k not in indexes:
                                lab2 = lab_i
                                lab1 = lab_j
                            else:
                                msg = '''\
Only 2 differentiations available for cubic force constants but all 3 active.
I am confused.'''
                                raise ValueError(msg)
                            print(fmt_err_ij.format(
                                dqty=lab_qty, qi=lab1, qj=lab2,
                                fi=fact*row[0], fj=fact*row[1], err=err))
                        data['d1'][ijk] = sum(row)/2.
                    elif nvals == 1:
                        data['d1'][ijk] = row[0]
        elif derord == 1:
            if do_d1:
                dord = 2
                lab_comp = pstruct_to_labels(qinfo, linearize=True)
                fmt_qty = 'D2(' + fmt_mode + ',' + fmt_mode + ')'
                for ij, row in data['d1'].items():
                    i, j = ij
                    wi = freq[i]
                    lab_i = fmt_mode.format(i+1)
                    wj = freq[j]
                    lab_j = fmt_mode.format(j+1)
                    lab_qty = fmt_qty.format(i+1, j+1)
                    nvals = len(row)
                    fact = phys_fact('mwq2q')**2/sqrt(wi*wj)
                    if nvals == 2:
                        comps = []
                        if qinfo.size == 1:  # quantity is a scalar
                            vals = sorted(row, key=abs)
                            err = abs((vals[0]-vals[-1])/vals[-1])
                            if (err > THRESH_NUMERR
                                    and abs(vals[-1]) > THRESH_NUMERR):
                                print(fmt_err_ij.format(
                                    dqty=lab_qty, qi=lab_j, qj=lab_i,
                                    fi=fact*row[0], fj=fact*row[1],
                                    err=err))
                            data['d1'][ij] = sum(vals)/2.
                        else:
                            comps = []
                            for ixyz in range(len(row[0])):
                                vals = sorted([item[ixyz] for item in row],
                                              key=abs)
                                err = abs((vals[0]-vals[-1])/vals[-1])
                                if (err > THRESH_NUMERR
                                        and abs(vals[-1]) > THRESH_NUMERR):
                                    print(fmt_err_ij_x.format(
                                        dqty=lab_qty,
                                        comp=lab_comp[ixyz].upper(),
                                        qi=lab_j, qj=lab_i,
                                        fi=fact*row[0][ixyz],
                                        fj=fact*row[1][ixyz],
                                        err=err))
                                comps.append(sum(vals)/2.)
                            data['d1'][ij] = np.array(comps)
                    elif nvals == 1:
                        data['d1'][ij] = row[0]

    return data


def build_2dx_num_diff(qlabel: QLabel,
                       derord: int,
                       qinfo: QBaseInfo,
                       file_pattern: str,
                       natoms: int,
                       indexes: Sequence[int],
                       der_step: float,
                       do_d1: bool,
                       do_d2: bool,
                       fref: DataFile,
                       error_on_missing: bool = True,
                       print_warn: bool = True):
    """Build derivatives from 2-steps numerical differentiation.

    Builds derivatives from 2-steps numerical differentiation.

    Parameter
    ---------
    qlabel
        Quantity label.
    derord
        Analytical derivative order, for Cartesian/normal conversion.
    qinfo
        Quantity information.
    file_pattern
        Pattern for the filename for the displaced geometries, as format
        specification.
    natoms
        Total number of atoms.
    der_step
        Size of the displacement step for the differentation.
    indexes
        Differentiation indexes to take into account.
    do_d1
        Do first numerical differentiation.
    do_d2
        Do second numerical differentiation.
    fref
        Datafile for the reference geometry, used for `do_d2 = True`.
    error_on_missing
        Raise an error if a differentiation point is missing.
    print_warn
        Print warning messages if missing files.
    """
    raise NotImplementedError('Cartesian derivatives NYI')


def get_dq_diff(qtag: int,
                derords: Sequence[int],
                qinfo: QBaseInfo,
                fref: DataFile,
                file_pattern: str,
                der_step: float,
                indexes: Sequence[int],
                lwmat: npt.NDArray | None,
                freq: npt.NDArray | None,
                output_fmt: str,
                output: tp.IO):
    """Get derivatives with respect to normal coordinates.

    Get derivatives with respect to normal coordinates of quantity
    `qtag`.

    Parameters
    ----------
    qtag
        Quantity label.
    derords
        Derivative orders of interest.
    qinfo
        Quantity information.
    fref
        Reference file, as DataFile.
    file_pattern
        Pattern for the filename for the displaced geometries, as format
        specification.
    der_step
        Size of the displacement step for the differentation.
    indexes
        Differentiation indexes to take into account.
    lwmat
        Conversion matrix Cartesian -> normal coordinates.
    freq
        Harmonic frequencies (for consistency check).
    output_fmt
        Output format.
    output
        Output file object.
    """
    if lwmat is None:
        raise ValueError('Missing eigenvector matrix')
    nvib = lwmat.shape[0]

    # Now check what is available analytically and not.
    # Available analytic data should be in fref, while derivatives are
    # constructed from displaced geometries.
    # We first find a step file to use to check data
    for i in indexes:
        fshift_ref = file_pattern.format(coord=i+1, dir=LABEL_PSTEP)
        if os.path.exists(fshift_ref):
            break
    else:
        msg = 'Could not find any displaced geoemtry file. Check pattern.'
        raise FileNotFoundError(msg)
    dshift_ref = DataFile(fshift_ref)

    # Now let us check what is available or not:
    # Build list of derivatives to do.
    # Status values:
    # 0: not done
    # 1: available analytically
    # 2: done numerically
    der_stat = {i: 0 for i in sorted(derords, reverse=True)}
    conv_to_nm = output_fmt.endswith('nm')
    ddata = {}
    for derord, derstat in der_stat.items():
        if not derstat:
            qkey = QLabel(quantity=qtag, derorder=derord, dercoord='X')
            qkey_lab = str(qkey)
            try:
                qdat = fref.get_data(qkey, error_noqty=False)
            except NotImplementedError:
                msg = f'''\
Data parser reported it did not support {qkey}.
Assuming that the derivative is not yet available.'''
                print(msg)
                qdat = {qkey_lab: None}
            if qdat[qkey_lab] is None:
                for i in range(derord-1, -1, -1):
                    qlabel = QLabel(quantity=qtag, derorder=i, dercoord='X')
                    try:
                        qdata1 = dshift_ref.get_data(qlabel, error_noqty=False)
                    except NotImplementedError:
                        msg = f'''\
Data parser reported it did not support {qlabel}.
Assuming that the derivative is not yet available.'''
                        print(msg)
                        qdata1 = {str(qlabel): None}
                    if qdata1[str(qlabel)] is not None:
                        ader = i
                        break
                else:
                    raise QuantityError(
                        f'No data relative to {qinfo.name} found')
                i = derord - ader
                if i > 2:
                    msg = f'''\
Unable to build derivatives of order {derord} for {qinfo.name}.
Analytic derivatives are only available up to order {ader}.'''
                    raise QuantityError(qinfo.name, msg)
                elif i == 2:
                    do_d2 = True
                    do_d1 = False
                    # Let us check if lower order also required and not
                    # available analytically
                    if derord - 1 in der_stat:
                        do_d1 = not der_stat[derord-1]
                elif i == 1:
                    do_d1 = True
                    do_d2 = False
                else:
                    raise ValueError(
                        'Unexpected value for analytic derivative order found.'
                        )
                lmat = lwmat if conv_to_nm else None
                data = build_2dq_num_diff(qlabel, ader, qinfo, file_pattern,
                                          nvib, indexes, der_step, do_d1,
                                          do_d2, fref, lmat, freq)
                if do_d2:
                    ddata[derord] = data['d2']
                    der_stat[derord] = 2
                    if do_d1:
                        ddata[derord-1] = data['d1']
                        der_stat[derord-1] = 2
            elif not conv_to_nm:
                ddata[derord] = np.array(qdat[qkey_lab].data)  # type: ignore
                der_stat[derord] = 1
            else:
                der_stat[derord] = 1
                if derord == 2:
                    if qtag == 1:
                        # We need to symmetrize matrix for energy 2nd deriv.
                        dat = np.einsum(
                            'ij,jk,lk', lwmat,
                            square_ltmat(qdat[qkey_lab].data),  # type: ignore
                            lwmat)
                        ddata[derord] = {}
                        for i in range(nvib):
                            for j in range(i+1):
                                ij = (i, j)
                                ddata[derord][ij] = dat[i][j]
                    else:
                        raise NotImplementedError(
                            'Unknown analytic second derivs.')
                elif derord == 1:
                    if qinfo.size == 1:
                        shape = (-1, )
                        formula = 'ij,j->i'
                    else:
                        if qinfo.size < 0:
                            raise NotImplementedError(
                                'Atom-dependent properties NYI')
                        shape = (-1, qinfo.dim)
                        formula = 'ij,jk->ik'
                    # Shutting up syntax analyzers below, since they cannot get
                    # it right.
                    dat = np.einsum(
                        formula, lmat,  # type: ignore
                        np.reshape(qdat[qkey_lab].data,  # type: ignore
                                   shape))  # type: ignore
                    ddata[derord] = {}
                    for i in range(nvib):
                        ddata[derord][i] = dat[i]
                elif derord == 0:
                    if qinfo.size == 1:
                        ddata[derord] = qdat[qkey_lab].data  # type: ignore
                    else:
                        ddata[derord] = np.array(
                            qdat[qkey_lab].data)  # type: ignore
    if output_fmt.lower() == 'indatanm':
        print_indatanm(qtag, qinfo, ddata, output)
    elif output_fmt.lower() == 'indatanm':
        print_indatax(qtag, qinfo, lwmat.shape[1]//3, ddata, der_stat, output)
    else:
        raise NotImplementedError(f'Unrecognized output format: {output_fmt}')


def get_dx_diff(qtag: int,
                derords: Sequence[int],
                qinfo: QBaseInfo,
                fref: DataFile,
                file_pattern: str,
                der_step: float,
                indexes: Sequence[int],
                output_fmt: str,
                output: tp.IO):
    """Get derivatives with respect to Cartesian coordinates.

    Get derivatives with respect to Cartesian coordinates of quantity
    `qtag`.

    Parameters
    ----------
    qtag
        Quantity label.
    derords
        Derivative orders of interest.
    qinfo
        Quantity information.
    fref
        Reference file, as DataFile.
    file_pattern
        Pattern for the filename for the displaced geometries.
    der_step
        Size of the displacement step for the differentation.
    indexes
        Differentiation indexes to take into account.
    output_fmt
        Output format.
    output
        Output file object.
    """
    raise NotImplementedError('Cartesian differentiation NYI.')


#
# MAIN
# ====
def main_diff(fname_ref: str,
              file_pattern: str,
              quantity: str,
              der_coord: str,
              der_step: float,
              opts: argparse.Namespace):
    """Manage the differentiation run.

    Manages the operations to run in differentiation mode.

    Parameters
    ----------
    fname_ref
        Filename for the reference geometry.
    file_pattern
        Pattern of the files with displaced geometries.
    quantity
        Quantity/ies of interest.
    der_coord
        Derivation coordinate.
    der_step
        Size of the displacement step for the differentation.
    opts
        User options.
    """
    # Parse quantity information
    n_steps = 1
    # if n_steps != 1:
    #     print('ERROR: Multi-steps differentiations/fitting NYI.')
    #     sys.exit(1)

    # First set the output for the data
    if opts.output:
        output = open(opts.output, 'w', encoding='utf-8')
    else:
        output = sys.stdout

    sec_header(sys.stdout, -1, 'DERIVEUR - Differentiation Mode')

    sec_header(sys.stdout, 1, 'Basic Parameters')
    try:
        qdata = parse_qty(quantity)
    except KeyError as err:
        print(f'ERROR: {str(err).strip("\'")}')
        sys.exit(1)

    patt_info = get_tmpl_fmt_from_file(file_pattern)
    if patt_info['multi']:
        raise NotImplementedError('Multi-steps differentiation NYI')

    # Extract basic information from reference file
    fref = DataFile(fname_ref)
    dref = fref.get_data(error_noqty=True, **QLABS_ATGEOM)
    if dref['atnum'] is None:
        print('ERROR: Could not parse atomic numbers')
        sys.exit(1)
    else:
        atnum_ref = dref['atnum'].data
    if dref['atcrd'] is None:
        print('ERROR: Could not parse atomic coordinates')
        sys.exit(1)
    else:
        atcrd_ref = dref['atcrd'].data
    if dref['atmas'] is None:
        print('ERROR: Could not parse atomic masses')
        sys.exit(1)
    else:
        atmas_ref = dref['atmas'].data
    atsmb = convert_labsymb(True, *atnum_ref)
    natoms = len(atsmb)

    # # Check if rotation/translation information present.
    # qdata = fref.get_data(error_noqty=False, rotmat=QLabel(quantity=92),
    #                                          trvec=QLabel(quantity=93))
    # if qdata['rotmat'] is None:
    #     drotmat = np.eye(3)
    # else:
    #     drotmat = np.array(qdata['rotmat'].data)
    # print(drotmat)

    # atcrd = np.array(dref['atcrd'].data)
    # print(atcrd)
    # print(np.einsum('ij,jk->ik', atcrd,
    # drotmat)+np.array(qdata['trvec'].data))

    # sys.exit(1)

    sec_header(sys.stdout, 1, 'Differentiation Coordinates')

    # Extract displacement coordinates if necessary
    if der_coord in ('Q', 'q'):
        if opts.Lmatfile is not None:
            if not os.path.exists(opts.Lmatfile):
                print(f'ERROR: File {opts.Lmatfile} does not exist.')
                sys.exit(1)
            dfile = DataFile(opts.Lmatfile)
            qdata = dfile.get_data(error_noqty=True, **QLABS_ATGEOM)
            if qdata['atnum'] is None:
                print('ERROR: Could not extract atomic numbers from '
                      + f'{opts.Lmatfile}')
                sys.exit(1)
            else:
                atnum_tmp = qdata['atnum'].data
            if qdata['atcrd'] is None:
                print('ERROR: Could not extract atomic coordinates from '
                      + f'{opts.Lmatfile}')
                sys.exit(1)
            else:
                atcrd_tmp = qdata['atcrd'].data
            if atnum_tmp != atnum_ref:
                print('ERROR: Inconsistency in the atomic numbers between '
                      + f'{opts.Lmatfile} and {fname_ref}.')
                print('       Cannot ensure that the structures are '
                      + 'consistent (superposition).')
                sys.exit(1)
            # Deactivating typing since error_noqty ensures the absence of data
            # is impossible...
            rotmat = superpose(atcrd_ref, atcrd_tmp, atmas_ref)['rmat']
            try:
                hessdat = dfile.get_hess_data(get_evec=False,
                                              get_lweigh='lwmat',
                                              norm_lweigh=True)
            except QuantityError as err:
                print('Failed to extract normal-modes coordinates')
                print(err)
                sys.exit(2)
            nvib = hessdat['lwmat'].shape[0]
            # Ignore type because Pylance does not understand type of rotmat.
            lwmat = np.einsum('ijk,kl -> ijl',
                              hessdat['lwmat'].reshape((nvib, -1, 3)),
                              rotmat)  # type: ignore
            if opts.no_nmorient:
                lwmat = np.reshape(lwmat, (nvib, -1))
            elif opts.force_nmorient:
                facts = opts.force_nmorient.split(',')
                if len(facts) != nvib:
                    print('ERROR: Number of coefficients for nmorient '
                          + 'different from number of normal modes')
                    sys.exit(1)
                lwmat = np.reshape(lwmat, (nvib, -1))
                for i, fact in enumerate(facts):
                    if float(fact) < 0:
                        lwmat[i, :] *= -1.0
            else:
                lwmat = orient_modes(np.reshape(lwmat, (nvib, -1)))
            raise NotImplementedError('LMatfile not yet tested.')
        else:
            try:
                hessdat = fref.get_hess_data(get_evec=False,
                                             get_lweigh='lwmat',
                                             norm_lweigh=True)
            except QuantityError as err:
                print('Failed to extract normal-modes coordinates')
                print(err)
                sys.exit(2)
            nvib = hessdat['lwmat'].shape[0]
            if opts.no_nmorient:
                lwmat = hessdat['lwmat']
            elif opts.force_nmorient:
                facts = opts.force_nmorient.split(',')
                if len(facts) != nvib:
                    print('ERROR: Number of coefficients for nmorient '
                          + 'different from number of normal modes')
                    sys.exit(1)
                lwmat = hessdat['lwmat']
                for i, fact in enumerate(facts):
                    if float(fact) < 0:
                        lwmat[i, :] *= -1.0
            else:
                lwmat = orient_modes(hessdat['lwmat'])
        nmax_indexes = nvib
        if der_coord == 'q':
            raise NotImplementedError('Support of reduced q NYI.')
    else:
        hessdat = {'lwmat': None, 'eval': None}
        lwmat = None
        nvib = 0
        nmax_indexes = natoms

    if opts.indexes is None:
        indexes = range(nmax_indexes)
    else:
        indexes = convert_range_spec(opts.indexes, py_index=True)
    if max(indexes) >= nmax_indexes:
        print('Incorrect index specification, outside valid range.')
        print(f'Maximum value of the indexes is {nmax_indexes}.')
        sys.exit(1)

    sec_header(sys.stdout, 1, 'Execution')
    # == Main loop
    if opts.refheader:
        sec_header(sys.stdout, 2, 'Header data')
        if opts.format.startswith('indata'):
            if der_coord not in ('Q', 'q'):
                print('ERROR: InDataNM/X only for Cart. displacements NYI.')
                sys.exit(9)
            if lwmat is None or hessdat['eval'] is None:
                print('ERROR: Uninitialized data to print header for InData')
                sys.exit(9)
            print_header_indata(output, atcrd_ref, atsmb, atmas_ref, lwmat,
                                hessdat['eval'], indexes)
        else:
            print_header(output, atcrd_ref, atsmb, hessdat['lwmat'])

    for qlab, derords in qdata.items():
        qinfo = pinfo(qlab)
        sec_header(sys.stdout, 2, f'Building Data for {qinfo.name.title()}')
        if der_coord in ('Q', 'q'):
            if n_steps == 1:
                try:
                    get_dq_diff(qlab, derords, qinfo, fref, patt_info['fmt'],
                                der_step, indexes, lwmat, hessdat['eval'],
                                opts.format, output)
                except QuantityError as err:
                    print('ERROR: Failed to build derivatives for '
                          + f'{qinfo.name}')
                    print(err)
                    sys.exit(2)
            else:
                print('Multi-steps derivatives NYI.')
                sys.exit(9)
        elif der_coord == 'X':
            if n_steps == 1:
                try:
                    get_dx_diff(qlab, derords, qinfo, fref, patt_info['fmt'],
                                der_step, indexes, opts.format, output)
                except QuantityError as err:
                    print('ERROR: Failed to build derivatives for '
                          + f'{qinfo.name}')
                    print(err)
                    sys.exit(2)
            else:
                print('Multi-steps derivatives NYI.')
                sys.exit(9)
        else:
            print(f'ERROR: Unsupported type of coordinates, {der_coord}')
            sys.exit(9)
