"""Program MIRAGE.

MIRAGE: Molecular Image Rendition for Analysis and Graphical Experiment

A simple program to extract molecular information from a file and
generate alternative input commands for rendering software to build
publication-quality images.
"""

import sys
import os
import argparse
import typing as tp

import numpy as np

from estampes.base import QLabel
from estampes.parser import DataFile
from estampes.data.physics import PHYSFACT
from estampes.tools.mol import list_bonds
from estampes.visual.povrender import build_box, set_cam_zpos, write_pov

try:
    from PySide6 import QtGui
    Qt_avail = True
    from estampes.visual.molui import MolWin
except ModuleNotFoundError:
    Qt_avail = False

VIEW_MODELS = {
    'sticks': ('stick', 'sticks'),
    'balls': ('ball', 'balls', 'ballsandsticks'),
    'spheres': ('sphere', 'spheres', 'vdW'),
}
VIB_MODELS = {
    'arrows': ('arrow', 'arrows'),
    'spheres': ('sphere', 'spheres')
}

VIEW_ALIAS = sorted([item.lower()
                     for items in VIEW_MODELS.values()
                     for item in items])

VIB_ALIAS = sorted([item.lower()
                    for items in VIB_MODELS.values()
                    for item in items])


def build_opts(parser: argparse.ArgumentParser):
    """Build commandline options.

    Builds commandline options inside input `parser`.

    Parameters
    ----------
    parser
        Parser to update with supported options.
    """
    parser.add_argument('infiles', help="File with coordinates", nargs='+')
    parser.add_argument('--bond-tol', type=float, default=1.2,
                        help='Tolerance for the identification of bond as' +
                        '--bond-tol * covalence_radii')
    msg = 'Display (Qt) instead of generating the Pov file.'
    parser.add_argument('-D', '--display', action='store_true',
                        default=Qt_avail, help=msg)
    msg = 'Merge multiple file in 1 rendition.  ' \
        + 'Each mol. has a different color'
    parser.add_argument('--merge', action='store_true', default=False,
                        help=msg)
    msg = f'''Representation model:
- balls and stick: {', '.join(VIEW_MODELS['balls'])}
- sticks: {', '.join(VIEW_MODELS['sticks'])}
- spheres: {', '.join(VIEW_MODELS['spheres'])}
'''
    parser.add_argument('-m', '--model', type=str.lower,
                        default='sticks', choices=VIEW_ALIAS, help=msg)
    parser.add_argument('--pov', action='store_true', default=not Qt_avail,
                        help='Scaling factor for all objects')
    parser.add_argument('-s', '--scale', type=float,
                        help='Scaling factor for all objects')
    parser.add_argument('--scale-atoms', '--scale-at', type=float,
                        help='Scaling factor for all objects')
    parser.add_argument('--scale-bonds', '--scale-bo', type=float,
                        help='Scaling factor for all bonds')
    msg = '''Vibrational mode to display.
NOTE: only available if input file contains eigenvectors.'''
    parser.add_argument('--vib', type=int, help=msg)
    msg = f'''Model to represent vibrations:
- spheres: {', '.join(VIB_MODELS['spheres'])}
- arrows: {', '.join(VIB_MODELS['arrows'])}
'''
    parser.add_argument('--vibmodel', type=str.lower, default='arrows',
                        choices=VIB_ALIAS, help=msg)
    # material?
    # scale


def parse_args(args: tp.Sequence[str]) -> argparse.Namespace:
    """Parse arguments.

    Parameters
    ----------
    args
        Commandline arguments

    Returns
    -------
    :obj:`argparse.Namespace`
        Object holding results as attributes.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    build_opts(parser)
    return parser.parse_args(args)


def extract_data(infile: tp.Sequence[str],
                 tol_bond: float,
                 read_vib: bool = False) -> tp.Dict[str, tp.Any]:
    """Extract data from one input file.

    Extracts atomic data and coordinates, as well as on vibrational
    modes if requested.

    Parameters
    ----------
    infile
        Input file.
    tol_bond
        Tolerance factors for the definition of bonds.
    read_vib
        Read vibrational-mode data.

    Returns
    -------
    dict
        Returns dictionary of data.
    """
    data = {}

    dkeys = {
        'atcrd': QLabel(quantity='AtCrd'),
        'atnum': QLabel(quantity='AtNum')
    }

    dfile = DataFile(infile)
    dobjs = dfile.get_data(**dkeys)
    data['atnum'] = dobjs['atnum'].data
    data['coord'] = np.array(dobjs['atcrd'].data)*PHYSFACT.bohr2ang
    data['bonds'] = list_bonds(data['atnum'], data['coord'], tol_bond)

    if read_vib:
        data['evec'] = dfile.get_hess_data(True, False)[0]

    return data


def main():
    """Run the main program."""
    opts = parse_args(sys.argv[1:])
    read_vib = opts.vib is not None

    # Set representation models
    for key, aliases in VIEW_MODELS.items():
        if opts.model in aliases:
            mol_model = key
            break
    if read_vib:
        for key, aliases in VIB_MODELS.items():
            if opts.vibmodel in aliases:
                vib_model = key
                break

    # Set the scaling factor for the modeling objects
    scale_at = opts.scale_atoms if opts.scale_atoms is not None \
        else opts.scale if opts.scale is not None else 1.0
    scale_bo = opts.scale_bonds if opts.scale_bonds is not None \
        else opts.scale if opts.scale is not None else 1.0

    if opts.display:
        if not Qt_avail:
            print('ERROR: Qt not available. Switching to POV mode.')
            do_pov = True
        else:
            do_pov = False
    else:
        do_pov = True
    num_mols = len(opts.infiles)
    merge = opts.merge and num_mols > 1
    if num_mols:
        mols_atcrd = []
        mols_atnum = []
        mols_bonds = []
        mols_evec = []

    if do_pov:
        if read_vib:
            print('Vibrational modes not yet available for rendering.')
            sys.exit(1)
        mols_zcams = 0
        for fname in opts.infiles:
            moldat = extract_data(fname, opts.bond_tol, read_vib)

            box_mol = build_box(moldat['atnum'], moldat['coord'], True)
            zcam = set_cam_zpos(box_mol)
            if not merge:
                new_file = os.path.splitext(fname)[0] + '.pov'
                write_pov(new_file, 1, moldat['atnum'], moldat['coord'],
                          moldat['bonds'], zcam, True, mol_model, scale_at,
                          scale_bo)
            else:
                mols_atcrd.append(moldat['coord'])
                mols_atnum.append(moldat['atnum'])
                mols_bonds.append(moldat['bonds'])
                if abs(zcam) > abs(mols_zcams):
                    mols_zcams = zcam
        if merge:
            fmt_fname = 'mols_merge_{:03d}.pov'
            i = 1
            while True:
                fpov = fmt_fname.format(i)
                if not os.path.exists(fpov):
                    break
                i = i + 1
                if i > 999:
                    print('Too many files generated. Time to clean up.')
                    sys.exit()

            write_pov(fpov, num_mols, mols_atnum, mols_atcrd, mols_bonds,
                      mols_zcams, True, mol_model, scale_at, scale_bo)
    else:
        if num_mols > 1 and read_vib:
            txt = 'ERROR: Cannot visualize normal modes with multiple '\
                + 'molecules.'
            print(txt)
            sys.exit(1)
        for fname in opts.infiles:
            moldat = extract_data(fname, opts.bond_tol, read_vib)

            if num_mols > 1:
                mols_atcrd.append(moldat['coord'])
                mols_atnum.append(moldat['atnum'])
                mols_bonds.append(moldat['bonds'])
                if read_vib:
                    mols_evec.append(moldat['evec'])
            else:
                mols_atcrd = moldat['coord']
                mols_atnum = moldat['atnum']
                mols_bonds = moldat['bonds']
                if read_vib:
                    mols_evec = moldat['evec']

        app = QtGui.QGuiApplication()
        rad_atom_as_bond = opts.model == 'sticks'
        molwin = MolWin(num_mols, mols_atnum, mols_atcrd, mols_bonds, True,
                        rad_atom_as_bond)
        if read_vib:
            if opts.vib <= 0 or opts.vib > mols_evec.shape[0]:
                print('ERROR: Incorrect normal mode.')
                sys.exit(2)
            molwin.add_vibmode(mols_evec[opts.vib].reshape(-1, 3),
                               repr=vib_model)
        molwin.show()
        sys.exit(app.exec())


if __name__ == '__main__':
    main()
