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
from estampes.data.visual import MODELS, MATERIALS
from estampes.tools.mol import list_bonds
from estampes.visual.povrender import POVBuilder

try:
    from PySide6 import QtGui
    Qt_avail = True
    from estampes.visual.molui import MolWin
except ModuleNotFoundError:
    Qt_avail = False

MOL_ALIAS = sorted([item.lower()
                    for items in MODELS['mol'].values()
                    for item in items['alias']])

VIB_ALIAS = sorted([item.lower()
                    for items in MODELS['vib'].values()
                    for item in items['alias']])

RENDERING = ('display', 'povray', 'interactive', 'D', 'pov', 'I')


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
    # msg = 'Display (Qt) instead of generating the Pov file.'
    # parser.add_argument('-D', '--display', action='store_true',
    #                     help=msg)
    msg = '''Material to be used (for now, only for rendering):
- glass: glassy aspect
- metal: metallic aspect
- plastic: plastic-like effect (default)
'''
    parser.add_argument('--material', choices=MATERIALS.keys(),
                        default='plastic', help=msg)
    msg = 'Merge multiple file in 1 rendition.  ' \
        + 'Each mol. has a different color'
    parser.add_argument('--merge', action='store_true', default=False,
                        help=msg)
    msg = f'''Representation model:
- balls and stick: {', '.join(MODELS['mol']['balls']['alias'])}
- sticks: {', '.join(MODELS['mol']['sticks']['alias'])}
- spheres: {', '.join(MODELS['mol']['spheres']['alias'])}
'''
    parser.add_argument('-m', '--model', type=str.lower,
                        default='sticks', choices=MOL_ALIAS, help=msg)
    parser.add_argument('--mol-mat', choices=MATERIALS.keys(),
                        help='Material for the molecule.')
    msg = '''Rendering model:
- display: 3D interactive display (synonym: interactive, I, D).
- povray: build a PovRay description file (synonym: pov).'''
    parser.add_argument('-o', '--output',
                        help='Output file.')
    parser.add_argument('-r', '--render', choices=RENDERING,
                        default=None, help=msg)
    parser.add_argument('-s', '--scale', type=float,
                        help='Scaling factor for all objects')
    parser.add_argument('--scale-atoms', '--scale-at', type=float,
                        help='Scaling factor for all atoms')
    parser.add_argument('--scale-bonds', '--scale-bo', type=float,
                        help='Scaling factor for all bonds')
    parser.add_argument('--silent', action='store_true', default=False,
                        help='Silent mode, also deactivate help messages.')
    msg = '''Vibrational mode to display.
NOTE: only available if input file contains eigenvectors.'''
    parser.add_argument('--vib', type=int, help=msg)
    parser.add_argument('--vib-mat', choices=MATERIALS.keys(),
                        help='Material for the vibration.')
    msg = f'''Model to represent vibrations:
- spheres: {', '.join(MODELS['vib']['spheres']['alias'])}
- arrows: {', '.join(MODELS['vib']['arrows']['alias'])}
'''
    parser.add_argument('--vib-model', type=str.lower, default='arrows',
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
    # Set mode
    if opts.render is None:
        opmode = 'display' if Qt_avail else 'povray'
    elif opts.render in ('display', 'interactive', 'D', 'I'):
        if not Qt_avail:
            print('ERROR: Qt not available. Switching to POV mode.')
            opmode = 'povray'
        else:
            opmode = 'display'
    else:
        opmode = 'povray'

    read_vib = opts.vib is not None

    # Set representation models
    for key, pars in MODELS['mol'].items():
        if opts.model in pars['alias']:
            mol_model = key
            break
    if read_vib:
        for key, pars in MODELS['vib'].items():
            if opts.vib_model in pars['alias']:
                vib_model = key
                break
    else:
        vib_model = None

    # Set material
    if opts.mol_mat is not None:
        mol_mat = opts.mol_mat
    else:
        mol_mat = opts.material
    if opts.vib_mat is not None:
        vib_mat = opts.vib_mat
    else:
        vib_mat = opts.material

    # Set the scaling factor for the modeling objects
    scale_at = opts.scale_atoms if opts.scale_atoms is not None \
        else opts.scale if opts.scale is not None else 1.0
    scale_bo = opts.scale_bonds if opts.scale_bonds is not None \
        else opts.scale if opts.scale is not None else 1.0

    num_mols = len(opts.infiles)
    merge = opts.merge and num_mols > 1
    if num_mols:
        mols_atcrd = []
        mols_atnum = []
        mols_bonds = []
        mols_evec = []

    if opmode == 'povray':
        if read_vib and num_mols > 1:
            print('Vibrational modes not yet available for rendering.')
            sys.exit(1)
        if merge:
            if opts.output is None:
                fmt_fname = 'mols_merge_{:03d}.pov'
                i = 1
                while True:
                    povfile = fmt_fname.format(i)
                    if not os.path.exists(povfile):
                        break
                    i = i + 1
                    if i > 999:
                        print('Too many files generated. Time to clean up.')
                        sys.exit()
            else:
                povfile = opts.output
            builder = POVBuilder(opts.infiles, load_vibs=read_vib,
                                 tol_bonds=opts.bond_tol)
            builder.write_pov(povname=povfile, merge_mols=True,
                              mol_repr=mol_model, vib_repr=vib_model,
                              mol_mater=mol_mat, vib_mater=vib_mat,
                              scale_atoms=scale_at, scale_bonds=scale_bo,
                              verbose=not opts.silent)
        else:
            if read_vib:
                vib = opts.vib - 1
            else:
                vib = None
            for fname in opts.infiles:
                builder = POVBuilder(fname, load_vibs=read_vib,
                                     tol_bonds=opts.bond_tol)
                if opts.output is None:
                    povfile = os.path.splitext(fname)[0] + '.pov'
                else:
                    povfile = opts.output
                builder.write_pov(povname=povfile, id_vib=vib,
                                  mol_repr=mol_model, vib_repr=vib_model,
                                  mol_mater=mol_mat, vib_mater=vib_mat,
                                  scale_atoms=scale_at, scale_bonds=scale_bo,
                                  verbose=not opts.silent)
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
            molwin.add_vibmode(mols_evec[opts.vib-1].reshape(-1, 3),
                               repr=vib_model)
        molwin.show()
        sys.exit(app.exec())


if __name__ == '__main__':
    main()
