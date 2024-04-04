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
from estampes.visual.molview import build_box, set_cam_zpos, write_pov


VIEW_MODELS = {
    'sticks': ('stick', 'sticks'),
    'balls': ('ball', 'balls', 'ballsandsticks'),
    'spheres': ('sphere', 'spheres', 'vdW'),
}

VIEW_ALIAS = sorted([item.lower()
                     for items in VIEW_MODELS.values()
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
    # msg = 'Display (Qt) instead of generating the Pov file.'
    # parser.add_argument('-D', '--display', action='store_true',
    #                     default=False,
    #                     help=msg)
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
    parser.add_argument('-s', '--scale', type=float,
                        help='Scaling factor for all objects')
    parser.add_argument('--scale-atoms', '--scale-at', type=float,
                        help='Scaling factor for all objects')
    parser.add_argument('--scale-bonds', '--scale-bo', type=float,
                        help='Scaling factor for all bonds')
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


def main():
    """Run the main program."""
    opts = parse_args(sys.argv[1:])

    # Set representation model
    for key, aliases in VIEW_MODELS.items():
        if opts.model in aliases:
            model = key
            break

    # Set the scaling factor for the modeling objects
    scale_at = opts.scale_atoms if opts.scale_atoms is not None \
        else opts.scale if opts.scale is not None else 1.0
    scale_bo = opts.scale_bonds if opts.scale_bonds is not None \
        else opts.scale if opts.scale is not None else 1.0

    dkeys = {
        'atcrd': QLabel(quantity='AtCrd'),
        'atnum': QLabel(quantity='AtNum')
    }

    # if opts.display:
    #     if not Qt_avail:
    #         print('ERROR: Qt not available. Switching to POV mode.')
    #         do_pov = True
    #     else:
    #         do_pov = False
    # else:
    #     do_pov = True
    do_pov = True
    merge = opts.merge and len(opts.infiles) > 1
    if merge:
        mols_atcrd = []
        mols_atnum = []
        mols_bonds = []
        mols_zcams = 0
        num_mols = len(opts.infiles)
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
    if do_pov:
        for fname in opts.infiles:
            dfile = DataFile(fname)
            moldat = dfile.get_data(**dkeys)
            atnum = moldat['atnum'].data
            atcrd = np.array(moldat['atcrd'].data)*PHYSFACT.bohr2ang

            bonds = list_bonds(atnum, atcrd, opts.bond_tol)
            box_mol = build_box(atnum, atcrd, True)
            zcam = set_cam_zpos(box_mol)
            if not merge:
                new_file = os.path.splitext(fname)[0] + '.pov'
                write_pov(new_file, 1, atnum, atcrd, bonds, zcam, True, model,
                          scale_at, scale_bo)
            else:
                mols_atcrd.append(atcrd)
                mols_atnum.append(atnum)
                mols_bonds.append(bonds)
                if abs(zcam) > abs(mols_zcams):
                    mols_zcams = zcam
        if merge:
            write_pov(fpov, num_mols, mols_atnum, mols_atcrd, mols_bonds,
                      mols_zcams, True, model, scale_at, scale_bo)
    else:
        print('Graphical representation NYI')
        # if len(opts.infiles) > 1:
        #     print('ERROR: Only 1 XYZ file supported in display mode.')
        #     sys.exit()
        # app = QtGui.QGuiApplication(sys.argv)
        # fname = opts.infiles[0]
        # mol = parse_xyz(fname)
        # bonds = list_bonds(mol['atoms'])
        # view = Window(mol['atoms'], bonds, True, False)
        # view.show()
        # sys.exit(app.exec_())


if __name__ == '__main__':
    main()
