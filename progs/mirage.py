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
from estampes.tools.char import parse_argval_options
from estampes.tools.mol import list_bonds
from estampes.visual.povrender import POVBuilder

try:
    from PySide6 import QtCore, QtWidgets
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
    msg = '''Show list of modes together with visualizer.'''
    parser.add_argument('--vibs', action='store_true', help=msg)
    parser.add_argument('--vib-mat', choices=MATERIALS.keys(),
                        help='Material for the vibration.')
    parser.add_argument('--vib-scale', type=float, default=1.0,
                        help='Scaling factor to the vibrations')
    msg = f'''Model to represent vibrations:
- spheres: {', '.join(MODELS['vib']['spheres']['alias'])}
- arrows: {', '.join(MODELS['vib']['arrows']['alias'])}
- midarrows: {', '.join(MODELS['vib']['midarrows']['alias'])}
- dualarrows: {', '.join(MODELS['vib']['dualarrows']['alias'])}
'''
    parser.add_argument('--vib-model', help=msg)
    # parser.add_argument('--vib-model', type=str.lower, default='arrows',
    #                     choices=VIB_ALIAS, help=msg)
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
        data['eval'], data['evec'] = dfile.get_hess_data(
            get_evec=False, get_eval=True, get_lweigh=True)
        # data['evec'] = dfile.get_hess_data(get_evec=True, get_eval=False,
        #                                    get_lweigh=False)

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

    # Check if reading vibrational modes requested and check that
    # options are not overlapping.
    if opts.vib is not None and opts.vibs:
        print('ERROR: Option --vib and --vibs are incompatible.')
        sys.exit(1)
    if opts.vibs and opmode != 'display':
        print('ERROR: --vibs only available for interactive viewer')
        sys.exit(1)
    read_vib = opts.vib is not None or opts.vibs

    # Set representation models
    for key, pars in MODELS['mol'].items():
        if opts.model in pars['alias']:
            mol_model = key
            break
    else:
        print('ERROR: Molecular representation model unsupported.')
        print('       This should not happen.')
        sys.exit(1)
    vib_opts = {'col0': None, 'col1': None}
    if read_vib:
        if opts.vib_model is None:
            vib_model = 'arrows'
            args = ()
            kwargs = {}
        else:
            vib_model, args, kwargs = parse_argval_options(
                opts.vib_model, '(', VIB_ALIAS, str.lower)
        for key, pars in MODELS['vib'].items():
            if vib_model in pars['alias']:
                vib_model = key
                break
        if kwargs:
            for key in ('color0', 'color1', 'color+', 'color', 'col1', 'col0',
                        'col+', 'col'):
                if key in kwargs:
                    vib_opts['col0'] = kwargs[key]
                    break
            for key in ('color2', 'color-', 'col2', 'col-'):
                if key in kwargs:
                    vib_opts['col1'] = kwargs[key]
                    break
        if vib_opts['col0'] is None and args:
            vib_opts['col0'] = args[0]
            if vib_opts['col1'] is None and len(args) >= 2:
                vib_opts['col1'] = args[1]
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
        mols_eval = []

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
                              mol_model=mol_model, vib_model=vib_model,
                              mol_mater=mol_mat, vib_mater=vib_mat,
                              scale_atoms=scale_at, scale_bonds=scale_bo,
                              verbose=not opts.silent, vib_opts=vib_opts)
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
                                  mol_model=mol_model, vib_model=vib_model,
                                  mol_mater=mol_mat, vib_mater=vib_mat,
                                  scale_atoms=scale_at, scale_bonds=scale_bo,
                                  verbose=not opts.silent, vib_opts=vib_opts)
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
                    mols_eval.append(moldat['eval'])
            else:
                mols_atcrd = moldat['coord']
                mols_atnum = moldat['atnum']
                mols_bonds = moldat['bonds']
                if read_vib:
                    mols_evec = moldat['evec']
                    mols_eval = moldat['eval']

        if not opts.vibs:
            app = QtWidgets.QApplication()
            molwin = MolWin(num_mols, mols_atnum, mols_atcrd, mols_bonds,
                            mol_model, mol_mat, True, fnames=opts.infiles,
                            skip_guide=opts.silent)
            if read_vib:
                if opts.vib <= 0 or opts.vib > mols_evec.shape[0]:
                    print('ERROR: Incorrect normal mode.')
                    sys.exit(2)
                molwin.add_vibmode(mols_evec[opts.vib-1].reshape(-1, 3),
                                   vib_model, vib_mat, color=vib_opts['col0'],
                                   color2=vib_opts['col1'])
            molwin.show()
            sys.exit(app.exec())
        else:
            app = QtWidgets.QApplication()
            win = QtWidgets.QWidget()
            layout = QtWidgets.QHBoxLayout()
            molwin = MolWin(num_mols, mols_atnum, mols_atcrd, mols_bonds,
                            mol_model, mol_mat, True, fnames=opts.infiles,
                            skip_guide=opts.silent)
            view3D = QtWidgets.QWidget.createWindowContainer(molwin)
            layout.addWidget(view3D, 4)

            table = QtWidgets.QTableWidget()
            table.setRowCount(len(mols_eval))
            table.setColumnCount(1)
            table.setHorizontalHeaderLabels(["Energy"])
            for i, freq in enumerate(mols_eval):
                item_freq = QtWidgets.QTableWidgetItem(
                    f'{freq:12.4f} cm⁻¹')
                item_freq.setTextAlignment(QtCore.Qt.AlignRight)
                table.setItem(i, 0, item_freq)
            table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            table.resizeColumnToContents(0)
            table.horizontalHeader().setSectionResizeMode(
                QtWidgets.QHeaderView.Stretch)
            table.cellClicked.connect(
                lambda row, col: molwin.upd_vibmode(
                    mols_evec[row].reshape(-1, 3), model=vib_model,
                    material=vib_mat, color=vib_opts['col0'],
                    color2=vib_opts['col1'])
            )
            layout.addWidget(table, 1)

            win.setLayout(layout)
            win.setGeometry(300, 200, 800, 600)
            win.show()
            sys.exit(app.exec())


if __name__ == '__main__':
    main()
