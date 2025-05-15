"""Program: SOAR.

SOAR: Structural Overlap Analysis and Representation.

Simple program to analyze differences in structures.
"""
import sys
import argparse
import typing as tp

import numpy as np

# import matplotlib.pyplot as plt

from estampes.base import QLabel
from estampes.parser import DataFile
from estampes.tools.math import superpose
from estampes.tools.mol import eckart_orient
from estampes.tools.vib import build_dusch_J, build_dusch_K
from estampes.visual.plotmat import plot_jmat, plot_kvec

try:
    from estampes.data.physics import PHYSFACT
    from estampes.data.visual import MOLCOLS
    from estampes.tools.mol import list_bonds
    # from PySide6 import QtWidgets
    from matplotlib.backends.backend_qtagg import FigureCanvas
    from matplotlib.backends.backend_qtagg import \
        NavigationToolbar2QT as NavigationToolbar
    from matplotlib.backends.qt_compat import QtWidgets
    from matplotlib.figure import Figure
    from estampes.visual.molui import MolWin
    QT_AVAIL = True
except ModuleNotFoundError:
    QT_AVAIL = False

PROGNAME = 'soar.py'


def parse_args(args: tp.Optional[tp.Sequence[str]] = None
               ) -> argparse.Namespace:
    """Define and check options given in arguments.

    Defines available options in the program, parses a list of arguments
    given in input, and returns the results as a populated namespace.

    Parameters
    ----------
    args
        List of options/arguments to parse.

    Returns
    -------
    obj:`argparse.Namespace`
        Returned results from the parsing.
    """
    parser = argparse.ArgumentParser(
        prog=PROGNAME,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('reffile',
                        help='File with data from the reference structure.')
    parser.add_argument('newfile',
                        help='File with data from the other structure.')
    if QT_AVAIL:
        msg = 'Do not show the overlapped structures.'
        parser.add_argument('--no-mols', action='store_false',
                            dest='show_mols', default=None, help=msg)
    msg = 'Do not plot the matrices and vectors.'
    parser.add_argument('--no-plots', action='store_false', dest='show_plots',
                        default=None, help=msg)
    msg = 'Show the overlapped structures.'
    if QT_AVAIL:
        parser.add_argument('--show-mols', action='store_true',
                            dest='show_mols', default=None, help=msg)
    msg = 'Plot the matrices and vectors (default).'
    parser.add_argument('--show-plots', action='store_true', dest='show_plots',
                        default=None, help=msg)
    msg = 'Superpose the structure before checking the overlap.'
    parser.add_argument('-s', '--superpose', action='store_true',
                        help=msg)

    if args is not None:
        opts = parser.parse_args(args)
    else:
        opts = parser.parse_args()

    return opts


def main():
    """Run the main program."""
    dopts = parse_args()

    # Set defaults
    if dopts.show_plots is None:
        dopts.show_plots = True
    if dopts.show_mols is None:
        dopts.show_mols = False
    show_ui = dopts.show_plots or dopts.show_mols

    dfile_ref = DataFile(dopts.reffile)
    dfile_new = DataFile(dopts.newfile)

    keys = {
        # 'E': QLabel(quantity=1),
        'mass': QLabel(quantity='AtMas'),
        'num': QLabel(quantity='AtNum'),
        'crd': QLabel(quantity='AtCrd')
    }
    Lmat_ref, freq_ref = dfile_ref.get_hess_data()
    Lmat_new, freq_new = dfile_new.get_hess_data()

    data_ref = dfile_ref.get_data(**keys)
    data_new = dfile_new.get_data(**keys)

    at_mass = np.array(data_ref['mass'].data)
    c_ref = np.array(data_ref['crd'].data)
    c_new = np.array(data_new['crd'].data)
    n_vib = Lmat_ref.shape[0]
    if dopts.superpose:
        # Put first molecule in Eckart orientation and rotate the eigenvectors
        c_ref, rotmat = eckart_orient(c_ref, at_mass, True).values()
        Lmat_ref = np.reshape(np.reshape(Lmat_ref, (n_vib, -1, 3)) @ rotmat,
                              (n_vib, -1))
        # Superpose second geometry on top of first
        rotmat, _, c_new = superpose(c_ref, c_new, at_mass, get_ctrans=True)
        Lmat_new = np.reshape(np.reshape(Lmat_new, (n_vib, -1, 3)) @ rotmat,
                              (n_vib, -1))

    jmat = build_dusch_J(Lmat_new, Lmat_ref)
    kvec = build_dusch_K(Lmat_ref, at_mass, c_new - c_ref)

    print('''\
 num. | Ref. freq. |  New freq.
 -----+------------+-----------''')
    i = 0
    for x, y in zip(freq_ref, freq_new):
        i += 1
        print(' {:4d} | {:10.4f} | {:10.4f}'.format(i, x, y))

    if show_ui:
        qapp = QtWidgets.QApplication.instance()
        if not qapp:
            qapp = QtWidgets.QApplication()
        qwin = QtWidgets.QWidget()
        main_layout = QtWidgets.QHBoxLayout(qwin)

        if dopts.show_plots:
            figsize = (10, 8)
            layout = QtWidgets.QVBoxLayout()
            fig = Figure(figsize=figsize)
            canvas = FigureCanvas(fig)
            layout.addWidget(NavigationToolbar(canvas))
            layout.addWidget(canvas)
            subp = fig.subplots(2, 1)
            fig.set_size_inches(figsize)
            _ = plot_jmat(jmat, subp[0])
            _ = plot_kvec(kvec, subp[1])
            # fig.colorbar(plot)
            main_layout.addLayout(layout, 50)

        if dopts.show_mols:
            mols_atcrd = [c_ref*PHYSFACT.bohr2ang, c_new*PHYSFACT.bohr2ang]
            mols_atnum = [data_ref['num'].data, data_new['num'].data]
            mols_bonds = [list_bonds(mols_atnum[0], mols_atcrd[0], 1.2),
                          list_bonds(mols_atnum[1], mols_atcrd[1], 1.2)]
            mols_cols = MOLCOLS[:2]
            molwin = MolWin(2, mols_atnum, mols_atcrd, mols_bonds,
                            'sticks', 'plastic', True, mols_cols,
                            skip_guide=True)
            view3D = QtWidgets.QWidget.createWindowContainer(molwin)
            main_layout.addWidget(view3D, 50)
            # Force some minimum size of the window if only 3D view, otherwise
            # it may end up being minimized.
            if not dopts.show_plots:
                qwin.setGeometry(300, 200, 800, 600)

        qwin.show()

        sys.exit(qapp.exec())


if __name__ == '__main__':
    main()
