#!/usr/bin/env python3

"""
    ESTAMPES: Experimental Support Toolbox for the Analysis, Modelling,
              Plotting and Elucidation of Spectra.

This is the main program for data parser, which also acts as an
  illustration of the toolbox.
"""

# flake8: noqa: F402
# deactivate "F402" regarding modules being at top.
#   We need some tests based on the availability of Qt on machines.

import sys
import os
import argparse
import typing as tp

import numpy as np
import matplotlib.pyplot as plt

# Make PySide optional
try:
    from PySide6 import QtCore, QtGui
    from PySide6.Qt3DCore import Qt3DCore
    from PySide6.Qt3DRender import Qt3DRender
    from PySide6.Qt3DExtras import Qt3DExtras
    QtYes = True
except ModuleNotFoundError:
    QtYes = False

# Deactivate molecular visualization if Qt not available
has_molview = QtYes

from estampes.base import QuantityError, TypeAtCrdM, TypeAtLabM, \
    TypeBondsM, TypeColor, QLabel
from estampes.base.spectrum import VSPC2DATA, Spectrum
from estampes.parser import DataFile
from estampes.data.physics import PHYSFACT
from estampes.tools.atom import convert_labsymb
from estampes.tools.mol import list_bonds
if has_molview:
    from estampes.visual.molview import Molecule, MOLCOLS
from estampes.visual.plotmat import plot_jmat, plot_cmat, plot_kvec
from estampes.visual.plotspec import format_label, plot_spec_2D

FCHT_QTIES = {
   'jmat': {
        'JMat': QLabel(quantity='fcdat', descriptor='JMat')
    },
    'fulljmat': {
        'JFul': QLabel(quantity='fcdat', descriptor='JMatF')
    },
    'cmat': {
        'CMat': QLabel(quantity='fcdat', descriptor='CMat')
    },
    'kvec': {
        'KVec': QLabel(quantity='fcdat', descriptor='KVec')
    },
    'spec': {
        'Spec': QLabel(quantity='fcdat', descriptor='Spec'),
        'Pars': QLabel(quantity='fcdat', descriptor='SpcPar'),
    }
}
if has_molview:
    FCHT_QTIES['mols'] = {
        'atnum': QLabel(quantity='atnum'),
        'IniS': QLabel(quantity='fcdat', descriptor='GeomIS'),
        'FinS': QLabel(quantity='fcdat', descriptor='GeomFS'),
        'MidS': QLabel(quantity='fcdat', descriptor='GeomMS'),
        'ExtG': QLabel(quantity='fcdat', descriptor='ExGeom')
    }


if QtYes and has_molview:
    class MolWin(Qt3DExtras.Qt3DWindow):
        """Qt3D Window for the visualization of molecule(s)
    
        Attributes
        ----------
        nmols
            Number of molecules stored in `atlabs`, `atcrds` and `bonds`.
        atlabs
            Atomic labels.
            If `nmols>1`, list of lists.
        atcrds
            3-tuples with atomic coordinates, in Ang.
            If `nmols>1`, list of lists.
        bonds
            2-tuples listing connected atoms.
            If `nmols>1`, list of lists.
        col_bond_as_atom
            If true, bonds are colored based on the connected atoms
        rad_atom_as_bond
            If true, atomic radii are set equal to the bonds (tubes).
        molcols
            If not `None`, color of the each molecule.
        """
        def __init__(self, nmols: int,
                    atlabs: TypeAtLabM,
                    atcrds: TypeAtCrdM,
                    bonds: TypeBondsM,
                    col_bond_as_atom: bool = False,
                    rad_atom_as_bond: bool = False,
                    molcols: tp.Optional[TypeColor] = None):
            super(MolWin, self).__init__()
    
            # Camera
            self.camera().lens().setPerspectiveProjection(45, 16 / 9, 0.1,
                                                          1000)
            self.camera().setPosition(QtGui.QVector3D(0, 1, 40))
            self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))
    
            # For camera controls
            self.rootEntity = Qt3DCore.QEntity()
            if nmols == 1:
                self.mol = Molecule(atlabs, atcrds, bonds, col_bond_as_atom,
                                    rad_atom_as_bond, molcols, self.rootEntity)
                self.mol.addMouse(self.camera)
            else:
                self.mols = []
                for i in range(nmols):
                    self.mols.append(Molecule(atlabs[i], atcrds[i], bonds[i],
                                              col_bond_as_atom,
                                              rad_atom_as_bond, molcols[i],
                                              self.rootEntity))
                    self.mols[-1].addMouse(self.camera)
            self.camController = \
                Qt3DExtras.QOrbitCameraController(self.rootEntity)
            self.camController.setLinearSpeed(50)
            self.camController.setLookSpeed(180)
            self.camController.setCamera(self.camera())
            self.obj_light = Qt3DCore.QEntity(self.camera())
            self.camLight = Qt3DRender.QPointLight(self.obj_light)
            self.cam_tvec = Qt3DCore.QTransform()
            self.cam_tvec.setTranslation(QtGui.QVector3D(0, 50, 100))
            self.obj_light.addComponent(self.camLight)
            self.obj_light.addComponent(self.cam_tvec)
            # self.camLight.setIntensity(100)
            # for mol in mols:
            self.setRootEntity(self.rootEntity)
    
        def mousePressEvent(self, mouseEvent):
            if (mouseEvent.button() == QtCore.Qt.RightButton):
                self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))
            # print(mouseEvent.x())
            # self.camera().setViewCenter(QtGui.QVector3D(mouseEvent.x(),
            #                                             mouseEvent.y(), 0))
            super(MolWin, self).mousePressEvent(mouseEvent)


def parse_args(args: tp.Sequence[str]) -> argparse.Namespace:
    """Parses arguments.

    Parses commandline arguments

    Parameters
    ----------
    args
        Commandline arguments

    Returns
    -------
    :obj:`argparse.Namespace`
        Object holding results as attributes
    """
    msg_desc = '''A simple script to parse files supported by ESTAMPES.
This is also used as a proof-of-concept, implementation test of the library.'''
    msg_epilog = ''''''
    if not has_molview:
        msg_epilog += '''
WARNING: The Qt library (PySide6) was not found.  \
Molecular visualization is not available.'''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=msg_desc, epilog=msg_epilog)
    # Basic
    psubs = parser.add_subparsers(help='Graphical interface (NYI)')
    parser.set_defaults(mode='gui')

    # MolView
    if has_molview:
        pmol = psubs.add_parser('molview', aliases=['mol'],
                                help='Molecular viewer')
        pmol.add_argument('datafile',
                          help='Data file.')
        pmol.set_defaults(mode='mol')

    # Vibrational spectroscopy
    pvib = psubs.add_parser('vibrational', aliases=['vib', 'l717'],
                            help='Vibrational spectroscopy')
    pvib.add_argument('datafile', help='Data file.')
    pvib.add_argument('-b', '--broaden',
                      choices=('gaussian', 'lorentzian', 'stick'),
                      default='stick',
                      help='Broadening function')
    pvib.add_argument('-l', '--level', choices=('H', 'A'),
                      default='A',
                      help='Vibrational level of theory: H(arm), A(nharm)')
    pvib.add_argument('-o', '--output',
                      help='Output file.')
    pvib.add_argument('-s', '--spec', choices=VSPC2DATA.keys(),
                      default='IR',
                      help='Type of spectroscopy')
    pvib.add_argument('--xmin', type=float,
                      default=0.0,
                      help='Lower bound')
    pvib.add_argument('--xmax', type=float,
                      default=4000.0,
                      help='Upper bound')
    pvib.add_argument('--xres', type=float,
                      default=5.0,
                      help='Resolution')
    pvib.add_argument('-w', '--hwhm', type=float,
                      default=5.0,
                      help='Half-width at half-maximum')
    pvib.set_defaults(mode='vib')

    # Vibrationally-resolved electronic spectroscopy
    pvel = psubs.add_parser('vibronic', aliases=['FCHT', 'l718'],
                            help='Vibronic specroscopy')
    pvel.add_argument('-o', '--output',
                      help='Output file.')
    fmt = 'Quantity to show.  Possible values:\n{}'
    pvel.add_argument('-q', '--quantity', type=str.lower,
                      choices=[item.lower() for item in FCHT_QTIES],
                      default='mols' if 'mols' in FCHT_QTIES else 'jmat',
                      help=fmt.format(', '.join([item.lower()
                                                 for item in FCHT_QTIES])))
    pvel.add_argument('-t', '--title',
                      help='Title for the figure.')
    pvel.add_argument('datafile',
                      help='Data file.')
    pvel.set_defaults(mode='vel')
    # Vibrationally-resolved electronic spectroscopy
    pspc = psubs.add_parser('spectra', aliases=['spec', 'spc'],
                            help='Printing and comparison of spectra')
    pspc.add_argument('optfile', nargs='?',
                      help='Option file (INI style).')
    # pspc.add_argument('-o', '--output',
    #                   help='Output file.')
    pspc.add_argument('-c', '--colors', action='append',
                      help='Spectral colors.')
    msg = '''\
Colors of the spectra.  By default, it follows the order of input files.
It is possible to change the order by putting a number followed by ":".
Ex. '3:Test' means that the label 'Test' is for the 3rd file (start at 1).
'r'/'e'/'0' refers to the reference data.
'''
    pspc.add_argument('-i', '--inpfile', action='append',
                      help='Input data file.')
    msg = '''\
Labels for the legend.  By default, it follows the order of input files.
It is possible to change the order by putting a number followed by ":".
Ex. '3:Test' means that the label 'Test' is for the 3rd file (start at 1).
'r'/'e'/'0' refers to the reference data.
'''
    pspc.add_argument('-l', '--label', action='append',
                      help=msg)
    pspc.add_argument('-r', '--refdata',
                      help='Reference spectrum file.')
    fmt = '''Type of spectra:
{}
Default for vibrational spectroscopy: IR
Default for vibronic spectroscopy is extracted from the log.'''
    spc = ('auto', 'IR', 'VCD', 'Raman', 'ROA', 'OPA', 'OPE', 'ECD', 'CPL',
           'RR', 'RROA')
    pspc.add_argument('-t', '--type', choices=spc, default='auto',
                      help=fmt.format(', '.join(spc)))
    pspc.set_defaults(mode='spc')
    # Vibrationally-resolved electronic spectroscopy
    pgui = psubs.add_parser('gui', aliases=['GUI', 'main'],
                            help='General interface')
    pgui.add_argument('datafile', nargs='?', action='append',
                      help='Data file.')
    pgui.set_defaults(mode='gui')

    return parser.parse_args(args)


def mode_molview(dfile: DataFile):
    """Molview Mode.

    Main function managing molecule viewer.

    Parameters
    ----------
    dfile
        `ep.DataFile` object.
    """
    if not QtYes:
        print('Missing PySide.  Cannot display.')
        sys.exit(1)
    dkeys = {
        'atcrd': QLabel(quantity='atcrd', descriptor='last'),
        'atnum': QLabel(quantity='atnum')
    }
    dobjs = dfile.get_data(**dkeys)
    atlab = convert_labsymb(True, *dobjs['atnum'].data)
    atcrd = np.array(dobjs['atcrd'].data)*PHYSFACT.bohr2ang
    bonds = list_bonds(atlab, atcrd, 1.2)

    app = QtGui.QGuiApplication(sys.argv)
    view = MolWin(1, atlab, atcrd, bonds, True, False)
    view.show()
    sys.exit(app.exec_())


def mode_vibronic(dfile: DataFile,
                  qty: str,
                  **kwargs: tp.Dict[str, tp.Any]):
    """Vibrational-spectroscopy Mode.

    Main function managing vibratonal spectroscopy.

    Parameters
    ----------
    dfile
        `ep.DataFile` object.
    qty
        Quantity of interest.
    kwargs
        Keyword-type arguments.
    """
    # For geometries, some data may not be available
    error_noqty = qty != 'mols'
    dkeys = FCHT_QTIES[qty]
    try:
        dobjs = dfile.get_data(**dkeys.values(), error_noqty=error_noqty)
    except IndexError:
        print('Data not available in file.')
        sys.exit()
    except QuantityError:
        print('Quantity not supported. Someone was lazy...')
        sys.exit(1)
    if qty == 'mols':
        if dobjs['IniS'] is None:
            print('Data not available in file.')
            sys.exit()
        # Check with geometry to use for the final state
        # If extrapolated geometry available (VH, VG), uses it
        # If intermediate state defined, RR, so use it
        # Otherwise, use standard final state definition.
        if dobjs['ExtG'] is not None:
            fs = 'ExtG'
        elif dobjs['MidS'] is not None:
            fs = 'MidS'
        else:
            fs = 'FinS'
        if not dobjs[fs].data:
            print('ERROR: Something went wrong, final-state geom. missing.')
            sys.exit()
        atlabs = []
        atcrds = []
        bonds = []
        molcols = []
        i = 0
        for sta in ('IniS', fs):
            atlabs.append(convert_labsymb(True, *dobjs['atnum'].data))
            atcrds.append(np.array(dobjs[sta].data)*PHYSFACT.bohr2ang)
            bonds.append(list_bonds(atlabs[-1], atcrds[-1], 1.2))
            molcols.append(MOLCOLS[i])
            i += 1
        app = QtGui.QGuiApplication(sys.argv)
        view = MolWin(2, atlabs, atcrds, bonds, True, True, molcols)
        view.show()
        sys.exit(app.exec_())
    else:
        if qty == 'jmat' and not dobjs['JMat'].data:
            # First check that it is not a reduced-dimensionality case
            #   with only the full matrix printed and not the reddim one.
            qkey = FCHT_QTIES['fulljmat']
            key = qkey.keys()[0]
            dobj2 = dfile.get_data(**qkey, error_noqty=error_noqty)
            if dobj2[key].data:
                print('''J is missing. Only the full matrix is available.
If you want to display that one, use "fulljmat" instead.''')
            else:
                print('J is the identity matrix.')
        else:
            figsize = (10, 8)
            fig, subp = plt.subplots(1, 1, tight_layout=True)
            fig.set_size_inches(figsize)
            if qty == 'jmat':
                mat = np.array(dobjs['JMat'].data)
                plot = plot_jmat(mat, subp)
                fig.colorbar(plot)
            elif qty == 'fulljmat':
                mat = np.array(dobjs['JFul'].data)
                plot = plot_jmat(mat, subp)
                fig.colorbar(plot)
            elif qty == 'cmat':
                mat = np.array(dobjs['CMat'].data)
                norm, plot = plot_cmat(mat, subp)
                print(f'Normalization factor: {norm:15.6e}')
                fig.colorbar(plot)
            elif qty == 'kvec':
                mat = np.array(dobjs['KVec'].data)
                plot = plot_kvec(mat, subp)
            elif qty == 'spec':
                if 'y1' in dobjs['Pars'].extra_fields():
                    leg = dobjs['Pars'].extra_fields()
                else:
                    leg = None
                stick = dobjs['Pars'].get('func').lower() == 'stick'
                _ = plot_spec_2D(dobjs['Spec'].extra_fields(), subp,
                                 legends=leg, is_stick=stick)
            if kwargs['title'] is not None:
                subp.set_title(kwargs['title'])
            plt.show()


def mode_vibspec(dfile: DataFile,
                 **kwargs: tp.Dict[str, tp.Any]):
    """Vibrational-spectroscopy Mode.

    Main function managing vibratonal spectroscopy.

    Parameters
    ----------
    dfile
        DataFile object.
    kwargs
        Keyword-type arguments.
    """

    # For geometries, some data may not be available
    data = Spectrum(dfile, kwargs['spec'], kwargs['level'], None)
    stick = kwargs['broaden'] == 'stick'
    if not stick:
        data.set_broadening(kwargs['hwhm'], kwargs['broaden'], 'default',
                            kwargs['xres'], kwargs['xmin'], kwargs['xmax'])
    outfile = kwargs['output'] or False
    if outfile:
        ext = os.path.splitext(outfile)[1][1:].lower()
        if ext in ('csv', 'xy', 'txt'):
            save_img = False
            fmt = '{:12.5f}, {:15.6e}\n'
            with open(outfile, 'w', encoding='utf-8') as fobj:
                for x, y in zip(data.xaxis, data.yaxis):
                    fobj.write(fmt.format(x, y))
        else:
            save_img = True
    figsize = (10, 8)
    fig, subp = plt.subplots(1, 1, tight_layout=True)
    fig.set_size_inches(figsize)
    xunit = format_label(data.xunit)
    yunit = format_label(data.yunit)
    _ = plot_spec_2D({'x': data.xaxis, 'y': data.yaxis}, subp,
                     xlabel=xunit, ylabel=yunit, is_stick=stick)
    if save_img:
        plt.savefig(outfile, bbox_inches='tight')
    # if kwargs['title'] is not None:
    #     subp.set_title(kwargs['title'])
    plt.show()


def mode_spectra(optfile: tp.Optional[str] = None):
    """
    """


def main() -> tp.NoReturn:
    """Main function.
    """
    args = parse_args(sys.argv[1:])
    if args.mode == 'gui':
        print('ERROR: Not yet available')
        sys.exit(1)
    elif args.mode == 'mol':
        fname = args.datafile
        if not os.path.exists(fname):
            print(f'ERROR: File "{fname}"" not found.')
            sys.exit(2)
        dfile = DataFile(fname)
        mode_molview(dfile)
    elif args.mode == 'vel':
        fname = args.datafile
        if not os.path.exists(fname):
            print(f'ERROR: File "{fname}"" not found.')
            sys.exit(2)
        dfile = DataFile(fname)
        mode_vibronic(dfile, args.quantity, **args.__dict__)
    elif args.mode == 'vib':
        fname = args.datafile
        if not os.path.exists(fname):
            print(f'ERROR: File "{fname}"" not found.')
            sys.exit(2)
        dfile = DataFile(fname)
        mode_vibspec(dfile, **args.__dict__)
    elif args.mode == 'spc':
        if not args.inpfile and not args.optfile:
            print('ERROR: Missing files or option file.')
            sys.exit(2)
        elif args.inpfile and args.optfile:
            msg = 'ERROR: Option file and single files cannot be treated' \
                + ' together'
            print(msg)
            sys.exit(2)
        else:
            msg = '''\
NYI: The spectra-centric mode is not yet implemented.
     For simple plotting, check module-specific modes, like 'l718', 'l717'...
'''
            print(msg)
            sys.exit(1)
    else:
        print('ERROR: Unsupported feature "{}"'.format(args.mode))
        sys.exit(1)


if __name__ == '__main__':
    main()
