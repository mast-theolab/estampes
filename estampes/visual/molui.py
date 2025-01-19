"""Build user interfaces for molecular visualization.

Provides classes to build simple or more advanced interfaces to
visualize molecules and associated properties.
"""
import typing as tp

from PySide6 import QtCore, QtGui
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base.types import Type1Vib, TypeAtCrdM, TypeAtLabM, TypeBondsM, \
    TypeColor
from estampes.visual.molview import Molecule
from estampes.visual.vibview import VibMode


class MolWin(Qt3DExtras.Qt3DWindow):
    """The Molecule Window class (Qt3D).

    Generates a Qt3D window instance for the visualization of
    molecule(s).
    """

    click_mol = QtCore.Signal(list)

    def __init__(self, nmols: int,
                 atlabs: TypeAtLabM,
                 atcrds: TypeAtCrdM,
                 bonds: TypeBondsM,
                 col_bond_as_atom: bool = False,
                 rad_atom_as_bond: bool = False,
                 molcols: tp.Optional[TypeColor] = None):
        """Build an instance of MolWin.

        Builds an instance of MolWin.

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
        super(MolWin, self).__init__()

        self.__vib = None
        self.__mol = None
        self.__nmols = 0

        # Camera
        self.camera().lens().setPerspectiveProjection(45, 16 / 9, 0.1, 1000)
        self.camera().setPosition(QtGui.QVector3D(0, 1, 40))
        self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))

        # For camera controls
        self.rootEntity = Qt3DCore.QEntity()
        Molecule.reset_counter()
        if nmols == 1:
            self.__nmols = 1
            self.__mol = Molecule(atlabs, atcrds, bonds, col_bond_as_atom,
                                  rad_atom_as_bond, molcols, self.rootEntity)
            self.__mol.addMouse(self.camera)
            self.__mol.click_molatom.connect(self.atom_clicked)
        else:
            self.__nmols = nmols
            self.__mol = []
            for i in range(nmols):
                if molcols is None:
                    molcol = None
                else:
                    molcol = molcols[i]
                self.__mol.append(Molecule(atlabs[i], atcrds[i], bonds[i],
                                           col_bond_as_atom, rad_atom_as_bond,
                                           molcol, self.rootEntity))
                self.__mol[-1].addMouse(self.camera)
                self.__mol[-1].click_molatom.connect(self.atom_clicked)
        self.camController = Qt3DExtras.QOrbitCameraController(self.rootEntity)
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

    def add_vibmode(self, vib_mode: Type1Vib,
                    repr: tp.Optional[str] = None,
                    color: tp.Optional[tp.Any] = None,
                    color2: tp.Optional[tp.Any] = None,
                    mol: int = -1):
        """Add vibrational mode.

        Creates a VibMode class and populate.

        Parmaters
        ---------
        vib_mode
            List of atomic displacements for a given mode.
        repr
            Representation of the vibration.
        color
            Color of the object representing the displacement.
        color2
            Secondary color of the object representing the displacement.
        mol
            Set target molecule if several of interest.
            Default: last included molecule.
        """
        if self.__vib is not None:
            del self.__vib
        if self.__nmols > 1:
            try:
                atcrd = self.__mol[mol].at_crd
            except IndexError as err:
                raise ValueError('Molecule does not exist.') from err
        else:
            atcrd = self.__mol.at_crd
        self.__vib = VibMode(atcrd, vib_mode, vizmode=repr,
                             color=color, color2=color2,
                             rootEntity=self.rootEntity)

    def upd_vibmode(self, vib_mode: Type1Vib):
        if self.__vib is None:
            self.__add_vibmode(vib_mode)
        else:
            self.__vib.update_vib(vib_mode)

    @QtCore.Slot(list)
    def atom_clicked(self, molatom: tp.Sequence[int]):
        """Capture the information on the clicked atom.

        Captures the information from the clicked atom on one of the
        molecules.

        Parameters
        ----------
        molatom
            Molar atom, as a list of integers.
        """
        self.click_mol.emit(molatom)

    def mousePressEvent(self, mouseEvent):
        """Define the event for mouse pressing."""
        if mouseEvent.button() == QtCore.Qt.RightButton:
            self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))
        # print(mouseEvent.x())
        # self.camera().setViewCenter(QtGui.QVector3D(mouseEvent.x(),
        #                                             mouseEvent.y(), 0))
        super(MolWin, self).mousePressEvent(mouseEvent)
