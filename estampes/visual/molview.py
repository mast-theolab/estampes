"""Module related to molecular visualization for ESTAMPES.

This module provides basic methods for the molecular visualization.
"""

from math import sqrt
import typing as tp

import numpy as np

from PySide6 import QtCore, QtGui
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base.types import TypeAtCrd, TypeAtLab, \
    TypeBonds, TypeColor
from estampes.base.errors import ArgumentError
from estampes.data.atom import atomic_data
from estampes.data.visual import BONDDATA, RAD_VIS_SCL
from estampes.tools.math import vrotate_3D


# ==============
# Module Classes
# ==============

class Molecule(Qt3DCore.QEntity):
    """The Molecule class for the visualization.

    Defines a molecule for the visualization.

    Methods
    -------
    update_geom(at_lab, at_crd, bonds, render=True)
        Updates geometry information, and renders the new molecule
    set_display_setting(col_bond_as_atom, rad_atom_as_bond, molcol)
        Sets display settings for the molecule, but does not re-render.
    update_render
        Renders the molecule, with up-to-date internal data
    addMouse(cam)
        Add mouse support.
    """

    click_molatom = QtCore.Signal(list)
    __mol_count = 0

    def __init__(self,
                 at_lab: TypeAtLab,
                 at_crd: TypeAtCrd,
                 bonds: TypeBonds,
                 col_bond_as_atom: bool = False,
                 rad_atom_as_bond: bool = False,
                 molcol: tp.Optional[TypeColor] = None,
                 rootEntity: tp.Optional[Qt3DCore.QEntity] = None):
        """Initialize a Molecule instance.

        Initializes an instance of the Molecule class.

        Parameters
        ----------
        at_lab
            Atomic labels.
        at_crd
            3-tuples with atomic coordinates, in Ang.
        bonds
            2-tuples listing connected atoms.
        col_bond_as_atom
            If true, bonds are colored based on the connected atoms
        rad_atom_as_bond
            If true, atomic radii are set equal to the bonds (tubes).
        molcol
            If not None, the whole molecule is set with the given color.
        rootEntity
            Qt root entity to connect the new `Molecule` QEntity.
        """
        super().__init__(rootEntity)

        type(self).__mol_count += 1
        self.__mol_id = type(self).__mol_count

        # The initialization below should be properly handled through
        # functions.  Just making sure everything declared.
        self.__atlist = self.__atdata = self.__at_mat = self.__at_rad = None
        self.__cam = None
        self.__bo_obj = self.__bo_mesh = self.__bo_trro = None
        self.__at_obj = self.__at_pick = self.__at_mesh = self.__at_tvec = None

        self.set_display_settings(
            col_bond_as_atom=col_bond_as_atom,
            rad_atom_as_bond=rad_atom_as_bond,
            molcol=molcol)
        self.update_geom(at_lab, at_crd, bonds)

    @property
    def at_crd(self):
        """Get atomic coordinates.

        Returns atomic coordinates in the Cartesian space.
        """
        return self.__atcrd

    def update_geom(self,
                    at_lab: TypeAtLab,
                    at_crd: TypeAtCrd,
                    bonds: TypeBonds,
                    render: bool = True):
        """Update geometry information.

        Updates atomic labels and coordinates.

        Parameters
        ----------
        at_lab
            Atomic labels, as string.
        at_crd
            Atomic coordinates, as (N, 3) Numpy array.
        bonds
            List of bonds as `(atom1, atom2)`.
        render
            If True, the molecule is re-rendered.

        Raises
        ------
        IndexError
            Inconsistency size between labels and coordinates.
        """
        if len(at_lab) != len(at_crd):
            raise IndexError('Coordinates do not match atomic labels.')
        self.__atlab = at_lab
        self.__atcrd = at_crd
        self.__bonds = bonds
        self.__upd_atdat()
        if render:
            self.update_render()

    def set_display_settings(self, *,
                             col_bond_as_atom: bool = False,
                             rad_atom_as_bond: bool = False,
                             molcol: tp.Optional[TypeColor] = None):
        """Set display settings for the molecule.

        Sets color information and rendering.

        Parameters
        ----------
        col_bond_as_atom
            Each bond half uses the color of the connected atom.
        rad_atom_as_bond
            Atoms are rendered with the same radii as bonds (tubes).
        molcol
            If not None, use the given color for the whole molecule.

        The internal databases are updated.

        Raises
        ------
        ArgumentError
            Incorrect argument types.
        """
        self.__material = Qt3DExtras.QDiffuseSpecularMaterial
        self.__optview = {
            'col_bond': col_bond_as_atom,
            'rad_atom': rad_atom_as_bond,
        }
        if molcol is None:
            _molcol = None
        elif isinstance(molcol, str):
            if not molcol.startswith('#'):
                raise ArgumentError('Wrong color for molecule.')
            _molcol = [int(molcol[i:i+2], 16) for i in range(1, 6, 2)]
        else:
            _molcol = list(molcol)
        self.__optview['molcol'] = _molcol
        if not col_bond_as_atom:
            self.__bo_mat = self.__material(self)
            self.__bo_mat.setAmbient(QtGui.QColor(*BONDDATA['rgb']))
        else:
            self.__bo_mat = None

        self.__bo_rad = BONDDATA['rvis']*RAD_VIS_SCL

    def update_render(self):
        """Update rendering of the molecules.

        Updates the rendering of the molecules with the current display
          settings.
        """
        self.__build_bonds()
        self.__build_atoms()

    def addMouse(self, cam: Qt3DRender.QCamera):
        """Add mouse support.

        Adds mouse support (clicks) on the molecule object.

        Parameters
        ----------
        cam
            A Qt3Render.QCamera method.
        """
        self.__cam = cam
        for at in self.__at_pick:
            at.pressed.connect(self.__clickAtom)

    def __upd_atdat(self):
        """Update internal atomic data information.

        Builds arrays with unique atomic data:
        - list of unique atoms by alphabetical order
        - list of atomic radii
        - list of atomic textures/materials

        Raises
        ------
        ArgumentError
            Error in colors.
        """
        self.__atlist = sorted(set(self.__atlab))
        self.__atdata = atomic_data(*self.__atlist)
        self.__at_mat = {}
        self.__at_rad = {}

        __molcol = self.__optview['molcol']
        for atom in self.__atlist:
            if __molcol is None:
                r, g, b = self.__atdata[atom]['rgb']
            else:
                r, g, b = __molcol
            self.__at_mat[atom] = self.__material(self)
            self.__at_mat[atom].setAmbient(QtGui.QColor(r, g, b))
            if self.__optview['rad_atom']:
                rval = self.__bo_rad
            else:
                rval = self.__atdata[atom]['rvis']*RAD_VIS_SCL
            self.__at_rad[atom] = rval

    def __build_bonds(self):
        """Build bonds objects and associated properties.

        Builds bonds objects and data in the molecule.
        """
        self.__bo_obj = []
        self.__bo_mesh = []
        self.__bo_trro = []
        for bond in self.__bonds:
            iat1, iat2 = bond
            xyzat1 = np.array(self.__atcrd[iat1])
            xyzat2 = np.array(self.__atcrd[iat2])
            xyzmid = (xyzat1+xyzat2)/2.
            # 1st half of the bond
            # --------------------
            # Initialization of the objects
            bo_obj = Qt3DCore.QEntity(self)
            bo_mesh = Qt3DExtras.QCylinderMesh()
            bo_trro = Qt3DCore.QTransform()
            # Operations
            bo_mesh.setRadius(self.__bo_rad)
            delta = xyzmid - xyzat1
            bo_len = sqrt(np.dot(delta, delta))
            bo_mesh.setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            bo_trro.setRotation(QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            bo_trro.setTranslation(QtGui.QVector3D(*(xyzat1+delta/2)))
            bo_obj.addComponent(bo_mesh)
            bo_obj.addComponent(bo_trro)
            if self.__bo_mat is None:
                bo_obj.addComponent(self.__at_mat[self.__atlab[iat1]])
            else:
                bo_obj.addComponent(self.__bo_mat)
            # Update DB
            self.__bo_obj.append(bo_obj)
            self.__bo_trro.append(bo_trro)
            self.__bo_mesh.append(bo_mesh)

            # 2nd half of the bond
            # --------------------
            # Initialization of the objects
            bo_obj = Qt3DCore.QEntity(self)
            bo_mesh = Qt3DExtras.QCylinderMesh()
            bo_trro = Qt3DCore.QTransform()
            # Operations
            bo_mesh.setRadius(self.__bo_rad)
            delta = xyzat2 - xyzmid
            bo_len = sqrt(np.dot(delta, delta))
            bo_mesh.setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            bo_trro.setRotation(QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            bo_trro.setTranslation(QtGui.QVector3D(*(xyzmid+delta/2)))
            bo_obj.addComponent(bo_mesh)
            bo_obj.addComponent(bo_trro)
            if self.__bo_mat is None:
                bo_obj.addComponent(self.__at_mat[self.__atlab[iat2]])
            else:
                bo_obj.addComponent(self.__bo_mat)
            # Update DB
            self.__bo_obj.append(bo_obj)
            self.__bo_trro.append(bo_trro)
            self.__bo_mesh.append(bo_mesh)

    def __build_atoms(self):
        """Build atom objects and associated properties.

        Builds atoms objects and data in the molecule.
        """
        self.__at_obj = []
        self.__at_pick = []
        self.__at_mesh = []
        self.__at_tvec = []
        for atlab, atxyz in zip(self.__atlab, self.__atcrd):
            # Initialization
            at_obj = Qt3DCore.QEntity(self)
            at_mesh = Qt3DExtras.QSphereMesh()
            at_tvec = Qt3DCore.QTransform()
            # Operations
            at_mesh.setRadius(self.__at_rad[atlab])
            at_tvec.setTranslation(QtGui.QVector3D(*atxyz))
            at_obj.addComponent(at_mesh)
            at_obj.addComponent(at_tvec)
            at_obj.addComponent(self.__at_mat[atlab])
            at_pick = Qt3DRender.QObjectPicker(at_obj)
            at_obj.addComponent(at_pick)
            # Update DB
            self.__at_obj.append(at_obj)
            self.__at_pick.append(at_pick)
            self.__at_mesh.append(at_mesh)
            self.__at_tvec.append(at_tvec)

    def __clickAtom(self, clickEvent: Qt3DRender.QPickEvent):
        """Click Atom events.

        Adds events in case an atom is clicked.

        Parameters
        ----------
        clickEvent
            Qt QPickEvent.
        """
        if clickEvent.button() == Qt3DRender.QPickEvent.RightButton:
            # abs_pos: absolute position of the clicked point
            abs_pos = clickEvent.worldIntersection()
            # loc_pos: local position of the clicked point in the object
            loc_pos = clickEvent.localIntersection()
            # Subtracting them give us the origin of the sphere
            self.__cam().setViewCenter(abs_pos-loc_pos)
        elif clickEvent.button() == Qt3DRender.QPickEvent.LeftButton:
            self.click_molatom.emit(
                [self.__mol_id,
                 self.__at_obj.index(clickEvent.entity())+1])

    @classmethod
    def reset_counter(cls):
        """Reset molecule instance counter."""
        cls.__mol_count = 0
