"""Module related to molecular visualization for ESTAMPES.

This module provides basic methods for the molecular visualization.
"""

from math import sqrt
import typing as tp

import numpy as np
import numpy.typing as npt

from PySide6 import QtCore, QtGui
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base.types import TypeAtCrd, TypeAtLab, \
    TypeBonds, TypeColor
from estampes.base.errors import ArgumentError
from estampes.data.atom import atomic_data
from estampes.data.visual import BONDDATA, MATERIALS, MODELS, RAD_VIS_SCL
from estampes.tools.math import vrotate_3D


# ==============
# Module Classes
# ==============

class Molecule(Qt3DCore.QEntity):
    """The Molecule class for the visualization.

    Defines a molecule for the visualization.
    """

    click_molatom = QtCore.Signal(list)
    __mol_count = 0

    def __init__(self,
                 at_lab: TypeAtLab,
                 at_crd: TypeAtCrd,
                 bonds: TypeBonds,
                 model: tp.Optional[str] = None,
                 material: tp.Optional[str] = None,
                 col_bond_as_atom: tp.Optional[bool] = None,
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
        model
            Representation model of the molecule.
        material
            Material simulated for the molecule.
        col_bond_as_atom
            If true, bonds are colored based on the connected atoms
        molcol
            If not None, the whole molecule is set with the given color.
        rootEntity
            Qt root entity to connect the new `Molecule` QEntity.
        """
        super().__init__(rootEntity)

        type(self).__mol_count += 1
        self.__mol_id = type(self).__mol_count

        # Initialize dictionary of options for visualization
        self.__optview = {
            'model': None,
            'material': None,
            'col_bond': None,
            'molcol': None
        }
        # The initialization below should be properly handled through
        # functions.  Just making sure everything declared.
        self.__atlist = self.__atdata = self.__at_mat = self.__at_rad = None
        self.__cam = None
        self.__bo_obj = self.__bo_mesh = self.__bo_trro = None
        self.__at_obj = self.__at_pick = self.__at_mesh = self.__at_tvec = None
        # Atomic displacement for animation
        self.__at_displ = self.__at_displ_step = None
        self.vibrate = None

        self.set_display_settings(
            model=model,
            material=material,
            col_bond_as_atom=col_bond_as_atom,
            molcol=molcol)
        self.update_geom(at_lab, at_crd, bonds)

    def get_at_crd(self, actual: bool = False) -> TypeAtCrd:
        """Return atomic coordinates.

        Returns atomic coordinates in the Cartesian space.

        Parameters
        ----------
        actual
            Print the actual coordinates, including atomic shift if set.
        """
        if actual:
            if self.__at_displ_step is not None:
                return self.__atcrd + self.__at_displ_step*self.__at_displ
            return self.__atcrd
        return self.__atcrd

    at_crd = property(get_at_crd)

    @property
    def at_lab(self) -> TypeAtLab:
        """Return atomic labels.

        Returns atomic labels.
        """
        return self.__atlab

    def viz_options(self, key: tp.Optional[str] = None) -> tp.Optional[tp.Any]:
        """Return one or more visual-centric options

        Returns information specific to a given key or the whole array.
        """
        if key is None:
            return self.__optview
        else:
            if key in ('molcol', 'color'):
                return self.__optview.get('molcol')
            elif key in ('model', 'representation'):
                return self.__optview.get('model')
            elif key in ('material',):
                return self.__optview.get('material')
            elif key in ('col_bond', 'col_bond_as_atom'):
                return self.__optview.get('col_bond')
            else:
                return None

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
                             model: tp.Optional[str] = None,
                             material: tp.Optional[str] = None,
                             col_bond_as_atom: tp.Optional[bool] = None,
                             molcol: tp.Optional[TypeColor] = None):
        """Set display settings for the molecule.

        Sets color information and rendering.

        Parameters
        ----------
        model
            Representation model of the molecule.
        material
            Material simulated for the molecule.
        col_bond_as_atom
            Each bond half uses the color of the connected atom.
        molcol
            If not None, use the given color for the whole molecule.

        The internal databases are updated.

        Raises
        ------
        ArgumentError
            Incorrect argument types.
        """
        # Check options with specific choices
        if model is None:
            model_key = 'balls'
        else:
            for key, val in MODELS['mol'].items():
                if model in val['alias']:
                    model_key = key
                    break
            else:
                raise ArgumentError('model', 'Unrecognized molecular model')
        self.__optview['model'] = model_key

        if material is None:
            material_key = 'plastic'
        else:
            if material not in MATERIALS:
                raise ArgumentError('material', 'Unrecognized material')
            else:
                material_key = material
        self.__optview['material'] = material_key

        self.__opacity = 255
        if material_key == 'glass':
            self.__opacity = 180

        if molcol is None:
            _molcol = None
            # If there is no overall color molecule, then check how to handle
            # bond colors.
            self.__optview['col_bond'] = col_bond_as_atom
        elif isinstance(molcol, str):
            if not molcol.startswith('#'):
                raise ArgumentError('Wrong color for molecule.')
            _molcol = [int(molcol[i:i+2], 16) for i in range(1, 6, 2)]
        else:
            _molcol = list(molcol)
        self.__optview['molcol'] = _molcol
        if not col_bond_as_atom:
            self.__bo_mat = self.__set_material(BONDDATA['rgb'],
                                                self.__opacity)
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

        for atom in self.__atlist:
            if self.__optview['molcol'] is None:
                r, g, b, a = *self.__atdata[atom]['rgb'], \
                              self.__opacity
            else:
                r, g, b, a = *self.__optview['molcol'], self.__opacity
            self.__at_mat[atom] = self.__set_material([r, g, b], a)
            if self.__optview['model'] == 'balls':
                rval = self.__atdata[atom]['rvis']*RAD_VIS_SCL
            elif self.__optview['model'] == 'sticks':
                rval = self.__bo_rad
            elif self.__optview['model'] == 'spheres':
                rval = self.__atdata[atom]['rvdw']
            else:
                raise KeyError('Unsupported model for the radius setup')
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

    def animate(self,
                displ: tp.Optional[npt.ArrayLike] = None,
                vibrate: bool = True,
                amplitude: float = 1.0,
                duration: int = 1000):
        """Move atoms and correct bonds along chosen displacement.

        Creates an animation of the atoms displacement.

        Parameters
        ----------
        displ
            Displacement array.  If None, use existing one.
        vibrate
            Consider symmetric displacements.
        amplitude
            Maximum amplitude of the displacement.
            With vibrate, the same value is used in both directions.
        duration
            Duration of a 1 period (in ms).
        """
        if self.vibrate is not None:
            if self.vibrate.state() == QtCore.QAbstractAnimation.Running:
                self.vibrate.pause()
                return
            elif (self.vibrate.state() == QtCore.QAbstractAnimation.Paused
                    and displ is None):
                self.vibrate.resume()
                return
        if displ is not None:
            displ = np.asarray(displ)
            shape = displ.shape
            if len(shape) == 1:
                self.__at_displ.reshape(-1, 3)
            elif len(shape) != 2:
                raise ArgumentError(
                    'displ',
                    'Too many dimension in displacement parameters')
            else:
                self.__at_displ = displ.copy()
            if self.__at_displ.shape[0] != len(self.__at_obj):
                raise IndexError(
                    'Number of displacements inconsistent with molecule')

        self.vibrate = QtCore.QPropertyAnimation(self, b"shift_atoms")
        self.vibrate.setStartValue(0)
        if vibrate:
            self.vibrate.setKeyValueAt(.25, -abs(amplitude))
            self.vibrate.setKeyValueAt(.75, abs(amplitude))
            self.vibrate.setEndValue(0)
            self.vibrate.setDuration(abs(int(duration)))
            self.vibrate.setLoopCount(-1)
        else:
            self.vibrate.setEndValue(1)
            self.vibrate.setDuration(abs(int(duration)))
        self.vibrate.start()

    @property
    def animation_status(self) -> str:
        """Return the animation status."""
        if self.vibrate is None:
            status = 'unknown'
        else:
            if self.vibrate.state() == QtCore.QAbstractAnimation.Running:
                status = 'running'
            elif self.vibrate.state() == QtCore.QAbstractAnimation.Paused:
                status = 'paused'
            elif self.vibrate.state() == QtCore.QAbstractAnimation.Stopped:
                status = 'stopped'
            else:
                raise ValueError('Unknown animation state.')
        return status

    @animation_status.setter
    def animation_status(self, status: str):
        """Set the animation status."""
        if self.vibrate is not None:
            if status.lower() == 'running':
                if self.vibrate.state() == QtCore.QAbstractAnimation.Paused:
                    self.vibrate.resume()
                elif self.vibrate.state() == QtCore.QAbstractAnimation.Stopped:
                    self.vibrate.start()
            elif status.lower() == 'stopped':
                self.vibrate.stop()
            elif status.lower() == 'paused':
                self.vibrate.pause()
            elif status.lower() == 'reset':
                self.vibrate.stop()
                self.vibrate.setCurrentTime(0)

    def redraw(self):
        """Simply redraw the molecule with the default positions."""
        for atvec, atxyz in zip(self.__at_tvec, self.__atcrd):
            atvec.setTranslation(QtGui.QVector3D(*atxyz))

        for i, bond in enumerate(self.__bonds):
            iat1, iat2 = bond
            xyzat1 = self.__atcrd[iat1]
            xyzat2 = self.__atcrd[iat2]
            xyzmid = (xyzat1+xyzat2)/2.
            # First half of bond
            delta = xyzmid - xyzat1
            bo_len = sqrt(np.dot(delta, delta))
            self.__bo_mesh[2*i].setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            self.__bo_trro[2*i].setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            self.__bo_trro[2*i].setTranslation(
                QtGui.QVector3D(*(xyzat1+delta/2)))
            # Second half of bond
            delta = xyzat2 - xyzmid
            bo_len = sqrt(np.dot(delta, delta))
            self.__bo_mesh[2*i+1].setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            self.__bo_trro[2*i+1].setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            self.__bo_trro[2*i+1].setTranslation(
                QtGui.QVector3D(*(xyzmid+delta/2)))

    def __shift_atoms(self, step: float):
        """Shift atoms by a given step.

        Shifts atoms by a given step along the displacement vectors
        stored internall and corrects the bonds accordingly.

        Parameters
        ----------
        step
            Value of the step.

        Notes
        -----
        * The bonds are not broken by the displacement, they are simply
          corrected.  To break it, it is necessary to provide new
          geometries and require the system to rebuild them.
        """
        if self.__at_displ is None:
            return
        self.__at_displ_step = step
        new_crds = []
        for atvec, atxyz, displ in zip(self.__at_tvec, self.__atcrd,
                                       self.__at_displ):
            new_crds.append(atxyz+step*displ)
            atvec.setTranslation(QtGui.QVector3D(*new_crds[-1]))

        for i, bond in enumerate(self.__bonds):
            iat1, iat2 = bond
            xyzat1 = new_crds[iat1]
            xyzat2 = new_crds[iat2]
            xyzmid = (xyzat1+xyzat2)/2.
            # First half of bond
            delta = xyzmid - xyzat1
            bo_len = sqrt(np.dot(delta, delta))
            self.__bo_mesh[2*i].setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            self.__bo_trro[2*i].setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            self.__bo_trro[2*i].setTranslation(
                QtGui.QVector3D(*(xyzat1+delta/2)))
            # Second half of bond
            delta = xyzat2 - xyzmid
            bo_len = sqrt(np.dot(delta, delta))
            self.__bo_mesh[2*i+1].setLength(bo_len)
            bo_rot = vrotate_3D(np.array([0, 1, 0]), delta/bo_len)
            rot3x3 = QtGui.QMatrix3x3(np.array(bo_rot).reshape(9).tolist())
            self.__bo_trro[2*i+1].setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            self.__bo_trro[2*i+1].setTranslation(
                QtGui.QVector3D(*(xyzmid+delta/2)))

        self.step_changed.emit()

    def __get_shift(self):
        """Returns object"""
        return self.__at_displ_step

    step_changed = QtCore.Signal()
    shift_atoms = QtCore.Property(float, __get_shift, __shift_atoms,
                                  notify=step_changed)

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
            self.__cam.setViewCenter(abs_pos-loc_pos)
        elif clickEvent.button() == Qt3DRender.QPickEvent.LeftButton:
            self.click_molatom.emit(
                [self.__mol_id,
                 self.__at_obj.index(clickEvent.entity())+1])

    def __set_material(self,
                       color: tp.Sequence[int],
                       opacity: int = 255
                       ) -> Qt3DExtras.QDiffuseSpecularMaterial:
        """Set material for a given object.

        Returns the material parameters for a given color and opacity,
        given the default material set internally.

        Parameters
        ----------
        color
            color parameters as a RGB sequence.
        opacity
            Opacity (alpha).
        """
        if self.__optview['material'] == 'plastic':
            matter = Qt3DExtras.QDiffuseSpecularMaterial(self)
            matter.setShininess(100)
            matter.setAmbient(QtGui.QColor(*color))
        elif self.__optview['material'] == 'metal':
            matter = Qt3DExtras.QDiffuseSpecularMaterial(self)
            matter.setSpecular(QtGui.QColor(*[min(1.5*i, 255) for i in color]))
            matter.setDiffuse(QtGui.QColor(100, 100, 100))
            matter.setAmbient(QtGui.QColor(*color))
            matter.setShininess(2000)
        elif self.__optview['material'] == 'glass':
            matter = Qt3DExtras.QDiffuseSpecularMaterial(self)
            matter.setDiffuse(QtGui.QColor(255, 255, 255, 150))
            matter.setAlphaBlendingEnabled(True)
            matter.setShininess(300)
            matter.setAmbient(QtGui.QColor(*color, opacity))
        else:
            raise KeyError('Unspecified material')
        return matter

    @classmethod
    def reset_counter(cls):
        """Reset molecule instance counter."""
        cls.__mol_count = 0
