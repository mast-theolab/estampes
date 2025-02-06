"""Module related to mode/vibration visualization for ESTAMPES.

This modules provides low-level elements to represent modes.
"""

import os
import typing as tp

import numpy as np

from PySide6 import QtCore, QtGui
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base import Type1Vib, TypeAtCrd, TypeColor
from estampes.base.errors import ArgumentError
from estampes.data.colors import to_rgb_list
from estampes.data.visual import PATH_OBJ3D, MODELS, VIBCOLS
from estampes.tools.math import vrotate_3D


class VibMode(Qt3DCore.QEntity):
    """The VibMode class for  visualization.

    Builds a representation of a molecular vibration.

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

    def __init__(self,
                 at_crd: TypeAtCrd,
                 vib_mode: Type1Vib,
                 vizmode: tp.Optional[str] = None,
                 color: tp.Optional[tp.Any] = None,
                 color2: tp.Optional[tp.Any] = None,
                 rootEntity: tp.Optional[Qt3DCore.QEntity] = None):
        """Initialize a VibMode instance.

        Initializes an instance of the VibMode class.

        Parameters
        ----------
        at_crd
            List of atomic centers, as origins for the vibrations, in Ang.
        vib_mode
            List of atomic displacements for a given mode.
        vizmode
            Visualization mode.
            arrow: displacements are represented as arrow.
        color
            Color of the object representing the displacement.
        color2
            Secondary color of the object representing the displacement.
        rootEntity
            Qt root entity to connect the new `VibMode` QEntity.
        """
        super().__init__(rootEntity)

        self.__atcrd = None
        self.__vmode = None
        # The initialization below should be properly handled through
        # functions.  Just making sure everything declared.
        self.__vib_objs = self.__vib_mesh = self.__vib_trro = None

        self.set_display_settings(vizmode=vizmode, color=color, color2=color2)
        self.update_vib(vib_mode, at_crd)

    def update_vib(self,
                   vib_mode: Type1Vib,
                   at_crd: tp.Optional[TypeAtCrd] = None,
                   render: bool = True):
        """Update geometry information.

        Updates atomic labels and coordinates.

        Parameters
        ----------
        vib_mode
            List of atomic displacements for a given mode.
        at_crd
            List of atomic centers, as origins for the vibrations, in Ang.
        render
            If True, the vibration is re-rendered.
        """
        self.__vmode = vib_mode
        if at_crd is not None:
            self.__atcrd = at_crd
        elif self.__atcrd is None:
            raise ArgumentError('Missing coordinates of atomic centers.')
        if render:
            self.update_render()

    def set_display_settings(self, *,
                             vizmode: tp.Optional[str] = None,
                             color: tp.Optional[tp.Any] = None,
                             scale_mode: float = 1.0,
                             color2: tp.Optional[tp.Any] = None):
        """Set display settings for the vibration.

        Sets color information and rendering.

        Parameters
        ----------
        vizmode
            Visualization mode.
            See [[update_renderer]] for details.
        color
            Primary color of the vibration.
        scale_mode
            Scaling factor to be applied to the displacements.
        color2
            Secondary color of the vibration (typically: negative sign).

        The internal databases are updated.

        Raises
        ------
        ArgumentError
            Incorrect argument types.
        """
        self.__material = Qt3DExtras.QDiffuseSpecularMaterial

        self.__optview = {}
        if vizmode is None:
            _mode = 'arrows'
        for key, info in MODELS['vib'].items():
            if vizmode.lower() in info['alias']:
                _mode = key
                break
        else:
            raise ArgumentError('Unrecognized visualization mode')
        self.__optview['mode'] = _mode

        _rgb0 = None
        _rgb1 = None
        if color is None:
            if _mode in ('arrows', 'midarrows'):
                # _rgb0 = (54, 117, 188)
                _rgb0 = VIBCOLS['arrow']
            elif _mode == 'dualarrows':
                # _rgb0 = (28, 214, 46)
                _rgb0 = VIBCOLS['arrow+']
        else:
            try:
                _rgb0 = to_rgb_list(color)
            except (ValueError, ArgumentError) as err:
                raise ArgumentError('vibcol', 'Incorrect color specification'
                                    ) from err
        if _mode in ('dualarrows', ):
            if color2 is None:
                # _rgb1 = (242, 46, 46)
                _rgb1 = VIBCOLS['arrow-']
            else:
                try:
                    _rgb1 = to_rgb_list(color2)
                except (ValueError, ArgumentError) as err:
                    raise ArgumentError('vibcol',
                                        'Incorrect color specification'
                                        ) from err

        self.__optview['col'] = _rgb0
        self.__vibmat = self.__material(self)
        self.__vibmat.setAmbient(QtGui.QColor(*_rgb0))
        if _rgb1 is not None:
            self.__optview['col1'] = _rgb1
            self.__vibmat1 = self.__material(self)
            self.__vibmat1.setAmbient(QtGui.QColor(*_rgb1))

    def update_render(self):
        """Update rendering of the vibration.

        Updates the rendering of the vibration with the current display
          settings.
        """
        if self.__optview['mode'] == 'arrows':
            self.__build_arrows()
        elif self.__optview['mode'] == 'midarrows':
            self.__build_centered_arrows()
        elif self.__optview['mode'] == 'dualarrows':
            self.__build_dual_arrows()
        elif self.__optview['mode'] == 'spheres':
            self.__build_spheres()

    def __build_arrows(self):
        """Build list of arrows showing the vibration.

        Builds a list of arrows representing the atomic displacements
        for a given vibration.
        """
        self.__vib_objs = []
        self.__vib_trro = []
        self.__vib_mesh = []
        shaft_r = 0.03
        shaft_l = 3
        head_r = 0.06
        head_l = 0.1
        for xyz, dxyz in zip(self.__atcrd, self.__vmode):
            # Initialization
            shaft_obj = Qt3DCore.QEntity(self)
            shaft_mesh = Qt3DExtras.QCylinderMesh()
            shaft_trro = Qt3DCore.QTransform()
            head_obj = Qt3DCore.QEntity(self)
            head_mesh = Qt3DExtras.QConeMesh()
            head_trro = Qt3DCore.QTransform()
            # Operations
            ldxyz = np.linalg.norm(dxyz)
            shaft_mesh.setRadius(shaft_r)
            shaft_mesh.setLength(ldxyz*shaft_l)
            rot = vrotate_3D(np.array([0, 1, 0]), dxyz/ldxyz)
            rot3x3 = QtGui.QMatrix3x3(np.array(rot).reshape(9).tolist())
            shaft_trro.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            shaft_trro.setTranslation(QtGui.QVector3D(*(xyz + dxyz*shaft_l/2)))
            shaft_obj.addComponent(shaft_mesh)
            shaft_obj.addComponent(shaft_trro)
            shaft_obj.addComponent(self.__vibmat)
            head_mesh.setBottomRadius(head_r)
            head_mesh.setLength(head_l)
            head_trro.setRotation(QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            head_trro.setTranslation(QtGui.QVector3D(*(xyz + dxyz*shaft_l)))
            head_obj.addComponent(head_mesh)
            head_obj.addComponent(head_trro)
            head_obj.addComponent(self.__vibmat)
            # # Update DB
            self.__vib_objs.append((shaft_obj, head_obj))
            self.__vib_mesh.append((shaft_mesh, head_mesh))
            self.__vib_trro.append((shaft_trro, head_trro))

    def __build_centered_arrows(self):
        """Build list of centered arrows showing the vibration.

        Builds a list of centered arrows representing the atomic
        displacements for a given vibration.
        In this representation, the arrow passes through the atom 
        """
        self.__vib_objs = []
        self.__vib_trro = []
        self.__vib_mesh = []
        shaft_r = 0.03
        shaft_l = 4
        head_r = 0.06
        head_l = 0.1
        for xyz, dxyz in zip(self.__atcrd, self.__vmode):
            # Initialization
            shaft_obj = Qt3DCore.QEntity(self)
            shaft_mesh = Qt3DExtras.QCylinderMesh()
            shaft_trro = Qt3DCore.QTransform()
            head_obj = Qt3DCore.QEntity(self)
            head_mesh = Qt3DExtras.QConeMesh()
            head_trro = Qt3DCore.QTransform()
            # Operations
            ldxyz = np.linalg.norm(dxyz)
            shaft_mesh.setRadius(shaft_r)
            shaft_mesh.setLength(ldxyz*shaft_l)
            rot = vrotate_3D(np.array([0, 1, 0]), dxyz/ldxyz)
            rot3x3 = QtGui.QMatrix3x3(np.array(rot).reshape(9).tolist())
            shaft_trro.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            shaft_trro.setTranslation(QtGui.QVector3D(*xyz))
            shaft_obj.addComponent(shaft_mesh)
            shaft_obj.addComponent(shaft_trro)
            shaft_obj.addComponent(self.__vibmat)
            head_mesh.setBottomRadius(head_r)
            head_mesh.setLength(head_l)
            head_trro.setRotation(QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            head_trro.setTranslation(QtGui.QVector3D(*(xyz + dxyz*shaft_l/2)))
            head_obj.addComponent(head_mesh)
            head_obj.addComponent(head_trro)
            head_obj.addComponent(self.__vibmat)
            # # Update DB
            self.__vib_objs.append((shaft_obj, head_obj))
            self.__vib_mesh.append((shaft_mesh, head_mesh))
            self.__vib_trro.append((shaft_trro, head_trro))

    def __build_dual_arrows(self):
        """Build list of centered arrows showing the vibration.

        Builds a list of centered arrows representing the atomic
        displacements for a given vibration.
        In this representation, the arrow passes through the atom 
        """
        self.__vib_objs = []
        self.__vib_trro = []
        self.__vib_mesh = []
        shaft_r = 0.03
        shaft_l = 2
        head_r = 0.06
        head_l = 0.1
        for xyz, dxyz in zip(self.__atcrd, self.__vmode):
            # Initialization
            shaft_obj1 = Qt3DCore.QEntity(self)
            shaft_obj2 = Qt3DCore.QEntity(self)
            shaft_mesh = Qt3DExtras.QCylinderMesh()
            shaft_trro1 = Qt3DCore.QTransform()
            shaft_trro2 = Qt3DCore.QTransform()
            head_obj1 = Qt3DCore.QEntity(self)
            head_obj2 = Qt3DCore.QEntity(self)
            head_mesh = Qt3DExtras.QConeMesh()
            head_trro1 = Qt3DCore.QTransform()
            head_trro2 = Qt3DCore.QTransform()
            # Operations
            ldxyz = np.linalg.norm(dxyz)
            shaft_mesh.setRadius(shaft_r)
            shaft_mesh.setLength(ldxyz*shaft_l)
            rot1 = vrotate_3D(np.array([0, 1, 0]), dxyz/ldxyz)
            rot2 = vrotate_3D(np.array([0, -1, 0]), dxyz/ldxyz)
            rot3x3 = QtGui.QMatrix3x3(np.array(rot1).reshape(9).tolist())
            shaft_trro1.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            shaft_trro1.setTranslation(
                QtGui.QVector3D(*(xyz + dxyz*shaft_l/2)))
            shaft_obj1.addComponent(shaft_mesh)
            shaft_obj1.addComponent(shaft_trro1)
            shaft_obj1.addComponent(self.__vibmat)
            shaft_trro2.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            shaft_trro2.setTranslation(
                QtGui.QVector3D(*(xyz - dxyz*shaft_l/2)))
            shaft_obj2.addComponent(shaft_mesh)
            shaft_obj2.addComponent(shaft_trro2)
            shaft_obj2.addComponent(self.__vibmat1)
            head_mesh.setBottomRadius(head_r)
            head_mesh.setLength(head_l)
            head_trro1.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            head_trro1.setTranslation(QtGui.QVector3D(*(xyz + dxyz*shaft_l)))
            head_obj1.addComponent(head_mesh)
            head_obj1.addComponent(head_trro1)
            head_obj1.addComponent(self.__vibmat)
            head_trro2.setRotation(QtGui.QQuaternion.fromRotationMatrix(
                QtGui.QMatrix3x3(np.array(rot2).reshape(9).tolist())
            ))
            head_trro2.setTranslation(QtGui.QVector3D(*(xyz - dxyz*shaft_l)))
            head_obj2.addComponent(head_mesh)
            head_obj2.addComponent(head_trro2)
            head_obj2.addComponent(self.__vibmat1)
            # # Update DB
            self.__vib_objs.append((shaft_obj1, head_obj1, shaft_obj2,
                                    head_obj2))
            self.__vib_mesh.append((shaft_mesh, head_mesh))
            self.__vib_trro.append((shaft_trro1, head_trro1, shaft_trro2,
                                    head_trro2))

    def __build_spheres(self, scale_vib: float = 1.0):
        """Build atom objects and associated properties.

        Builds atoms objects and data in the molecule.
        """
        self.__vib_objs = []
        self.__vib_trro = []
        self.__vib_mesh = []
        path_obj = QtCore.QUrl.fromLocalFile(
            os.path.join(PATH_OBJ3D, 'vib2cols.obj'))
        for xyz, dxyz in zip(self.__atcrd, self.__vmode):
            # Initialization
            obj = Qt3DCore.QEntity(self)
            mesh = Qt3DRender.QSceneLoader()
            trro = Qt3DCore.QTransform()
            # Operations
            ldxyz = np.linalg.norm(dxyz)
            mesh.setSource(path_obj)
            rot = vrotate_3D(np.array([0, 1, 0]), dxyz/ldxyz)
            rot3x3 = QtGui.QMatrix3x3(np.array(rot).reshape(9).tolist())
            trro.setScale(ldxyz*.5)
            trro.setRotation(
                QtGui.QQuaternion.fromRotationMatrix(rot3x3))
            trro.setTranslation(QtGui.QVector3D(*xyz))
            obj.addComponent(mesh)
            obj.addComponent(trro)
            # Update DB
            self.__vib_objs.append(obj)
            self.__vib_mesh.append(mesh)
            self.__vib_trro.append(trro)
