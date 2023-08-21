"""Module related to molecular visualization for ESTAMPES.

This module provides basic methods for the molecular visualization.

Attributes
----------
BONDDATA
    Bond data for visualization.
MOLCOLS
    List of molecular colors with a sufficiently good contrast.
RAD_VIS_SCL
    Default scaling factors of radius for visualization.
"""

from math import sqrt, inf, tan, radians
import typing as tp

import numpy as np

from PySide6 import QtCore, QtGui
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base import TypeAtCrd, TypeAtCrdM, TypeAtLab, TypeAtLabM, \
    TypeBonds, TypeBondsM, TypeColor, ArgumentError
from estampes.data.atom import atomic_data
from estampes.tools.math import vrotate_3D


# ================
# Module Constants
# ================

BONDDATA = {
    'rvis': 0.15,
    'rgb': (200, 200, 200)
}

MOLCOLS = [
    (0x1F, 0x77, 0xB4),
    (0xFF, 0x7F, 0x0C),
    (0x2C, 0xA0, 0x2C),
    (0xD6, 0x27, 0x28),
    (0x94, 0x67, 0xBD),
    (0x8C, 0x56, 0x4B),
    (0xE3, 0x77, 0xC2),
    (0x7F, 0x7F, 0x7F),
    (0xBC, 0xBD, 0x22),
    (0x17, 0xBE, 0xCF),
]

RAD_VIS_SCL = .6


# ==============
# Module Classes
# ==============

class Molecule(Qt3DCore.QEntity):
    """Molecule class for the visualization.

    Attributes
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
        super(Molecule, self).__init__(rootEntity)

        self.__class__.__mol_count += 1
        self.__mol_id = self.__class__.__mol_count
        self.set_display_settings(
            col_bond_as_atom=col_bond_as_atom,
            rad_atom_as_bond=rad_atom_as_bond,
            molcol=molcol)
        self.update_geom(at_lab, at_crd, bonds)

    def update_geom(self,
                    at_lab: TypeAtLab,
                    at_crd: TypeAtCrd,
                    bonds: TypeBonds,
                    render: bool = True):
        """Updates geometry information.

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
        """Sets display settings for the molecule.

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
            self.__bo_mat.setAmbient(QtGui.QColor(BONDDATA['rgb']))
        else:
            self.__bo_mat = None

        self.__bo_rad = BONDDATA['rvis']*RAD_VIS_SCL

    def update_render(self):
        """Updates rendering of the molecules.

        Updates the rendering of the molecules with the current display
          settings.
        """
        self.__build_bonds()
        self.__build_atoms()

    def addMouse(self, cam: Qt3DRender.QCamera):
        """Adds mouse support.

        Adds mouse support (clicks) on the molecule object.

        Parameters
        ----------
        cam
            A Qt3Render.QCamera method.
        """
        self._cam = cam
        for at in self.__at_pick:
            at.pressed.connect(self.__clickAtom)

    def __upd_atdat(self):
        """Updates internal atomic data information.

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
        """Builds bonds objects and associated properties.

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
        """Builds atoms objects and associated properties.

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
        if clickEvent.button() == QtCore.Qt.RightButton:
            # abs_pos: absolute position of the clicked point
            abs_pos = clickEvent.worldIntersection()
            # loc_pos: local position of the clicked point in the object
            loc_pos = clickEvent.localIntersection()
            # Subtracting them give us the origin of the sphere
            self._cam().setViewCenter(abs_pos-loc_pos)
        elif clickEvent.button() == QtCore.Qt.LeftButton:
            self.click_molatom.emit(
                [self.__mol_id,
                 self.__at_obj.index(clickEvent.entity())+1])

    @classmethod
    def reset_counter(cls):
        Molecule.__mol_count = 0


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
    click_mol = QtCore.Signal(list)

    def __init__(self, nmols: int,
                 atlabs: TypeAtLabM,
                 atcrds: TypeAtCrdM,
                 bonds: TypeBondsM,
                 col_bond_as_atom: bool = False,
                 rad_atom_as_bond: bool = False,
                 molcols: tp.Optional[TypeColor] = None):
        super(MolWin, self).__init__()

        # Camera
        self.camera().lens().setPerspectiveProjection(45, 16 / 9, 0.1, 1000)
        self.camera().setPosition(QtGui.QVector3D(0, 1, 40))
        self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))

        # For camera controls
        self.rootEntity = Qt3DCore.QEntity()
        Molecule.reset_counter()
        if nmols == 1:
            self.mol = Molecule(atlabs, atcrds, bonds, col_bond_as_atom,
                                rad_atom_as_bond, molcols, self.rootEntity)
            self.mol.addMouse(self.camera)
            self.mol.click_molatom.connect(self.atom_clicked)
        else:
            self.mols = []
            for i in range(nmols):
                if molcols is None:
                    molcol = None
                else:
                    molcol = molcols[i]
                self.mols.append(Molecule(atlabs[i], atcrds[i], bonds[i],
                                          col_bond_as_atom, rad_atom_as_bond,
                                          molcol, self.rootEntity))
                self.mols[-1].addMouse(self.camera)
                self.mols[-1].click_molatom.connect(self.atom_clicked)
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

    @QtCore.Slot(list)
    def atom_clicked(self, molatom: tp.Sequence[int]):
        """Captures the information on the clicked atom.

        Captures the information from the clicked atom on one of the
        molecules.

        Parameters
        ----------
        molatom
            Molar atom, as a list of integers.
        """
        self.click_mol.emit(molatom)

    def mousePressEvent(self, mouseEvent):
        if mouseEvent.button() == QtCore.Qt.RightButton:
            self.camera().setViewCenter(QtGui.QVector3D(0, 0, 0))
        # print(mouseEvent.x())
        # self.camera().setViewCenter(QtGui.QVector3D(mouseEvent.x(),
        #                                             mouseEvent.y(), 0))
        super(MolWin, self).mousePressEvent(mouseEvent)


# ==============
# Module Methods
# ==============

def build_box(at_lab: TypeAtLab,
              at_crd: TypeAtCrd,
              rad_atom_as_bond: bool = False) -> tp.Dict[str, float]:
    """Builds the outer box containing the molecule.

    Builds the outer box containing the whole molecule and returns its
    bounds:

    * xmin: lower X value
    * xmax: upper X value
    * ymin: lower Y value
    * ymax: upper Y value
    * zmin: lower Z value
    * zmax: upper Z value

    Arguments
    ---------
    at_lab
        Atomic labels
    at_crd
        3-tuples with atomic coordinates, in Ang.
    rad_atom_as_bond
        If true, atomic radii are set equal to the bonds (tubes).

    Returns
    -------
    dict
        Dictionary with positions of `xmin`, `xmax`, `ymin`, `ymax`,
        `zmin`, `zmax`.
    """
    xmin, ymin, zmin, xmax, ymax, zmax = +inf, +inf, +inf, -inf, -inf, -inf
    atrad = rad_atom_as_bond and BONDDATA['rvis']*RAD_VIS_SCL or None
    if atrad is None:
        atdat = atomic_data(*sorted(set(at_lab)))
    for at, xyz in zip(at_lab, at_crd):
        x, y, z = xyz
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        if y > ymax:
            ymax = y
        if z < zmin:
            zmin = z
        if z > zmax:
            zmax = z
    if atrad is None:
        _rad = atdat[at]
    else:
        _rad = atrad
    xmin = xmin - _rad
    xmax = xmax + _rad
    ymin = ymin - _rad
    ymax = ymax + _rad
    zmin = zmin - _rad
    zmax = zmax + _rad

    return {
        'xmin': xmin, 'ymin': ymin, 'zmin': zmin,
        'xmax': xmax, 'ymax': ymax, 'zmax': zmax
    }


def set_cam_zpos(box_mol: tp.Dict[str, float],
                 xangle: tp.Union[float, int] = 68,
                 yangle: tp.Union[float, int] = 53) -> float:
    """Sets Z position of camera in POV-Ray.

    Sets the Z position of the camera in POV-Ray, assuming that the axes
    are the following:

    .. code-block:: text

        y ^   ^ z
          |  /
          | /
          -------------> x

    References
    ----------
    http://www.povray.org/documentation/view/3.7.0/246/

    Parameters
    ----------
    box_mol
        Outer box containing the molecule.
    xangle
        Full angle along x of the camera (in degrees).
    yangle
        Full angle along y of the camera (in degrees).

    Returns
    -------
    float
        Position of the camera along -Z.
    """
    Zx = max(abs(box_mol['xmin'])/tan(radians(xangle/2.)),
             abs(box_mol['xmax'])/tan(radians(xangle/2.)))
    Zy = max(abs(box_mol['ymin'])/tan(radians(yangle/2.)),
             abs(box_mol['ymax'])/tan(radians(yangle/2.)))
    rval = box_mol['zmin'] - max(Zx, Zy)
    scale = 100.
    # round to 2 decimal digits (slightly overestimating the value)
    return (int(rval*scale)-1)/scale


def write_pov(fname: str,
              nmols: int,
              at_lab: TypeAtLabM,
              at_crd: TypeAtCrdM,
              bonds: TypeBondsM,
              zcam: float = -8.0,
              col_bond_as_atom: bool = False,
              rad_atom_as_bond: bool = False):
    """Writes a Pov-Ray file.

    Builds and writes a Pov-Ray file.
    If `nmols` > 1, `at_lab`, `at_crd`, `bonds` are lists, with each
    item corresponding to a molecule.

    Parameters
    ----------
    nmols
        Number of molecules stored in `at_lab`, `at_crd` and `bonds`.
    at_lab
        Atomic labels.
    at_crd
        3-tuples with atomic coordinates, in Ang.
    bonds
        2-tuples listing connected atoms.
    zcam
        Position of the camera along -Z.
    col_bond_as_atom
        If true, bonds are colored based on the connected atoms
    rad_atom_as_bond
        If true, atomic radii are set equal to the bonds (tubes).
    """
    bo_rad = BONDDATA['rvis']*RAD_VIS_SCL
    if rad_atom_as_bond:
        rad_at = 1.
        rad_bo = 1.
    else:
        rad_at = 1.
        rad_bo = 1.
    if nmols == 1:
        list_atoms = sorted(set(at_lab))
    else:
        list_atoms = sorted(set([item for mol in at_lab for item in mol]))
    atdat = atomic_data(*list_atoms)
    header = '''\
#declare dist = 1.;
#declare scl_bond = {scl_bo:.2f};
#declare scl_rat = {scl_at:.2f};
#declare Trans = 0.0;

global_settings {{
    ambient_light rgb <0.200, 0.200, 0.200>
    max_trace_level 15
}}

camera {{
    location <   0.00000,   0.00100, {zcam:10.5f}>*dist
    look_at  <   0.00000,   0.00000,    0.00000>
}}

light_source {{
    <10.000, 10.000, -20.000>*dist
    color rgb <1.000, 1.000, 1.000>
    fade_distance 25.000
    fade_power 0
    parallel
    point_at <0.000, 0.000, 0.000>
}}

light_source {{
    <-10.000, -15.000, -15.000>*dist
    color rgb <1.000, 1.000, 1.000>*.3
    fade_distance 30
    fade_power 0
    parallel
    point_at <-2.000, -2.000, 0.000>
}}

#declare Aspect = finish {{
    ambient .8
    diffuse 1
    specular 1
    roughness .005
    metallic 0.5
}}

'''
    fmt_col_at = '''\
#declare col_at_{at:2s} = pigment {{\
 rgbt < {c[0]:3d}, {c[1]:3d}, {c[2]:3d}, Trans >/255. \
}}
'''
    fmt_tex_at = '''\
#declare tex_at_{at:2s} = texture {{
    pigment {{ col_at_{at:2s} }}
    finish {{ Aspect }}
}}
'''
    fmt_tex_bo = '''\
#declare tex_bo_{at:2s} = texture {{
    pigment {{ col_{ref} }}
    finish {{ Aspect }}
}}
'''
    fmt_col_mol = '''\
#declare col_mol{id:02d} = pigment {{\
 rgbt < {c[0]:3d}, {c[1]:3d}, {c[2]:3d}, Trans >/255. \
}}
'''
    fmt_tex_mol = '''\
#declare tex_mol{id:02d} = texture {{
    pigment {{ col_mol{id:02d} }}
    finish {{ Aspect }}
}}
'''
    fmt_rad_at = '''\
#declare r_at_{at:2s} = {r:.6f}*scl_rat;
'''
    if nmols == 1:
        fmt_obj_bo = '''\
    cylinder {{ // Bond {at1}({id1})- {at2}({id2})
        <{xyz1[0]:14.6f},{xyz1[1]:14.6f},{xyz1[2]:14.6f}>,
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        r_bond
        texture {{ tex_bo_{at1} }}
    }}
    cylinder {{ // Bond {at1}({id1}) -{at2}({id2})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        <{xyz2[0]:14.6f},{xyz2[1]:14.6f},{xyz2[2]:14.6f}>,
        r_bond
        texture {{ tex_bo_{at2} }}
    }}
'''
        fmt_obj_at = '''\
    sphere {{ // Atom {at}({id})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>, r_at_{at}
        texture {{ tex_at_{at} }}
    }}
'''
    else:
        fmt_obj_bo = '''\
    cylinder {{ // Bond {at1}({id1})- {at2}({id2})
        <{xyz1[0]:14.6f},{xyz1[1]:14.6f},{xyz1[2]:14.6f}>,
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        r_bond
        texture {{ tex_mol{idmol:02d} }}
    }}
    cylinder {{ // Bond {at1}({id1}) -{at2}({id2})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        <{xyz2[0]:14.6f},{xyz2[1]:14.6f},{xyz2[2]:14.6f}>,
        r_bond
        texture {{ tex_mol{idmol:02d} }}
    }}
'''
        fmt_obj_at = '''\
    sphere {{ // Atom {at}({id})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>, r_at_{at}
        texture {{ tex_mol{idmol:02d} }}
    }}
'''

    with open(fname, 'w') as fobj:
        fobj.write(header.format(scl_bo=rad_bo, scl_at=rad_at, zcam=zcam))
        # Set atoms colors
        if nmols == 1:
            fobj.write('\n// ATOMS COLORS AND TEXTURES\n')
            for atom in list_atoms:
                fobj.write(fmt_col_at.format(at=atom, c=atdat[atom]['rgb']))
            fobj.write('\n')
            # Set atoms texture
            for atom in list_atoms:
                fobj.write(fmt_tex_at.format(at=atom))
            # Set bonds colors/textures
            fobj.write('\n// BONDS COLORS AND TEXTURES\n')
            if not col_bond_as_atom:
                fmt = '''\
#declare col_bond = pigment {{
    rgbt < {c[0]:3d}, {c[1]:3d}, {c[2]:3d}, Trans >/255.
}}
'''
                fobj.write(fmt.format(c=BONDDATA['rgb']))
                fmt = 'bond'
            else:
                fmt = 'at_{at}'
            for atom in list_atoms:
                fobj.write(fmt_tex_bo.format(at=atom, ref=fmt.format(at=atom)))
        else:
            fobj.write('\n// MOLECULES COLORS AND TEXTURES\n')
            # Set Molecule colors
            for i in range(nmols):
                fobj.write(fmt_col_mol.format(id=i+1, c=MOLCOLS[i]))
            # Set Molecule textures
            for i in range(nmols):
                fobj.write(fmt_tex_mol.format(id=i+1))
        # Defines bonds/atoms radi
        fobj.write('\n// ATOMS AND BONDS RADII\n')
        fmt = '#declare r_bond  = {:.6f}*scl_bond;\n'
        fobj.write(fmt.format(bo_rad))
        for atom in list_atoms:
            if rad_atom_as_bond:
                rval = bo_rad
            else:
                rval = atdat[atom]['rvis']*RAD_VIS_SCL
            fobj.write(fmt_rad_at.format(at=atom, r=rval))
        # Molecule specification
        if nmols == 1:
            fobj.write('\n// MOLECULE DEFINITION\n\n')
            fobj.write('#declare molecule = union {\n')
            # -- First build bonds
            for bond in bonds:
                iat1, iat2 = bond
                xyz1 = at_crd[iat1]
                xyz2 = at_crd[iat2]
                xyzmid = (xyz1 + xyz2)/2.
                fobj.write(fmt_obj_bo.format(
                    xyz=xyzmid,
                    at1=at_lab[iat1], id1=iat1+1, xyz1=xyz1,
                    at2=at_lab[iat2], id2=iat2+1, xyz2=xyz2))
            # -- Next build atoms
            for iat in range(len(at_lab)):
                fobj.write(fmt_obj_at.format(
                    at=at_lab[iat], id=iat+1, xyz=at_crd[iat]))
            # -- Close and add molecules
            fobj.write('''}

object {
    molecule
}
''')

        else:
            fobj.write('\n// MOLECULES DEFINITION\n')
            for imol in range(nmols):
                fobj.write('\n#declare mol{:02d} = union {{\n'.format(imol+1))
                # -- First build bonds
                for bond in bonds[imol]:
                    iat1, iat2 = bond
                    xyz1 = at_crd[imol][iat1]
                    xyz2 = at_crd[imol][iat2]
                    xyzmid = (xyz1 + xyz2)/2.
                    fobj.write(fmt_obj_bo.format(
                        xyz=xyzmid,
                        at1=at_lab[imol][iat1], id1=iat1+1, xyz1=xyz1,
                        at2=at_lab[imol][iat2], id2=iat2+1, xyz2=xyz2,
                        idmol=imol+1))
                # -- Next build atoms
                for iat in range(len(at_lab[imol])):
                    fobj.write(fmt_obj_at.format(
                        at=at_lab[imol][iat], id=iat+1,
                        xyz=at_crd[imol][iat], idmol=imol+1))
                # -- Close
                fobj.write('}\n')
            # Add block to visualize molecules
            fobj.write('union {\n')
            for i in range(nmols):
                fobj.write('    object {{ mol{:02d} }}\n'.format(i+1))
            fobj.write('}\n')
