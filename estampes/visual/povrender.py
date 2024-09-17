"""Persistence of Vision Raytracing: Rendering tools

This module provides functions and instruments to build input files for
POV-Ray (http://www.povray.org/).
"""

import os
import re
import typing as tp
from math import inf, tan, radians

import numpy as np

from estampes.base import QLabel
from estampes.base.types import Type1Vib, TypeAtCrd, TypeAtCrdM, TypeAtLab, \
    TypeAtLabM, TypeBondsM, TypeVibsM
from estampes.base.errors import ArgumentError, QuantityError
from estampes.parser import DataFile
from estampes.data.atom import atomic_data
from estampes.data.physics import PHYSFACT
from estampes.data.visual import BONDDATA, MOLCOLS, RAD_VIS_SCL
from estampes.tools.atom import convert_labsymb
from estampes.tools.math import vrotate_3D
from estampes.tools.mol import list_bonds


OBJ_2COL_SPHERE = '''
#declare vib_sphere = union {{
    difference {{
        sphere {{ <0, 0, 0>, .5 }}
        box {{ <-1, 0, -1>, <1, -1, 1> }}
        {mat}(255, 188, 0)
    }}
    difference {{
        sphere {{ <0, 0, 0>, .5 }}
        box {{ <-1, 0, -1>, <1, 1, 1> }}
        {mat}(0, 125, 255)
    }}
}}
'''

OBJ_ARROW = '''
#macro Build_Arrow(l_shaft)
    union {{
        cylinder {{<0, 0, 0>, <0, l_shaft, 0>, .025}}
        cone {{
            <0, l_shaft, 0>, 0.06
            <0, l_shaft+.08, 0>, 0
        }}
        {mat}(0, 125, 255)
    }}
#end
'''

MAT_PLASTIC = '''
#macro Mater_Plastic(r, g, b)
    material {
        texture {
            finish {
                ambient .8
                diffuse 1
                specular 1
                roughness .005
                metallic 0.5
            }
            pigment { rgb <r, g, b>/255. }
        }
    }
#end
'''

MAT_GLASS = '''
#macro Mater_Glass(r, g, b)
    material {
        texture {
            finish {
                specular 0.6
                roughness .001
                ambient 0.2
                diffuse 0.2
                emission 0.05
                reflection {
                    0.2, 1.0
                    fresnel on
                }
                conserve_energy
            }
            pigment {
                rgb <r, g, b>/255.
                filter 0.9
            }
        }
        interior {
            ior 1.5
            caustics 1.0
            fade_distance 1
            fade_power 2
        }
    }
#end
'''

MAT_METAL = '''
#macro Mater_Metal(r, g, b)
    material {
        texture {
            finish {
                metallic
                ambient 0.1
                diffuse 0.75
                brilliance 6.0
                reflection 0.3
                phong 0.8
                phong_size 120
            }
            pigment {
                rgb <r, g, b>/255.
            }
        }
    }
#end
'''

CHOICES_MATER = {
    'glass': 'Mater_Glass',
    'plastic': 'Mater_Plastic',
    'metal': 'Mater_Metal'
}


class POVBuilder():
    """POV-Ray input builder.

    Builds and manages a POV-ray input file.
    """
    def __init__(self,
                 infiles: tp.Optional[tp.Union[tp.List[str], str]] = None,
                 *,
                 load_data: bool = True,
                 load_vibs: bool = False,
                 povfile: tp.Optional[str] = None,
                 at_crds: tp.Optional[TypeAtCrdM] = None,
                 at_labs: tp.Optional[TypeAtLabM] = None,
                 v_modes: tp.Optional[TypeVibsM] = None,
                 nmols: int = 0,
                 in_au: bool = True,
                 **params):
        """Initialize instance with core informations.

        Initializes internal data for the POV-Ray file builder.
        The constructor can be called in three modes:

        1. no input, the data will be provided later.
        2. providing data files, that will be parsed internally.
        3. providing directly the data.

        NOTE: The modes are incompatible and would raise an error.

        The latter can be useful to limit memory usage and avoid
        duplicating information.

        Parameters
        ----------
        infiles
            List of input files, containing data to load.
        load_data
            Load data during the instantiation (only for mode 2).
        load_vibs
            Load vibrational data (only for mode 2).
            NOTE: If vibrations are absent, this does not raise an error.
        povfile
            Output file to store the POV-Ray scene description.
            The file can be given as template, supporting the following
            Python replacement fields:

            - `{mol}`: replaced by molecule number.
            - `{vib}`: replaced by vibration number

            Format specifications can be indicated (ex: `{vib:02d}`).
        at_crds
            List containing the atomic coordinates.
        at_labs
            List of atomic labels.
        v_modes
            List of vibrational modes.
        nmols
            Number of molecules (only for mode 3).
            `nmols` must be set for mode 3.
        in_au
            Atomic coordinates are in atomic units.
        params
            Additional parameters, transmitted to:

            - load_data
        """
        if infiles is not None and nmols > 0:
            raise ArgumentError(
                'POVBuilder should be initialized with a single mode.')

        self.__load_vibs = load_vibs
        self.__nmols = 0
        self.__atcrd = []
        self.__atlab = []
        self.__bonds = []
        self.__vmodes = []
        self.__infiles = []
        self.__loaded = []
        if infiles is not None:
            if isinstance(infiles, (tuple, list)):
                if load_data:
                    for infile in infiles:
                        self.load_data(infile=infile, **params)
                else:
                    self.__nmols = len(infiles)
                    self.__infiles = infiles[:]
                    self.__atcrd = [None for _ in range(self.__nmols)]
                    self.__atlab = [None for _ in range(self.__nmols)]
                    self.__vmodes = [None for _ in range(self.__nmols)]
                    self.__loaded = [False for _ in range(self.__nmols)]
            else:
                if load_data:
                    self.load_data(infile=infiles, **params)
                else:
                    self.__infiles = [infiles]
                    self.__atcrd = [None]
                    self.__atlab = [None]
                    self.__vmodes = [None]
                    self.__loaded = [False]
        else:
            if at_crds is None or at_labs is None:
                raise ArgumentError(
                    'Input files and molecule specifications both missing.\n'
                    + 'Nothing to do.')
            if nmols > 1:
                self.__nmols = nmols
                self.__atcrd = np.array(at_crds)
                if self.__atcrd.ndim != 3:
                    raise ArgumentError('Expected array of coordinates')
                if self.__atcrd.shape[0] < nmols:
                    raise ArgumentError(f'Expected {nmols} structures.')
                self.__atlab = np.array(at_labs)
                if self.__atlab.ndim == 1:
                    # We expand atlab to be the same for each quantity
                    self.__atlab = np.tile(self.__atlab, (nmols, 1))
                elif self.__atlab.ndim == 2:
                    if self.__atlab.shape[0] < nmols:
                        raise ArgumentError(f'Expected {nmols} atomic labels.')
                if v_modes is not None:
                    self.__vmodes = np.array(v_modes)
                    if self.__vmodes.ndim < 3:
                        raise ArgumentError('Expected array of vib. modes')
                    if self.__vmodes.ndim == 4:
                        raise ArgumentError(
                            'Mode specifications should be a 1D vector')
                    if self.__vmodes.ndim > 4:
                        raise ArgumentError(
                            'Too many dimensions for vib. mode specifications')
                    if self.__vmodes.shape[0] < nmols:
                        raise ArgumentError(
                            f'Expected {nmols} definitions of modes')
                else:
                    self.__vmodes = [None for _ in range(nmols)]
                self.__loaded = [True for _ in range(nmols)]
                if in_au:
                    self.__atcrd *= PHYSFACT.bohr2ang
            elif nmols in (0, 1):
                self.__nmols = 1
                self.__atcrd = [np.array(at_crds)]
                if self.__atcrd[0].ndim != 2:
                    raise ArgumentError('Expected coordinates as 2D array')
                self.__atlab = [np.array(at_labs)]
                if self.__atlab[0].ndim != 1:
                    raise ArgumentError('Expected atomic numbers as 1D array')
                if v_modes is not None:
                    self.__vmodes = [np.array(v_modes)]
                    if self.__vmodes[0].ndim != 2:
                        raise ArgumentError(
                            'Expected vibrational modes as 2D array')
                else:
                    self.__vmodes = [None]
                self.__loaded = [True]
                if in_au:
                    self.__atcrd[0] *= PHYSFACT.bohr2ang
            else:
                raise ArgumentError(
                    'Negative number of arguments not allowed.')
            _tol_bonds = params.get('tol_bonds', 1.2)
            self.__bonds = []
            for imol in range(nmols):
                self.__bonds.append(
                    list_bonds(self.__atlab[imol],
                               self.__atcrd[imol], _tol_bonds))

        if povfile is not None:
            self.__povfile = povfile
        else:
            self.__povfile = None

    def load_data(self, *,
                  infile: tp.Optional[tp.Union[str, int]] = None,
                  reset: bool = True,
                  force_reload: bool = False,
                  **params):
        """Load data from `infile` or files already in internal DB.

        Loads data from `infile` into internal databases or from files
        stored internally if no file is provided.  If the file is an
        integer, the data from the corresponding file in the internal
        DB are loaded.

        If a file is given and `reset` is True, the previous content is
        discarded.  If no file is given, `reset` is ignored.

        Parameters
        ----------
        infile
            Input file or index in internal DB.
        reset
            Reset internal database and fill it with content of file.
            Ignored if no file is provided.
        force_reload
            If file has already been loaded, its content is reloaded.
        params
            Additional parameters. Supported:

            - tol_bonds: tolerance pour bond definition.
        """
        dkeys = {
            'natoms': QLabel(quantity='natoms'),
            'atnum': QLabel(quantity='atnum'),
            'atcrd': QLabel(quantity='atcrd'),
        }
        _tol_bonds = params.get('tol_bonds', 1.2)
        extfile = isinstance(infile, str)
        if extfile:
            ifile = 0
            if reset:
                del self.__atcrd[:]
                del self.__atlab[:]
                del self.__vmodes[:]
                del self.__infiles[:]
                del self.__loaded[:]
                del self.__bonds[:]
                self.__nmols = 0
            if infile in self.__infiles:
                if force_reload:
                    ifile = self.__infiles.index(infile)
                else:
                    return
            dfile = DataFile(infile)
            qdata = dfile.get_data(**dkeys)
            if ifile == 0:
                self.__infiles.append(infile)
                self.__nmols += 1
                self.__loaded.append(True)
                self.__atcrd.append(
                    np.array(qdata['atcrd'].data)*PHYSFACT.bohr2ang)
                self.__atlab.append(
                    convert_labsymb(True, *qdata['atnum'].data))
                self.__bonds.append(
                    list_bonds(self.__atlab[-1], self.__atcrd[-1], _tol_bonds))
                if self.__load_vibs:
                    try:
                        self.__vmodes.append(
                            dfile.get_hess_data(True, False)[0])
                    except QuantityError:
                        self.__vmodes.append(None)
            else:
                self.__loaded[ifile] = True
                self.__atcrd[ifile] = \
                    np.array(qdata['atcrd'].data)*PHYSFACT.bohr2ang
                self.__atlab[ifile] = \
                    convert_labsymb(True, *qdata['atnum'].data)
                self.__bonds[ifile] = \
                    list_bonds(self.__atlab[ifile],
                               self.__atcrd[ifile], _tol_bonds)
                if self.__load_vibs:
                    try:
                        self.__vmodes[ifile] = \
                            dfile.get_hess_data(True, False)[0]
                    except QuantityError:
                        self.__vmodes[ifile] = None
        elif isinstance(infile, int):
            # Check if content already loaded:
            if self.__loaded[infile] and not force_reload:
                return
            dfile = DataFile(self.__infiles[infile])
            qdata = dfile.get_data(**dkeys)
            self.__loaded[infile] = True
            self.__atcrd[infile] = \
                np.array(qdata['atcrd'].data)*PHYSFACT.bohr2ang
            self.__atlab[infile] = \
                convert_labsymb(True, *qdata['atnum'].data)
            self.__bonds[infile] = \
                list_bonds(self.__atlab[infile],
                           self.__atcrd[infile], _tol_bonds)
            if self.__load_vibs:
                try:
                    self.__vmodes[infile] = \
                        dfile.get_hess_data(True, False)[0]
                except QuantityError:
                    self.__vmodes[infile] = None
        else:
            for ifile, fname in self.__infiles:
                # Check if content already loaded:
                if self.__loaded[ifile] and not force_reload:
                    continue
                dfile = DataFile(fname)
                qdata = dfile.get_data(**dkeys)
                self.__loaded[ifile] = True
                self.__atcrd[ifile] = \
                    np.array(qdata['atcrd'].data)*PHYSFACT.bohr2ang
                self.__atlab[ifile] = \
                    convert_labsymb(True, *qdata['atnum'].data)
                self.__bonds[ifile] = \
                    list_bonds(self.__atlab[ifile],
                               self.__atcrd[ifile], _tol_bonds)
                if self.__load_vibs:
                    try:
                        self.__vmodes[ifile] = \
                            dfile.get_hess_data(True, False)[0]
                    except QuantityError:
                        self.__vmodes[ifile] = None

    def write_pov(self, *, povname: tp.Optional[str] = None,
                  id_mol: tp.Optional[int] = None,
                  id_vib: tp.Optional[int] = None,
                  merge_mols: bool = True,
                  mol_repr: str = 'balls',
                  mol_mater: str = 'plastic',
                  vib_repr: str = 'arrows',
                  vib_mater: str = 'plastic',
                  verbose: bool = True,
                  **params):
        """Write POV-Ray file.

        Writes POV-Ray description file.

        Parameters
        ----------
        povname
            POV-Ray filename/template.  Overrides internal one.
        id_mol
            Identifier of the molecule (starting at 0).
            If not provided, all molecules are plotted.
        id_vib
            Identifier of the vibration of interest (starting at 0).
            If not provided, no vibration is shown.
        merge_mols
            Merge molecules in a single file.
        mol_repr
            Representation model for atoms.
        mol_mater
            Material to be used for the molecule.
        vib_repr
            Representation model for vibrations.
        vib_mater
            Material to be used for the vibration.
        verbose
            Prints details on operations and helps.
        params
            Additional parameters, transmitted to:

            - self.load_data
            - write_pov_head
            - write_pov_mol
            - write_pov_vib
        """
        fmt_write_file = 'Writing POV-Ray description scene in file: {file}'
        msg_cmds = '''
To run POVRay, you can use the following input (on Unix systems):
povray +W1600 +H1200 +A +UA -D {file}

+W: Set the width of the image, in pixels
+H: Set the height of the image, in pixels
+A: Enable anti-aliasing for smoother images
+UA: Activate background transparency
-D: deactivate the display of the image while it is generated.
    Remove the option to see the image being generated in real-time.
'''
        msg_cmds_win = '''
For Windows, if you run POVRay for Windows, the easiest way would be to use
the graphical interface if provided.
'''
        # Initial check on arguments validity/consistency
        if id_mol is not None:
            mol = id_mol
            try:
                if not self.__loaded[id_mol]:
                    self.load_data(infile=id_mol, **params)
            except IndexError as err:
                raise ArgumentError(f'Molecule id. {id_mol} does not exist.') \
                    from err
        else:
            if self.__nmols == 1:
                mol = 0
            if not all(self.__loaded):
                self.load_data()
        vib = None
        if id_vib is not None:
            if id_mol is None and self.__nmols > 1:
                raise ArgumentError(
                    'Vibrations for multiple molecules not yet available.')
            if self.__vmodes[mol] is None:
                if id_mol is not None:
                    msg = 'Vibrational modes not available for molecule ' \
                        + f'{id_mol}'
                    raise ArgumentError(msg)
                else:
                    raise ArgumentError(
                        'Vibrational modes not available for this molecule')
            try:
                vib = self.__vmodes[mol][id_vib].reshape(-1, 3)
            except IndexError as err:
                raise ArgumentError(
                    f'Cannot find vibration with index {id_vib}') from err
        # Build name template
        if povname:
            _povfile = povname
        elif self.__povfile:
            _povfile = self.__povfile
        else:
            if id_mol is None:
                if self.__nmols == 1:
                    _povfile = 'molecule'
                elif merge_mols:
                    _povfile = 'molecules'
                else:
                    _povfile = 'molecule_{mol}'
            else:
                if self.__nmols == 1:
                    _povfile = 'molecule'
                else:
                    _povfile = 'molecule_{mol}'
            if id_vib is not None:
                lmax = len(str(self.__vmodes[id_mol].shape[0]+1))
                _povfile += f'_vib{{vib:0{lmax}d}}'
            _povfile += '.pov'
        # First, consider all possible cases where only 1 molecule needs to be
        # printed
        if self.__nmols == 1 or id_mol is not None:
            box = build_box(self.__atlab[mol], self.__atcrd[mol], mol_repr)
            if vib is not None:
                ivib = id_vib+1
            else:
                ivib = 0
                _povfile = re.sub(r'\{vib:?[^}]*\}', 'NA', _povfile)
            _file = _povfile.format(mol=mol+1, vib=ivib)
            with open(_file, 'w', encoding='utf-8') as fobj:
                if verbose:
                    print(fmt_write_file.format(file=_file))
                zcam = set_cam_zpos(box)
                write_pov_head(fobj, zcam, **params)
                write_pov_mol(fobj, 1, self.__atlab[mol], self.__atcrd[mol],
                              self.__bonds[mol],
                              col_bond_as_atom=True,
                              representation=mol_repr, material=mol_mater,
                              show_mol=False,
                              **params)
                if vib is not None:
                    write_pov_vib(fobj, self.__atcrd[mol], vib,
                                  representation=vib_repr, material=vib_mater,
                                  **params)
        else:
            if merge_mols:
                box = build_box(self.__atlab[0], self.__atcrd[0], mol_repr)
                for i in range(1, self.__nmols):
                    box_ = build_box(self.__atlab[i], self.__atcrd[i],
                                     mol_repr)
                    if box['xmin'] > box_['xmin']:
                        box['xmin'] = box_['xmin']
                    if box['ymin'] > box_['ymin']:
                        box['ymin'] = box_['ymin']
                    if box['zmin'] > box_['zmin']:
                        box['zmin'] = box_['zmin']
                    if box['xmax'] < box_['xmax']:
                        box['xmax'] = box_['xmax']
                    if box['ymax'] < box_['ymax']:
                        box['ymax'] = box_['ymax']
                    if box['zmax'] < box_['zmax']:
                        box['zmax'] = box_['zmax']
                # Check if format template present
                if '{mol}' in _povfile:
                    _povfile = re.sub(r'\{mol:?[^}]*\}', 'all', _povfile)
                if '{vib}' in _povfile:
                    _povfile = re.sub(r'\{vib:?[^}]*\}', 'NA', _povfile)
                _file = _povfile
                with open(_file, 'w', encoding='utf-8') as fobj:
                    if verbose:
                        print(fmt_write_file.format(file=_file))
                    zcam = set_cam_zpos(box)
                    write_pov_head(fobj, zcam, **params)
                    write_pov_mol(fobj, self.__nmols, self.__atlab,
                                  self.__atcrd, self.__bonds,
                                  representation=mol_repr, material=mol_mater,
                                  **params)
            else:
                if '{mol' not in _povfile:
                    lmax = len(str(self.__nmols+1))
                    _povfile = f'_mol_{{mol:0{lmax}d}}'.join(
                        os.path.splitext(_povfile))
                if '{vib}' in _povfile:
                    _povfile = re.sub(r'\{vib:?[^}]*\}', 'N/A', _povfile)
                for i in range(self.__nmols):
                    box = build_box(self.__atlab[i], self.__atcrd[i], mol_repr)
                    _file = _povfile.format(mol=i+1)
                    with open(_file, 'w', encoding='utf-8') as fobj:
                        if verbose:
                            print(fmt_write_file.format(file=_file))
                        zcam = set_cam_zpos(box)
                        write_pov_head(fobj, zcam, **params)
                        write_pov_mol(fobj, 1, self.__atlab[i],
                                      self.__atcrd[i], self.__bonds[i],
                                      representation=mol_repr,
                                      material=mol_mater, **params)
        if verbose:
            print(msg_cmds.format(file=_file))
            if os.sys.platform == 'win32':
                print(msg_cmds_win)


def build_box(at_lab: TypeAtLab,
              at_crd: TypeAtCrd,
              at_repr: str = 'balls') -> tp.Dict[str, float]:
    """Build the outer box containing the molecule.

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
    at_repr
        Atomic representation to build properly the box.

    Returns
    -------
    dict
        Dictionary with positions of `xmin`, `xmax`, `ymin`, `ymax`,
        `zmin`, `zmax`.
    """
    xmin, ymin, zmin, xmax, ymax, zmax = +inf, +inf, +inf, -inf, -inf, -inf
    if at_repr != 'sticks':
        rkey = 'rdvw' if at_repr == 'spheres' else 'rvis'
        atdat = atomic_data(*sorted(set(at_lab)))
        _rad = max([x[rkey] for x in atdat.values()])
    else:
        _rad = BONDDATA['rvis']*RAD_VIS_SCL
    for x, y, z in at_crd:
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
    """Set Z position of camera in POV-Ray.

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


def write_pov_head(fobj: tp.TextIO,
                   zcam: float = -8.0,
                   incl_material: str = 'all',
                   **kwargs):
    """Write the header for POV-Ray input file.

    Writes the header in a POV-Ray input file.

    Parameters
    ----------
    fobj
        File object where the header will be written.
    zcam
        Position of the camera along -Z.
    incl_material
        Include material definition.
        Possible values: glass, metal, plastic, all.
    """
    try:
        if not fobj.writable():
            raise OSError('Cannot write on file.')
    except AttributeError() as err:
        raise OSError('The file must be opened before') from err

    _mater = incl_material.lower()
    if _mater not in set(CHOICES_MATER) | {'all'}:
        _mater = 'all'

    fobj.write(f'''\
global_settings {{
    ambient_light rgb <0.200, 0.200, 0.200>
    max_trace_level 15
}}

#declare dist = 1.;
#declare Trans = 0.0;

background {{ 1.0 }}

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
''')
    if _mater in ('plastic', 'all'):
        fobj.write(MAT_PLASTIC)
    if _mater in ('glass', 'all'):
        fobj.write(MAT_GLASS)
    if _mater in ('metal', 'all'):
        fobj.write(MAT_METAL)


def write_pov_mol(fobj: tp.TextIO,
                  nmols: int,
                  at_lab: TypeAtLabM,
                  at_crd: TypeAtCrdM,
                  bonds: TypeBondsM,
                  col_bond_as_atom: bool = False,
                  representation: str = 'balls',
                  scale_atoms: float = 1.0,
                  scale_bonds: float = 1.0,
                  material: str = 'plastic',
                  show_mol: bool = True,
                  **kwargs):
    """Write molecular specifications in a POV-Ray file.

    Builds and writes in an opened POV-Ray file the molecular
    specifications.
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
    col_bond_as_atom
        If true, bonds are colored based on the connected atoms
    representation
        Molecular representation: balls, spheres, stick...
    scale_atoms
        Scaling factor of the atoms.
    scale_bonds
        Scaling factor of the bonds.
    material
        Material to use for the molecular representation.
    show_mol
        Add block to display the molecule.
    """
    try:
        if not fobj.writable():
            raise OSError('Cannot write on file.')
    except AttributeError() as err:
        raise OSError('The file must be opened before') from err

    if representation not in ('balls', 'sticks', 'spheres'):
        raise ArgumentError('Unsupported molecular representation.')

    if material not in CHOICES_MATER:
        raise ArgumentError('Unsupported material definition.')

    show_bonds = representation != 'spheres'

    bo_rad = BONDDATA['rvis']*RAD_VIS_SCL
    rad_at = scale_atoms
    rad_bo = scale_bonds

    if nmols == 1:
        list_atoms = sorted(set(at_lab))
    else:
        list_atoms = sorted(set([item for mol in at_lab for item in mol]))
    atdat = atomic_data(*list_atoms)

    fmt_mat_at = f'''
#declare mat_at_{{at:2s}} = \
{CHOICES_MATER[material]}({{c[0]:3d}}, {{c[1]:3d}}, {{c[2]:3d}})
'''
    fmt_mat_bo = '''
#declare mat_bo_{at:2s} = material {{ mat_{ref} }}
'''
    fmt_mat_mol = f'''
#declare mat_mol{{id:02d}} = \
{CHOICES_MATER[material]}({{c[0]:3d}}, {{c[1]:3d}}, {{c[2]:3d}})
'''
    fmt_rad_at = '''
#declare r_at_{at:2s} = {r:.6f}*scl_rat;
'''
    if nmols == 1:
        fmt_obj_bo = '''\
    cylinder {{ // Bond {at1}({id1})- {at2}({id2})
        <{xyz1[0]:14.6f},{xyz1[1]:14.6f},{xyz1[2]:14.6f}>,
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        r_bond
        material {{ mat_bo_{at1} }}
    }}
    cylinder {{ // Bond {at1}({id1}) -{at2}({id2})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        <{xyz2[0]:14.6f},{xyz2[1]:14.6f},{xyz2[2]:14.6f}>,
        r_bond
        material {{ mat_bo_{at2} }}
    }}
'''
        fmt_obj_at = '''\
    sphere {{ // Atom {at}({id})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>, r_at_{at}
        material {{ mat_at_{at} }}
    }}
'''
    else:
        fmt_obj_bo = '''\
    cylinder {{ // Bond {at1}({id1})- {at2}({id2})
        <{xyz1[0]:14.6f},{xyz1[1]:14.6f},{xyz1[2]:14.6f}>,
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        r_bond
        material {{ mat_mol{idmol:02d} }}
    }}
    cylinder {{ // Bond {at1}({id1}) -{at2}({id2})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>,
        <{xyz2[0]:14.6f},{xyz2[1]:14.6f},{xyz2[2]:14.6f}>,
        r_bond
        material {{ mat_mol{idmol:02d} }}
    }}
'''
        fmt_obj_at = '''\
    sphere {{ // Atom {at}({id})
        <{xyz[0]:14.6f},{xyz[1]:14.6f},{xyz[2]:14.6f}>, r_at_{at}
        material {{ mat_mol{idmol:02d} }}
    }}
'''

    # Add declarations specific to molecular representation
    fobj.write(f'#declare scl_rat = {rad_at:.3f};\n')
    if show_bonds:
        fobj.write(f'#declare scl_bond = {rad_bo:.3f};\n')

    # Case for 1 molecules
    if nmols == 1:
        # Set atoms materials
        fobj.write('\n// MATERIALS FOR ATOMS\n')
        for atom in list_atoms:
            fobj.write(fmt_mat_at.format(at=atdat[atom]['symb'],
                                         c=atdat[atom]['rgb']))
        fobj.write('\n')
        # Set bonds materials
        if show_bonds:
            fobj.write('\n// BONDS COLORS AND TEXTURES\n')
            if not col_bond_as_atom:
                fmt = f'''
#declare mat_bond = \
{CHOICES_MATER[material]}({{c[0]:3d}}, {{c[1]:3d}}, {{c[2]:3d}})
'''
                fobj.write(fmt.format(c=BONDDATA['rgb']))
                fmt = 'bond'
            else:
                fmt = 'at_{at}'
            for atom in list_atoms:
                at = atdat[atom]['symb']
                fobj.write(fmt_mat_bo.format(at=at, ref=fmt.format(at=at)))
    else:
        fobj.write('\n// MOLECULES COLORS AND TEXTURES\n')
        # Set Molecule materials
        for i in range(nmols):
            fobj.write(fmt_mat_mol.format(id=i+1, c=MOLCOLS[i]))
    # Defines bonds/atoms radi
    if show_bonds:
        fobj.write('\n// ATOMIC AND BOND RADII\n')
        fmt = '#declare r_bond  = {:.6f}*scl_bond;\n'
        fobj.write(fmt.format(bo_rad))
    else:
        fobj.write('\n// ATOMIC RADII\n')
    for atom in list_atoms:
        if representation == 'sticks':
            rval = bo_rad
        elif representation == 'spheres':
            rval = atdat[atom]['rvdw']
        else:
            rval = atdat[atom]['rvis']*RAD_VIS_SCL
        fobj.write(fmt_rad_at.format(at=atdat[atom]['symb'], r=rval))
    # Molecule specification
    if nmols == 1:
        fobj.write('\n// MOLECULE DEFINITION\n\n')
        fobj.write('#declare molecule = union {\n')
        if show_bonds:
            # -- First build bonds
            for bond in bonds:
                iat1, iat2 = bond
                xyz1 = at_crd[iat1]
                xyz2 = at_crd[iat2]
                xyzmid = (xyz1 + xyz2)/2.
                fobj.write(fmt_obj_bo.format(
                    xyz=xyzmid,
                    at1=atdat[at_lab[iat1]]['symb'], id1=iat1+1, xyz1=xyz1,
                    at2=atdat[at_lab[iat2]]['symb'], id2=iat2+1,
                    xyz2=xyz2))
        # -- Next build atoms
        for iat, lab in enumerate(at_lab):
            fobj.write(fmt_obj_at.format(
                at=atdat[lab]['symb'], id=iat+1, xyz=at_crd[iat]))
        # -- Close and add molecules
        fobj.write('}\n')
        if show_mol:
            fobj.write('''
object {
    molecule
}
''')

    else:
        fobj.write('\n// MOLECULES DEFINITION\n')
        for imol in range(nmols):
            fobj.write(f'\n#declare mol{imol+1:02d} = union {{\n')
            # -- First build bonds
            if show_bonds:
                for bond in bonds[imol]:
                    iat1, iat2 = bond
                    xyz1 = at_crd[imol][iat1]
                    xyz2 = at_crd[imol][iat2]
                    xyzmid = (xyz1 + xyz2)/2.
                    fobj.write(fmt_obj_bo.format(
                        xyz=xyzmid,
                        at1=atdat[at_lab[imol][iat1]]['symb'],
                        id1=iat1+1, xyz1=xyz1,
                        at2=atdat[at_lab[imol][iat2]]['symb'],
                        id2=iat2+1, xyz2=xyz2,
                        idmol=imol+1))
            # -- Next build atoms
            for iat in range(len(at_lab[imol])):
                fobj.write(fmt_obj_at.format(
                    at=atdat[at_lab[imol][iat]]['symb'], id=iat+1,
                    xyz=at_crd[imol][iat], idmol=imol+1))
            # -- Close
            fobj.write('}\n')
        # Add block to visualize molecules
        fobj.write('union {\n')
        for i in range(nmols):
            fobj.write(f'    object {{ mol{i+1:02d} }}\n')
        fobj.write('}\n')


def write_pov_vib(fobj: tp.TextIO,
                  at_crd: TypeAtCrdM,
                  evec: Type1Vib,
                  representation: str = 'arrow',
                  scale_mode: float = 1.0,
                  material: str = 'glass',
                  show_obj: str = 'both',
                  **kwargs):
    """Write normal mode description in a POV-Ray file.

    Builds and writes in an opened POV-Ray file the specifications to
    represent a vibration.

    Parameters
    ----------
    at_crd
        Coordinates of the atomic centers, in Ang.
    evec
        Description of the vibrational displacement from each center.
    representation
        Representation of the displacement: arrows, spheres
    scale_mode
        Scaling factor to be applied to the displacements.
    material
        Material to use for the vibrational representation.
    show_obj
        Add commands to show object(s). Possible values:

        * both: show combined molecule and vibration
        * vib: show only vibration
        * no: do not show
    """
    try:
        if not fobj.writable():
            raise OSError('Cannot write on file.')
    except AttributeError() as err:
        raise OSError('The file must be opened before') from err

    if material not in CHOICES_MATER:
        raise ArgumentError('Unsupported material definition.')

    fobj.write(f'\n#declare scl_vib = {scale_mode:.3f};\n')

    if representation in ('arrow', 'arrows'):
        fobj.write(OBJ_ARROW.format(mat=CHOICES_MATER[material]))
        fmt_vib = '''\
    object {{
        Build_Arrow({lvib})
        scale scl_vib
        matrix < {rot[0][0]:9.6f}, {rot[0][1]:9.6f}, {rot[0][2]:9.6f},
                 {rot[1][0]:9.6f}, {rot[1][1]:9.6f}, {rot[1][2]:9.6f},
                 {rot[2][0]:9.6f}, {rot[2][1]:9.6f}, {rot[2][2]:9.6f},
                 {trans[0]:9.6f}, {trans[1]:9.6f}, {trans[2]:9.6f} >
    }}
'''

    elif representation in ('sphere', 'spheres'):
        fobj.write(OBJ_2COL_SPHERE.format(mat=CHOICES_MATER[material]))
        fmt_vib = '''\
    object {{
        vib_sphere
        scale {lvib}*scl_vib
        matrix < {rot[0][0]:9.6f}, {rot[0][1]:9.6f}, {rot[0][2]:9.6f},
                 {rot[1][0]:9.6f}, {rot[1][1]:9.6f}, {rot[1][2]:9.6f},
                 {rot[2][0]:9.6f}, {rot[2][1]:9.6f}, {rot[2][2]:9.6f},
                 {trans[0]:9.6f}, {trans[1]:9.6f}, {trans[2]:9.6f} >
    }}
'''
    else:
        raise ArgumentError('Unsupported representation of vibrations.')

    fobj.write('\n#declare vib = union{\n')
    for xyz, dxyz in zip(at_crd, evec):
        ldxyz = np.linalg.norm(dxyz)
        rot = np.transpose(vrotate_3D(np.array([0, 1, 0]), dxyz/ldxyz))
        fobj.write(fmt_vib.format(rot=rot.tolist(), trans=xyz.tolist(),
                                  lvib=ldxyz))
    fobj.write('}\n')

    # Conclude with options for showing
    if show_obj.lower() == 'both':
        fobj.write('''
union {
    object { molecule }
    object { vib }
}
''')
    elif show_obj.lower() == 'vib':
        fobj.write('''
object {
    vib
}
''')
