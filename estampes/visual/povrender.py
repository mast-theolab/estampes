"""Persistence of Vision Raytracing: Rendering tools

This module provides functions and instruments to build input files for
POV-Ray (http://www.povray.org/).
"""

from math import inf, tan, radians
import typing as tp

from estampes.base.types import TypeAtCrd, TypeAtCrdM, TypeAtLab, TypeAtLabM, \
    TypeBondsM
from estampes.base.errors import ArgumentError
from estampes.data.atom import atomic_data
from estampes.data.visual import BONDDATA, MOLCOLS, RAD_VIS_SCL


def build_box(at_lab: TypeAtLab,
              at_crd: TypeAtCrd,
              rad_atom_as_bond: bool = False) -> tp.Dict[str, float]:
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
        _rad = max(atdat, key=lambda x: x['rvis'])
    else:
        _rad = atrad
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


def write_pov(fname: str,
              nmols: int,
              at_lab: TypeAtLabM,
              at_crd: TypeAtCrdM,
              bonds: TypeBondsM,
              zcam: float = -8.0,
              col_bond_as_atom: bool = False,
              representation: str = 'balls',
              scale_atoms: float = 1.0,
              scale_bonds: float = 1.0):
    """Write a Pov-Ray file.

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
    representation
        Molecular representation: balls, spheres, stick...
    scale_atoms
        Scaling factor of the atoms.
    scale_bonds
        Scaling factor of the bonds.
    """
    if representation not in ('balls', 'sticks', 'spheres'):
        raise ArgumentError('Unsupported molecular representation.')

    show_bonds = representation != 'spheres'

    bo_rad = BONDDATA['rvis']*RAD_VIS_SCL
    rad_at = scale_atoms
    rad_bo = scale_bonds

    if nmols == 1:
        list_atoms = sorted(set(at_lab))
    else:
        list_atoms = sorted(set([item for mol in at_lab for item in mol]))
    atdat = atomic_data(*list_atoms)
    header = '#declare dist = 1.;\n'
    if show_bonds:
        header += '#declare scl_bond = {scl_bo:.2f};\n'
    header += '''\
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

    with open(fname, 'w', encoding="utf-8") as fobj:
        fobj.write(header.format(scl_bo=rad_bo, scl_at=rad_at, zcam=zcam))
        # Set atoms colors
        if nmols == 1:
            fobj.write('\n// ATOMS COLORS AND TEXTURES\n')
            for atom in list_atoms:
                fobj.write(fmt_col_at.format(at=atdat[atom]['symb'],
                                             c=atdat[atom]['rgb']))
            fobj.write('\n')
            # Set atoms texture
            for atom in list_atoms:
                fobj.write(fmt_tex_at.format(at=atdat[atom]['symb']))
            # Set bonds colors/textures
            if show_bonds:
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
                    at = atdat[atom]['symb']
                    fobj.write(fmt_tex_bo.format(at=at,
                                                 ref=fmt.format(at=at)))
        else:
            fobj.write('\n// MOLECULES COLORS AND TEXTURES\n')
            # Set Molecule colors
            for i in range(nmols):
                fobj.write(fmt_col_mol.format(id=i+1, c=MOLCOLS[i]))
            # Set Molecule textures
            for i in range(nmols):
                fobj.write(fmt_tex_mol.format(id=i+1))
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
            fobj.write('''}

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
