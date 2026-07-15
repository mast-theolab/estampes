"""Program DERIVEUR - Common core module.

The module contains global variable and common
"""
import glob
import os
import re
import typing as tp

from estampes.base import QLabel, ArgumentError

KWORD_DISP = '[DISP]'
KWORD_STEP = '[STEP]'
KWORD_NUM = '[NUM]'
KWORD_XYZ = '[XYZ]'
KWORD_DIR = '[DIR]'
KWORD_GEOM = '[GEOM]'
KWORD_DESC = '[DESC]'
KWORD_CHK = '[FULLCHK]'
HELP_TMPL_FNAME = f'''\
- {KWORD_DISP}: replaced with the full displacement information.
- {KWORD_STEP}: replaced with the step number (normally for multi-steps \
displacements).
- {KWORD_NUM}: coordinate number (e.g., mode or atom).
- {KWORD_XYZ}: Cartesian coordinate for the displacement, where relevant.
- {KWORD_DIR}: direction (P: positive, N: negative).'''
HELP_TMPL_MAIN = f'''\
- {KWORD_GEOM} will be replaced with the generated geometry.
- {KWORD_DESC} to print information on the step.'''
HELP_TMPL_OPT = f'''\
- {KWORD_CHK} generates a full chk filename based on template filename and \
step.'''
HELP_QTY = """\
Multiple quantities can be given, separated by commas.
Supported quantities are (keywords after ->, aliases separated by commas):
- energy -> Energy, V
- electric dipole -> ElecDip, ED, u
- magnetic dipole -> MagDip, MD, m
- static polarizability tensor -> PolTens, PT, a
- frequency-dependent polarizability -> FDAlpha, FDA, a(w)
- frequency-dependent optical rotation -> FDOptRot, FDOR, O(w)
- frequency-dependent dipole-quadrupole -> FDDipQuad, FDDQ, Q(w)
Each quantity can be preceded by "D" followed by the order(s) of derivations.
- "D0" is equivalent to reference property
- ranges can be indicated with a "-" sign (ex: D2-3u)
Note: the keywords are case-insensitive.
"""
HELP_FMT = """\
Output format (case-insensitive).  Possible values are:
- InDataX: Output to be read by Gaussian with DataSrc=InDataX.
- InDataNM (default): Output to be read by Gaussian with DataSrc=InDataNM.
"""
LIM_STEP = 1.0e-8
XYZ = ('X', 'Y', 'Z')

QLABS_ATGEOM = {
    'atnum': QLabel(quantity='atnum'),
    'atcrd': QLabel(quantity='atcrd'),
    'atmas': QLabel(quantity='atmas')
}
KEYS_COORDS = {
    'Q': ('Q', 'modes', 'normal', 'mode', 'qweigh'),
    'q': ('q', 'qred', 'reduced'),
    'X': ('x', 'cart', 'cartesian')
}
INFO_PM_STEPS = {'P': ('up', 1.0), 'M': ('down', -1.0)}
LABEL_PSTEP = 'P'
LABEL_MSTEP = 'M'
THRESH_NUMERR = 0.1
THRESH_SMALL = 1.0e-8

USE_FORTRAN_FMT = True

# Format specification for component inside full displacement specification.
FMT_FDISP = '{}'
FMT_FSTEP = '_{}'
FMT_FNUM = '{}'
FMT_FXYZ = '{}'
FMT_FDIR = '_{}'


class TypePattInfo(tp.TypedDict):
    """Type for Pattern Information data."""

    fmt: str
    multi: bool


#
# UTILITIES
# =========
def get_tmpl_fmt_from_file(pattern: str, basedir: str = '.'
                           ) -> TypePattInfo:
    """Get template format specification from filenames.

    Generates a template format specification to be used to build
    filenames containing displaced geometries based on the pattern
    specification and the files present in the root directory.

    Parameters
    ----------
    pattern
        Filename pattern.
    basedir
        Base directory from which files are searched.

    Returns
    -------
    dict
        Format specification.
    """
    # Check that basedir exists
    if not os.path.exists(os.path.abspath(basedir)):
        raise ArgumentError('basedir',
                            f'Base director "{basedir}" does not exist.')

    # We now transfor all specifications to wildcards for the search
    path_mask = os.path.join(
        basedir,
        pattern.replace(KWORD_DISP, '*').replace(KWORD_STEP, '*')
        .replace(KWORD_NUM, '*').replace(KWORD_XYZ, '*')
        .replace(KWORD_DIR, '*'))
    # Check we do not end up with doubled wildcards, which could be
    # interpreted incorrectly in the search:
    while '**' in path_mask:
        path_mask = path_mask.replace('**', '*')

    # Now get list of files
    list_files = glob.glob(path_mask)
    if not list_files:
        raise FileNotFoundError('Could not find any file matching the pattern')

    # Check that the files have a consistent size
    # We expect that, by definitions, the filenames have greater size than the
    # wildcard-based pattern.
    files_by_len = {}
    for file in list_files[1:]:
        if (length := len(file)) > len(path_mask):
            if length not in files_by_len:
                files_by_len[length] = []
            files_by_len[length].append(file)

    # 1 or 2 files could be the reference and template files mixed up inside.
    # More is very unlikely correct, the user should clean up their stuff.
    lfiles = sorted(files_by_len, reverse=True)
    if sum(lfiles[1:]) > 2:
        raise ValueError('Inconsistent file lengths; this is confusing.')

    # Now let us try to construct the template.
    # We take the first file name and try to guess its structure.
    pat_coord = r'(?P<coord>\d+)'
    pat_step = r'(?P<step>\d+)'
    pat_xyz = r'(?P<xyz>X|Y|Z)'
    pat_dir = rf'(?P<dir>{LABEL_PSTEP}|{LABEL_MSTEP})'
    pat_disp = r'(?P<disp>' \
        + FMT_FNUM.format(pat_coord.replace('coord', 'fcoord')) \
        + FMT_FXYZ.format(pat_xyz.replace('xyz', 'fxyz')) + '?' \
        + FMT_FDIR.format(pat_dir.replace('dir', 'fdir')) + '?' \
        + '(?:' + FMT_FSTEP.format(pat_step.replace('step', 'fstep')) + ')?)'

    file_pattern = re.compile(os.path.basename(pattern)
                              .replace(KWORD_DISP, pat_disp)
                              .replace(KWORD_STEP, pat_step)
                              .replace(KWORD_NUM, pat_coord)
                              .replace(KWORD_XYZ, pat_xyz)
                              .replace(KWORD_DIR, pat_dir), re.I)

    fname_ref = os.path.basename(files_by_len[lfiles[0]][0])

    res = file_pattern.match(fname_ref)
    if res is None:
        raise ValueError('Could not match file pattern with existing files.')

    fmt_coord = ''
    fmt_fcoord = ''
    fmt_xyz = ''
    fmt_fxyz = ''
    fmt_dir = ''
    fmt_fdir = ''
    fmt_step = ''
    fmt_fstep = ''
    fmt_disp = ''
    comps = res.groupdict()
    if not comps.get('coord') and not comps.get('fcoord'):
        raise ValueError('Missing coordinate specifications in filenames.')
    if comps.get('coord'):
        fmt_coord = f'{{coord:0{len(comps["coord"])}d}}'
    if comps.get('fcoord'):
        fmt_fcoord = FMT_FNUM.format(f'{{coord:0{len(comps["fcoord"])}d}}')
    if comps.get('xyz'):
        fmt_xyz = '{xyz:1s}'
    if comps.get('fxyz'):
        fmt_fxyz = FMT_FXYZ.format('{xyz:1s}')
    if comps.get('dir'):
        fmt_dir = '{dir:1s}'
    if comps.get('fdir'):
        fmt_fdir = FMT_FDIR.format('{dir:1s}')
    if comps.get('step'):
        fmt_step = f'{{step:0{len(comps["step"])}d}}'
    if comps.get('fstep'):
        fmt_fstep = FMT_FSTEP.format(f'{{step:0{len(comps["fstep"])}d}}')
    if comps.get('disp'):
        fmt_disp = fmt_fcoord + fmt_fxyz + fmt_fdir + fmt_fstep

    return {
        'fmt': pattern.replace(KWORD_DISP, fmt_disp)
                      .replace(KWORD_STEP, fmt_step)
                      .replace(KWORD_NUM, fmt_coord)
                      .replace(KWORD_XYZ, fmt_xyz)
                      .replace(KWORD_DIR, fmt_dir),
        'multi': bool(fmt_fstep)
    }


def gen_fname_from_tmpl(n_coords: int,
                        n_steps: int,
                        i_coord: int,
                        i_step: int,
                        label_xyz: str | None = None,
                        label_dir: str | None = None) -> dict[str, str]:
    """Generate file name from template and step specifications..

    Checks step differentiation parameters and generates the filename.

    Parameters
    ----------
    fname_tmpl
        Template filename from which the final name is generated.
    n_coords
        Total number of unique displacement coordinates.
    n_steps
        Number of steps along a direction (half if symmetric displ.).
    i_coord
        Index of the coordinate.
    i_step
        Index of the step.
    label_xyz
        Label for the XYZ coordinates.
    label_dir
        Label for the displacement direction.

    Returns
    -------
    dict
        The dictionary associated each key supported for template
        filenames with the corresponding text.
    """
    fmt_coord = f'{{num:0{len(str(n_coords))}d}}'
    fmt_step = f'{{idx:0{len(str(n_steps))}d}}'

    txt_coord = fmt_coord.format(num=i_coord)
    txt_step = fmt_step.format(num=i_step) if n_steps > 1 else ''
    txt_dir = f'{label_dir:1s}' if label_dir is not None else ''
    txt_xyz = f'{label_xyz:1s}' if label_xyz is not None else ''

    txt_fcoord = FMT_FNUM.format(txt_coord) if txt_coord else ''
    txt_fxyz = FMT_FXYZ.format(txt_xyz) if txt_xyz else ''
    txt_fstep = FMT_FSTEP.format(txt_step) if txt_step else ''
    txt_fdir = FMT_FDIR.format(txt_dir) if txt_dir else ''
    txt_disp = txt_fcoord + txt_fxyz + txt_fdir + txt_fstep

    return {
        KWORD_DISP: txt_disp,
        KWORD_STEP: txt_step,
        KWORD_NUM: txt_coord,
        KWORD_XYZ: txt_xyz,
        KWORD_DIR: txt_dir}


def check_fname_pattern(name: str, use_cart: bool) -> str:
    """Check filename pattern.

    Checks that filename pattern is valid and return a clean version
    if necessary.

    If the placeholders are missing, one is added at the end of the
    filename, before the extension.

    Parameters
    ----------
    name
        Template for the file names.
    use_cart
        Displacement will be done along Cartesian coordinates.

    Returns
    -------
    str
        Checked/fixed template filename.

    Notes
    -----
    * Redundant specifications are allowed.
    """
    txt = name
    # Check full displacement specification
    if KWORD_DISP not in name:
        if KWORD_DISP.lower() in name:
            print(f'{KWORD_DISP} found with an incorrect case, fixing.')
            pattern = f'\\{KWORD_DISP[:-1]}\\]'
            txt = re.sub(pattern, KWORD_DISP, txt, flags=re.I)

    # Check full step specification
    if KWORD_STEP not in name:
        if KWORD_STEP.lower() in name:
            print(f'{KWORD_STEP} found with an incorrect case, fixing.')
            pattern = f'\\{KWORD_STEP[:-1]}\\]'
            txt = re.sub(pattern, KWORD_STEP, txt, flags=re.I)

    # Check placeholder for index/number
    if KWORD_NUM not in name:
        if KWORD_NUM.lower() in name:
            print(f'{KWORD_NUM} found with an incorrect case, fixing.')
            pattern = f'\\{KWORD_NUM[:-1]}\\]'
            txt = re.sub(pattern, KWORD_NUM, txt, flags=re.I)

    # Check placeholder for direction specification
    if KWORD_DIR not in name:
        if KWORD_DIR.lower() in name:
            print(f'{KWORD_DIR} found with an incorrect case, fixing.')
            pattern = f'\\{KWORD_DIR[:-1]}\\]'
            txt = re.sub(pattern, KWORD_DIR, txt, flags=re.I)

    # Check placeholder for Cartesian coordinate (only for Cart. displacement)
    has_xyz = KWORD_XYZ in name
    if not has_xyz:
        if use_cart:
            if KWORD_XYZ.lower() in name:
                print(f'{KWORD_XYZ} found with an incorrect case, fixing.')
                pattern = f'\\{KWORD_XYZ[:-1]}\\]'
                txt = re.sub(pattern, KWORD_XYZ, txt, flags=re.I)
                has_xyz = True
    if has_xyz and not use_cart:
        print('Cartesian coordinate specification will be ignored.')

    return txt
