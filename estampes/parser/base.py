"""General module for the computational chemistry programs.

Provides general routines and wrappers for software-independent parsing
operations.

Notes
-----
* The submodules of parsers should provide a function get_data, which
  will be used in the main class.
"""

import os
import re
import typing as tp

from estampes.base import TypeQData, TypeQLab, TypeRSta
from estampes.parser import csv, xyz
from estampes.parser.gaussian import glog, fchk


def parse_qlabel(qlabel: str) -> TypeQLab:
    """Parses a quantity label and returns a list of items.

    Parses a label in the form:

    .. code-block:: text

        id|label[:[subopt][:[ord][:[coord][:[state(s)][:[level]]]]]

    with

    ``id|label``
        quantity identifier (see below) or textual label
    ``subopt``
        sub-option (e.g., gauge, specific components, sub-field)
    ``ord``
        derivative order (0 for reference value)
    ``coord``
        reference coordinates
    ``state(s)``
        reference electronic state or transition (as ``"i->j"``)
    ``level``
        Level of theory use to generate the quantity.
        ex: vibrational (harm, anharm), electronic...

    Parameters
    ----------
    qlabel
        Full quantity label to parse

    Returns
    -------
    list
        Parsed quantity label, as (`None` if undefined):

        1. quantity identifier or label
        2. quantity-specific options
        3. derivative order
        4. reference coordinates
        5. reference state or transition

            :"a": refers to all available electronic states
            :"c": refers to the current electronic state

        6. level of theory: vibrational (harm, anharm), electronic...

    Notes
    -----
    Below is the list of possible identifiers or labels

    *Possible identifiers or labels*

    =========  =====================================================
      Value      Description
    =========  =====================================================
           *Special labels*
    ----------------------------------------------------------------
    Title      Title of the job
    NAtoms     Number of atoms
    NVib       Number of vibrations
    AtMas      Atomic masses
    AtNum      Atomic numbers
    AtLab      Atomic labels
    MolSym     Symmetry
    Charge     Charge
    Multip     Multiplicity
    AtCrd      Coordinates
    Atoms      Atoms (can be numbers or labels)
    HessVec    Hessian eigenvectors
    HessDat    Hessian-related data: frequencies, red. mass...
    SWOpt      Software runtime options
    SWVer      Software version
    VTrans     Vibrational transitions
    VLevel     Vibrational energy levels
    DipStr     Dipole strength
    RotStr     Rotatory strength
    RamAct     Raman activity
    ROAAct     ROA activity
    AnySpc     Generic quantity for spectra (e.g., for CSV files)
    FCDat      Franck-Condon-related data
    VPTDat     VPT2-related data
    Intens     Spectral intensity
    *Basic properties*
    ----------------------------------------------------------------
        1      Energy
        2      Coordinates
    *Special electronic properties*
    ----------------------------------------------------------------
        50     Non-adiabatic couplings
    *Special quantities*
    ----------------------------------------------------------------
        91     Coriolis Couplings
        92     Rotation Matrix
        93     Transition vector
    *Static electric/mixed-field properties*
    ----------------------------------------------------------------
       101     Electric dipole
       102     Magnetic dipole
       103     Polarizability tensor
       104     Optical rotations
       105     Dipole-quadrupole polarizability
       106     Hyperpolarizability
       107     Quadrupole
    *Magnetic-field properties*
    ----------------------------------------------------------------
       201     Magnetic susceptibility       (NMR)
       202     Fake rotational g-Tensor      (EPR)
       203     NMR shielding tensors         (NMR)
       204     Spin-rotation tensors         (EPR)
       205     Anisotropic hyperfine tensors (EPR)
       206     Isotropic (Fermi) terms       (EPR)
       207     ESR g-tensor                  (EPR)
       208     Nuclear quadrupole tensors    (NMR)
       209     Isotropic Spin-Spin coupling  (NMR)
    *Dynamic (frequency-dependent) properties*
    ----------------------------------------------------------------
       301     Polarizability Alpha(-w,w)
       302     Optical rotations
       303     Polarizability Alpha(w,0)
       304     Dipole-quadrupole polarizability
       305     Hyperpolarizability Beta(-w,w,0)
       306     Hyperpolarizability Beta(w,w,-2w)
    *Vibrational transition moments of properties*
    ----------------------------------------------------------------
      1300     List of incident frequencies
      1301     Electric dipole-electric dipole tensor
      1302     Induced electric dipole-magnetic dipole tensor
      1303     Electric dipole-induced magnetic dipole tensor
      1304     Induced electric dipole-el. quadrupole tensor
      1305     Electric dipole-induced el. quadrupole tensor
    =========  =====================================================

    **Sub-options**

    =========  ========  ===========================================
     Option     Sub       Description
    =========  ========  ===========================================
     Intens    IR         Infrared intensity
        101    Len        Length-gauge dipole strength

               Vel        Velocity-gauge dipole strength
     DipStr    Len        Length-gauge dipole strength

               Vel        Velocity-gauge dipole strength
     HessDat   freq       Frequency of each vibration

               redmas     Reduced mass of each vibration
    =========  ========  ===========================================

    **Possible coordinates**

    +---------+----------------------------------------------------+
    |  Value  |   Description                                      |
    +=========+====================================================+
    |   X     | Derivatives with respect to Cartesian coordinates  |
    |   Q     | Derivatives with respect to normal coordinates     |
    |   I     | Derivatives with respect to internal coordinates   |
    |   QX    | Derivatives in ixed normal-Cartesian coordinates   |
    +---------+----------------------------------------------------+
    """
    nparts = 6  # Number of parts expected in full label
    # Parse qlabel and build full list (filling missing information)
    qlist = [item or None for item in
             qlabel.split(':')]  # type: tp.List[tp.Any]
    qlist.extend((nparts-len(qlist))*[None])

    qty_tag = 0  # type: tp.Union[int, str]
    qty_opt = None  # type: tp.Union[str, None]
    der_ord = 0  # type: int
    der_crd = None  # type: tp.Union[str, None]
    ref_sta = 'c'  # type: TypeRSta
    qty_lvl = None  # type: tp.Union[str, None]
    # Main label
    try:
        qty_tag = int(qlist[0])
    except ValueError:
        qty_tag = qlist[0].lower()
    # Label-specific keyword (sub-option)
    if qty_tag in ('ramact', 'roaact'):
        if qlist[1] is None or not qlist[1].strip():
            qty_opt = 'static'
        elif qlist[1].lower() in ('all', 'dynamic', 'static'):
            qty_opt = qlist[1].lower()
            if qty_opt == 'all':
                qty_opt = 'dynamic'
            if qty_tag == 'roaact' and qty_opt == 'static':
                raise ValueError('Static ROA not supported')
        else:
            try:
                val = float(qlist[1])
                qty_opt = qlist[1]
            except ValueError:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'intens':
        if qlist[1] is None or not qlist[1].strip():
            qty_opt = 'IR'
        else:
            key = qlist[1].upper()
            if key == 'IR':
                qty_opt = 'IR'
            else:
                raise ValueError('Unsupported spectral intensity')
    elif qty_tag == 'hessdat':
        if qlist[1] is None or not qlist[1].strip():
            qty_opt = 'freq'
        else:
            key = qlist[1].lower()
            if key in ('freq', 'frq', 'eval', 'eigval', 'w'):
                qty_opt = 'freq'
            elif key in ('redmass', 'rmass', 'redmas', 'u'):
                qty_opt = 'redmas'
            else:
                raise ValueError('Unsupported spectral intensity')
    elif qty_tag in ('atcrd', 2, 1):
        if qlist[1] is None:
            qty_opt = 'last'
        else:
            if qlist[1].lower() in ('all', 'last', 'first', 'scan'):
                qty_opt = qlist[1].lower()
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'anyspc':
        if qlist[1] is None:
            qty_opt = 'Spec'
        else:
            key = qlist[1].upper()
            if key in ('SPEC', 'SPECTRUM', 'SPECTRA'):
                qty_opt = 'Spec'
            elif key in ('SPCPAR', 'SPCPARS', 'SPCPARAMS', 'SPCPARAMETERS'):
                qty_opt = 'SpcPar'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'fcdat':
        if qlist[1] is None:
            qty_opt = 'JMat'
        else:
            key = qlist[1].upper()
            if key in ('J', 'JMAT'):
                qty_opt = 'JMat'
            elif key in ('JF', 'JMATF'):
                qty_opt = 'JMatF'
            elif key in ('K', 'KVEC'):
                qty_opt = 'KVec'
            elif key in ('A', 'AMAT', 'SRAMAT'):
                qty_opt = 'SRAMat'
            elif key in ('B', 'BVEC', 'SRBVEC'):
                qty_opt = 'SRBVec'
            elif key in ('C', 'CMAT', 'SRCMAT'):
                qty_opt = 'SRCMat'
            elif key in ('D', 'DVEC', 'SRDVEC'):
                qty_opt = 'SRDVec'
            elif key in ('E', 'EMAT', 'SREMAT'):
                qty_opt = 'SREMat'
            elif key in ('SPEC', 'SPECTRUM', 'SPECTRA'):
                qty_opt = 'Spec'
            elif key in ('CONV', 'CONVERGENCE', 'PROGRESSION'):
                qty_opt = 'Conv'
            elif key in ('ASSIGN', 'ASSIGNMENT'):
                qty_opt = 'Assign'
            elif key in ('REDDIM', 'RED-DIM', 'REDUCED', 'REDUCED-DIMENSION'):
                qty_opt = 'RedDim'
            elif key in ('IS', 'ISGEOM', 'GEOMIS'):
                qty_opt = 'GeomIS'
            elif key in ('FS', 'FSGEOM', 'GEOMFS'):
                qty_opt = 'GeomFS'
            elif key in ('MS', 'MSGEOM', 'GEOMMS'):
                qty_opt = 'GeomMS'
            elif key in ('EXTRGEOM', 'EXGEOM'):
                qty_opt = 'ExGeom'
            elif key in ('SIMINF', 'SIMINFO', 'INFO'):
                qty_opt = 'SimInf'
            elif key in ('SPCPAR', 'SPCPARS', 'SPCPARAMS', 'SPCPARAMETERS'):
                qty_opt = 'SpcPar'
            elif key in ('E(0-0)', 'E0-0', 'DE(0-0)', 'DE0-0', 'E00', 'DE00'):
                qty_opt = 'E(0-0)'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'vptdat':
        if qlist[1] is None:
            qty_opt = 'XMat'
        else:
            key = qlist[1].upper()
            if key in ('X', 'XMAT'):
                qty_opt = 'XMat'
            elif key == 'GMAT':
                qty_opt = 'GMat'
            elif key in ('CICOEF', 'VARCOEF'):
                qty_opt = 'CICoef'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag in ('vtrans', 'vlevel'):
        if qlist[1] is not None:
            key = qlist[1].upper()
            if key == 'SOS':
                qty_opt = 'SOS'
            elif key in ('RR', 'RROA'):
                qty_opt = 'RR'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag in (101, 'dipstr', 'rotstr'):
        if qlist[1] is None:
            if qty_tag == 'rotstr':
                qty_opt = 'vel'
            else:
                qty_opt = 'len'
        else:
            if qlist[1][:3].lower() in ('len', 'vel'):
                qty_opt = qlist[1][:3].lower()
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag in range(300, 310) or qty_tag in range(1300, 1310):
        if qlist[1] is None:
            qty_opt = 0
        else:
            try:
                val = int(qlist[1])
                if val >= 0:
                    qty_opt = val
                else:
                    fmt = 'Incorrect sup-opt for {}'
                    raise ValueError(fmt.format(qty_tag))
            except ValueError:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    else:
        qty_opt = qlist[1]
    # Derivative order
    if qlist[2] is None or not qlist[2].strip():
        der_ord = 0
    else:
        try:
            der_ord = int(qlist[2])
        except ValueError:
            raise ValueError('Incorrect derivative order')
    # Derivative coordinate
    if qlist[3] is None or not qlist[3].strip():
        der_crd = None
    else:
        der_crd = qlist[3].upper()
        if der_crd not in ['X', 'Q', 'QX', 'I']:
            raise ValueError('Incorrect derivatives coordinates')
    # State(s)
    if qlist[4] is not None:
        key = r'(?P<start>[aAcC]|\d+)?' \
            + r'(.*(?<=->)(?P<end>(?(start)([aAcC]?|\d*)|([aAcC]|\d+)))|)\s*$'
        res = re.match(key, qlist[4])
        if res is None:
            raise ValueError('Incorrect state specification')
        else:
            data = res.groupdict()
        if data['start'] == '':
            val1 = 'c'  # type: tp.Union[str, int]
        elif data['start'].lower() in ('c', 'a'):
            val1 = data['start'].lower()
        else:
            try:
                val1 = int(data['start'])
            except ValueError:
                raise ValueError('Incorrect reference state definition')
        if data['end'] is None:
            ref_sta = val1
        else:
            if data['end'] == '':
                val2 = 'c'  # type: tp.Union[str, int]
            elif data['end'].lower() in ('c', 'a'):
                val2 = data['end'].lower()
            else:
                try:
                    val2 = int(data['end'])
                except ValueError:
                    raise ValueError('Incorrect final state definition')
            if val2 == val1:
                raise ValueError('Equivalent states in transition')
            ref_sta = (val1, val2)
    # Level of theory
    if qlist[5] is None or not qlist[5].strip():
        if qty_tag in ('dipstr', 'rotstr'):
            if isinstance(ref_sta, tuple):
                qty_lvl = 'E'
            else:
                qty_lvl = 'H'
        elif qty_tag in ('ramact', 'roaact', 'vtrans', 'vlevel'):
            qty_lvl = 'H'
        elif qty_tag in ('fcdat', 'hessvec', 'hessval'):
            qty_lvl = 'H'
        elif qty_tag in ('vptdat', ):
            qty_lvl = 'A'
        elif isinstance(qty_tag, int) or qty_tag in ('atcrd', ):
            qty_lvl = 'E'
        else:
            qty_lvl = None
    else:
        lvl = qlist[5].strip().upper()
        if lvl not in ('E', 'H', 'A', 'VE'):
            raise ValueError('Incorrect level for {}'.format(qty_tag))
        else:
            qty_lvl = lvl
    # Return parse data
    return qty_tag, qty_opt, der_ord, der_crd, ref_sta, qty_lvl


def build_qlabel(qtag: tp.Union[str, int],
                 qopt: tp.Optional[str] = None,
                 dord: int = 0,
                 dcrd: tp.Optional[str] = None,
                 state: tp.Union[tp.Tuple[int], int, str] = 'c',
                 level: tp.Optional[str] = None
                 ) -> str:
    """Builds a `qlabel` from option elements.

    This is a very simple routine to quickly generate a *qlabel*.
    It should not be used without some control from `parse_qlabel`.
    The full qlabel is generated::

        qtag:[qopt]:[dord]:[dcrd]:[state]:[level].

    Parameters
    ----------
    qtag
        Quantity tag (see data.property).
    qopt
        Quantity sub-option (where relevant).
    dord
        Derivative order.
    dcrd
        Derivation coordinate.
    state
        Reference electronic state or transition.
    level
        Level of theory used to generate the quantity.

    Returns
    -------
    str
        Full `qlabel`.
    """
    a = str(qtag)
    b = qopt or ''
    c = str(dord)
    d = dcrd or ''
    if isinstance(state, tuple):
        e = '{}->{}'.format(*state)
    else:
        e = str(state)
    f = level or ''
    return '{}:{}:{}:{}:{}:{}'.format(a, b, c, d, e, f)


def reshape_dblock(dblock: tp.Sequence[tp.Any],
                   dshape: tp.Sequence[int]
                   ) -> tp.List[tp.Any]:
    """Reshapes a data block.

    Reshapes a data block, normally extracted from fchk,
    assuming a Fortran-like array was stored in memory.
    If the shape is smaller than the length of `dblock`, the function
    tries to add a dimension, checking the consistency of this new
    dimension.

    Parameters
    ----------
    dblock
        Data block to reshape.
    dshape
        Shape of the returned data block.

        `0`/`-1` can be used for automatic definition of a dimension,

    Returns
    -------
    list
        Reshaped data block.

    Raises
    ------
    ValueError
        Raised if shape inconsistent with the size of `dblock`.
    """
    lshape = len(dshape)
    lblock = len(dblock)
    sshape = sum(dshape)
    # Check consistency between shape and data block size
    sshape = 1
    nzero = []
    for i, dim in enumerate(dshape):
        if dim <= 0:
            nzero.append(i)
        else:
            sshape *= dim
    _dshape = list(dshape)
    if nzero:
        if len(nzero) > 1:
            raise ValueError('Shape inconsistent with data block size')
        _dshape[nzero[0]] = lblock // sshape
        if _dshape[nzero[0]]*sshape != lblock:
            raise ValueError('Shape inconsistent with data block size')
    elif sshape != lblock:
        # If shape inconsistent, check if last dim was implicit
        # If so, product of other dimensions must be multiple with the
        #   full size
        new_dim = lblock // sshape
        if sshape*new_dim == lblock:
            lshape += 1
        else:
            raise ValueError('Shape inconsistent with data block size')
    # Reshape data block
    if lshape == 1:
        data = dblock
    elif lshape == 2:
        dim1 = _dshape[0]
        data = []
        for i in range(0, lblock, dim1):
            data.append([dblock[i+j] for j in range(dim1)])
    elif lshape == 3:
        dim1 = _dshape[0]
        dim2 = _dshape[1]
        dim12 = dim1*dim2
        data = []
        for i in range(0, lblock, dim12):
            data.append([])
            for j in range(0, dim12, dim1):
                data[-1].append([dblock[i+j+k] for k in range(dim1)])
    else:
        raise NotImplementedError()
    return data

# ==============
# Module Classes
# ==============


class DataFile(object):
    """Main class for data files.

    Acts as a wrapper to the format-/program-specific wrappers.

    Parameters
    ----------
    filename : str
        Filename.
    filetype : str, optional
        Filetype.
        Supported: 'fchk', 'glog', 'xyz'
    """
    def __init__(self, filename: str,
                 filetype: tp.Optional[str] = None):
        if filetype is None:
            ftype = os.path.splitext(filename)[1][1:].lower()
        else:
            ftype = filetype.lower()
        if ftype in ('fchk', 'fch'):
            self._dfile = fchk.FChkIO(filename)
            self._parser = fchk
        elif ftype in ('glog', 'log', 'out'):
            self._dfile = glog.GLogIO(filename)
            self._parser = glog
        elif ftype in ('csv', 'txt'):
            self._dfile = csv.FileCSV(filename)
            self._parser = csv
        elif ftype == 'xyz':
            self._dfile = xyz.FileXYZ(filename)
            self._parser = xyz
        else:
            raise NotImplementedError('Unsupported filetype')

    @property
    def filename(self) -> str:
        """Name of the file used to extract data."""
        return self._dfile.filename

    @property
    def version(self) -> tp.Tuple[str, tp.Any]:
        """Version of the program used to generate the file.

        The version is a tuple with:

        # The program name
        # Program/Format-specific version information
        """
        return self._dfile.full_version

    def get_data(self,
                 *qlabels,
                 error_noqty: bool = True,
                 **keys4qlabels) -> TypeQData:
        return self._parser.get_data(self._dfile, *qlabels,
                                     error_noqty=error_noqty,
                                     **keys4qlabels)

    def get_hess_data(self,
                      get_evec: bool = True,
                      get_eval: bool = True,
                      pre_data: tp.Optional[TypeQData] = None
                      ) -> tp.Tuple[tp.Any]:
        return self._parser.get_hess_data(self._dfile, get_evec=get_evec,
                                          get_eval=get_eval, pre_data=pre_data)
