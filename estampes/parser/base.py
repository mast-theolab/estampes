"""General module for the computational chemistry programs.

Provides general routines and wrappers for software-independent parsing
  operations.

Attributes
----------

Methods
-------
parse_qlabel
    Parses a full quantity label into sub-items.
build_qlabel
    Builds a `qlabel` from option elements.
reshape_dblock
    Reshapes data block.

Classes
-------
DataFile
    Main class for data parsing.

Notes
-----
- The submodules of parsers should provide a function get_data, which
  will be used in the main class.
"""

import os
import re
import typing as tp
from estampes.parser.gaussian import fchk, glog
from estampes.parser import xyz
from estampes.base import TypeQLab, TypeRSta


def parse_qlabel(qlabel: str) -> TypeQLab:
    """Parses a quantity label and returns a list of items.

    Parses a label in the form:
        id|label[:[subopt][:[ord][:[coord][:[state(s)]]]]]
    with
        id|label
            quantity identifier (see below) or textual label
        subopt
            sub-option (e.g., gauge, specific components, sub-field)
        ord
            derivative order (0 for reference value)
        coord
            reference coordinates
        state(s)
            reference electronic state or transition (as `i->j`)

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
           `a`: refers to all available electronic states
           `c`: refers to the current electronic state

    Notes
    -----
    Possible identifiers or labels
    =========  ==================================================
      Value    Description
    =========  ==================================================
            *Special labels*
    -------------------------------------------------------------
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
     HessVal    Hessian eigenvalues
     SWOpt      Software runtime options
     SWVer      Software version
     DipStr     Dipole strength
     RotStr     Rotatory strength
     RamAct     Raman activity
     ROAAct     ROA activity
     FCDat      Franck-Condon-related data
     VPTDat     VPT2-related data
            *Basic properties*
    -------------------------------------------------------------
         1      Energy
         2      Coordinates
            *Special electronic properties*
    -------------------------------------------------------------
        50      Non-adiabatic couplings
            *Special quantities*
    -------------------------------------------------------------
        91      Coriolis Couplings
        92      Rotation Matrix
        93      Transition vector
            *Static electric/mixed-field properties*
    -------------------------------------------------------------
       101      Electric dipole
       102      Magnetic dipole
       103      Polarizability tensor
       104      Optical rotations
       105      Dipole-quadrupole polarizability
       106      Hyperpolarizability
       107      Quadrupole
            *Magnetic-field properties*
    -------------------------------------------------------------
       201      Magnetic susceptibility       (NMR)
       202      Fake rotational g-Tensor      (EPR)
       203      NMR shielding tensors         (NMR)
       204      Spin-rotation tensors         (EPR)
       205      Anisotropic hyperfine tensors (EPR)
       206      Isotropic (Fermi) terms       (EPR)
       207      ESR g-tensor                  (EPR)
       208      Nuclear quadrupole tensors    (NMR)
       209      Isotropic Spin-Spin coupling  (NMR)
            *Dynamic (frequency-dependent) properties*
    -------------------------------------------------------------
       301      Polarizability Alpha(-w,w)
       302      Optical rotations
       303      Polarizability Alpha(w,0)
       304      Dipole-quadrupole polarizability
       305      Hyperpolarizability Beta(-w,w,0)
       306      Hyperpolarizability Beta(w,w,-2w)
    =========  ==================================================

    Sub-Options:
    ========  ========  =========================================
     Option     Sub      Description
    ========  ========  =========================================
     DipStr      H      Harmonic dipole strength
                 A      Anharmonic dipole strength
    -------------------------------------------------------------
     RotStr      H      Harmonic rotatory strength
                 A      Anharmonic rotatory strength
    -------------------------------------------------------------
     FCDat      JMat     Duschinsky matrix
                JMatF    Duschinsky matrix (full, only for reddim)
                KVec     Shift vector
               SRAMat    Sharp and Rosenstock A matrix
               SRBVec    Sharp and Rosenstock B vector
               SRCMat    Sharp and Rosenstock C matrix
               SRDVec    Sharp and Rosenstock D vector
               SREMat    Sharp and Rosenstock E matrix
               GeomIS    Initial-state geometry
               GeomFS    Final-state geometry
               GeomMS    Intermediate geometry
               ExGeom    Extrapolated geometry
               SimInf    Simulation information
                Spec     Final spectrum/a
               SpcLeg    Spectrum legend
               BShape    Band-shape broadening
                Conv     Intensity convergence/spectrum progress
               Assign    Transition assignment data
    -------------------------------------------------------------
     VPTDat     XMat     Anharmonic X matrix
                GMat     Variational correction matrix
    -------------------------------------------------------------
     AtCrd      all      All present geometries (scan, opt)
                last     Only the last geometry if more present
               first     Only the first geometry
                scan     Extract scan-related data
    ========  ========  =========================================


    Possible coordinates:
    =========  ==================================================
      Value     Description
    =========  ==================================================
       X        Derivatives with respect to Cartesian coordinates
       Q        Derivatives with respect to normal coordinates
       I        Derivatives with respect to internal coordinates
       QX       Derivatives in ixed normal-Cartesian coordinates
    =========  ==================================================
    """
    nparts = 5  # Number of parts expected in full label
    # Parse qlabel and build full list (filling missing information)
    qlist = [item or None for item in
             qlabel.split(':')]  # type: tp.List[tp.Any]
    qlist.extend((nparts-len(qlist))*[None])

    qty_tag = 0  # type: tp.Union[int, str]
    qty_opt = None  # type: tp.Union[str, None]
    der_ord = 0  # type: int
    der_crd = None  # type: tp.Union[str, None]
    ref_sta = 'c'  # type: TypeRSta
    # Main label
    try:
        qty_tag = int(qlist[0])
    except ValueError:
        qty_tag = qlist[0].lower()
    # Label-specific keyword (sub-option)
    if qty_tag in ('dipstr', 'rotstr'):
        if qlist[1] is None:
            qty_opt = 'H'
        else:
            if qlist[1][0].upper() in ('H', 'A'):
                qty_opt = qlist[1][0].upper()
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'atcrd':
        if qlist[1] is None:
            qty_opt = 'last'
        else:
            if qlist[1].lower() in ('all', 'last', 'first', 'scan'):
                qty_opt = qlist[1].lower()
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'fcdat':
        if qlist[1] is None:
            qty_opt = 'JMat'
        else:
            val = qlist[1].upper()
            if val in ('J', 'JMAT'):
                qty_opt = 'JMat'
            elif val in ('JF', 'JMATF'):
                qty_opt = 'JMatF'
            elif val in ('K', 'KVEC'):
                qty_opt = 'KVec'
            elif val in ('A', 'AMAT', 'SRAMAT'):
                qty_opt = 'SRAMat'
            elif val in ('B', 'BVEC', 'SRBVEC'):
                qty_opt = 'SRBVec'
            elif val in ('C', 'CMAT', 'SRCMAT'):
                qty_opt = 'SRCMat'
            elif val in ('D', 'DVEC', 'SRDVEC'):
                qty_opt = 'SRDVec'
            elif val in ('E', 'EMAT', 'SREMAT'):
                qty_opt = 'SREMat'
            elif val in ('SPEC', 'SPECTRUM', 'SPECTRA'):
                qty_opt = 'Spec'
            elif val in ('CONV', 'CONVERGENCE', 'PROGRESSION'):
                qty_opt = 'Conv'
            elif val in ('ASSIGN', 'ASSIGNMENT'):
                qty_opt = 'Assign'
            elif val in ('IS', 'ISGEOM', 'GEOMIS'):
                qty_opt = 'GeomIS'
            elif val in ('FS', 'FSGEOM', 'GEOMFS'):
                qty_opt = 'GeomFS'
            elif val in ('MS', 'MSGEOM', 'GEOMMS'):
                qty_opt = 'GeomMS'
            elif val in ('EXTRGEOM', 'EXGEOM'):
                qty_opt = 'ExGeom'
            elif val in ('SIMINF', 'SIMINFO', 'INFO'):
                qty_opt = 'SimInf'
            elif val in ('SPCLEG', 'LEGEND', 'SPCLEGEND', 'SPECLEGEND'):
                qty_opt = 'SpcLeg'
            elif val in ('BSHAPE', 'BANDSHAPE', 'BROADENING'):
                qty_opt = 'BShape'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 'vptdat':
        if qlist[1] is None:
            qty_opt = 'XMat'
        else:
            val = qlist.upper()
            if val in ('X', 'XMAT'):
                qty_opt = 'XMat'
            elif val == 'GMAT':
                qty_opt = 'GMat'
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag == 101:
        if qlist[1] is None:
            qty_opt = 'len'
        else:
            if qlist[1][:3].lower() in ('len', 'vel'):
                qty_opt = qlist[1][:3].lower()
            else:
                raise ValueError('Incorrect sup-opt for {}'.format(qty_tag))
    elif qty_tag in range(300, 310):
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
    if qlist[2] is None:
        der_ord = 0
    else:
        try:
            der_ord = int(qlist[2])
        except ValueError:
            raise ValueError('Incorrect derivative order')
    # Derivative coordinate
    if qlist[3] is None:
        der_crd = None
    else:
        der_crd = qlist[3].upper()
        if der_crd not in ['X', 'Q', 'QX', 'I']:
            raise ValueError('Incorrect derivatives coordinates')
    # State(s)
    if qlist[4] is not None:
        key = r'(?P<start>[aAcC]|\d+)?' \
            + r'(.*(?<=->)(?P<end>(?(start)([aAcC]?|\d*)|([aAcC]|\d+)))|)'
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
                val2 = data['start'].lower()
            else:
                try:
                    val1 = int(data['end'])
                except ValueError:
                    raise ValueError('Incorrect reference state definition')
            if val2 == val1:
                raise ValueError('Equivalent states in transition')
            ref_sta = (val1, val2)
    # Return parse data
    return qty_tag, qty_opt, der_ord, der_crd, ref_sta


def build_qlabel(qtag: tp.Union[str, int],
                 qopt: tp.Optional[str] = None,
                 dord: int = 0,
                 dcrd: tp.Optional[str] = None,
                 state: tp.Union[tp.Tuple[int], int, str] = 'c'
                 ) -> str:
    """Builds a `qlabel` from option elements.

    This is a very simple routine to quickly generate a *qlabel*.
    It should not be used without some control from `parse_qlabel`.
    The full qlabel is generated: qtag:[qopt]:[dord]:[dcrd]:[state].

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
    return '{}:{}:{}:{}:{}'.format(a, b, c, d, e)


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

    Attributes
    ----------
    version : tuple
        1. Program used to generate file/exchange format
        2. Program/Format-specific version information

    Methods
    -------
    get_data
        Extract data corresponding to the provided labels.
    """
    def __init__(self, filename: str,
                 filetype: tp.Optional[str] = None) -> None:
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
        elif ftype == 'xyz':
            self._dfile = xyz.FileXYZ(filename)
            self._parser = xyz
        else:
            raise NotImplementedError('Unsupported filetype')

    @property
    def version(self) -> tp.Tuple[str, tp.Any]:
        """Returns the program version

        The version is a tuple with:
        * The program name
        * The version (software-dependent)
        """
        return self._dfile.full_version

    def get_data(self,
                 *qlabels,
                 error_noqty: bool = True):
        return self._parser.get_data(self._dfile, *qlabels,
                                     error_noqty=error_noqty)

    def get_hess_data(self,
                      natoms: int,
                      get_evec: bool = True,
                      get_eval: bool = True,
                      mweigh: bool = True,
                      hessvec: tp.Optional[tp.List[float]] = None,
                      hessval: tp.Optional[tp.List[float]] = None,
                      atmass: tp.Optional[tp.List[float]] = None,
                      fccart: tp.Optional[tp.List[float]] = None
                      ) -> tp.Tuple[tp.Any]:
        return self._parser.get_hess_data(
            natoms=natoms, get_evec=get_evec, get_eval=get_eval, mweigh=mweigh,
            dfobj=self._dfile, hessvec=hessvec, hessval=hessval, atmass=atmass,
            fccart=fccart)
