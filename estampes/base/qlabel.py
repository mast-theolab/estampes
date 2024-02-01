"""Define Quantity Label (QLabel) class.

Defines the QLabel class, used to communicate with the parser(s) to
retrieve the information of interest.
"""
import re
import typing as tp

from estampes.base.errors import ArgumentError
from estampes.base.types import TypeRSta


class QLabel():
    """The QLabel class to store quantity information.

    Stores and manages information on a quantity.
    """

    def __init__(self, qstring: tp.Optional[str] = None,
                 *,
                 quantity: tp.Optional[tp.Union[str, int]] = None,
                 descriptor: tp.Optional[str] = None,
                 derorder: tp.Optional[tp.Union[str, int]] = None,
                 dercoord: tp.Optional[str] = None,
                 refstate: tp.Optional[TypeRSta] = None,
                 level: tp.Optional[str] = None) -> None:
        """Build QLabel instance.

        Builds components of QLabel instance.

        Parameters
        ----------
        qstring
            Full quantity label as string (using all format).
        quantity
            General label, the type of quantity.  Can be an integer.
        descriptor
            Quantity descriptor: component, category, gauge...
        derorder
            Derivative order, if relevant (0: reference value).
        dercoord
            Reference coordinates system for the derivatives.
        refstate
            Reference electronic state or transition of the quantity.
        level
            Level of theory used to generate the quantity.
            Ex: vibrational (harm, anharm), electronic

        Notes
        -----
        `qstring` is incompatible with the quantity-specific parameters..
        If it is specified, the rest is simply ignored.
        """
        qkeys = ('tag', 'opt', 'ord', 'crd', 'sta', 'lvl')
        qdata = {key: None for key in qkeys}
        if qstring is not None:
            for i, item in enumerate(qstring.split(':')):
                if i > len(qkeys):
                    raise ArgumentError('Too many parameters in qstring')
                qdata[qkeys[i]] = item
        else:
            qdata['tag'] = quantity
            qdata['opt'] = descriptor
            qdata['ord'] = derorder
            qdata['crd'] = dercoord
            qdata['sta'] = refstate
            qdata['lvl'] = level

        if qdata['tag'] is None:
            raise ArgumentError('Missing quantity information to build QLabel')
        self.__qtype = None  # Quantity type
        self.__qdesc = None  # Quantity sub-description
        self.__dord = None  # Derivative order
        self.__dcrd = None  # Derivation coordinates
        self.__rsta = None  # Reference state(s) (incl. transition)
        self.__qlvl = None  # Level of theory used to generate quantity

        self.set(quantity=qdata['tag'], descriptor=qdata['opt'],
                 derorder=qdata['ord'], dercoord=qdata['crd'],
                 refstate=qdata['sta'], level=qdata['lvl'])

    @property
    def label(self) -> tp.Optional[str]:
        """Return quantity label."""
        return self.__qtype

    @property
    def kind(self) -> tp.Optional[tp.Union[str, int]]:
        """Return quantity descriptor (sub-reference to label)."""
        return self.__qdesc

    @property
    def derord(self) -> tp.Optional[int]:
        """Return derivative order."""
        return self.__dord

    @property
    def dercrd(self) -> tp.Optional[str]:
        """Return derivation coordinates."""
        return self.__dcrd

    @property
    def rstate(self) -> tp.Optional[TypeRSta]:
        """Return reference state or transition."""
        return self.__rsta

    @property
    def level(self) -> tp.Optional[str]:
        """Return reference level."""
        return self.__qlvl

    def __repr__(self):
        """Return the string representation of the QLabel object."""
        return self.build_qstring()

    def set(self, *, quantity: tp.Optional[tp.Union[str, int]] = None,
            descriptor: tp.Optional[str] = None,
            derorder: tp.Optional[tp.Union[str, int]] = None,
            dercoord: tp.Optional[str] = None,
            refstate: tp.Optional[TypeRSta] = None,
            level: tp.Optional[str] = None) -> None:
        """Set QLabel information.

        Sets the QLabel information by calling specific procedures after
        having checked data consistency.
        """
        if quantity is not None:
            self.reset()
            self.__set_quantity(quantity)

        if descriptor is not None or quantity is not None:
            if self.__qtype is None:
                raise ArgumentError('Quantity type is not set. Cannot proceed')
            self.__set_descriptor(descriptor)

        if ((derorder is not None or dercoord is not None)
                or quantity is not None):
            if self.__qtype is None:
                raise ArgumentError('Quantity type is not set. Cannot proceed')
            self.__set_deriv_info(derorder, dercoord)

        if refstate is not None or quantity is not None:
            if self.__qtype is None:
                raise ArgumentError('Quantity type is not set. Cannot proceed')
            self.__set_state_trans(refstate)

        if level is not None or quantity is not None:
            if self.__qtype is None:
                raise ArgumentError('Quantity type is not set. Cannot proceed')
            self.__set_level_info(level)

    def build_qstring(self) -> str:
        """Build qlabel string.

        Builds the old-style qlabel.
        """
        tags = []
        tags.append(str(self.__qtype))
        tags.append(self.__qdesc or '')
        tags.append(str(self.__dord) if self.__dord is not None else '')
        tags.append(self.__dcrd or '')
        if isinstance(self.__rsta, tuple):
            tags.append(f'{self.__rsta[0]}->{self.__rsta[1]}')
        else:
            tags.append(str(self.__rsta))
        tags.append(self.__qlvl or '')
        return ':'.join(tags)

    def reset(self) -> None:
        """Reset all quantity information in QLabel."""
        self.__qtype = None  # Quantity type
        self.__qdesc = None  # Quantity sub-description
        self.__dord = None  # Derivative order
        self.__dcrd = None  # Derivation coordinates
        self.__rsta = None  # Reference state(s) (incl. transition)
        self.__qlvl = None  # Level of theory used to generate quantity

    def __set_quantity(self,
                       quantity: tp.Union[str, int]) -> None:
        """Set quantity type.

        Supported quantity types:

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

        Parameters
        ----------
        quantity
            General label, the type of quantity.  Can be an integer.

        Notes
        -----
        * The system does very little testing to facilitate extensions
          and uses the beyond standard supported methods.
        """
        try:
            self.__qtype = int(quantity)
        except ValueError:
            self.__qtype = quantity.lower()

    def __set_descriptor(self, descriptor: tp.Optional[str]) -> None:
        """Set quantity descriptor.

        Sets the quantity descriptor, for instance a specific version of
        a quantity (e.g., gauge), a component for meta-quantities like
        Franck-Condon calculations data...

        Recognized sub-options:

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

        Parameters
        ----------
        descriptor
            Quantity descriptor: component, category, gauge...

        Notes
        -----
        * For unsupported quantity types, the function simply stores
          the descriptor to facilitate integration of new quantities.
        """
        # Label-specific keyword (sub-option)
        if self.__qtype in ('ramact', 'roaact'):
            if descriptor is None or not descriptor.strip():
                self.__qdesc = 'static'
            elif descriptor.lower() in ('all', 'dynamic', 'static'):
                self.__qdesc = descriptor.lower()
                if self.__qdesc == 'all':
                    self.__qdesc = 'dynamic'
                if self.__qtype == 'roaact' and self.__qdesc == 'static':
                    raise ArgumentError('Static ROA not supported')
            else:
                try:
                    val = float(descriptor)
                    self.__qdesc = descriptor
                except ValueError as err:
                    msg = f'Incorrect sup-opt for {self.__qtype}'
                    raise ArgumentError(msg) from err
        elif self.__qtype == 'intens':
            if descriptor is None or not descriptor.strip():
                self.__qdesc = 'IR'
            else:
                key = descriptor.upper()
                if key == 'IR':
                    self.__qdesc = 'IR'
                else:
                    raise ArgumentError('Unsupported spectral intensity')
        elif self.__qtype == 'hessdat':
            if descriptor is None or not descriptor.strip():
                self.__qdesc = 'freq'
            else:
                key = descriptor.lower()
                if key in ('freq', 'frq', 'eval', 'eigval', 'w'):
                    self.__qdesc = 'freq'
                elif key in ('redmass', 'rmass', 'redmas', 'u'):
                    self.__qdesc = 'redmas'
                else:
                    raise ArgumentError('Unsupported option for hessdat')
        elif self.__qtype in ('atcrd', 2, 1):
            if descriptor is None:
                self.__qdesc = 'last'
            else:
                if descriptor.lower() in ('all', 'last', 'first', 'scan'):
                    self.__qdesc = descriptor.lower()
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif self.__qtype == 'anyspc':
            if descriptor is None:
                self.__qdesc = 'Spec'
            else:
                key = descriptor.upper()
                if key in ('SPEC', 'SPECTRUM', 'SPECTRA'):
                    self.__qdesc = 'Spec'
                elif key in ('SPCPAR', 'SPCPARS', 'SPCPARAMS',
                             'SPCPARAMETERS'):
                    self.__qdesc = 'SpcPar'
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif self.__qtype == 'fcdat':
            if descriptor is None:
                self.__qdesc = 'JMat'
            else:
                key = descriptor.upper()
                if key in ('J', 'JMAT'):
                    self.__qdesc = 'JMat'
                elif key in ('JF', 'JMATF'):
                    self.__qdesc = 'JMatF'
                elif key in ('K', 'KVEC'):
                    self.__qdesc = 'KVec'
                elif key in ('A', 'AMAT', 'SRAMAT'):
                    self.__qdesc = 'SRAMat'
                elif key in ('B', 'BVEC', 'SRBVEC'):
                    self.__qdesc = 'SRBVec'
                elif key in ('C', 'CMAT', 'SRCMAT'):
                    self.__qdesc = 'SRCMat'
                elif key in ('D', 'DVEC', 'SRDVEC'):
                    self.__qdesc = 'SRDVec'
                elif key in ('E', 'EMAT', 'SREMAT'):
                    self.__qdesc = 'SREMat'
                elif key in ('SPEC', 'SPECTRUM', 'SPECTRA'):
                    self.__qdesc = 'Spec'
                elif key in ('CONV', 'CONVERGENCE', 'PROGRESSION'):
                    self.__qdesc = 'Conv'
                elif key in ('ASSIGN', 'ASSIGNMENT'):
                    self.__qdesc = 'Assign'
                elif key in ('REDDIM', 'RED-DIM', 'REDUCED',
                             'REDUCED-DIMENSION'):
                    self.__qdesc = 'RedDim'
                elif key in ('IS', 'ISGEOM', 'GEOMIS'):
                    self.__qdesc = 'GeomIS'
                elif key in ('FS', 'FSGEOM', 'GEOMFS'):
                    self.__qdesc = 'GeomFS'
                elif key in ('MS', 'MSGEOM', 'GEOMMS'):
                    self.__qdesc = 'GeomMS'
                elif key in ('EXTRGEOM', 'EXGEOM'):
                    self.__qdesc = 'ExGeom'
                elif key in ('SIMINF', 'SIMINFO', 'INFO'):
                    self.__qdesc = 'SimInf'
                elif key in ('SPCPAR', 'SPCPARS', 'SPCPARAMS',
                             'SPCPARAMETERS'):
                    self.__qdesc = 'SpcPar'
                elif key in ('E(0-0)', 'E0-0', 'DE(0-0)', 'DE0-0', 'E00',
                             'DE00'):
                    self.__qdesc = 'E(0-0)'
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif self.__qtype == 'vptdat':
            if descriptor is None:
                self.__qdesc = 'XMat'
            else:
                key = descriptor.upper()
                if key in ('X', 'XMAT'):
                    self.__qdesc = 'XMat'
                elif key == 'GMAT':
                    self.__qdesc = 'GMat'
                elif key in ('CICOEF', 'VARCOEF'):
                    self.__qdesc = 'CICoef'
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif self.__qtype in ('vtrans', 'vlevel'):
            if descriptor is not None:
                key = descriptor.upper()
                if key == 'SOS':
                    self.__qdesc = 'SOS'
                elif key in ('RR', 'RROA'):
                    self.__qdesc = 'RR'
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif self.__qtype in (101, 'dipstr', 'rotstr'):
            if descriptor is None:
                if self.__qtype == 'rotstr':
                    self.__qdesc = 'vel'
                else:
                    self.__qdesc = 'len'
            else:
                if descriptor[:3].lower() in ('len', 'vel'):
                    self.__qdesc = descriptor[:3].lower()
                else:
                    raise ArgumentError(
                        f'Incorrect sup-opt for {self.__qtype}')
        elif (self.__qtype in range(300, 310)
              or self.__qtype in range(1300, 1310)):
            if descriptor is None:
                self.__qdesc = 0
            else:
                try:
                    val = int(descriptor)
                    if val >= 0:
                        self.__qdesc = val
                    else:
                        fmt = 'Incorrect sup-opt for {}'
                        raise ArgumentError(fmt.format(self.__qtype))
                except ValueError as err:
                    msg = f'Incorrect sup-opt for {self.__qtype}'
                    raise ArgumentError(msg) from err
        else:
            self.__qdesc = descriptor

    def __set_deriv_info(self, derorder: tp.Optional[tp.Union[str, int]],
                         dercoord: tp.Optional[str]) -> None:
        """Set derivation information.

        Sets the derivation information:

        * Derivation order
        * Derivation coordinates

        Recognized derivation coordinates:

        +---------+----------------------------------------------------+
        |  Value  |   Description                                      |
        +=========+====================================================+
        |   X     | Derivatives with respect to Cartesian coordinates  |
        |   Q     | Derivatives with respect to normal coordinates     |
        |   I     | Derivatives with respect to internal coordinates   |
        |   QX    | Derivatives in ixed normal-Cartesian coordinates   |
        +---------+----------------------------------------------------+

        Parameters
        ----------
        derorder
            Derivative order, if relevant (0: reference value).
        dercoord
            Reference coordinates system for the derivatives.
        """
        # Derivative order
        if derorder is None:
            if isinstance(self.__qtype, int):
                self.__dord = None
            else:
                self.__dord = 0
        elif isinstance(derorder, int):
            self.__dord = derorder
        else:
            if not derorder.strip():
                if isinstance(self.__qtype, int):
                    self.__dord = None
                else:
                    self.__dord = 0
            else:
                try:
                    self.__dord = int(derorder)
                except ValueError as err:
                    raise ArgumentError('Incorrect derivative order') from err
        # Derivative coordinate
        if dercoord is None or not dercoord.strip():
            self.__dcrd = None
        else:
            self.__dcrd = dercoord.upper()
            if self.__dcrd not in ['X', 'Q', 'QX', 'I']:
                raise ValueError('Incorrect derivatives coordinates')

    def __set_state_trans(self, refstate: tp.Optional[TypeRSta]) -> None:
        """Set reference electronic state or transition.

        Sets the reference electronic state or transition for which the
        quantity was computed.

        Parameters
        ----------
        refstate
            Reference electronic state or transition of the quantity.

        Notes
        -----
        Special values:

            :"a": refers to all available electronic states
            :"c": refers to the current electronic state
        """
        # State(s)
        key = r'(?P<start>[aAcC]|\d+)?' \
            + r'(.*(?<=->)(?P<end>(?(start)([aAcC]?|\d*)|([aAcC]|\d+)))|)\s*$'
        if refstate is not None:
            res = re.match(key, refstate)
            if res is None:
                raise ArgumentError('Incorrect state specification')
            data = res.groupdict()
            if data['start'] == '':
                val1 = 'c'  # type: tp.Union[str, int]
            elif data['start'].lower() in ('c', 'a'):
                val1 = data['start'].lower()
            else:
                try:
                    val1 = int(data['start'])
                except ValueError as err:
                    msg = 'Incorrect reference state definition'
                    raise ArgumentError(msg) from err
            if data['end'] is None:
                self.__rsta = val1
            else:
                if data['end'] == '':
                    val2 = 'c'  # type: tp.Union[str, int]
                elif data['end'].lower() in ('c', 'a'):
                    val2 = data['end'].lower()
                else:
                    try:
                        val2 = int(data['end'])
                    except ValueError as err:
                        raise ValueError('Incorrect final state definition') \
                            from err
                if val2 == val1:
                    raise ValueError('Equivalent states in transition')
                self.__rsta = (val1, val2)
        else:
            if isinstance(self.__qtype, int):
                self.__rsta = 'c'
            else:
                self.__rsta = None

    def __set_level_info(self, level: tp.Optional[str]) -> None:
        """Set the reference level of theory.

        Sets the reference level of theory for the quantity.
        Typical values are: E(lectronic), H(armonic), A(nharmonic)

        Parameters
        ----------
        level
            Level of theory used to generate the quantity.
        """
        # Level of theory
        if level is None or not level.strip():
            if self.__qtype in ('dipstr', 'rotstr'):
                if isinstance(self.__rsta, tuple):
                    self.__qlvl = 'E'
                else:
                    self.__qlvl = 'H'
            elif self.__qtype in ('ramact', 'roaact', 'vtrans', 'vlevel'):
                self.__qlvl = 'H'
            elif self.__qtype in ('fcdat', 'hessvec', 'hessdat'):
                self.__qlvl = 'H'
            elif self.__qtype in ('vptdat', ):
                self.__qlvl = 'A'
            elif isinstance(self.__qtype, int) or self.__qtype in ('atcrd', ):
                self.__qlvl = 'E'
            else:
                self.__qlvl = None
        else:
            lvl = level.strip().upper()
            if lvl not in ('E', 'H', 'A', 'VE'):
                raise ArgumentError(f'Incorrect level for {self.__qtype}')
            self.__qlvl = lvl


TypeQInfo = tp.Dict[str, QLabel]
