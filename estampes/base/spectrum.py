"""Provide basic tools to manipulate spectrum-related data.

A basic module providing the main class for manipulating spectral data.
"""

from math import ceil
import typing as tp

from estampes.parser import DataFile
from estampes.base import QLabel, QData, TypeColor
from estampes.tools.spec import broaden, convert_y


# =================
# Module Attributes
# =================

# SPEC2DATA provide some equivalence and default information
# name: full spectrosocpy name
# unit: default unit for the broadening
# qlabels for the spectroscopies

VSPC2DATA = {
    'IR': {
        'name': 'Infrared',
        'unit': 'I:/M/cm',
        'H': {'freq': QLabel(quantity='vlevel', level='H'),
              # 'int': QLabel(quantity='dipstr', level='H'),
              'int': QLabel(quantity='intens', descriptor='IR', level='H'),
              'assign': QLabel(quantity='vtrans', level='H')},
        'A': {'freq': QLabel(quantity='vlevel', level='A'),
              # 'int': QLabel(quantity='dipstr', level='A'),
              'int': QLabel(quantity='intens', descriptor='IR', level='A'),
              'assign': QLabel(quantity='vtrans', level='A')},
        'DS': 'Dipole strength',
        'II': 'Integrated intensity'
    },
    'VCD': {
        'name': 'Vibrational Circular Dichroism',
        'unit': 'I:/M/cm',
        'H': {'freq': QLabel(quantity='vlevel', level='H'),
              'int': QLabel(quantity='rotstr', level='H'),
              'assign': QLabel(quantity='vtrans', level='H')},
        'A': {'freq': QLabel(quantity='vlevel', level='A'),
              'int': QLabel(quantity='rotstr', level='A'),
              'assign': QLabel(quantity='vtrans', level='A')},
        'RS': 'Rotatory strength'
    },
    'RS0': {
        'name': 'Raman Scattering',
        'unit': 'I:cm3/mol/sr',
        'H': {'freq': QLabel(quantity='vlevel', level='H'),
              'int': QLabel(quantity='ramact', descriptor='static', level='H'),
              'assign': QLabel(quantity='vtrans', level='H')},
        'A': {'freq': QLabel(quantity='vlevel', level='A'),
              'int': QLabel(quantity='ramact', descriptor='static', level='A'),
              'assign': QLabel(quantity='vtrans', level='A')},
        'RA': 'Raman activity'
    },
    'RS': {
        'name': 'Raman Scattering',
        'unit': 'I:cm3/mol/sr',
        'H': {'freq': QLabel(quantity='vlevel', level='H'),
              'int': QLabel(quantity='ramact', descriptor='dynamic',
                            level='H'),
              'assign': QLabel(quantity='vtrans', level='H')},
        'A': {'freq': QLabel(quantity='vlevel', level='A'),
              'int': QLabel(quantity='ramact', descriptor='dynamic',
                            level='A'),
              'assign': QLabel(quantity='vtrans', level='A')},
        'RA': 'Raman activity'
    },
    'ROA': {
        'name': 'Raman Optical Activity',
        'unit': 'I:cm3/mol/sr',
        'H': {'freq': QLabel(quantity='vlevel', level='H'),
              'int': QLabel(quantity='roaact', descriptor='dynamic',
                            level='H'),
              'assign': QLabel(quantity='vtrans', level='H')},
        'A': {'freq': QLabel(quantity='vlevel', level='A'),
              'int': QLabel(quantity='roaact', descriptor='dynamic',
                            level='A'),
              'assign': QLabel(quantity='vtrans', level='A')},
        'ROA': 'Raman optical activity'
    }
}

ESPC2DATA = {
    'OPA': {
        'name': 'One-Photon Absorption',
        'unit': 'I:/M/cm',
        'E': {'ener': QLabel(quantity=1, refstate=(0, 'a')),
              'int': QLabel(quantity='dipstr', refstate=(0, 'a'))},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')},
        'DS': 'Dipole strength'
    },
    'OPE': {
        'name': 'One-Photon Emission',
        'unit': 'I:uJ/mol',
        'E': {},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')}
    },
    'ECD': {
        'name': 'Electronic Circular Dichroism',
        'unit': 'I:/M/cm',
        'E': {'ener': QLabel(quantity=1, refstate=(0, 'a')),
              'int': QLabel(quantity='rotstr', refstate=(0, 'a'))},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')},
        'RS': 'Rotatory strength'
    },
    'CPL': {
        'name': 'Circularly Polarized Luminescence',
        'unit': 'I:uJ/mol',
        'E': {},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')}
    },
    'RR': {
        'name': 'Resonance Raman',
        'unit': None,
        'E': {},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')}
    },
    'RROA': {
        'name': 'Resonance Raman Optical Activity',
        'unit': None,
        'E': {},
        'H': {'spc': QLabel(quantity='fcdat', descriptor='Spec'),
              'par': QLabel(quantity='fcdat', descriptor='SpcPar'),
              'info': QLabel(quantity='fcdat', descriptor='SimInf')}
    },
}

SPEC2DATA = {**VSPC2DATA, **ESPC2DATA}


# ==============
# Module Classes
# ==============

class Spectrum():
    """Main class to manipulate spectra.

    Provides the basic instruments to access and manipulate spectral
    data.
    Only one spectroscopy can be stored.  Supported:

    :"OPA": One-Photon Absorption
    :"OPE": One-Photon Emission
    :"ECD": Electronic Circular Dichroism
    :"CPL": Circularly Polarized Luminescence
    :"RR": Resonance Raman
    :"RROA": Resonance Raman Optical Activity
    :"IR": Infrared
    :"VCD": Vibrational Circular Dichroism
    :"RS0": Raman Scattering (static)
    :"RS": Raman Scattering (dynamic)
    :"ROA": Raman Optical Activity

    Levels of theory:

    :"E[le[ctronic]]": Pure electronic level (only electronic trans.)
    :"H[arm]": Harmonic approximation of nuclear vibrations
    :"A[nh[arm]]": Anharmonic representation of nuclear vibrations

    To avoid any ambiguity, spectroscopies and levels of theory **must**
      be provided.

    Raises
    ------
    KeyError
        Unrecognized spectroscopy
    IndexError
        Mismatch in normal modes
    """

    def __init__(self, fname: tp.Union[str, DataFile],
                 specabbr: str,
                 level: str,
                 ylabel: tp.Optional[str] = None,
                 load_data: bool = True,
                 ftype: tp.Optional[str] = None,
                 **params: tp.Dict[str, tp.Any]):
        """Initialize Spectrum instance.

        Initializes an instance of the class Spectrum.

        Parameters
        ----------
        fname
            Filename from which data are to be extracted.
        specabbr
            Spectroscopy name (abbreviated version).
        level
            Level of theory:
            - `Anharm`, `Anh`, `A`: anharmonic
            - `Harm, `H`: harmonic (including Franck-Condon)
            - `Electronic`, `Ele`, `E`: pure electronic
        ylabel
            If multiple spectra may be present, label for the Y axis.
            Ex: TD OP spectra w/o temperature, Raman scattering setups.
        load_data
            If True, load spectral data in memory
        ftype
            File type (sent to `DataFile`).
            If `ftype` is `CSV`, `specabbr` and `level` are ignored.
        params
            Spectroscopy-specific parameters:

            Raman/ROA

                :"incfrq": incident frequency
                :"setup": experimental setup (e.g., SCP(180))
        """
        # Initialization internal parameters
        if isinstance(fname, DataFile):
            self.__dfile = fname
        else:
            self.__dfile = DataFile(fname, ftype)
        self.__spec = specabbr.upper()
        _level = level.upper()
        if _level in ('H', 'HARM', 'HARMONIC'):
            self.__theory = 'H'
        elif _level in ('A', 'ANH', 'ANHARM', 'ANHARMONIC'):
            self.__theory = 'A'
        elif _level in ('E', 'ELE', 'ELEC', 'ELECTRONIC'):
            self.__theory = 'E'
        else:
            raise KeyError('Unrecognized level of theory: '+level)
        self.__ylab = ylabel
        # Spectroscopy-specific parameters
        self.__params = {}
        if params:
            for key, val in params.items():
                if key.lower() in ('incfrq', 'incfreq'):
                    self.__params['incfrq'] = val
                elif key.lower() == 'setup':
                    self.__params['setup'] = val
        # Initialize main data array
        self.__xaxis = [None, None]
        self.__yaxis = [None, None]
        self.__xlabel = [None, None]
        self.__ylabel = [None, None]
        self.__xunit = [None, None]
        self.__yunit = [None, None]
        self.__broad = [{'func': None, 'hwhm': None},
                        {'func': None, 'hwhm': None}]
        self.__idref = 0
        self.__info = None
        self.__broad_ok = None
        self.__ytags = None
        # Initialize basic data
        if load_data:
            self.load_data(ylabel)
        # Create aliases
        self.reset()
        self.__linecol = None
        self.__linesty = None
        self.__linewdt = None

    def load_data(self, ylabel: tp.Optional[str] = None):
        """Load data from data file.

        Parameters
        ----------
        ylabel
            If multiple spectra may be present, label for the Y axis.
            Ex: TD OP spectra w/o temperature, Raman scattering setups.

        Raises
        ------
        KeyError
            `ylabel` not found.
        """
        def dobj_to_ydat(dobj: QData, ykey: str,
                         yidx: tp.Optional[int] = None) -> tp.Any:
            if yidx is not None:
                return dobj.get(ykey)[yidx]
            else:
                return dobj.get(ykey)

        if self.__spec not in SPEC2DATA:
            raise KeyError('Unsupported spectroscopy: '+self.__spec)
        if self.__theory not in SPEC2DATA[self.__spec]:
            msg = 'Level of theory not available for spectroscopy.'
            raise KeyError(msg)
        if ylabel is None:
            ylabel = self.__ylab
        if self.__dfile.version[0] == 'CSV':
            qkeys = {
                'spc': QLabel(quantity='anyspc', descriptor='Spec'),
                'par': QLabel(quantity='anyspc', descriptor='SpcPar')
            }
        else:
            qkeys = SPEC2DATA[self.__spec][self.__theory]
            self.__fullspec = SPEC2DATA[self.__spec]['name']
            if not qkeys:
                raise NotImplementedError('Keywords not available')
        dobjs = self.__dfile.get_data(**qkeys)
        if 'spc' in qkeys:
            if isinstance(dobjs['par'].get('x'), list):
                # Processing if multiple blocks:
                # TODO: Add check that parameters and spectra are consistent
                #       on Y as well
                nblocks = len(dobjs['par'].get('x'))
                if len(dobjs['spc'].get('x')) != nblocks:
                    msg = 'Inconsistency between spectral ranges and ' \
                        + 'parameters.'
                    raise KeyError(msg)
                # X axis should always be the same
                # TODO: Add proper check that X data are consistent.
                self.__xaxis[0] = dobjs['spc'].get('x')[0]
                self.__xlabel[0] = dobjs['par'].get('x')[0]
                # For Y, generate key information to retrieve data
                ykeys = (key for key in dobjs['spc'].extra_fields()
                         if key.startswith('y'))
                nidx = len(ykeys[0]) - 1
                yfmt = f'y{{:0{nidx:d}d}}'
                offset = 0
                __yindexes = {}
                for bloc in range(nblocks):
                    ylabels = [item for item in ykeys
                               if (item[0] == 'y' and
                                   len(dobjs['par'].get(item)) > bloc)]
                    for yax in ylabels:
                        i = int(yax[1:])
                        y = yfmt.format(i+offset)
                        __yindexes[y] = {'ykey': yax, 'yidx': bloc}
                    offset += len(ylabels)
            else:
                self.__xaxis[0] = dobjs['spc'].get('x')
                self.__xlabel[0] = dobjs['par'].get('x')
                __yindexes = {key: {'ykey': key}
                              for key in dobjs['spc'].extra_fields()
                              if key.startswith('y')}
            self.__xunit[0] = dobjs['par']['unitx']
            self.__ytags = {key: dobj_to_ydat(dobjs['par'], **val)
                            for key, val in __yindexes.items()}
            if ylabel is None:
                _ydat = __yindexes[__yindexes.keys()[0]]
            else:
                try:
                    _ydat = __yindexes[ylabel.lower()]
                except KeyError as err:
                    raise KeyError('Non-existent Y label') from err
            self.__yaxis[0] = dobj_to_ydat(dobjs['spc'], **_ydat)
            self.__ylabel[0] = dobj_to_ydat(dobjs['par'], **_ydat)
            self.__yunit[0] = dobjs['par'].get('unity')
            self.__label = self.__yunit[0]
            self.__broad[0]['func'] = dobjs['par'].get('func')
            self.__broad[0]['hwhm'] = dobjs['par'].get('hwhm')
            self.__broad_ok = self.__broad[0]['func'] == 'stick'
            if 'info' in qkeys:
                self.__info = dobjs['info'].extra_fields()
            else:
                self.__info = None
        elif 'freq' in qkeys:
            modes = []
            for key in dobjs['freq'].data:
                if isinstance(key, int):
                    modes.append(key)
            modes.sort()
            self.__xaxis[0] = []
            self.__yaxis[0] = []
            if self.__spec in ('RS', 'ROA'):
                incfrq = self.__params.get('incfrq')
                if incfrq is None:
                    for key in dobjs['int'].data:
                        try:
                            _ = float(key)
                            incfrq = key
                            break
                        except ValueError:
                            continue
                    if incfrq is None:
                        raise KeyError('No incident frequency data found.')
                    else:
                        self.__params['incfrq'] = incfrq
                else:
                    if incfrq not in dobjs['int'].dobjs:
                        vals = [item for item in dobjs['int'].data
                                if item.replace('.', '', 1).isdigit()]
                        fmt = '''No incident frequency matches the value {}
Available: {}'''
                        raise KeyError(fmt.format(incfrq, ', '.join(vals)))
                setup = self.__params.get('setup')
                if setup is None:
                    for key in dobjs['int'].data[incfrq].keys():
                        if key in ('SCP(180)', 'SCP(180)u'):
                            setup = key
                            break
                    if setup is None:
                        setup = dobjs['int'].data[incfrq].keys()[0]
                    self.__params['setup'] = setup
                else:
                    if setup not in dobjs['int'].data[incfrq].keys():
                        vals = dobjs['int'].data[incfrq].keys()
                        fmt = '''Setup "{}" not found
Available: {}'''
                        raise KeyError(fmt.format(setup, vals))
                ydata = dobjs['int'].data[incfrq][setup]
            else:
                ydata = dobjs['int'].data
            for mode in modes:
                try:
                    if (isinstance(dobjs['freq'].data[mode], float) and
                            isinstance(ydata[mode], float)):
                        self.__xaxis[0].append(dobjs['freq'].data[mode])
                        self.__yaxis[0].append(ydata[mode])
                except KeyError as err:
                    fmt = 'Inconsistency in the list of states between ' \
                        + 'energies and intensities.\nState {} was not ' \
                        + 'found in one of the lists.'
                    raise IndexError(fmt.format(mode)) from err
            self.__xlabel[0] = 'Wavenumbers'
            self.__xunit[0] = dobjs['freq'].unit
            self.__yunit[0] = dobjs['int'].unit
            self.__ylabel[0] = \
                SPEC2DATA[self.__spec][self.__yunit[0].split(':')[0]]
            self.__label = self.__yunit[0]
            self.__broad[0]['func'] = 'stick'
            self.__broad[0]['hwhm'] = None
            self.__broad_ok = self.__broad[0]['func'] == 'stick'
        elif 'ener' in qkeys:
            self.__xaxis[0] = []
            self.__yaxis[0] = []
            for i, ener in enumerate(dobjs['ener'].data):
                self.__xaxis[0].append(ener)
                self.__yaxis[0].append(dobjs['int'].data[i])
            self.__xlabel[0] = 'Energy'
            self.__xunit[0] = dobjs['ener'].unit
            self.__yunit[0] = dobjs['int'].unit
            self.__ylabel[0] = \
                SPEC2DATA[self.__spec][self.__yunit[0].split(':')[0]]
            self.__label = self.__yunit[0]
            self.__broad[0]['func'] = 'stick'
            self.__broad[0]['hwhm'] = None
            self.__broad_ok = self.__broad[0]['func'] == 'stick'
        else:
            raise NotImplementedError()

    def get_yaxes(self) -> tp.Dict[str, str]:
        """Return all available Y axes in data file."""
        if self.__ytags is None:
            self.load_data()
        return self.__ytags

    def reset(self):
        """Reset default axes to original axes.

        Resets the pointers/aliases to the original data.
        Note that this is only possible if the data have not been
        overwritten.
        """
        self.__idref = 0
        self.__xaxis[1] = None
        self.__yaxis[1] = None
        self.__xlabel[1] = None
        self.__ylabel[1] = None

    def overwrite_axis(self, data: tp.Sequence[tp.Union[float, int]],
                       axis: str = 'x'):
        """Overwrite one of the original axes.

        Note that this should only be used in special cases.

        Parameters
        ----------
        data
            New axis data.
        axis
            Axis to overwrite.

        Raises
        ------
        TypeError
            Cannot convert to list.
        IndexError
            Inconsistency in length size between the axes.
        """
        if axis.lower() == 'x':
            self.__xaxis[0] = list(data)
        elif axis.lower() == 'y':
            self.__yaxis[0] = list(data)
        if len(self.__xaxis[0]) != len(self.__yaxis[0]):
            raise IndexError('Inconsistency in the size of the axes')
        self.reset()

    def get_xaxis(self, origin: bool = False) -> tp.List[float]:
        """Show the X axis.

        Parameter
        ---------
        origin
            Use original X axis instead of current one.
        """
        if not origin:
            return self.__xaxis[self.__idref]
        else:
            return self.__xaxis[0]
    xaxis = property(get_xaxis)

    def get_yaxis(self, origin: bool = False) -> tp.List[float]:
        """Show the Y axis.

        Parameters
        ----------
        ylabel
            Y axis of interest (string label).
        origin
            Use original Y axis instead of current one.

        Raises
        ------
        KeyError
            Unrecognized label.
        """
        if not origin:
            return self.__yaxis[self.__idref]
        else:
            return self.__yaxis[0]
    yaxis = property(get_yaxis)

    def get_xunit(self, origin: bool = False) -> str:
        """Get the X unit.

        Parameter
        ---------
        origin
            Use unit for original X axis instead of current one.
        """
        i = self.__idref if not origin else 0
        return '{} / {}'.format(self.__xlabel[i],
                                self.__xunit[i].split(':')[-1])
    xunit = property(get_xunit)

    def get_yunit(self, origin: bool = False) -> str:
        """Get the Y unit.

        Parameter
        ---------
        origin
            Use unit for original Y axis instead of current one.
        """
        i = self.__idref if not origin else 0
        return '{} / {}'.format(self.__ylabel[i],
                                self.__yunit[i].split(':')[-1])
    yunit = property(get_yunit)

    @property
    def label(self) -> str:
        """Get or set the label for the spectrum."""
        return self.__label

    @label.setter
    def label(self, val: str):
        self.__label = val

    def get_broadening(self, info: tp.Optional[str] = None,
                       origin: bool = False
                       ) -> tp.Union[tp.Dict[str, tp.Union[str, float]], str,
                                     float]:
        """Show the broadening data.

        Parameters
        ----------
        info
            Information to return: func, hwhm
            If `None`, return all datas as a dictionary.
        origin
            Use original broadening data instead of current one.

        Raises
        ------
        KeyError
            Unrecognized label.
        """
        if info is None:
            if not origin:
                res = self.__broad[self.__idref]
            else:
                res = self.__broad[0]
        else:
            if info.lower() in ('f', 'func', 'function'):
                if not origin:
                    res = self.__broad[self.__idref]['func']
                else:
                    res = self.__broad[0]['func']
            elif info.lower() in ('hw', 'hwhm'):
                if not origin:
                    res = self.__broad[self.__idref]['hwhm']
                else:
                    res = self.__broad[0]['hwhm']
            else:
                raise KeyError('Unrecognized broadening parameter.')
        return res

    def set_broadening(self, hwhm: tp.Optional[float] = None,
                       func: tp.Optional[str] = None,
                       yunit: str = 'normalized',
                       xres: float = 10.0,
                       xmin: tp.Optional[float] = None,
                       xmax: tp.Optional[float] = None,
                       origin: bool = False):
        """Set the broadening parameters.

        Parameters
        ----------
        hwhm
            Half-width at half maximum, in unit of X.
        func
            Broadening function.
        yunit
            Unit for the Y axis, compatible with `convert_y`, except:
            * `n`, `norm`, `normalize` for normalized intensities.
            * `default`: use SPEC2DATA[...]['unit'].
        xres
            Resolution for the generation of the X axis.
        xmin
            Lower bound of the X axis.
        xmax
            Upper bound of the X axis.
        origin
            Use original broadening data instead of current one.

        Notes
        -----
        It is possible to override the original data only if the
        function is `None`.

        Raises
        ------
        ValueError
            Incorrect parameters for the broadening or X resolution.
        """
        # Deactivate Pylint being annoying with xaxis for the bounds
        # pylint: disable=unsubscriptable-object
        if not origin:
            _ids = 1
        else:
            _ids = 0
        if hwhm is not None:
            self.__broad[_ids]['hwhm'] = float(hwhm)
        if func is not None:
            if func.lower() in ('g', 'gau', 'gaussian'):
                _func = 'gaussian'
            elif func.lower() in ('l', 'lor', 'lorentzian'):
                _func = 'lorentzian'
            else:
                raise ValueError('Broadening function not supported.')
            self.__broad[_ids]['func'] = _func
        elif not origin and self.__broad[_ids]['func'] is None:
            if self.__spec in ('RR', 'RROA', 'IR', 'VCD', 'RS', 'RS0', 'ROA'):
                self.__broad[_ids]['func'] = 'lorentzian'
            else:
                self.__broad[_ids]['func'] = 'gaussian'
        if not origin:
            if self.__broad[_ids]['hwhm'] is None:
                raise ValueError('HWHM not set for the broadening.')
            if xmin is None:
                _xmin = self.__xaxis[0][0] - 10*self.__broad[_ids]['hwhm']
            else:
                _xmin = xmin
            if xmax is None:
                _xmax = self.__xaxis[0][-1] + 10*self.__broad[_ids]['hwhm']
            else:
                _xmax = xmax
            if xres <= 0.0:
                raise ValueError('Wrong value for `xres`.')
            if _xmax < _xmin:
                _xmin, _xmax = _xmax, _xmin
            npoints = int(ceil((_xmax-_xmin)/xres))
            self.__xaxis[_ids] = [_xmin + i*xres for i in range(npoints)]
            if yunit.lower() in ('n', 'norm', 'normalized'):
                _unit_dest = SPEC2DATA[self.__spec]['unit']
                _yunit = True
            else:
                _yunit = False
                if yunit.lower() == 'default':
                    _unit_dest = SPEC2DATA[self.__spec]['unit']
                else:
                    _unit_dest = yunit
            d = {}
            if 'incfrq' in self.__params:
                d['incfrq'] = float(self.__params['incfrq'])
            _yfac, _xfun = convert_y(self.__spec, _unit_dest, self.__yunit[0],
                                     **d)
            self.__yaxis[_ids] = broaden(self.__xaxis[0], self.__yaxis[0],
                                         self.__xaxis[_ids],
                                         self.__broad[_ids]['func'],
                                         self.__broad[_ids]['hwhm'],
                                         _yfac, _xfun, _yunit, False)
            self.__idref = 1
            self.__xunit[self.__idref] = self.__xunit[0]
            self.__xlabel[self.__idref] = self.__xlabel[0]
            self.__yunit[self.__idref] = _unit_dest
            self.__ylabel[self.__idref] = self.__ylabel[0]

    @property
    def hwhm(self) -> tp.Optional[float]:
        """Get or set the HWHM for the broadening."""
        return self.get_broadening(info='hwhm')

    @hwhm.setter
    def hwhm(self, val: float):
        self.set_broadening(hwhm=val)

    @property
    def func(self) -> tp.Optional[str]:
        """Get or set the broadening function."""
        return self.get_broadening(info='func')

    @func.setter
    def func(self, val: str):
        self.set_broadening(func=val)

    def set_display(self, color: tp.Optional[TypeColor] = None,
                    linestyle: tp.Optional[str] = None,
                    linewidth: tp.Optional[float] = None):
        """Set the display parameters.

        Sets options related to display.

        Parameters
        ----------
        color
            Color of the curve.
        linestyle
            Linestyle of the curve.
        linewidth
            Linewidth of the curve.

        Raises
        ------
        TypeError
            Wrong type for input data.
        """
        if color is not None:
            self.__linecol = color
        if linestyle is not None:
            self.__linesty = linestyle
        if linewidth is not None:
            try:
                self.__linewdt = float(linewidth)
            except ValueError as err:
                raise TypeError('Incorrect type for linewidth') from err

    @property
    def linecolor(self) -> tp.Optional[TypeColor]:
        """Line color."""
        return self.__linecol

    @linecolor.setter
    def linecolor(self, val: TypeColor):
        self.set_display(color=val)

    @property
    def linestyle(self) -> tp.Optional[str]:
        """Line style."""
        return self.__linesty

    @linestyle.setter
    def linestyle(self, val: str):
        self.set_display(linestyle=val)

    @property
    def linewidth(self) -> tp.Optional[str]:
        """Line width."""
        return self.__linewdt

    @linewidth.setter
    def linewidth(self, val: str):
        self.set_display(linewidth=val)
