"""Module providing basic tools to manipulate spectrum-related data

A basic module providing the main class for manipulating spectral data.

Attributes
----------

Methods
-------

Classes
-------
Spectrum
"""

from math import ceil
import typing as tp

from estampes.base import TypeColor
import estampes.parser as ep
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
        'H': {'freq': ep.build_qlabel('vlevel', 'H'),
              'int': ep.build_qlabel('dipstr', 'H'),
              'assign': ep.build_qlabel('vtrans', 'H')},
        'A': {'freq': ep.build_qlabel('vlevel', 'A'),
              'int': ep.build_qlabel('dipstr', 'A'),
              'assign': ep.build_qlabel('vtrans', 'A')},
        'DS': 'Dipole strength',
        'II': 'Integrated intensity'
    },
    'VCD': {
        'name': 'Vibrational Circular Dichroism',
        'unit': 'I:/M/cm',
        'H': {'freq': ep.build_qlabel('vlevel', 'H'),
              'int': ep.build_qlabel('rotstr', 'H'),
              'assign': ep.build_qlabel('vtrans', 'H')},
        'A': {'freq': ep.build_qlabel('vlevel', 'A'),
              'int': ep.build_qlabel('rotstr', 'A'),
              'assign': ep.build_qlabel('vtrans', 'A')},
        'RS': 'Rotatory strength'
    },
    'RS': {
        'name': 'Raman Scattering',
        'unit': None,
        'H': {},
        'A': {}
    },
    'ROA': {
        'name': 'Raman Optical Activity',
        'unit': None,
        'H': {},
        'A': {}
    }
}

ESPC2DATA = {
    'OPA': {
        'name': 'One-Photon Absorption',
        'unit': 'I:/M/cm',
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
    },
    'OPE': {
        'name': 'One-Photon Emission',
        'unit': 'I:uJ/mol',
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
    },
    'ECD': {
        'name': 'Electronic Circular Dichroism',
        'unit': 'I:/M/cm',
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
    },
    'CPL': {
        'name': 'Circularly Polarized Luminescence',
        'unit': 'I:uJ/mol',
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
    },
    'RR': {
        'name': 'Resonance Raman',
        'unit': None,
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
    },
    'RROA': {
        'name': 'Resonance Raman Optical Activity',
        'unit': None,
        'E': {},
        'H': {'spc': ep.build_qlabel('fcdat', 'Spec'),
              'par': ep.build_qlabel('fcdat', 'SpcPar'),
              'info': ep.build_qlabel('fcdat', 'SimInf')}
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
    * `OPA`: One-Photon Absorption
    * `OPE`: One-Photon Emission
    * `ECD`: Electronic Circular Dichroism
    * `CPL`: Circularly Polarized Luminescence
    * `RR`: Resonance Raman
    * `RROA`: Resonance Raman Optical Activity
    * `IR`: Infrared
    * `VCD`: Vibrational Circular Dichroism
    * `RS`: Raman Scattering
    * `ROA`: Raman Optical Activity
    Levels of theory:
    * `E[le[ctronic]]`: Pure electronic level (only electronic trans.)
    * `H[arm]`: Harmonic approximation of nuclear vibrations
    * `A[nh[arm]]`: Anharmonic representation of nuclear vibrations
    To avoid any ambiguity, spectroscopies and levels of theory *must*
      be provided.

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
        File type (sent to `ep.DataFile`).
        If `ftype` is `CSV`, `specabbr` and `level` are ignored.

    Attributes
    ----------
    xaxis
        Shows reference X axis.
    yaxis
        Shows reference Y axis.
    label
        Gets/sets label of the spectrum.
    hwhm
        Gets/sets HWHM for the broadening.
    func
        Gets/sets broadening function.
    linecolor
        Gets/sets the line color.
    linestyle
        Gets/sets the line style.
    linewidth
        Gets/sets the line width.

    Methods
    -------
    load_data(ylabel)
        Loads data from data file.
    reset()
        Resets reference axes to original axes.
    get_ytags()
        Returns all available Y axes in data file (will load if needed).
    get_xaxis(origin)
        Shows X axis.
    get_yaxis(origin)
        Shows Y axis.
    get_broadening(info, origin)
        Shows broadening data.
    set_broadening(hwhm, func, yunit, origin)
        Sets broadening data.
    set_display(color, linestyle, linewidth)
        Sets display parameters.

    Raises
    ------
    IndexError
        Unrecognized spectroscopy
    """
    def __init__(self, fname: str,
                 specabbr: str,
                 level: str,
                 ylabel: tp.Optional[str] = None,
                 load_data: bool = True,
                 ftype: tp.Optional[str] = None):
        # Initialization internal parameters
        self.__dfile = ep.DataFile(fname, ftype)
        self.__spec = specabbr.upper()
        _level = level.upper()
        if _level in ('H', 'HARM', 'HARMONIC'):
            self.__theory = 'H'
        elif _level in ('A', 'ANH', 'ANHARM', 'ANHARMONIC'):
            self.__theory = 'A'
        elif _level in ('E', 'ELE', 'ELEC', 'ELECTRONIC'):
            self.__theory = 'E'
        else:
            raise IndexError('Unrecognized level of theory: '+level)
        self.__ylab = ylabel
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
            if ftype is None or ftype.lower() != 'csv':
                self.load_data(ylabel)
            else:
                raise NotImplementedError('CSV File not yet supported.')
        # Create aliases
        self.reset()
        self.__linecol = None
        self.__linesty = None
        self.__linewdt = None

    def load_data(self, ylabel: tp.Optional[str] = None) -> tp.NoReturn:
        """Loads data from data file.

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
        if self.__spec not in SPEC2DATA:
            raise IndexError('Unsupported spectroscopy: '+self.__spec)
        if self.__theory not in SPEC2DATA[self.__spec]:
            msg = 'Level of theory not available for spectroscopy.'
            raise IndexError(msg)
        if ylabel is None:
            ylabel = self.__ylab
        if self.__dfile.version[0] == 'CSV':
            qkeys = {
                'spc': ep.build_qlabel('anyspc', 'Spec'),
                'par': ep.build_qlabel('anyspc', 'SpcPar')
            }
        else:
            qkeys = SPEC2DATA[self.__spec][self.__theory]
            self.__fullspec = SPEC2DATA[self.__spec]['name']
            if not qkeys:
                raise NotImplementedError('Keywords not available')
        data = self.__dfile.get_data(*qkeys.values())
        if 'spc' in qkeys:
            self.__xaxis[0] = data[qkeys['spc']]['x']
            self.__xlabel[0] = data[qkeys['par']]['x']
            self.__xunit[0] = data[qkeys['par']]['unitx']
            _yindexes = [y for y in data[qkeys['spc']] if y.startswith('y')]
            _yindexes.sort()
            self.__ytags = {y: data[qkeys['par']] for y in _yindexes}
            if ylabel is None:
                _ylab = _yindexes[0]
            else:
                try:
                    _id = [y.lower() for y in _yindexes].index(ylabel.lower())
                except ValueError:
                    raise KeyError('Non-existent Y label')
                _ylab = _yindexes[_id]
            self.__yaxis[0] = data[qkeys['spc']][_ylab]
            self.__ylabel[0] = data[qkeys['par']][_ylab]
            self.__yunit[0] = data[qkeys['par']]['unity']
            self.__label = self.__yunit[0]
            self.__broad[0]['func'] = data[qkeys['par']]['func']
            self.__broad[0]['hwhm'] = data[qkeys['par']]['hwhm']
            self.__broad_ok = self.__broad[0]['func'] == 'stick'
            if 'info' in qkeys:
                self.__info = data[qkeys['info']]
            else:
                self.__info = None
        elif 'freq' in qkeys:
            modes = []
            for key in data[qkeys['freq']]:
                if isinstance(key, int):
                    modes.append(key)
            modes.sort()
            self.__xaxis[0] = []
            self.__yaxis[0] = []
            for mode in modes:
                if isinstance(data[qkeys['freq']][mode], float) and \
                    isinstance(data[qkeys['int']][mode], float):
                    self.__xaxis[0].append(data[qkeys['freq']][mode])
                    self.__yaxis[0].append(data[qkeys['int']][mode])
            self.__xlabel[0] = 'Wavenumbers'
            self.__xunit[0] = data[qkeys['freq']]['unit']
            self.__yunit[0] = data[qkeys['int']]['unit']
            self.__ylabel[0] = \
                SPEC2DATA[self.__spec][self.__yunit[0].split(':')[0]]
            self.__label = self.__yunit[0]
            self.__broad[0]['func'] = 'stick'
            self.__broad[0]['hwhm'] = None
            self.__broad_ok = self.__broad[0]['func'] == 'stick'
        else:
            raise NotImplementedError()

    def get_yaxes(self) -> tp.Dict[str, str]:
        """Returns all available Y axes in data file."""
        if self.__ytags is None:
            self.load_data()
        return self.__ytags

    def reset(self):
        """Resets default axes to original axes.

        Resets the pointers/aliases to the original data.
        Note that this is only possible if the data have not been
          overwritten.
        """
        self.__idref = 0
        self.__xaxis[1] = None
        self.__yaxis[1] = None
        self.__xlabel[1] = None
        self.__ylabel[1] = None

    def get_xaxis(self, origin: bool = False) -> tp.List[float]:
        """Shows X axis.

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
        """Shows Y axis.

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

    @property
    def label(self) -> str:
        """Gets or sets label for the spectrum."""
        return self.__label

    @label.setter
    def label(self, val: str) -> tp.NoReturn:
        self.__label = val

    def get_broadening(self, info: tp.Optional[str] = None,
                       origin: bool = False
                       ) -> tp.Union[tp.Dict[str, tp.Union[str, float]], str,
                                     float]:
        """Shows broadening data.

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
                       origin: bool = False) -> tp.NoReturn:
        """Sets broadening parameters.

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
            if self.__spec in ('RR', 'RROA', 'IR', 'VCD', 'RS', 'ROA'):
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
            _yfac, _xpow = convert_y(self.__spec, _unit_dest, self.__yunit[0])
            self.__yaxis[_ids] = broaden(self.__xaxis[0], self.__yaxis[0],
                                         self.__xaxis[_ids],
                                         self.__broad[_ids]['func'],
                                         self.__broad[_ids]['hwhm'],
                                         _yfac, _xpow, _yunit, False)
            self.__idref = 1

    @property
    def hwhm(self) -> tp.Optional[float]:
        """Gets or sets the HWHM for the broadening."""
        return self.get_broadening(info='hwhm')

    @hwhm.setter
    def hwhm(self, val: float) -> tp.NoReturn:
        self.set_broadening(hwhm=val)

    @property
    def func(self) -> tp.Optional[str]:
        """Gets or sets the broadening function."""
        return self.get_broadening(info='func')

    @func.setter
    def func(self, val: str) -> tp.NoReturn:
        self.set_broadening(func=val)

    def set_display(self, color: tp.Optional[TypeColor] = None,
                    linestyle: tp.Optional[str] = None,
                    linewidth: tp.Optional[float] = None) -> tp.NoReturn:
        """Sets display parameters.

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
            except ValueError:
                raise TypeError('Incorrect type for linewidth')

    @property
    def linecolor(self) -> tp.Optional[TypeColor]:
        """Gets or sets the line color."""
        return self.__linecol

    @linecolor.setter
    def linecolor(self, val: TypeColor) -> tp.NoReturn:
        self.set_display(color=val)

    @property
    def linestyle(self) -> tp.Optional[str]:
        """Gets or sets the line style."""
        return self.__linesty

    @linestyle.setter
    def linestyle(self, val: str) -> tp.NoReturn:
        self.set_display(linestyle=val)

    @property
    def linewidth(self) -> tp.Optional[str]:
        """Gets or sets the line width."""
        return self.__linewdt

    @linewidth.setter
    def linewidth(self, val: str) -> tp.NoReturn:
        self.set_display(linewidth=val)
