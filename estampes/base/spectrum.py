"""Provide basic tools to manipulate spectrum-related data.

A basic module providing the main class for manipulating spectral data.
"""

from math import ceil
import typing as tp

from estampes.base.qdata import QData
from estampes.parser import DataFile
from estampes.base import QLabel, TypeColor
from estampes.base.errors import ArgumentError
from estampes.tools.char import unit_to_tex
from estampes.tools.spec import broaden, convert_y
from estampes.tools.vib import get_vib_trans, weigh_trans_progress


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

    The class can take one or more filenames, and in the latter case
    will combine them using the weights (or equal fractions if absent).

    Parameters
    ----------
    fname
        Filename from which data are to be extracted.
    specabbr
        Spectroscopy name (abbreviated version).
    level
        Level of theory:

        * anharmonic: ``"Anharm"``, ``"Anh"``, ``"A"``
        * harmonic (including Franck-Condon): ``"Harm"``, ``"H"``
        * pure electronic: ``"Electronic"``, ``"Ele"``, ``"E"``
    ylabel
        If multiple spectra may be present, label for the Y axis.
        Ex: TD OP spectra w/o temperature, Raman scattering setups.
    load_data
        If True, load spectral data in memory
    ftype
        File type (sent to `DataFile`).
        If `ftype` is `CSV`, `specabbr` and `level` are ignored.
    weigh_dsets
        Weigh data sets included in spectrum.
        The values can be provided as real numbers or with the keywords:

        * use electronic energies to set Boltzmann population
          weights: ``"bz"`` 
        * use harmonic ZPVE to set Boltzmann population weights:
          ``"bz_h"``
        * use anharmonic ZPVE to set Boltzmann population weights:
          ``"bz_a"``
    weigh_vtrans
        Weigh pure vibrational transitions, using one of the
        following methods:

        * use the harmonic probability of excitation of each mode,
          weighed by the quanta difference: ``"pni_dn"``
        * same as `"pni_dn"` but each transition is assumed to be a
          fundamental: ``"pni_d1"``
    params
        Spectroscopy-specific parameters:

        :"temp":
            Temperature for population estimates.
            By default, the temperature in first job is used.
        :"incfrq":
            (Raman/ROA) incident frequency.
        :"setup":
            (Raman/ROA) experimental setup (e.g., SCP(180))

    Notes
    -----
    If `incfrq` is given as integer (no decimal point), the
    incident frequencies are searched in data files truncated to
    the decimal point.  If `incfrq` is provided as float, the
    number are truncated to the lower integer.

    Raises
    ------
    KeyError
        Unrecognized spectroscopy.
    IndexError
        Mismatch in normal modes.
    """

    def __init__(self,
                 fname: tp.Union[str, DataFile,
                                 tp.Sequence[tp.Union[str, DataFile]]],
                 specabbr: str,
                 level: str,
                 ylabel: tp.Optional[str] = None,
                 load_data: bool = True,
                 ftype: tp.Optional[str] = None,
                 weigh_dbsets: tp.Optional[
                     tp.Union[str,
                              tp.Sequence[
                                tp.Union[str, float]]]] = None,
                 weigh_vtrans: tp.Optional[str] = None,
                 **params: tp.Dict[str, tp.Any]):
        self.__num_dfiles = 1
        # Initialization internal parameters
        # -- Build data files
        if isinstance(fname, DataFile):
            self.__dfiles = [fname]
        elif isinstance(fname, str):
            self.__dfiles = [DataFile(fname, ftype)]
        elif isinstance(fname, (tuple, list)):
            self.__dfiles = []
            self.__num_dfiles = len(fname)
            for name in fname:
                if isinstance(name, DataFile):
                    self.__dfiles.append(name)
                elif isinstance(name, str):
                    self.__dfiles.append(DataFile(fname, ftype))
                else:
                    raise TypeError(
                        f'Unsupported filename specification for {name}')
        else:
            raise TypeError(f'Unsupported filename specification for {fname}')
        # -- Set spectroscopy type
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
                    if isinstance(val, (int, float)):
                        self.__params['incfrq'] = str(int(val))
                    elif val is not None:
                        self.__params['incfrq'] = val
                elif key.lower() == 'setup':
                    if val is not None:
                        self.__params['setup'] = val
                elif key.lower() == 'temp':
                    self.__params['temp'] = val
        # Initialize main data array
        self.__xaxes = []
        self.__yaxes = []
        self.__xlabels = []
        self.__ylabels = []
        self.__xunits = []
        self.__yunits = []
        self.__broad_infos = []
        self.__spec_infos = []
        self.__xaxis = [None, None]
        self.__yaxis = [None, None]
        self.__xlabel = [None, None]
        self.__ylabel = [None, None]
        self.__xunit = [None, None]
        self.__yunit = [None, None]
        self.__broad_info = [None, None]
        self.__extras = []
        self.__spec_info = None
        self.__ytags = []
        self.__idref = 0
        self.__loaded = None
        self.__label = None
                # -- Set transition weights
        self.__weigh_vtrans = None
        if self.__spec in VSPC2DATA or self.__spec in ('RR', 'RROA'):
            if weigh_vtrans is not None:
                if weigh_vtrans.lower() in ('pni_dn', 'pni_d1'):
                    self.__weigh_vtrans = weigh_vtrans.lower()
                else:
                    raise ArgumentError(
                        'weigh_vtrans',
                        'Unrecognized weighing system for vibrational '
                        + 'transitions')
                if 'temp' not in self.__params:
                    self.__params['temp'] = 298.15
        # Initialize basic data
        if load_data:
            self.load_data(ylabel)
        # -- Set weights for multiple files
        self.__weights = None
        self.__calc_weight = None
        self.set_weights(weigh_dbsets)
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
        if self.__spec not in SPEC2DATA:
            raise KeyError('Unsupported spectroscopy: '+self.__spec)
        if self.__theory not in SPEC2DATA[self.__spec]:
            msg = 'Level of theory not available for spectroscopy.'
            raise KeyError(msg)
        if ylabel is None:
            ylabel = self.__ylab
        for dfile in self.__dfiles:
            if dfile.version[0] == 'CSV':
                qkeys = {
                    'spc': QLabel(quantity='anyspc', descriptor='Spec'),
                    'par': QLabel(quantity='anyspc', descriptor='SpcPar')
                }
            else:
                qkeys = SPEC2DATA[self.__spec][self.__theory]
                if not qkeys:
                    raise NotImplementedError('Keywords not available')
            dobjs = dfile.get_data(**qkeys)
            self.__xaxes.append([None, None])
            self.__yaxes.append([None, None])
            self.__xlabels.append([None, None])
            self.__ylabels.append([None, None])
            self.__xunits.append([None, None])
            self.__yunits.append([None, None])
            self.__spec_infos.append(None)
            self.__ytags.append(None)
            self.__extras.append(None)
            self.__broad_infos.append([{'func': None, 'hwhm': None},
                                       {'func': None, 'hwhm': None}])
            if 'spc' in qkeys:
                (self.__xaxes[-1][0], self.__xlabels[-1][0],
                 self.__xunits[-1][0], self.__yaxes[-1][0],
                 self.__ylabels[-1][0], self.__yunits[-1][0],
                 self.__broad_infos[-1][0], self.__ytags[-1],
                 self.__spec_infos[-1]
                 ) = _get_data_spectrum(dobjs, ylabel)
            elif 'freq' in qkeys:
                (self.__xaxes[-1][0], self.__xlabels[-1][0],
                 self.__xunits[-1][0], self.__yaxes[-1][0],
                 self.__ylabels[-1][0], self.__yunits[-1][0],
                 self.__broad_infos[-1][0], self.__extras[-1]
                 ) = _get_data_v_trans(dobjs, self.__spec, self.__params)
            elif 'ener' in qkeys:
                (self.__xaxes[-1][0], self.__xlabels[-1][0],
                 self.__xunits[-1][0], self.__yaxes[-1][0],
                 self.__ylabels[-1][0], self.__yunits[-1][0],
                 self.__broad_infos[-1][0]
                 ) = _get_data_e_trans(dobjs, self.__spec)
            else:
                raise NotImplementedError()
            if self.__label is None:
                self.__label = self.__yunits[-1][0]
            # Check if necessary to apply weighs to vib. transitions
            if (self.__weigh_vtrans is not None
                    and self.__extras[-1] is not None):
                for i, trans in enumerate(self.__extras[-1]['assign']):
                    if self.__weigh_vtrans == 'pni_dn':
                        dtrans = get_vib_trans(trans)
                    elif self.__weigh_vtrans == 'pni_d1':
                        dtrans = {i+1: 1}
                    else:
                        raise NotImplementedError(
                            'Unsupported weigh_vtrans model')
                    scale = weigh_trans_progress(
                        dtrans, self.__xaxes[-1][0], self.__params['temp'], -1)
                    self.__yaxes[-1][0][i] *= scale
        if self.__num_dfiles == 1:
            self.__xaxis = self.__xaxes[0]
            self.__yaxis = self.__yaxes[0]
            self.__xlabel = self.__xlabels[0]
            self.__ylabel = self.__ylabels[0]
            self.__xunit = self.__xunits[0]
            self.__yunit = self.__yunits[0]
            self.__broad_info = self.__broad_infos[0]

        self.__loaded = True

    def get_yaxes(self, dataset: int = 0
                  ) -> tp.Union[tp.Dict[str, str], tp.List[tp.Dict[str, str]]]:
        """Return all available Y axes in data file(s).

        If multiple datasets are stored, dataset can be provided to
        choose the dataset to use, starting at 1.
        If `dataset=0`, all yaxes are provides, as a list of tags.

        Parameters
        ----------
        dataset
            Dataset index (starting at 1).

        Returns
        -------
        dict, list
            Dictionary of the Y axes or list of dictionaries.

        Raises
        ------
        ArgumentError
            Incorrect dataset.
        """
        if not self.__loaded:
            self.load_data()
        if self.__num_dfiles == 1:
            return self.__ytags[0]
        else:
            if dataset == 0:
                return self.__ytags
            elif 1 < dataset <= len(self.__num_dfiles):
                return self.__ytags[dataset-1]
            else:
                raise ArgumentError('dataset',
                                    f'Dataset {dataset} not available.')

    def reset(self):
        """Reset default axes to original axes.

        Resets the pointers/aliases to the original data.
        Note that this is only possible if the data have not been
        overwritten.
        """
        self.__idref = 0
        for i in range(self.__num_dfiles):
            self.__xaxes[i][1] = None
            self.__yaxes[i][1] = None
            self.__xlabels[i][1] = None
            self.__ylabels[i][1] = None
        if self.__num_dfiles > 1:
            self.__xaxis = [None, None]
            self.__yaxis = [None, None]
            self.__xlabel = [None, None]
            self.__ylabel = [None, None]
            self.__xunit = [None, None]
            self.__yunit = [None, None]
            self.__broad_info = [None, None]
            self.__spec_info = None

    def overwrite_axis(self, data: tp.Sequence[tp.Union[float, int]],
                       axis: str = 'x', dataset: int = 1):
        """Overwrite one of the original axes.

        Note that this should only be used in special cases.
        If multiple datasets are stored, only one dataset is overridden.

        Parameters
        ----------
        data
            New axis data.
        axis
            Axis to overwrite.
        dataset
            Dataset index (starting at 1).

        Raises
        ------
        ArgumentError
            Incorrect value for dataset.
        IndexError
            Inconsistency in length size between the axes.
        TypeError
            Cannot convert to list.
        """
        if 0 < dataset <= self.__num_dfiles:
            if axis.lower() == 'x':
                self.__xaxes[dataset-1][0] = list(data)
            elif axis.lower() == 'y':
                self.__yaxes[dataset-1][0] = list(data)
            if (len(self.__xaxes[dataset-1][0])
                    != len(self.__yaxes[dataset-1][0])):
                raise IndexError('Inconsistency in the size of the axes')
        self.reset()

    def get_xaxis(self, origin: bool = False,
                  dataset: int = 0) -> tp.List[float]:
        """Show the X axis.

        Returns the X axis.
        For multiple datasets, the default is to provide the X axis from
        their combination.  If `dataset` is provided with a positive
        value, the X axis for the specific dataset is returned.

        Parameter
        ---------
        origin
            Use original X axis instead of current one.
        dataset
            Dataset index (starting at 1).  0 for the combined value.
        """
        if self.__num_dfiles == 1:
            if not origin:
                return self.__xaxis[self.__idref]
            else:
                return self.__xaxis[0]
        else:
            iaxis = 0 if origin else self.__idref
            if dataset == 0:
                if self.__xaxis[iaxis] is None:
                    self.__combine_datasets(iaxis)
                return self.__xaxis[iaxis]
            elif 0 < dataset <= self.__num_dfiles:
                return self.__xaxes[dataset-1][iaxis]
            else:
                raise ArgumentError('dataset',
                                    f'Dataset {dataset} not available.')

    xaxis = property(get_xaxis)

    def get_yaxis(self, origin: bool = False,
                  dataset: int = 0) -> tp.List[float]:
        """Show the Y axis.

        Returns the Y axis.
        For multiple datasets, the default is to provide the Y axis from
        their combination.  If `dataset` is provided with a positive
        value, the Y axis for the specific dataset is returned.

        Parameters
        ----------
        origin
            Use original Y axis instead of current one.
        dataset
            Dataset index (starting at 1).  0 for the combined value.

        Raises
        ------
        KeyError
            Unrecognized label.
        """
        if self.__num_dfiles == 1:
            if not origin:
                return self.__yaxis[self.__idref]
            else:
                return self.__yaxis[0]
        else:
            iaxis = 0 if origin else self.__idref
            if dataset == 0:
                if self.__yaxis[iaxis] is None:
                    self.__combine_datasets(iaxis)
                return self.__yaxis[iaxis]
            elif 0 < dataset <= self.__num_dfiles:
                return self.__yaxes[dataset-1][iaxis]
            else:
                raise ArgumentError('dataset',
                                    f'Dataset {dataset} not available.')

    yaxis = property(get_yaxis)

    def get_xunit(self,
                  tex_format: bool = True,
                  only_dot: bool = True,
                  use_cdot: bool = True,
                  origin: bool = False,
                  dataset: int = 0) -> str:
        """Get the X unit.

        Returns the X unit.
        For multiple datasets, the default is to provide the X unit from
        their combination.  If `dataset` is provided with a positive
        value, the X unit for the specific dataset is returned.

        Parameter
        ---------
        tex_format
            Use LaTeX/TeX formatting.
        only_dot
            If True, slashes are converted to dot and the exponent corrected.
            NOTE: only used for `tex_format = True`.
        use_cdot
            If True, cdot is used as "multiplier", otherwise the simple dot.
            NOTE: only used for `tex_format = True`.
        origin
            Use unit for original X axis instead of current one.
        dataset
            Dataset index (starting at 1).  0 for the combined spectrum.
        """
        iaxis = self.__idref if not origin else 0
        if self.__num_dfiles > 1:
            if self.__xaxis[iaxis] is None:
                self.__combine_datasets(iaxis)
        if self.__num_dfiles == 1 or dataset == 0:
            unit = self.__xunit
            label = self.__xlabel
        elif 0 < dataset <= self.__num_dfiles:
            unit = self.__xunits[dataset-1]
            label = self.__xlabel[dataset-1]
        else:
            raise ArgumentError('dataset',
                                f'Dataset {dataset} not available.')
        if unit[iaxis] is not None:
            if tex_format:
                unit = unit_to_tex(unit[iaxis].split(':')[-1], only_dot,
                                   use_cdot)
            else:
                unit = unit[iaxis].split(':')[-1]
            return f'{label[iaxis]} / {unit}'
        else:
            return None
    xunit = property(get_xunit)

    def get_yunit(self,
                  tex_format: bool = True,
                  only_dot: bool = True,
                  use_cdot: bool = True,
                  origin: bool = False,
                  dataset: int = 0) -> str:
        """Get the Y unit.

        Returns the Y unit.
        For multiple datasets, the default is to provide the Y unit from
        their combination.  If `dataset` is provided with a positive
        value, the Y unit for the specific dataset is returned.

        Parameter
        ---------
        tex_format
            Use LaTeX/TeX formatting.
        only_dot
            If True, slashes are converted to dot and the exponent corrected.
            NOTE: only used for `tex_format = True`.
        use_cdot
            If True, cdot is used as "multiplier", otherwise the simple dot.
            NOTE: only used for `tex_format = True`.
        origin
            Use unit for original Y axis instead of current one.
        dataset
            Dataset index (starting at 1).  0 for the combined spectrum.
        """
        iaxis = self.__idref if not origin else 0
        if self.__num_dfiles > 1:
            if self.__yaxis[iaxis] is None:
                self.__combine_datasets(iaxis)
        if self.__num_dfiles == 1 or dataset == 0:
            unit = self.__yunit
            label = self.__ylabel
        elif 0 < dataset <= self.__num_dfiles:
            unit = self.__yunits[dataset-1]
            label = self.__ylabels[dataset-1]
        else:
            raise ArgumentError('dataset',
                                f'Dataset {dataset} not available.')
        if unit[iaxis] is not None:
            if tex_format:
                unit = unit_to_tex(unit[iaxis].split(':')[-1], only_dot,
                                   use_cdot)
            else:
                unit = unit[iaxis].split(':')[-1]
            return f'{label[iaxis]} / {unit}'
        else:
            return None
    yunit = property(get_yunit)

    @property
    def label(self) -> str:
        """Get or set the label for the spectrum."""
        return self.__label

    @label.setter
    def label(self, val: str):
        self.__label = val

    def get_info(self, dataset: int = 0) -> tp.Optional[str]:
        """Get spectral information.

        Prints spectral information if available.

        Parameters
        ----------
        dataset
            Dataset index (starting at 1).  0 for the combined spectrum.

        Returns
        -------
        str
            Spectral information.
        """
        if dataset == 0 or self.__num_dfiles == 1:
            res = self.__spec_info
        elif 0 < dataset <= self.__num_dfiles:
            res = self.__spec_infos[dataset-1]
        else:
            raise KeyError('Unrecognized broadening parameter.')
        return res

    def get_broadening(self, info: tp.Optional[str] = None,
                       origin: bool = False,
                       dataset: int = 0
                       ) -> tp.Union[tp.Dict[str, tp.Union[str, float]], str,
                                     float]:
        """Get the broadening data.

        Returns the broadening information.
        For multiple datasets, the default is to provide the information
        for their combination.  If `dataset` is provided with a positive
        value, the broadening information for the specific dataset is
        returned.

        Parameters
        ----------
        info
            Information to return: func, hwhm
            If `None`, return all datas as a dictionary.
        origin
            Use original broadening data instead of current one.
        dataset
            Dataset index (starting at 1).  0 for the combined spectrum.

        Raises
        ------
        KeyError
            Unrecognized label.
        """
        iaxis = self.__idref if not origin else 0
        if self.__num_dfiles > 1:
            if self.__xaxis[iaxis] is None:
                self.__combine_datasets(iaxis)
        if self.__num_dfiles == 1 or dataset == 0:
            broad_info = self.__broad_info
        elif 0 < dataset <= self.__num_dfiles:
            broad_info = self.__broad_infos[dataset-1]
        else:
            raise ArgumentError('dataset',
                                f'Dataset {dataset} not available.')
        if info is None:
            res = broad_info[iaxis]
        else:
            if info.lower() in ('f', 'func', 'function'):
                res = broad_info[iaxis]['func']
            elif info.lower() in ('hw', 'hwhm'):
                res = broad_info[iaxis]['hwhm']
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
        """Set the broadening parameters and broaden the bands.

        Sets the broadening parameters to apply the broadening.
        For multiple datasets, the broadening is updated for each
        dataset.

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
        iaxis = 0 if origin else 1
        if func is None:
            bfun = None
        elif func.lower() in ('g', 'gau', 'gaussian'):
            bfun = 'gaussian'
        elif func.lower() in ('l', 'lor', 'lorentzian'):
            bfun = 'lorentzian'
        else:
            raise ValueError('Broadening function not supported.')
        if bfun is not None:
            for ifile in range(self.__num_dfiles):
                self.__broad_infos[ifile][iaxis]['hwhm'] = bfun
        else:
            if self.__num_dfiles == 1:
                bfun = self.__broad_infos[0][iaxis]['func']
            else:
                funcs = filter(lambda x: x[iaxis]['func'] is not None,
                               self.__broad_infos)
                if any(funcs):
                    bfun = funcs[0]
                    for fun in funcs[1:]:
                        if fun != bfun:
                            raise ValueError('Inconsistent functions')
        # we need to set consistently xmin and xmax for the broadening, so we
        # handle first bhw
        if hwhm is None:
            if self.__num_dfiles == 1:
                bhw = self.__broad_infos[0][iaxis]['hwhm']
            else:
                hwhms = filter(lambda x: x[iaxis]['hwhm'] is not None,
                               self.__broad_infos)
                if any(hwhms):
                    bhw = hwhms[0]
                    for hw in hwhms[1:]:
                        # WARNING: test is weak but precision should low.
                        if hw != bhw:
                            raise ValueError('Inconsistent HWHMs')
                else:
                    bhw = None
        else:
            bhw = float(hwhm)
            for ifile in range(self.__num_dfiles):
                self.__broad_infos[ifile][iaxis]['hwhm'] = bhw
        # Now the basic broadening parameters are set.
        # Check if we can proceed:
        if origin:
            return
        if bhw is None:
            raise ValueError('HWHM not set for the broadening.')
        # Reference axis switched to 1
        self.__idref = 1
        # Set X axis
        # -- set the lower bound, computing the minimum if not provided.
        if xmin is not None:
            x_min = xmin
        else:
            x_min = min((axis[iaxis] for axis in self.__xaxes)) - 10*bhw
        # -- do the same for the maximum
        if xmax is not None:
            x_max = xmax
        else:
            x_max = max((axis[iaxis] for axis in self.__xaxes)) + 10*bhw
        # -- check if inconsistent bounds
        if x_max < x_min:
            x_min, x_max = x_max, x_min
        # -- Set discretization points.
        if xres <= 0.0:
            raise ValueError('Wrong value for `xres`.')
        npoints = int(ceil((x_max-x_min)/xres))
        # -- Build reference X axis
        x_axis = [x_min + i*xres for i in range(npoints)]
        # Set parameters for Y
        # -- Set unit
        if yunit.lower() in ('n', 'norm', 'normalized'):
            y_unit_dest = SPEC2DATA[self.__spec]['unit']
            y_has_unit = True
        else:
            y_has_unit = False
            if yunit.lower() == 'default':
                y_unit_dest = SPEC2DATA[self.__spec]['unit']
            else:
                y_unit_dest = yunit
        # -- Set parameters
        y_params = {}
        if 'incfrq' in self.__params:
            y_params['incfrq'] = float(self.__params['incfrq'])
        # Now build the axes
        for ifile in range(self.__num_dfiles):
            self.__xaxes[ifile][iaxis] = x_axis[:]
            y_fac, x_fun, y_corr = convert_y(self.__spec, y_unit_dest,
                                             self.__yunits[ifile][0],
                                             **y_params)
            self.__yaxes[ifile][iaxis] = broaden(
                self.__xaxes[ifile][0], self.__yaxes[ifile][0],
                self.__xaxes[ifile][iaxis], bfun, bhw, y_fac, x_fun,
                y_corr, y_has_unit, False)
            self.__xunits[ifile][self.__idref] = self.__xunits[ifile][0]
            self.__xlabels[ifile][self.__idref] = self.__xlabels[ifile][0]
            self.__yunits[ifile][self.__idref] = y_unit_dest
            self.__ylabels[ifile][self.__idref] = self.__ylabels[ifile][0]

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

    def set_weights(self,
                    weights: tp.Optional[
                        tp.Union[str,
                                 tp.Sequence[
                                     tp.Union[str, float]]]] = None):
        """Set weights for multiple spectra.

        Sets weights to apply to each dataset.
        If only one data set is provided, the weight is set to 1.

        Parameters
        ----------
        weights
            Weights to be applied to each data set included in spectrum.
            See constructor for details.
        """
        self.__calc_weight = False
        if self.__num_dfiles == 1:
            self.__weights = 1.0
        else:
            if weights is None:
                self.__weights = \
                    [1.0/self.__num_dfiles for _ in range(self.__num_dfiles)]
            elif isinstance(weights, str):
                weight = weights.lower()
                if weight in ('bz', 'bz_e', 'bz_h', 'bz_a'):
                    self.__weights = ['bz_e' if weight == 'bz' else weight
                                      for _ in range(self.__num_dfiles)]
                else:
                    raise ArgumentError(
                        'weights', f'Unrecognized type of weights: {weights}')
                self.__calc_weight = 'bz_e' if weight == 'bz' else weight
            elif isinstance(weights, (list, tuple)):
                if len(weights) != self.__num_dfiles:
                    raise ArgumentError(
                        'weights',
                        'Number of weights do not match number of data files.'
                    )
                calc_w = None
                self.__weights = []
                for weight in weights:
                    if isinstance(weight, str):
                        if (item := weight.lower()) in ('bz', 'bz_e', 'bz_h',
                                                        'bz_a'):
                            if item == 'bz':
                                item = 'bz_e'
                            if calc_w is None:
                                calc_w = item
                                self.__calc_weight = calc_w
                            elif item != calc_w:
                                raise ArgumentError(
                                    'weights',
                                    'Inconsistent type of weights: only 1 '
                                    + 'type of computed weight allowed')
                            self.__weights.append(item)
                        else:
                            raise ArgumentError(
                                'weights',
                                f'Unrecognized type of weights: {weights}')
                    elif isinstance(weight, (float, int)):
                        self.__weights.append(float(weight))
                    else:
                        raise ArgumentError(
                            'weights',
                            f'Unrecognized weight parameter: {weight}')

        if self.__calc_weight:
            raise NotImplementedError('Boltzmann abundances NYI.')

        # To ensure that the calculations are correct, we need to reset
        # all computed elements.
        self.reset()

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

    def __combine_datasets(self, iaxis: int = 0):
        """Combine datasets into main axes.

        Combines multiple datasets and generate the information for
        the combined result.

        Parameters
        ----------
        iaxis
            Axis index (0 for origin).
        """
        if self.__num_dfiles > 1:
            self.__xaxis[iaxis], self.__yaxis[iaxis] = _combine_axes(
                self.__xaxes, self.__yaxes, self.__weights,
                self.__broad_infos, iaxis)
            (self.__broad_info[iaxis], self.__xunit[iaxis],
             self.__yunit[iaxis], self.__xlabel[iaxis], self.__ylabel[iaxis],
             self.__spec_info) = _set_combined_params(
                 self.__broad_infos, self.__spec_infos, self.__xunits,
                 self.__yunits, self.__xlabels, self.__ylabels)


def _get_data_v_trans(dobj: QData,
                      spec: str,
                      params: tp.Dict[str, tp.Any]
                      ) -> tp.Tuple[tp.List[float], str, str,
                                    tp.List[float], str, str,
                                    tp.Dict[str, tp.Any],
                                    tp.Dict[str, tp.Any]]:
    """Get data related to vibrational transitions from data object.

    Given a data object, extract information from vibrational transitions.

    Parameters
    ----------
    dobj
        Data object.
    spec
        Type of spectroscopy.
    params
        Spectroscopic parameters.
        May be modified on return to complete missing information.

    Returns
    -------
    list
        xaxis: List of X values.
    str
        xlabel: Label of X data.
    str
        xunit: Unit of X data.
    list
        yaxis: List of Y values.
    str
        ylabel: Label of Y data.
    str
        yunit: Unit of Y data.
    dict
        broadening: Broadening parameters.
    dict
        extra information, e.g., assignment data
    """
    modes = []
    for key in dobj['freq'].data:
        if isinstance(key, int):
            modes.append(key)
    modes.sort()
    xaxis = []
    yaxis = []
    extra = {'assign': []}
    if spec in ('RS', 'ROA'):
        incfrq = params.get('incfrq')
        if incfrq is None:
            for key in dobj['int'].data:
                try:
                    _ = float(key)
                    incfrq = key
                    break
                except ValueError:
                    continue
            if incfrq is None:
                raise KeyError('No incident frequency data found.')
            else:
                params['incfrq'] = incfrq
        else:
            if '.' not in incfrq:
                # truncate the incident frequencies stored in data.
                incfrqs = {item.split('.')[0]: item
                           for item in dobj['int'].data
                           if item.replace('.', '', 1).isdigit()}
            else:
                incfrqs = {item: item
                           for item in dobj['int'].data
                           if item.replace('.', '', 1).isdigit()}
            if incfrq not in incfrqs:
                vals = [item for item in dobj['int'].data
                        if item.replace('.', '', 1).isdigit()]
                raise KeyError(
                    f'''No incident frequency matches the value {incfrq}
Available: {', '.join(vals)}''')
            incfrq = incfrqs[incfrq]
            params['incfrq'] = incfrq
        setup = params.get('setup')
        if setup is None:
            for key in dobj['int'].data[incfrq].keys():
                if key in ('SCP(180)', 'SCP(180)u'):
                    setup = key
                    break
            if setup is None:
                setup = dobj['int'].data[incfrq].keys()[0]
            params['setup'] = setup
        else:
            if setup not in dobj['int'].data[incfrq].keys():
                raise KeyError(f'''Setup "{setup}" not found
Available: {dobj['int'].data[incfrq].keys()}''')
        ydata = dobj['int'].data[incfrq][setup]
    else:
        ydata = dobj['int'].data
    for mode in modes:
        try:
            if (isinstance(dobj['freq'].data[mode], float) and
                    isinstance(ydata[mode], float)):
                xaxis.append(dobj['freq'].data[mode])
                yaxis.append(ydata[mode])
                extra['assign'].append(dobj['assign'].data[mode])
        except KeyError as err:
            raise IndexError(f'''\
Inconsistency in the list of states between energies and intensities.
State {mode} was not found in one of the lists.''') from err
    xlabel = 'Wavenumbers'
    xunit = dobj['freq'].unit
    yunit = dobj['int'].unit
    ylabel = SPEC2DATA[spec][yunit.split(':')[0]]
    broadening = {'func': 'stick', 'hwhm': None}

    return xaxis, xlabel, xunit, yaxis, ylabel, yunit, broadening, extra


def _get_data_e_trans(dobj: QData,
                      spec: str
                      ) -> tp.Tuple[tp.List[float], str, str,
                                    tp.List[float], str, str,
                                    tp.Dict[str, tp.Any]]:
    """Get data related to electronic transitions from data object.

    Given a data object, extract information from electronic transitions.

    Parameters
    ----------
    dobj
        Data object.
    spec
        Type of spectroscopy.

    Returns
    -------
    list
        xaxis: List of X values.
    str
        xlabel: Label of X data.
    str
        xunit: Unit of X data.
    list
        yaxis: List of Y values.
    str
        ylabel: Label of Y data.
    str
        yunit: Unit of Y data.
    dict
        broadening: Broadening parameters.
    """
    xaxis = []
    yaxis = []
    for i, ener in enumerate(dobj['ener'].data):
        xaxis.append(ener)
        yaxis.append(dobj['int'].data[i])
    xlabel = 'Energy'
    xunit = dobj['ener'].unit
    yunit = dobj['int'].unit
    ylabel = SPEC2DATA[spec][yunit.split(':')]
    broadening = {'func': 'stick', 'hwhm': None}

    return xaxis, xlabel, xunit, yaxis, ylabel, yunit, broadening


def _get_data_spectrum(dobj: QData,
                       ylabel: tp.Optional[str] = None
                       ) -> tp.Tuple[tp.List[float], str, str,
                                     tp.List[float], str, str,
                                     tp.Dict[str, tp.Any],
                                     tp.Dict[str, tp.Any],
                                     tp.Optional[tp.Dict[str, tp.Any]]]:
    """Get data related to general spectral calc. from data object.

    Given a data object, extract information from spectral calculations,
    like vibronic calculations, where single transitions cannot be store
    explicitly.

    Parameters
    ----------
    dobj
        Data object.
    ylabel
        If multiple spectra may be present, label for the Y axis.
        Ex: TD OP spectra w/o temperature, Raman scattering setups.

    Returns
    -------
    list
        xaxis: List of X values.
    str
        xlabel: Label of X data.
    str
        xunit: Unit of X data.
    list
        yaxis: List of Y values.
    str
        ylabel: Label of Y data.
    str
        yunit: Unit of Y data.
    dict
        broadening: Broadening parameters.
    dict
        ytags: Labels of Y axes.
    dict
        info: Information on spectral simulation.
    """
    def dobj_to_ydat(dobj: QData, ykey: str,
                     yidx: tp.Optional[int] = None) -> tp.Any:
        if yidx is not None:
            return dobj.data.get(ykey)[yidx]
        else:
            return dobj.data.get(ykey)

    if isinstance(dobj['par'].data['x'], list):
        # Processing if multiple blocks:
        # TODO: Add check that parameters and spectra are consistent
        #       on Y as well
        nblocks = len(dobj['par'].data['x'])
        if len(dobj['spc'].data['x']) != nblocks:
            raise KeyError(
                'Inconsistency between spectral ranges and parameters.')
        # X axis should always be the same
        # TODO: Add proper check that X data are consistent.
        xaxis = dobj['spc'].data['x'][0]
        xlabel = dobj['par'].data['x'][0]
        # For Y, generate key information to retrieve data
        ykeys = (key for key in dobj['spc'].data if key.startswith('y'))
        nidx = len(ykeys[0]) - 1
        yfmt = f'y{{:0{nidx:d}d}}'
        offset = 0
        yindexes = {}
        for bloc in range(nblocks):
            ylabels = [item for item in ykeys
                       if (item[0] == 'y' and
                           len(dobj['par'].data[item]) > bloc)]
            for yax in ylabels:
                i = int(yax[1:])
                y = yfmt.format(i+offset)
                yindexes[y] = {'ykey': yax, 'yidx': bloc}
            offset += len(ylabels)
    else:
        xaxis = dobj['spc'].data['x']
        xlabel = dobj['par'].data['x']
        yindexes = {key: {'ykey': key} for key in dobj['spc'].data
                    if key.startswith('y')}
    xunit = dobj['par'].get('unitx')
    ytags = {key: dobj_to_ydat(dobj['par'], **val)
             for key, val in yindexes.items()}
    if ylabel is None:
        ydat = yindexes[list(yindexes.keys())[0]]
    else:
        try:
            ydat = yindexes[ylabel.lower()]
        except KeyError as err:
            raise KeyError('Non-existent Y label') from err
    yaxis = dobj_to_ydat(dobj['spc'], **ydat)
    ylabel = dobj_to_ydat(dobj['par'], **ydat)
    yunit = dobj['par'].get('unity')
    broadening = {'func': dobj['par'].get('func'),
                  'hwhm': dobj['par'].get('hwhm')}
    if 'info' in dobj:
        info = dobj['info'].extra_fields()
    else:
        info = None

    return xaxis, xlabel, xunit, yaxis, ylabel, yunit, broadening, ytags, info


def _combine_axes(xaxes: tp.Sequence[tp.Sequence[float]],
                  yaxes: tp.Sequence[tp.Sequence[float]],
                  weights: tp.Sequence[float],
                  broad_infos: tp.Sequence[tp.Sequence[tp.Dict[str, tp.Any]]],
                  axis: int = 0
                  ) -> tp.Tuple[tp.List[float], tp.List[float]]:
    """Combine multiple axes and return the resulting X and Y axes.

    Takes a list of X and Y axes and return their combinations with the
    Y axes scaled by their respective weights.
    The X axes can be interpolated if necessary.

    Parameters
    ----------
    xaxes
        List of X axes, as `[dataset][axis_index]`.
    yaxes
        List of Y axes, as `[dataset][axis_index]`.
    weights
        Weights to apply to the Y values of each dataset.
    broad_infos
        Broadening information for each dataset, as `[dataset][axis_index]`.
    axis
        Index of the axis of interest.

    Returns
    -------
    list
        Final X axis.
    list
        Final Y axis.

    Notes
    -----
    * The 0-th axes are assumed to contain "raw" data.  Interpolation is
      activated if datasets are found already broadened.
    """
    axes = [[], []]
    n_xaxes = len(xaxes)
    n_yaxes = len(yaxes)
    n_weights = len(weights)
    n_b_info = len(broad_infos)
    if not (n_xaxes == n_yaxes == n_weights == n_b_info):
        raise ArgumentError('arguments', 'Arguments sizes do not match.')
    if axis > min(n_xaxes, n_yaxes, n_weights, n_b_info) or axis < 0:
        raise IndexError('Wrong axis to combine')
    num_already_broad = len(tuple(filter(lambda x: x[axis]['func'] != 'stick',
                                         broad_infos)))
    if axis == 0:
        if 0 < num_already_broad < n_b_info:
            raise NotImplementedError('Interpolation not yet implemented')
        elif num_already_broad == n_b_info:
            # 1st check: X axes have same length
            l_xaxis = len(xaxes[0][axis])
            if not all((len(xaxis[axis]) == l_xaxis for xaxis in xaxes)):
                raise IndexError('Inconsistency in the X axes')
            # 2nd check, we control that all X have the same boundaries
            # and the same dx
            xmin = [xaxis[axis][0] for xaxis in xaxes]
            xmax = [xaxis[axis][-1] for xaxis in xaxes]
            step = [xaxis[axis][1]-xaxis[axis][0] for xaxis in xaxes]
            tol = 1.0e-4
            if not (abs(max(xmin) - min(xmin)) < tol
                    and abs(max(xmax) - min(xmax)) < tol
                    and abs(max(step) - min(step)) < tol):
                raise NotImplementedError('Interpolation not yet implemented')
            # Now let us combine spectra
            axes[0] = xaxes[0][axis][:]
            axes[1] = [0.0 for _ in range(l_xaxis)]
            for i, yaxis in enumerate(yaxes):
                if len(yaxis[axis]) != l_xaxis:
                    raise IndexError(
                        f'Y axis from dataset n. {i+1} inconsistent in size.')
                for iy, y in enumerate(yaxis[axis]):
                    axes[1][iy] += weights[i]*y
        else:
            axes[0] = [val for xaxis in xaxes for val in xaxis[axis]]
            axes[1] = [weight*val
                       for yaxis, weight in zip(yaxes, weights)
                       for val in yaxis[axis]]
            axes.sort(key=lambda x: x[0])
    else:
        if num_already_broad == 0:
            axes[0] = [val for xaxis in xaxes for val in xaxis[axis]]
            axes[1] = [weight*val
                       for yaxis, weight in zip(yaxes, weights)
                       for val in yaxis[axis]]
            axes.sort(key=lambda x: x[0])
        elif num_already_broad == n_b_info:
            # We assume that the processed spectra are consistent if broadened
            # So we take the first X axis range as reference.
            # As a basic check, we check all dimensions are equal
            l_xaxis = len(xaxes[0][axis])
            if not all((len(xaxis[axis]) == l_xaxis for xaxis in xaxes)):
                raise IndexError('Inconsistency in the X axes')
            axes[0] = xaxes[0][axis][:]
            axes[1] = [0.0 for _ in range(l_xaxis)]
            for i, yaxis in enumerate(yaxes):
                if len(yaxis[axis]) != l_xaxis:
                    raise IndexError(
                        f'Y axis from dataset n. {i+1} inconsistent in size.')
                for iy, y in enumerate(yaxis[axis]):
                    axes[1][iy] += weights[i]*y
        else:
            raise ArgumentError(
                'axes',
                'Mixed broadened/unbroadened data not supported')

    return axes[0], axes[1]


def _set_combined_params(broad_infos: tp.Sequence[tp.Sequence[
                                          tp.Dict[str, tp.Any]]],
                         spec_infos: tp.Optional[tp.Dict[str, tp.Any]],
                         xunits: tp.Sequence[tp.Sequence[str]],
                         yunits: tp.Sequence[tp.Sequence[str]],
                         xlabels: tp.Sequence[tp.Sequence[str]],
                         ylabels: tp.Sequence[tp.Sequence[str]],
                         axis: int = 0
                         ) -> tp.Tuple[
                             tp.Sequence[tp.Dict[str, tp.Any]],
                             str, str, str, str,
                             tp.Optional[tp.Dict[str, tp.Any]]]:
    """Set parameters for the combined spectra.

    Set the broadening information, units and labels for the combined spectra.

    Parameters
    ----------
    broad_infos
        Broadening information for each dataset, as `[dataset][axis_index]`.
    spec_infos
        Spectral information for each dataset, as `[dataset]`.
    xunits
        Unit of X values for each dataset, as `[dataset][axis_index]`.
    yunits
        Unit of Y values for each dataset, as `[dataset][axis_index]`.
    xlabels
        Label of X values for each dataset, as `[dataset][axis_index]`.
    ylabels
        Label of Y values for each dataset, as `[dataset][axis_index]`.
    """
    n_b_info = len(broad_infos)
    n_xunits = len(xunits)
    n_yunits = len(yunits)
    n_xlabels = len(xlabels)
    n_ylabels = len(ylabels)
    if not (n_b_info == n_xunits == n_yunits == n_xlabels == n_ylabels):
        raise ArgumentError('arguments', 'Arguments sizes do not match.')
    # initialize
    broad_info = broad_infos[0][axis]
    xunit = xunits[0][axis]
    yunit = yunits[0][axis]
    xlabel = xlabels[0][axis]
    ylabel = ylabels[0][axis]
    # Now check that the datasets are consistent
    for item in broad_infos[1:]:
        if item != broad_info:
            broad_info['func'] = None
            broad_info['hwhm'] = None
            break
    for item in xunits[1:]:
        if item != xunit:
            xunit = None
            break
    for item in yunits[1:]:
        if item != yunit:
            yunit = None
            break
    for item in xlabels[1:]:
        if item != xlabel:
            xlabel = None
            break
    for item in ylabels[1:]:
        if item != ylabel:
            ylabel = None
            break
    # Spectrum information only available in some specific cases.
    # We handle this case separately
    if all((item is not None for item in spec_infos)):
        spec_info = spec_infos[0]
        for item in spec_infos[1:]:
            if item != spec_info:
                spec_info = None
                break
    else:
        spec_info = None

    return broad_info, xunit, yunit, xlabel, ylabel, spec_info
