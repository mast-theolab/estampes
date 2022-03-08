"""Toolbox for spectra

Module providing tools to manipulate data relative to spectra.

Methods
-------
broaden
    Broadens a stick spectra with a given broadening function.
convert_y
    Returns parameters for the broadening to convert Y units.
"""

from math import ceil, log, pi
import typing as tp

from estampes.base import ArgumentError
from estampes.tools.char import convert_expr
from estampes.tools.math import f_gauss, f_lorentz
from estampes.data.physics import PHYSFACT, PHYSCNST, phys_fact


# ==============
# Module Methods
# ==============

def broaden(xval: tp.Sequence[float],
            yval: tp.Sequence[float],
            xaxis: tp.Sequence[float],
            funcname: str,
            hwhm: float,
            yfactor: float = 1.0,
            xfunc: tp.Optional[tp.Callable[[float], float]] = None,
            ynorm: bool = False,
            truncate: bool = False) -> tp.List[float]:
    """Broadens a stick spectra with a given broadening function.

    Performs the broadening operation given a list of transitions stored in
      two separate arrays, `xval` and `yval`.
    Returns the broadened spectrum for each point in `xaxis`.
    `yaxis` is computed as:
    yax[i] = sum_j yf*yv[j]*xfunc(xax[i])*f(xax[i]-xv[j])
    * yax: `yaxis`
    * jf: `yfactor`
    * yv: `yval`
    * xax: `xaxis`
    * xpow: `xpower`
    * f: `function`
    * xv: `xval`

    Parameters
    ----------
    xval
        List of X values.
    yval
        List of Y values.
    xaxis
        X axis of the broadened spectrum.
    funcname
        Broadening function: gaussian, lorentzian.
    hwhm
        Half-width at half-maximum (in unit of X).
    yfactor
        Scaling factor for the Y axis.
    xfunc
        Transformation of X to be included into the final Y value.
    ynorm
        Normalize Y (in this case, `yfactor` is ignored).
    truncate
        Truncate the broadening (for speed up).
        Note: `xaxis` must be linear and regular (delta x constant).

    Returns
    -------
    list
        Y axis with the broadened bands.

    Raises
    ------
    IndexError
        `xaxis` is empty.
    """
    supported_func = {'gaussian': f_gauss, 'lorentzian': f_lorentz}
    # Transform broadening function to be elemental
    func = supported_func[funcname]
    npoints = len(xaxis)
    if npoints == 0:
        raise IndexError('New X axis is empty')
    yaxis = [0.0 for _ in range(npoints)]
    if not truncate:
        for i in range(npoints):
            for x0, y0 in zip(xval, yval):
                yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
    else:
        # hwtimes: Function is computer over an interval of 2*hwhm*hwtimes
        # 3.5*hwhm represents 99.7% of the Gaussian area
        # For lorentzian, the same coverage requires 212*hwhm, truncated to 15.
        hwtimes = {'gaussian': 3.5, 'lorentzian': 15}[func]
        xmin = xaxis[0]
        dx = abs(xaxis[-1] - xaxis[0])/(npoints-1)
        nhwpoints = int(ceil(hwhm*hwtimes/dx))
        if xmin > xaxis[-1]:
            # X axis is reversed
            xmin = xaxis[-1]
            for x0, y0 in zip(xval, yval):
                i0 = -int((x0 - xmin)/dx) - 1
                imin = max(i0 - nhwpoints, -npoints)
                imax = min(i0 + nhwpoints, -1)
                for i in range(imin, imax+1):
                    yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
        else:
            for x0, y0 in zip(xval, yval):
                i0 = int((x0 - xmin)/dx)
                imin = max(i0 - nhwpoints, 0)
                imax = min(i0 + nhwpoints, npoints-1)
                for i in range(imin, imax+1):
                    yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
    if xfunc is not None:
        for i in range(npoints):
            yaxis[i] *= xfunc(xaxis[i])
    if ynorm:
        ymax = max(yaxis, key=abs)
        for i in range(npoints):
            yaxis[i] /= ymax
    return yaxis


def convert_y(specabbr: str,
              unit_to: str,
              unit_from: str,
              **subopts: tp.Dict[str, tp.Any]
              ) -> tp.Tuple[float, tp.Optional[tp.Callable[[float], float]]]:
    """Returns parameters for the conversion of Y units.

    Returns the scaling factor and the power of X needed to convert the
      Y unit from `unit_from` to `unit_to`.
    The units are strings composed of two parts, separated by colon:
    - a description: DS, RS, RA, I, II...
    - a unit: km/mol, cgs...
    A default unit may be available, in which case, the unit can be
      provided as:
    `desc:`

    Parameters
    ----------
    specabbr
        Spectroscopy name (abbreviated version).
    unit_to
        Final Y unit.
    unit_from
        Original Y unit.
    subopts
        Sub-options, dependent on quantity:
        - 'incfreq': incident frequency for Raman, ROA (in cm^-1)
                     if not present or None, use def. 532 nm.

    Returns
    -------
    float
        Scaling factor on the Y unit.
    xfunc
        Function to transform X in calculation of Y value.

    Raises
    ------
    ArgumentError
        Error in sub-options.
    IndexError
        Unrecognized spectroscopy or unit.
    NotImplementedError
        Spectroscopy or unit not yet available.

    Notes
    -----
    * Recognized unit descriptions:
        * `I` - intensity
        * `II` - integrated intensity
        * `DS` - dipole strength
        * `RS` - rotatory strength
    * Units are case-sensitive to avoid ambiguities.
      "/" and "." must be indicated correctly.
        * `M-1` or `/M` -> dm^3/mol
        * `M-1.cm-1` or `/M/cm` -> dm^3/mol/cm
        * `M-1.cm-2` or `/M/cm2` -> dm^3/mol/cm^2
    """
    msgNYI = f'For {specabbr}: conversion from "{unit_from}" to "{unit_to}"' +\
        ' not yet implemented.'
    msgNA = f'Conversion from "{unit_from}" to "{unit_to}" not recognized ' +\
        f'for {specabbr}'
    msgVal = 'Wrong value for option {subkey} in spectroscopy {specabbr}'
    UNIT_EPS = {  # Epsilon: Molar absorption coefficient
        '/M/cm': ('M-1.cm-1', '/M/cm', 'dm3.mol-1.cm-1', 'dm3/mol/cm',
                   'L.mol-1.cm-1', 'L/mol/cm', '')
    }
    UNIT_II = {  # Integrated intensity
        '/M/cm2': ('M-1.cm-2', '/M/cm2', 'dm3.mol-1.cm-2', 'dm3/mol/cm2',
                   'L.mol-1.cm-2', 'L/mol/cm2')
    }
    UNIT_DS = {  # Dipole strength
        'esu2.cm2': ('statC2.cm2', 'esu2.cm2')
    }
    UNIT_RS = {  # Rotatory strength
        'esu2.cm2': ('statA2.cm2', 'esu2.cm2')
    }
    # Initial setup
    _spec = specabbr.upper()
    res = unit_to.split(':')
    dfact = 1.0
    sfact = 1.0
    if len(res) == 1:
        _dest_unit = ''
    else:
        dunit = res[1].split()
        if len(dunit) == 1:
            _dest_unit = res[1].strip().replace('^', '')
            dfact = 1.0
        else:
            _dest_unit = dunit[1].strip().replace('^', '')
            dfact = eval(convert_expr(dunit[0], None))
    _dest_type = res[0].upper()
    res = unit_from.split(':')
    if len(res) == 1:
        _src_unit = ''
    else:
        sunit = res[1].split()
        if len(sunit) == 1:
            _src_unit = res[1].strip().replace('^', '')
            sfact = 1.0
        else:
            _src_unit = sunit[1].strip().replace('^', '')
            sfact = eval(convert_expr(sunit[0], None))
    _src_type = res[0].upper()
    yfactor = 1.0
    xfunc = None
    # Test
    if _spec == 'IR':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'DS':
                    if _src_unit in UNIT_DS['esu2.cm2']:
                        yfactor = 8*pi**3*PHYSCNST.avogadro * 1.0e-7 / \
                            (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'VCD':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'RS':
                    if _src_unit in UNIT_RS['esu2.cm2']:
                        yfactor = 32*pi**3*PHYSCNST.avogadro * 1.0e-7 / \
                            (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'OPA':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'II':
                    if _src_unit in UNIT_II['/M/cm2']:
                        yfactor = 1.0
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    else:
        raise NotImplementedError(msgNYI)

    return yfactor*sfact/dfact, xfunc
