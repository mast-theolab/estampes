"""Module for spectrum plotting.

This module provides simple function to display spectra.
"""

import re
import typing as tp

import numpy as np
import matplotlib.axes as mpl_axes

from estampes.base import ArgumentError


# ==============
# Module Classes
# ==============

class SpecLayout(object):
    """Class for the spectrum layout.

    Class managing information for the spectrum layout.

    Parameters
    ----------
    xleft
        Leftmost bound (X axis).
    xright
        Rightmost bound (X axis).
    ytop
        Top bound (Y axis).
    ybottom
        Bottom bound (Y axis).
    xscale
        Scale for the X axis.
        Acceptable values: linear
    yscale
        Scale for the Y axis.
        Acceptable values: linear
    xlabel
        Label for the X axis.
    ylabel
        Label for the Y axis.
    title
        Title for the spectrum layout.
    legpos
        Legend position.
    legcol
        Number of columns in the legend.
    plottag
        Tag (panel label) for the plot.

    Notes
    -----
    * `xlabel` and `ylabel` are assumed unset if a `None` value is given
      to the constructor (default).  To force the value to be considered
      `None`, the labels should be set through the setters instead.
    """
    def __init__(self,
                 xleft: tp.Optional[tp.Union[int, float, str]] = None,
                 xright: tp.Optional[tp.Union[int, float, str]] = None,
                 ytop: tp.Optional[tp.Union[int, float, str]] = None,
                 ybottom: tp.Optional[tp.Union[int, float, str]] = None,
                 xscale: tp.Optional[str] = None,
                 yscale: tp.Optional[str] = None,
                 xlabel: tp.Optional[str] = None,
                 ylabel: tp.Optional[str] = None,
                 title: tp.Optional[str] = None,
                 legpos: tp.Optional[str] = None,
                 legcol: tp.Optional[tp.Union[int, str]] = None,
                 plottag: tp.Optional[str] = None):
        self.__legpos = False
        self.__legcol = False
        self.__plottag = False
        self.__ylabel_set = False
        self.xleft = xleft
        self.xright = xright
        self.ytop = ytop
        self.ybottom = ybottom
        self.xscale = xscale
        self.yscale = yscale
        self.xlabel = xlabel
        if xlabel is None:
            self.__xlabel_set = False
        self.ylabel = ylabel
        if ylabel is None:
            self.__ylabel_set = False
        self.title = title
        self.legend(pos=legpos, ncols=legcol)
        self.def_panel(plottag)

    @property
    def xleft(self) -> tp.Optional[float]:
        """Leftmost bound (X axis), `None` if unset."""
        return self.__xleft

    @xleft.setter
    def xleft(self, val: tp.Optional[tp.Any]):
        if val is None or val.lower() == 'auto':
            self.__xleft = None
        else:
            self.__xleft = float(val)

    @property
    def xright(self) -> tp.Optional[float]:
        """Rightmost bound (X axis), `None` if unset."""
        return self.__xright

    @xright.setter
    def xright(self, val: tp.Optional[tp.Any]):
        if val is None or val.lower() == 'auto':
            self.__xright = None
        else:
            self.__xright = float(val)

    @property
    def xmin(self) -> tp.Optional[float]:
        """Lower bound (X axis), `None` if no bounds set."""
        if self.__xleft is None and self.__xright is None:
            return None
        else:
            return min([item for item in [self.__xleft, self.__xright]
                        if item is not None])

    @property
    def xmax(self) -> tp.Optional[float]:
        """Upper bound (X axis), `None` if no bounds set."""
        if self.__xleft is None and self.__xright is None:
            return None
        else:
            return max([item for item in [self.__xleft, self.__xright]
                        if item is not None])

    @property
    def ytop(self) -> tp.Optional[float]:
        """Top bound (Y axis), `None` if unset."""
        return self.__ytop

    @ytop.setter
    def ytop(self, val: tp.Optional[tp.Any]):
        if val is None or val.lower() == 'auto':
            self.__ytop = None
        else:
            self.__ytop = float(val)

    @property
    def ybottom(self) -> tp.Optional[float]:
        """Bottom bound (Y axis), `None` if unset."""
        return self.__ybottom

    @ybottom.setter
    def ybottom(self, val: tp.Optional[tp.Any]):
        if val is None or val.lower() == 'auto':
            self.__ybottom = None
        else:
            self.__ybottom = float(val)

    @property
    def ymin(self) -> tp.Optional[float]:
        """Lower value of the Y axis, `None` if no bounds set."""
        if self.__ytop is None and self.__ybottom is None:
            return None
        else:
            return min([item for item in [self.__ytop, self.__ybottom]
                        if item is not None])

    @property
    def ymax(self) -> tp.Optional[float]:
        """Upper value of the Y axis, `None` if no bounds set."""
        if self.__ytop is None and self.__ybottom is None:
            return None
        else:
            return max([item for item in [self.__ytop, self.__ybottom]
                        if item is not None])

    @property
    def xlabel(self) -> tp.Optional[str]:
        """Label for the X axis."""
        return self.__xlabel

    @xlabel.setter
    def xlabel(self, val: tp.Optional[str]):
        self.__xlabel = val
        self.__xlabel_set = True

    @property
    def ylabel(self) -> tp.Optional[str]:
        """Label for the X axis."""
        return self.__ylabel

    @ylabel.setter
    def ylabel(self, val: tp.Optional[str]):
        self.__ylabel = val
        self.__ylabel_set = True

    @property
    def xscale(self) -> tp.Dict[str, tp.Union[str, int]]:
        """Dictionary with the information on the scale used on X axis."""
        return self.__xscale

    @xscale.setter
    def xscale(self, val: tp.Optional[str]):
        if val is None:
            self.__xscale = {'value': 'linear'}
        else:
            key = val.lower()
            if key in ('lin', 'linear'):
                self.__xscale = {'value': 'linear'}
            elif key in ('log', 'log10'):
                self.__xscale = {'value': 'log', 'base': 10}
            elif key in ('ln', 'log2'):
                self.__xscale = {'value': 'log', 'base': 2}
            else:
                raise IndexError(f'Unrecognized scale: {val}')

    @property
    def yscale(self) -> tp.Dict[str, tp.Union[str, int]]:
        """Dictionary with the information on the scale used on Y axis."""
        return self.__yscale

    @yscale.setter
    def yscale(self, val: tp.Optional[str]):
        if val is None:
            self.__yscale = {'value': 'linear'}
        else:
            key = val.lower()
            if key in ('lin', 'linear'):
                self.__yscale = {'value': 'linear'}
            elif key in ('log', 'log10'):
                self.__yscale = {'value': 'log', 'base': 10}
            elif key in ('ln', 'log2'):
                self.__yscale = {'value': 'log', 'base': 2}
            else:
                raise IndexError(f'Unrecognized scale: {val}')

    @property
    def title(self) -> str:
        """Title of the spectrum."""
        return self.__title

    @title.setter
    def title(self, val: tp.Optional[str]):
        self.__title = val

    def is_label_set(self,
                     what: str = 'y'
                     ) -> tp.Union[bool, tp.Tuple[bool, bool]]:
        """Return if an axis label is set or not.

        Returns if an axis label is set or not.
        Accepted values are: 'y', 'x', 'b'/'both'.
        """
        if what.lower() == 'y':
            return self.__ylabel_set
        elif what.lower() == 'x':
            return self.__xlabel_set
        elif what.lower() in ('b', 'both'):
            return self.__xlabel_set, self.__ylabel_set
        else:
            raise ArgumentError('Unknown label specification to set.')

    def set_xbounds(self, *xaxes: tp.List[float], desc: bool = False):
        """Sets X bounds from a list of X axes.

        Given a list of X values, sets the bounds.
        By default, the bounds are set so that the the lower bound is on
        the left.  This can be reversed with `desc` = True.

        Parameters
        ----------
        xaxes
            List of X axes.
        desc
            Sets X bounds so that the highest value is on the left.
        """
        xmin = []
        xmax = []
        for xaxis in xaxes:
            xmin.append(min(xaxis))
            xmax.append(max(xaxis))
        if not desc:
            self.__xleft = min(xmin)
            self.__xright = max(xmax)
        else:
            self.__xright = min(xmin)
            self.__xleft = max(xmax)

    def set_ybounds(self, *yaxes: tp.List[float], desc: bool = False):
        """Sets Y bounds from a list of Y axes.

        Given a list of Y values, sets the bounds.
        By default, the bounds are set so that the the lower bound is at
        the bottom.  This can be reversed with `desc` = True.

        Parameters
        ----------
        yaxes
            List of Y axes.
        desc
            Sets Y bounds so that the highest value is at the bottom.
        """
        ymin = []
        ymax = []
        for yaxis in yaxes:
            ymin.append(min(yaxis))
            ymax.append(max(yaxis))
        if not desc:
            self.__ybottom = min(ymin)
            self.__ytop = max(ymax)
        else:
            self.__ytop = min(ymin)
            self.__ybottom = max(ymax)

    def legend(self, pos: tp.Optional[tp.Union[str, int, bool]] = None,
               ncols: tp.Optional[int] = None):
        """Sets the position and options of the legend.

        Sets the position and some parameters of the legend.

        Parameters
        ----------
        pos
            Position of the legend.
        ncols
            Number of columns.

        Raises
        ------
        IndexError
            Unrecognized position label.
        ValueError
            Unsupported number of columns.
        """
        LEGPOS = ['best',
                  'upper right',
                  'upper left',
                  'lower left',
                  'lower right',
                  'right',
                  'center left',
                  'center right',
                  'lower center',
                  'upper center',
                  'center']
        # Position definition
        if pos is not None:
            if isinstance(pos, int):
                if pos == -1:
                    self.__legpos = False
                else:
                    self.__legpos = LEGPOS[pos]
            elif isinstance(pos, str):
                if pos.lower() == 'no':
                    self.__legpos = False
                elif pos.lower() in LEGPOS:
                    self.__legpos = pos.lower()
                else:
                    raise IndexError(f'Unrecognized position: {pos}')
            elif not pos:
                self.__legpos = False
            else:
                raise IndexError('Unrecognized position definition.')
        # Number of columns
        if ncols is not None:
            self.__legcol = int(ncols)
            if self.__legcol < 1:
                raise ValueError('Error in number of columns.')

    def set_plot(self, canvas: mpl_axes.Axes):
        """Sets the spectrum parameters on a Matplotlib canvas.

        Sets the parameters for the spectrum layout on a Matplotlib
        canvas.

        Parameters
        ----------
        canvas
            Matplotlib plotting frame.
        """
        # X axis
        pars = {}
        if self.__xleft is not None:
            pars['left'] = self.__xleft
        if self.__xright is not None:
            pars['right'] = self.__xright
        if pars:
            canvas.set_xlim(**pars)
        canvas.set_xscale(**self.__xscale)
        if self.__xlabel is not None:
            canvas.set_xlabel(self.__xlabel)
        # Y axis
        pars = {}
        if self.__ybottom is not None:
            pars['bottom'] = self.__ybottom
        if self.__ytop is not None:
            pars['top'] = self.__ytop
        if pars:
            canvas.set_ylim(**pars)
        canvas.set_yscale(**self.__yscale)
        if self.__ylabel is not None:
            canvas.set_ylabel(self.__ylabel)
        # Legend
        if self.__legpos:
            pars = {'loc': self.__legpos}
            if self.__legcol is not None:
                pars['ncol'] = self.__legcol
            canvas.legend(**pars)
        # Title
        if self.__title is not None:
            canvas.set_title(self.__title)
        # Panel tag
        if self.__plottag:
            canvas.text(**self.__plottag, transform=canvas.transAxes)

    def def_panel(self, desc: tp.Optional[str] = None, **kwargs):
        """Define the parameters for a panel label on a canvas.

        Sets the parameters necessary to include a Places a label `text`
        in one of the corner of a canvas.

        Parameters
        ----------
        text
            Text to print with position as `text @ position`.
        kwargs
            Matplotlib Text properties (advanced usage).
            This can be used to override any properties, including
            the text.

        Raises
        ------
        ArgumentError
            Position specification in `desc` not recognized.
        """
        if desc is not None:
            self.__plottag = {
                'verticalalignment': 'center',
                'fontweight': 'semibold',
                'fontsize': 'x-large'
            }
            res = desc.split('@')
            self.__plottag['s'] = res[0].strip()
            if len(res) == 1:
                pos = 'top left'
            else:
                pos = res[1].strip().lower()
            if 'right' in pos:
                self.__plottag['x'] = 0.95
                self.__plottag['horizontalalignment'] = 'right'
                pos = pos.replace('right', '')
            else:
                self.__plottag['x'] = 0.05
                self.__plottag['horizontalalignment'] = 'left'
                pos = pos.replace('left', '')
            if 'bottom' in pos or 'lower' in pos:
                self.__plottag['y'] = 0.05
                pos = pos.replace('bottom', '')
                pos = pos.replace('lower', '')
            else:
                self.__plottag['y'] = 0.95
                pos = pos.replace('top', '')
                pos = pos.replace('upper', '')
            if pos.strip():
                raise ArgumentError('position in "desc"')
            if kwargs:
                for key, arg in kwargs.items():
                    self.__plottag[key] = arg


# ==============
# Module Methods
# ==============

def format_label(label: str, full: bool = False) -> str:
    """Format label.

    Parameters
    ----------
    label
        Original label
    full
        Full label, as 'label / unit'
    """
    if full:
        text, unit = label.rsplit('/', maxsplit=1)
    else:
        text, unit = False, label
    unit = re.sub(r'\^?([-+]?\d+)', r'$^{\1}$', unit)
    unit = re.sub(r'/(\w+)', r'\1$^{-1}$', unit)
    if text:
        res = f'{text} / {unit}'
    else:
        res = unit
    return res


def plot_spec_2D(axes: tp.Dict[str, tp.Sequence[float]],
                 canvas: mpl_axes.Axes,
                 xlabel: tp.Optional[str] = None,
                 ylabel: tp.Optional[str] = None,
                 legends: tp.Optional[tp.Dict[str, str]] = None,
                 colors: tp.Optional[tp.Dict[str, str]] = None,
                 *,
                 xleft: tp.Optional[tp.Union[int, float]] = None,
                 xright: tp.Optional[tp.Union[int, float]] = None,
                 ydown: tp.Optional[tp.Union[int, float]] = None,
                 yup: tp.Optional[tp.Union[int, float]] = None,
                 is_stick: bool = False,
                 add_y0: bool = False) -> tp.Dict[str, float]:
    """Plot 2D spectra.

    Plots one or more 2D spectra in the provided Matplotlib Axes
    `canvas`.
    If `legends` and/or `colors are given, they must have the same keys
    as `axes` (for the y axis).
    `xlabel` provides the label for the X axis.  If unset, `legends[x]`
    is used.  If the key does not exist, a generic label is used.
    `ylabel` provides the label for the Y axis.  If unset, `legends[I]`
    is used.  If the key does not exist, a generic label is used.

    Parameters
    ----------
    axes
        Dictionary of lists of data.  It should contain:
        * one 'x' key
        * one or more 'y'/'yX' keys for the y axis/axes (X=1+ digits.)
    canvas
        Matplotlib plotting frame.
    xlabel
        Label for the X axis.
        If `None`, the function tries to generate it from `Legends`.
    ylabel
        Label for the Y axis.
        If `None`, the function tries to generate it from `Legends`.
    legends
        Legends to display for each spectrum.
        If `None`, the legend is not displayed.
    colors
        Colors as hex string (ex: '#aaaaaa').
        If `None`, the colors are automatically generated.
    xleft
        Left bound on the X axis.
    xright
        Right bound on the X axis.
    ydown
        Down/bottom bound on the Y axis.
    yup
        Up/top on the Y axis.
    is_stick
        If `True`, the spectrum/a must be printed as sticks.
    add_y0
        Add 'y=0' horizontal line.

    Raises
    ------
    IndexError
        Missing quantity
    """
    # Build equivalence table for case-insensitive search in `axes`
    # -------------------------------------------------------------
    # We assume that all dictionaries have consistent cases
    keys = {}
    for key in axes:
        keys[key.lower()] = key
    if 'x' not in keys:
        raise IndexError('Missing X axis')

    # Build axis labels
    # -----------------
    if xlabel is not None:
        _xlab = xlabel.strip()
    elif legends is not None and 'x' in legends:
        _xlab = legends['x'].strip()
    else:
        _xlab = 'Energy / a.u.'
    if ylabel is not None:
        _ylab = ylabel.strip()
    elif legends is not None and 'I' in legends:
        _ylab = legends['I'].strip()
    else:
        _ylab = 'Intensity / a.u.'

    # Add plots
    # ---------
    # Put y=0 line in the bottom for visibility
    if add_y0:
        canvas.axhline(y=0., c='.5')
    _xax = axes[keys['x']]
    i = 0
    for key in keys:
        data = {}
        if key[0] == 'y':
            if legends is not None and keys[key] in legends:
                data['label'] = legends[keys[key]]
            if colors is not None and keys[key] in colors:
                data['color'] = colors[keys[key]]
            else:
                data['color'] = 'C{:d}'.format(i)
                i += 1
            if is_stick:
                zeros = np.zeros(len(axes[keys[key]]))
                canvas.vlines(_xax, zeros, axes[keys[key]], **data)
            else:
                canvas.plot(_xax, axes[keys[key]], **data)

    # Set axes bounds
    # ---------------
    data = {}
    if xleft is not None:
        data['left'] = xleft
    if xright is not None:
        data['right'] = xright
    if data:
        canvas.set_xlim(**data)
    data = {}
    if ydown is not None:
        data['bottom'] = ydown
    if yup is not None:
        data['top'] = yup
    if data:
        canvas.set_ylim(**data)
    # Get actual axes bounds
    # ----------------------
    xleft, xright = canvas.get_xlim()
    ydown, yup = canvas.get_ylim()
    bounds = {'xleft': xleft, 'xright': xright, 'ydown': ydown, 'yup': yup}
    # Add legend
    # ----------
    if legends is not None:
        canvas.legend(loc='best')
    # Add labels
    # ----------
    canvas.set_xlabel(_xlab)
    canvas.set_ylabel(_ylab)

    return bounds
