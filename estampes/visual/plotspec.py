"""Module for spectrum plotting.

This module provides simple function to display spectra.

Methods
-------
plot_spec_2D
"""

import typing as tp

import numpy as np
import matplotlib as mpl


# ==============
# Module Methods
# ==============

def plot_spec_2D(axes: tp.Dict[str, tp.Sequence[float]],
                 canvas: mpl.axes.Axes,
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
    """Plots 2D spectra.

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
    elif legends is not None and keys['x'] in legends:
        _xlab = legends[keys['x']].strip()
    else:
        _xlab = 'Energy / a.u.'
    if ylabel is not None:
        _ylab = ylabel.strip()
    elif legends is not None and keys['i'] in legends:
        _ylab = legends[keys['i']].strip()
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
