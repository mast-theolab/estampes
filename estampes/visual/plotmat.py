"""Module for matrix/vector plotting.

This module provides simple functions to build matrix/vector
  representations for Matplotlib.

Methods
-------
plot_cmat
    Plots a red-blue heatmap of a C-like matrix.
plot_jmat
    Plots a black-white heatmap of a J-like matrix.
plot_kvec
    Plots a horizontal bar chart to represent a K-like vector.
"""

import typing as tp

import numpy as np
import matplotlib as mpl

from estampes.base import ArgumentError


# ==============
# Module Methods
# ==============

def plot_jmat(mat: np.ndarray,
              canvas: mpl.axes.Axes,
              norm_mode: str = 'byrow') -> mpl.image.AxesImage:
    """Plots a black-white heatmap of a J-like matrix.

    Plots the squared elements of a matrix `mat` in the Matplotlib
      Axes `canvas`.
    Possible modes of normalization:
    none
        no normalization.
    byrow
        normalized by row.
    bycol
        normalized by col.
    highest
        sets highest squared element in the matrix to 1.

    Parameters
    ----------
    mat
        Matrix to plot.
    canvas
        Matplotlib plotting frame.
    norm_mode
        Normalization mode.

    Returns
    -------
    mpl.image.AxesImage
        Image of the matrix.

    Raises
    ------
    ArgumentError
        Error in arguments.
    """
    _mat = mat**2
    if norm_mode.lower() == 'byrow':
        norm = np.sum(_mat, axis=1)
        for i in range(np.shape(_mat)[0]):
            _mat[i, :] /= norm[i]
    elif norm_mode.lower() == 'bycol':
        norm = np.sum(_mat, axis=0)
        for i in range(np.shape(_mat)[1]):
            _mat[:, i] /= norm[i]
    elif norm_mode.lower() == 'highest':
        _mat /= _mat.max()
    elif norm_mode.lower() != 'none':

        raise ArgumentError('norm_mode')
    plot = canvas.matshow(_mat, cmap=mpl.cm.gray_r, vmin=0.0, vmax=1.0)
    canvas.set_xlabel('Final-state modes')
    canvas.set_ylabel('Initial-state modes')
    return plot


def plot_cmat(mat: np.ndarray,
              canvas: mpl.axes.Axes) -> tp.Union[float, mpl.image.AxesImage]:
    """Plots a red-blue heatmap of a C-like matrix.

    Plots a normalized matrix has a "hot-cold"-like matrix, normalizing
      the highest number in absolute value to 1 in Matplotlib Axes
      `canvas`.

    Parameters
    ----------
    mat
        Matrix to plot.
    canvas
        Matplotlib plotting frame.

    Returns
    -------
    float
        Normalization factor.
    mpl.image.AxesImage
        Image of the matrix.
    """
    _mat = mat.copy()
    norm = np.max(np.abs(_mat))
    _mat /= norm
    plot = canvas.matshow(_mat, cmap=mpl.cm.seismic_r, vmin=-1.0, vmax=1.0)
    canvas.set_xlabel('Final-state modes')
    canvas.set_ylabel('Final-state modes')

    return norm, plot


def plot_kvec(vec: np.ndarray,
              canvas: mpl.axes.Axes) -> mpl.container.BarContainer:
    """Plots a horizontal bar chart to represent a K-like vector.

    Plots a shift vector as a horizontal bars.

    Parameters
    ----------
    vec
        Vector to plot.
    canvas
        Matplotlib plotting frame.

    Returns
    -------
    mpl.container.BarContainer
        Image of the vector.
    """
    llo0 = '\N{SUBSCRIPT ZERO}'
    lloe = '\N{LATIN SUBSCRIPT SMALL LETTER E}'
    ldot = '\N{MIDDLE DOT}'
    xlab_kvec = f'Displacement / m{lloe}{ldot}a{llo0}'
    avec = np.abs(vec)
    num = len(vec)
    pos = np.arange(num) + 1
    plot = canvas.barh(pos, avec, align='center')
    for i in range(num):
        if vec[i] < 0:
            plot[i].set_color('red')
        else:
            plot[i].set_color('#73BFFF')
    canvas.set_xlabel(xlab_kvec)
    canvas.set_ylabel('Initial-state modes')
    canvas.set_ylim(bottom=0)

    return plot
