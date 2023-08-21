"""Module for matrix/vector plotting.

This module provides simple functions to build matrix/vector
representations for Matplotlib.
"""

import typing as tp

import numpy as np
import matplotlib as mpl
import matplotlib.axes as mpl_axes
import matplotlib.image as mpl_img
import matplotlib.container as mpl_cont
from matplotlib.ticker import FuncFormatter

from estampes.base import ArgumentError


# ==============
# Module Methods
# ==============

def plot_jmat(mat: np.ndarray,
              canvas: mpl_axes.Axes,
              norm_mode: str = 'byrow',
              show_grid: bool = True,
              *,
              top_down: bool = False,
              invert_modes: bool = False
              ) -> mpl_img.AxesImage:
    """Plots a black-white heatmap of a J-like matrix.

    Plots the squared elements of a matrix `mat` in the Matplotlib
    Axes `canvas`.

    The possible modes of normalization are::
    
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
    show_grid
        Show light grid to hel find the modes.
    top_down
        Orders the modes so that the first one is at the top.
    invert_modes
        Invert ordering of normal modes (== decreasing order).

    Returns
    -------
    mpl_img.AxesImage
        Image of the matrix.

    Raises
    ------
    ArgumentError
        Error in arguments.
    """
    def coords(x: float, y: float) -> str:
        num_rows, num_cols = np.shape(_mat)
        col = int(x+.5)
        row = int(y+.5)
        if col >= 0 and col < num_cols and row >= 0 and row < num_rows:
            z = _mat[row, col]
            fmt = 'i={x:1.0f}, k={y:1.0f}, J^2(i,k)={z:1.4f}'
        else:
            z = 0.0
            fmt = 'i={x:1.0f}, k={y:1.0f}'
        return fmt.format(x=x+1, y=y+1, z=z)

    def vmode(x, pos):
        return '{:d}'.format(int(x+1))

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
    order = top_down and 'upper' or 'lower'
    if invert_modes:
        _mat = _mat[::-1, ::-1]
    plot = canvas.matshow(_mat, cmap=mpl.cm.gray_r, vmin=0.0, vmax=1.0,
                          origin=order)
    if show_grid:
        canvas.grid(color='.9')
    canvas.tick_params(top=True, bottom=True, labeltop=False, labelbottom=True)
    canvas.set_xlabel('Final-state modes')
    canvas.set_ylabel('Initial-state modes')
    # Change tick labels to correct numbering (Python starts at 0)
    canvas.xaxis.set_major_formatter(FuncFormatter(vmode))
    canvas.yaxis.set_major_formatter(FuncFormatter(vmode))
    canvas.format_coord = coords

    return plot


def plot_cmat(mat: np.ndarray,
              canvas: mpl.axes.Axes,
              show_grid: bool = True,
              *,
              top_down: bool = False,
              invert_modes: bool = False
              ) -> tp.Tuple[float, mpl_img.AxesImage]:
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
    show_grid
        Show light grid to hel find the modes.
    top_down
        Orders the modes so that the first one is at the top.
    invert_modes
        Invert ordering of normal modes (== decreasing order).

    Returns
    -------
    float
        Normalization factor.
    mpl_img.AxesImage
        Image of the matrix.
    """
    def coords(x: float, y: float) -> str:
        num_rows, num_cols = np.shape(_mat)
        col = int(x+.5)
        row = int(y+.5)
        if col >= 0 and col < num_cols and row >= 0 and row < num_rows:
            z = _mat[row, col]
            fmt = 'i={x:1.0f}, k={y:1.0f}, C(i,k)={z:1.4f}'
        else:
            z = 0.0
            fmt = 'i={x:1.0f}, k={y:1.0f}'
        return fmt.format(x=x+1, y=y+1, z=z)

    _mat = mat.copy()
    norm = np.max(np.abs(_mat))
    _mat /= norm
    order = top_down and 'upper' or 'lower'
    if invert_modes:
        _mat = _mat[::-1, ::-1]
    plot = canvas.matshow(_mat, cmap=mpl.cm.seismic_r, vmin=-1.0, vmax=1.0,
                          origin=order)
    if show_grid:
        canvas.grid(color='.9')
    canvas.tick_params(top=True, bottom=True, labeltop=False, labelbottom=True)
    canvas.set_xlabel('Final-state modes')
    canvas.set_ylabel('Final-state modes')
    # Change tick labels to correct numbering (Python starts at 0)
    def vmode(x, pos): return '{:d}'.format(int(x+1))
    canvas.xaxis.set_major_formatter(FuncFormatter(vmode))
    canvas.yaxis.set_major_formatter(FuncFormatter(vmode))
    canvas.format_coord = coords

    return norm, plot


def plot_kvec(vec: np.ndarray,
              canvas: mpl.axes.Axes,
              show_grid: bool = True,
              *,
              top_down: bool = False,
              invert_modes: bool = False
              ) -> mpl_cont.BarContainer:
    """Plots a horizontal bar chart to represent a K-like vector.

    Plots a shift vector as a horizontal bars.

    Parameters
    ----------
    vec
        Vector to plot.
    canvas
        Matplotlib plotting frame.
    show_grid
        Show light grid to hel find the modes.
    top_down
        Orders the modes so that the first one is at the top.
    invert_modes
        Invert ordering of normal modes (== decreasing order).

    Returns
    -------
    mpl.container.BarContainer
        Image of the vector.
    """
    def coords(x: float, y: float) -> str:
        num_rows = len(_vec)
        row = int(y+.5) + 1
        if row > 0 and row <= num_rows:
            z = _vec[row-1]
            fmt = 'i={y:1.0f}, K(i)={z:.4f}'
        else:
            z = 0.0
            fmt = 'i={y:1.0f}'
        return fmt.format(y=y+1, z=z)

    def vmode(x, pos):
        return '{:d}'.format(int(x+1))

    llo0 = '\N{SUBSCRIPT ZERO}'
    lloe = '\N{LATIN SUBSCRIPT SMALL LETTER E}'
    ldot = '\N{MIDDLE DOT}'
    xlab_kvec = f'Displacement / m{lloe}{ldot}a{llo0}'
    if invert_modes:
        _vec = vec[::-1]
    else:
        _vec = vec
    avec = np.abs(_vec)
    num = len(_vec)
    pos = np.arange(num)
    plot = canvas.barh(pos, avec, align='center')
    for i in range(num):
        if _vec[i] < 0:
            plot[i].set_color('red')
        else:
            plot[i].set_color('#73BFFF')
    if show_grid:
        canvas.grid(color='.9')
    canvas.set_xlabel(xlab_kvec)
    canvas.set_ylabel('Initial-state modes')
    if top_down:
        canvas.set_ylim((num, -.5))
    else:
        canvas.set_ylim((-.5, num))
    canvas.yaxis.set_major_formatter(FuncFormatter(vmode))
    canvas.format_coord = coords

    return plot
