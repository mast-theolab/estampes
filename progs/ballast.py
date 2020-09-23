"""
    BALLAST: Builder Assistant to Lay out, Label and Arrange Spectra
                                Together

This is a simple program to combine and display spectra together.
"""


import sys
import os
import argparse
import typing as tp
import configparser as cfg

import numpy as np
import matplotlib.pyplot as plt

from estampes.base.spectrum import Spectrum
from estampes.visual.plotspec import SpecLayout


def build_opts(parser: argparse.ArgumentParser) -> tp.NoReturn:
    """Builds commandline options.

    Builds commandline options inside input `parser`.

    Parameters
    ----------
    parser
        Parser to update.
    """
    parser.add_argument('optfile', nargs='?',
                        help='Option file (INI style).')
    # parser.add_argument('-o', '--output',
    #                   help='Output file.')
    parser.add_argument('-c', '--colors', action='append',
                        help='Spectral colors.')
    msg = '''\
Colors of the spectra.  By default, it follows the order of input files.
It is possible to change the order by putting a number followed by ":".
Ex. '3:Test' means that the label 'Test' is for the 3rd file (start at 1).
'r'/'e'/'0' refers to the reference data.
'''
    parser.add_argument('-i', '--inpfile', action='append',
                        help='Input data file.')
    msg = '''\
Labels for the legend.  By default, it follows the order of input files.
It is possible to change the order by putting a number followed by ":".
Ex. '3:Test' means that the label 'Test' is for the 3rd file (start at 1).
'r'/'e'/'0' refers to the reference data.
'''
    parser.add_argument('-l', '--label', action='append',
                        help=msg)
    parser.add_argument('-r', '--refdata',
                        help='Reference spectrum file.')


def parse_args(args: tp.Sequence[str]) -> argparse.Namespace:
    """Parses arguments.

    Parses commandline arguments

    Parameters
    ----------
    args
        Commandline arguments

    Returns
    -------
    :obj:`argparse.Namespace`
        Object holding results as attributes
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    build_opts(parser)
    return parser.parse_args(args)


def parse_inifile(fname) -> tp.Tuple[tp.Dict[str, tp.Any],
                                     SpecLayout,
                                     tp.Dict[str, tp.Any]]:
    """Parses INI file.

    Parses a INI configuration file.

    Parameters
    ----------
    fname
        Filename

    Returns
    -------
    figdat
        Figure data
    spcdat
        Spectrum layout parameters
    curves
        Curves data

    Raises
    ------
    FileNotFoundError
        INI file or input file missing
    ValueError
        Incorrect parameter
    """
    if not os.path.exists(fname):
        raise FileNotFoundError('Missing INI file')
    opts = cfg.ConfigParser()
    opts.read(fname)
    secs = {key.lower(): key for key in opts.sections()}
    if 'layout' in secs:
        optlay = opts[secs['layout']]
        _title = optlay.get('Title', None)
        _xleft = optlay.get('XLeft', None) or optlay.get('XMin', None)
        _xright = optlay.get('XRight', None) or optlay.get('XMax', None)
        _ytop = optlay.get('YTop', None) or optlay.get('YMax', None)
        _ybottom = optlay.get('YBottom', None) or optlay.get('YMin', None)
        _xscale = optlay.get('XScale', 'linear')
        _yscale = optlay.get('YScale', 'linear')
        _xlabel = optlay.get('XLabel', None)
        _ylabel = optlay.get('YLabel', None)
        spcdat = SpecLayout(_xleft, _xright, _ytop, _ybottom, _xscale, _yscale,
                            _xlabel, _ylabel, _title)
        _legpos = optlay.get('Legend', 'best')
        _legcol = optlay.get('Legend_Cols', 1)
        spcdat.legend(_legpos, _legcol)
    else:
        spcdat = SpecLayout()

    figdat = {
        'geom': None,
    }
    if 'figure' in secs:
        optfig = opts[secs['figure']]
        _geom = optfig.get('Geometry', None)
        if _geom is not None:
            _res = _geom.replace('(', '').replace(')', '').split(',')
        if len(_res) == 1:
            raise ValueError('Incorrect value for geometry.')
        else:
            _geom = (float(_res[0]), float(_res[1]))
        figdat['geom'] = _geom

    curves = {}
    for sec in opts.sections():
        if sec.strip().lower().startswith('curve'):
            label, key = [item.strip() for item in sec.split(':', maxsplit=1)]
            if label.lower() != 'curve':
                print(sec, 'will be ignored as a curve definition.')
                continue  # This is not a right keyword, ignore.
            curve = opts[sec]
            if not curve.getboolean('show', fallback=True):
                continue
            curves[key] = {}
            if 'file' not in curve:
                print(f'WARNING: Missing file for "{sec}". Ignoring.')
                continue
            elif not os.path.exists(curve['file']):
                fmt = 'ERROR: File "{}" not found in "{}".'
                print(fmt.format(curve['file'], sec))
            spc = curve.get('spectroscopy', fallback=None)
            lvl = curve.get('level', fallback=None)
            if spc is None or lvl is None:
                raise ValueError('Spectroscopy not defined')
            yid = curve.get('yaxis', None)
            if yid is not None:
                yid = 'y' + yid
            curves[key]['data'] = Spectrum(curve['file'], spc, lvl, yid)
            if curve.getboolean('broaden', fallback=False):
                func = curve.get('function', None)
                hwhm = curve.get('HWHM', None)
                if hwhm is not None:
                    hwhm = float(hwhm)
                curves[key]['data'].set_broadening(hwhm, func, 'default')
            vizdata = {}
            for item in ('color', 'linestyle', 'linewidth'):
                if curve.get(item, False):
                    vizdata[item] = curve.get(item, False)
            if vizdata:
                curves[key]['data'].set_display(**vizdata)
            if curve.get('label', None) is not None:
                curves[key]['data'].label = curve.get('label')
            curves[key]['xshift'] = curve.getfloat('xshift', fallback=None)
            curves[key]['xscale'] = curve.getfloat('xscale', fallback=None)
            res = curve.get('yshift', None)
            if res is not None:
                try:
                    val = float(res)
                except ValueError:
                    if res.lower() in ('base', 'baseline'):
                        val = 'base'
                    else:
                        raise ValueError('Unsupported value for YShift')
            else:
                val = None
            curves[key]['yshift'] = val
            curves[key]['yscale'] = curve.getfloat('yscale', fallback=None)
            curves[key]['ynorm'] = curve.getboolean('normalize',
                                                    fallback=False)

    return figdat, spcdat, curves


def main() -> tp.NoReturn:
    """Main function.
    """
    args = parse_args(sys.argv[1:])
    if not args.inpfile and not args.optfile:
        print('ERROR: Missing files or option file.')
        sys.exit(2)
    elif args.inpfile and args.optfile:
        msg = 'ERROR: Option file and single files cannot be treated' \
            + ' together'
        print(msg)
        sys.exit(2)
    elif args.inpfile:
        print('ERROR: Files in input not yet supported')
    else:
        figdata, spcdata, curves = parse_inifile(args.optfile)

    fig, subp = plt.subplots(1, 1)
    if figdata['geom'] is not None:
        fig.set_size_inches(figdata['geom'])
    for key in curves:
        xaxis = np.array(curves[key]['data'].xaxis)
        if curves[key]['xshift'] is not None:
            xaxis += curves[key]['xshift']
        if curves[key]['xscale'] is not None:
            xaxis *= curves[key]['xscale']
        yaxis = np.array(curves[key]['data'].yaxis)
        if curves[key]['yshift'] is not None:
            if curves[key]['yshift'] == 'base':
                yshift = - np.min(yaxis)
            else:
                yshift = curves[key]['yshift']
            yaxis += yshift
        if curves[key]['ynorm']:
            yaxis /= np.max(np.abs(yaxis))
        elif curves[key]['yscale'] is not None:
            yaxis *= curves[key]['yscale']
        stick = curves[key]['data'].get_broadening('func') == 'stick'
        data = {}
        if curves[key]['data'].label is not None:
            data['label'] = curves[key]['data'].label
        if curves[key]['data'].linecolor is not None:
            data['color'] = curves[key]['data'].linecolor
        if curves[key]['data'].linewidth is not None:
            data['linewidth'] = curves[key]['data'].linewidth
        if stick:
            zeros = np.zeros(len(yaxis))
            subp.vlines(xaxis, zeros, yaxis, **data)
        else:
            if curves[key]['data'].linestyle is not None:
                data['linestyle'] = curves[key]['data'].linestyle
            subp.plot(xaxis, yaxis, **data)
        subp.legend()
        spcdata.set_plot(subp)
    plt.show()


if __name__ == '__main__':
    main()
