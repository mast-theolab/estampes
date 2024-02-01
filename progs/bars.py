"""Program: BARS.

BARS: Benchmark on Anharmonicity - Results Survey

Simple tool to merge, compare and plot results from a set of anharmonic
  calculations.
"""

import sys
import os
import argparse
from collections.abc import Sequence
import typing as tp

import numpy as np

from estampes.base import QLabel, ArgumentError, ParseKeyError, QuantityError
from estampes.parser import DataFile

try:
    import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # from matplotlib.ticker import FuncFormatter
    from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    rc('mathtext', fontset='stixsans')
    MPL_LOADED = True
    # mpl.rcParams['backend'] = 'SVG'
    # mpl.rcParams['backend.qt4'] = 'PyQt4'
    plt.rcParams['font.family'] = 'sans-serif'
except ImportError:
    MPL_LOADED = False


PROGNAME = 'bars'
HELP_OPTF = '''Option file with information on input files for the benchmark.
Each line contains the following columns (separated by ";"):
1. molecule name
2. label for the functional/electronic structure method
3. label for the basis set
4. label for the anharmonic level : A (anharmonic), H (harmonic), B (both), \
    X (auto, the default if blank field)
5. name of the file with data
6. (optional) tag to override automatic definition
7. (optional) associated color (Matplotlib-compatible)'''
HELP_MOL = 'Molecule of interest. \
"all" to treat all molecules in `optfile` together, \
"each" to treat them separately.'
HELP_XAXD = '''Type of X axis:
- linear: based on the quantity (e.g., energy) value
- range: results are gathered by range (bounds separated by commas)
- index: the values are equidistant, the x axis being the index
NOTE: In each case, the reference if the value from the `refdata` file'''
HELP_SORT = '''Sort the electronic structure methods:
- file: respect the order given in the file
- alpha: reorder by sorting the methods in alphabetic order'''

TypeMol = list[str]
TypeESM = list[str]
TypeGBS = list[str]
TypeVib = list[str]
TypeInF = list[str]
TypeTag = list[tp.Optional[str]]
TypeCol = list[tp.Optional[str]]
TypeOptDat = tuple[TypeMol, TypeESM, TypeGBS, TypeVib, TypeInF, TypeTag,
                   TypeCol]
# Type for calculated values (for a given method)
# [level][mode] = value
TypeCalc = dict[str, dict[int, float]]
# Type for labels: [label][level][mode] = value
TypeMeth = dict[str, TypeCalc]
# Type for molecules: [mol][method][level][mode] = value
TypeMols = dict[str, TypeMeth]
# Type for statistics: [method][level][mode] = value
TypeDiff = dict[str, dict[int, dict[str, float]]]
# Type for reference data: [mol][mode] = value
TypeRDat = dict[str, dict[int, float]]


def parse_args(args: Sequence[str]) -> argparse.Namespace:
    """Define and check options given in arguments.

    Defines available options in the program, parses a list of arguments
    given in input, and returns the results as a populated namespace.

    Parameters
    ----------
    args
        List of options/arguments to parse

    Returns
    -------
    obj:`argparse.Namespace`
        Returned results from the parsing
    """
    parser = argparse.ArgumentParser(
        prog=PROGNAME,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile',
                        type=argparse.FileType('r'),
                        help=HELP_OPTF)
    parser.add_argument('reffile',
                        type=argparse.FileType('r'),
                        help='File with reference (e.g., experimental) data')
    parser.add_argument('--fontsize',
                        type=int,
                        help='Font size for the plot(s)')
    txt = 'Level of theory of interest (harm/anharm/both/mixed)'
    parser.add_argument('--level',
                        type=str.lower, default='both',
                        choices=('a', 'anh', 'anharm', 'h', 'harm', 'b',
                                 'both', 'm', 'mixed'),
                        help=txt)
    parser.add_argument('-m', '--mol',
                        help=HELP_MOL)
    parser.add_argument('-o', '--output',
                        help='Output image file')
    parser.add_argument('-s', '--sort',
                        choices=('file', 'alpha'), default='file',
                        help=HELP_SORT)
    psub = parser.add_subparsers(help='Operating modes')
    pmean = psub.add_parser('mean',
                            help='Computes averages')
    pmean.set_defaults(mode='mean')
    pmean.add_argument('-c', '--compare', default=None,
                       help='File with methods to compare.')
    pdist = psub.add_parser('dist',
                            help='Shows distributions')
    pdist.set_defaults(mode='dist')
    pdist.add_argument('--absval',
                       action='store_true', default=False,
                       help='Use absolute values')
    pdist.add_argument('--linewidth', '--lw',
                       type=float,
                       help='Line width of the curves')
    pdist.add_argument('-x', '--xscale',
                       choices=('linear', 'range', 'index'), default='linear',
                       help=HELP_XAXD)
    opts = parser.parse_args(args)
    if not hasattr(opts, 'mode'):
        setattr(opts, 'mode', 'mean')
    if opts.mol is None:
        if opts.mode == 'mean':
            opts.mol = 'all'
        else:
            opts.mol = 'each'

    return opts


def parse_optfile(fobj: tp.IO) -> TypeOptDat:
    """Parse option file.

    Parses Option file and returns the data as series of list:
    1. List of molecule names
    2. List of electronic structure methods/functionals
    3. List of basis sets
    4. List of levels
    5. List of input filenames
    6. List of labels
    7. List of color codes

    Parameters
    ----------
    fobj
        File object from opened open file

    Returns
    -------
    list
        List of molecule names (one per line)
    list
        List of electronic structure methods/functionals
    list
        List of Gaussian basis sets
    list
        List of input filenames
    list
        List of color codes

    Raises
    ------
    IndexError
        Incorrect construct of one row
    """
    csep = ';'
    lmol = []
    lmet = []
    lbas = []
    lvib = []
    lnam = []
    ltag = []
    lcol = []
    for line in fobj:
        if not line.lstrip().startswith('#'):
            cols = line.split(csep)
            lcols = len(cols)
            if lcols == 4:  # Old minimal input
                id_nam = 3
                lvib.append('X')
                ltag.append(None)
                lcol.append(None)
            else:
                id_nam = 4
                if lcols == 7:
                    lcol.append(cols[6].strip())
                    ltag.append(cols[5].strip())
                elif lcols == 6:
                    lcol.append(None)
                    ltag.append(cols[5].strip())
                elif lcols == 5:
                    lcol.append(None)
                    ltag.append(None)
                else:
                    raise IndexError('At least 4 columns expected in row')
                vtag = (cols[3].strip().upper() or 'X')[0]
                if vtag not in ('A', 'H', 'B', 'X'):
                    raise IndexError('Incorrect label for vibrational level')
                lvib.append(vtag)
            lnam.append(cols[id_nam].strip())
            lbas.append(cols[2].strip())
            lmet.append(cols[1].strip())
            lmol.append(cols[0].strip())

    return lmol, lmet, lbas, lvib, lnam, ltag, lcol


def parse_reffile(fobj: tp.IO,
                  mol: str) -> TypeRDat:
    """Parse reference data file.

    Parses the reference data file and returns a dictionary, as:
    molecule -> mode -> value

    Parameters
    ----------
    fobj
        File object from opened reference data file
    mol
        Molecule of interest ("all" or "each" to include all)

    Returns
    -------
    dict
        Values as [mol][mode][value]

    Raises
    ------
    ValueError
        Incorrect format
    IndexError
        Molecule not found
    """
    data = {}
    label = ''
    for line in fobj:
        cols = line.strip().split()
        if len(cols) == 1:
            label = cols[0].rstrip(':')
            if mol.lower() in ('all', 'each') or label == mol:
                data[label] = {}
        elif label in data:
            try:
                mode = int(cols[0])
                val = float(cols[1])
            except ValueError as e:
                fmt = f'Error: Incorrect format\nCheck line: {line}'
                raise ValueError(fmt) from e
            data[label][mode] = val
    if not data:
        raise IndexError('No reference data for the molecule(s).')

    return data


def parse_compfile(fobj: tp.IO) -> list[tuple[tp.Optional[str],
                                              tp.Optional[str]]]:
    """Parse file with methods to compare.

    Accepted structure:
    * method ; basis || method ; basis
    * tag || tag

    Parameters
    ----------
    fobj
        File object from opened open file

    Returns
    -------
    List
        List of tuples of methods to compare.

    Raises
    ------
    IndexError
        Incorrect field construction.
    """
    csep = ';'
    gsep = '||'
    methods = []
    nline = 0
    for line in fobj:
        nline += 1
        if not line.lstrip().startswith('#'):
            groups = line.split(gsep)
            meths = [None, None]
            if len(groups) > 2:
                raise IndexError(
                    f'Only 2 methods/labels per line (l. {nline})')
            for i, group in enumerate(groups):
                if group.strip():
                    cols = [item.strip() for item in group.split(csep)]
                    if len(cols) > 2:
                        raise IndexError('Too many fields on line', nline)
                    elif len(cols) == 2:
                        if cols[0] == '' or cols[1] == '':
                            raise IndexError('Empty field on line', nline)
                        meths[i] = cols[0] + '/' + cols[1]
                    elif len(cols) == 1:
                        if cols[0] == '':
                            raise IndexError('Empty field on line', nline)
                        meths[i] = cols[0]
            methods.append(tuple(meths))

    return methods


def get_energies(dfname: str,
                 levels: Sequence[str]) -> TypeCalc:
    """Get vibrational energies from data file.

    Gets vibrational fundamental energies from a data file.

    Parameters
    ----------
    dfname
        Data filename.
    levels
        Sequence with the levels of theory of interest

    Returns
    -------
    dict
        Vibrational energies, sorted by level and mode index.
        H: harmonic data (if present in `levels`)
        A: anharmonic data (if present in `levels`)

    Raises
    ------
    ArgumentError
        Unsupported level.
    IndexError
        Quantity not found
    """
    dfile = DataFile(dfname)
    dkeys = {}
    for level in levels:
        if level[-1] not in ('H', 'A'):
            raise ArgumentError('Unsupported level of theory.')
    indata = {}
    dkeys = {}
    for item in levels:
        level = item[-1]
        qkeys = {
            'assign': QLabel(quantity='vtrans', level=level),
            'freq': QLabel(quantity='vlevel', level=level),
        }
    # keys = [item for key in dkeys for item in dkeys[key].values()]
        # Extract relevant data
        try:
            indata[level] = dfile.get_data(**qkeys)
            dkeys[level] = qkeys.copy()
        except ParseKeyError as e:
            if item[0] != '?':
                msg = 'Vibrational data missing' + str(e)
                raise IndexError(msg) from None
    # Check if we have the states are variational or harmonic/DVPT2
    if 'A' in dkeys:
        state1_nq = indata['A']['assign'].data[1][1][0][1]
        if state1_nq == 0:
            raise IndexError('Variational state notation found.')

    data = {}
    for level in dkeys:
        data[level] = {}
        for i, val in indata[level]['assign'].data.items():
            fsta = val[1]
            if len(fsta) == 1 and fsta[0][1] == 1:
                mode = fsta[0][0]
                data[level][mode] = indata[level]['freq'].data[i]
    return data


def print_data(mol: str,
               labels: list[str],
               dcalc: TypeMeth,
               dref: dict[int, float],
               skip_on_ref: bool = True) -> TypeDiff:
    """Print data.

    Prints data in CSV-compliant files.
    Returns the signed differences, in the form of dictionary
    [method][level][mode] = value

    Parameters
    ----------
    mol
        Molecule name (label)
    labels
        "Dictionary" of all chemical available methods, in proper order
    dcalc
        Computed data as `[method][mode][level] = data`
    dref
        Reference data as `[mode] = data`
    skip_on_ref
        Skip mode if reference data unavailable

    Returns
    -------
    dict
        Signed difference with respect to reference data for each mode,
        as `[method][level][mode] = value`
    """
    prec = 0
    fmt_fp = f'{{:{prec+7}.{prec}f}}'
    blank = (prec+7)*' '
    fmt_func = f'{{:{2*(prec+7)+3}s}}'  # 3 for ' ; '
    fmt_base = ' {:4d}'
    fmt_ref = f'{{:{prec+7}s}}'
    txt_head = ' Mode ; ' + fmt_ref.format('Ref.')
    txt_head2 = '      ; ' + fmt_ref.format(' ')
    fmt_delta = ' {:4s} ; ' + fmt_ref.format(' ')
    diffs = {}
    nmax = 0
    with open(f'res_{mol}.csv', 'w', encoding='utf-8') as fobj:
        txt0 = txt_head
        txt1 = txt_head2
        for key in labels:
            if key not in dcalc:  # pass if method not supported
                continue
            # Below: 3 for ' ; '
            fmt_func = f'{{:{len(dcalc[key].keys())*(prec+7)+3}s}}'
            txt0 += ' ; ' + fmt_func.format(key)
            diffs[key] = {}
            for lvl in dcalc[key]:
                if max(dcalc[key][lvl]) > nmax:
                    nmax = max(dcalc[key][lvl])
                txt1 += ' ; ' + fmt_ref.format(lvl)
                diffs[key][lvl] = {}
        txt0 += '\n'
        txt1 += '\n'
        fobj.write(txt0)
        fobj.write(txt1)
        if max(dref) > nmax:
            nmax = max(dref)
        # nmax is the highest mode, we use it to check all modes
        for i in range(1, nmax+1):
            res = []
            if i in dref:
                res.append(fmt_fp.format(dref[i]))
                dodiff = True
            else:
                dodiff = False
                if skip_on_ref:
                    continue
                res.append(blank)
            for key in labels:
                if key not in dcalc:  # pass if method not supported
                    continue
                for lvl in dcalc[key]:
                    if i not in dcalc[key][lvl]:
                        res.append(blank)
                    else:
                        res.append(fmt_fp.format(dcalc[key][lvl][i]))
                        if dodiff:
                            diffs[key][lvl][i] = dcalc[key][lvl][i] - dref[i]
            if (''.join(res)).strip():
                fobj.write(fmt_base.format(i) + ' ; '.join(res) + '\n')

        # Print statistics
        res1 = []
        res2 = []
        for key in labels:
            if key not in dcalc:  # pass if method not supported
                continue
            for lvl in diffs[key]:
                amax = 0.0
                dmue = 0.0
                nvib = 0
                for mode in diffs[key][lvl]:
                    val = abs(diffs[key][lvl][mode])
                    nvib += 1
                    dmue += val
                    if val > amax:
                        amax = val
                res1.append(fmt_fp.format(amax))
                if nvib > 0:
                    res2.append(fmt_fp.format(dmue/nvib))
                else:
                    res2.append(blank)
        fobj.write(fmt_delta.format('AMAX') + ' ; '.join(res1) + '\n')
        fobj.write(fmt_delta.format('MUE') + ' ; '.join(res2) + '\n')

        return diffs


def build_mean_data(op_mol: str,
                    labels: list[str],
                    levels: Sequence[str],
                    dcalc: dict[str, TypeMeth],
                    dref: TypeRDat,
                    do_plot: bool,
                    compfile: tp.Optional[str] = None,
                    imgfile: tp.Optional[str] = None):
    """Build (and display) average statistical data.

    Builds average statistical data and displays if requested.
    Also handles call to the function in charge of printing the data.

    Parameters
    ----------
    op_mol
        Option for the molecule
        - all: treat all molecules together
        - each: treat each molecule separately (NYI)
        - name of the molecule of interest
    labels
        "Dictionary" of all labels, in proper order.
    levels
        Sequence with the levels of theory of interest.
    dcalc
        Computed data as `[mol][method][mode][level] = data`
    dref
        Reference data as `[mol][mode] = data`
    do_plot
        Plot results
    compfile
        Comparison file
    imgfile
        Image file (only saved if not None)

    Raises
    ------
    NotImplementedError
        Treated separately each molecule is not yet possible
    KeyError
        Mismatch between molecules in `dref` and `dcalc`
    """
    def autolabel(ax, rects, color=None):
        """Attach a text label above each bar, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{height:.0f}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points", fontsize=6,
                        color=color,
                        ha='center', va='bottom', rotation=90)

    if op_mol.lower() == 'each':
        raise NotImplementedError('Separate treatment of molecules NYI')

    nmols = len(set(dcalc) & set(dref))
    if nmols == 0:
        raise KeyError('No overlap between `dref` and `dcalc`')

    full_amax = {}
    full_dmue = {}
    full_nvib = {}
    for key in labels:
        full_amax[key] = {}
        full_dmue[key] = {}
        full_nvib[key] = {}
        for lvl in levels:
            full_amax[key][lvl] = 0.0
            full_dmue[key][lvl] = 0.0
            full_nvib[key][lvl] = 0
    for mol in dcalc:
        if mol not in dref:
            print(f'Molecule "{mol}" not available in experiment')
        diffs = print_data(mol, labels, dcalc[mol], dref[mol])
        for key, diff in diffs.items():
            for lvl in levels:
                for val in (abs(item) for item in diff[lvl].values()):
                    full_nvib[key][lvl] += 1
                    full_dmue[key][lvl] += val
                    if val > full_amax[key][lvl]:
                        full_amax[key][lvl] = val

    compmeths = parse_compfile(open(compfile, 'r', encoding="utf-8")) \
        if compfile else []
    if do_plot:
        c_MAE = {'H': '#1317ff', 'A': '#bd2626'}
        c_MAX = {'H': '#6baff0', 'A': '#ff86db'}

        figid = 1
        # width2 = 0.4

        ind = np.arange(len(labels))
        if len(levels) == 2:
            width = 0.18
            if compmeths:
                figsize = (7, 12)
                ind = np.arange(len(compmeths))
                dmue_h = [[], []]
                dmue_a = [[], []]
                amax_h = [[], []]
                amax_a = [[], []]
                ddmue_h = []
                ddmue_a = []
                damax_h = []
                damax_a = []
                rlabels = []
                llabels = []
                err_max = 0.0
                derr_max = 0.0
                for key1, key2 in compmeths:
                    if key1 is not None:
                        v1 = full_dmue[key1]['H'] / full_nvib[key1]['H']
                        v2 = full_dmue[key1]['A'] / full_nvib[key1]['A']
                        v3 = full_amax[key1]['H']
                        v4 = full_amax[key1]['A']
                        vmax = max(v1, v2, v3, v4)
                        if vmax > err_max:
                            err_max = vmax
                        dmue_h[0].append(v1)
                        dmue_a[0].append(v2)
                        amax_h[0].append(v3)
                        amax_a[0].append(v4)
                        rlabels.append(key1)
                    else:
                        dmue_h[0].append(0.0)
                        dmue_a[0].append(0.0)
                        amax_h[0].append(0.0)
                        amax_a[0].append(0.0)
                        rlabels.append('')
                    if key2 is not None:
                        v1 = full_dmue[key2]['H'] / full_nvib[key2]['H']
                        v2 = full_dmue[key2]['A'] / full_nvib[key2]['A']
                        v3 = full_amax[key2]['H']
                        v4 = full_amax[key2]['A']
                        vmax = max(v1, v2, v3, v4)
                        if vmax > err_max:
                            err_max = vmax
                        dmue_h[1].append(-v1)
                        dmue_a[1].append(-v2)
                        amax_h[1].append(-v3)
                        amax_a[1].append(-v4)
                        llabels.append(key2)
                    else:
                        dmue_h[1].append(0.0)
                        dmue_a[1].append(0.0)
                        amax_h[1].append(0.0)
                        amax_a[1].append(0.0)
                        llabels.append('')
                    if key1 is not None and key2 is not None:
                        ddmue_h.append(dmue_h[0][-1] + dmue_h[1][-1])
                        ddmue_a.append(dmue_a[0][-1] + dmue_a[1][-1])
                        damax_h.append(amax_h[0][-1] + amax_h[1][-1])
                        damax_a.append(amax_a[0][-1] + amax_a[1][-1])
                        vmax = max([ddmue_h[-1], ddmue_a[-1], damax_h[-1],
                                    damax_h[-1]], key=abs)
                        if vmax > derr_max:
                            derr_max = vmax
                    else:
                        ddmue_h.append(0.0)
                        ddmue_a.append(0.0)
                        damax_h.append(0.0)
                        damax_a.append(0.0)
                derr_max = max(derr_max, 10)
                err_max += 10
            else:
                figsize = (11, 7)
                dmue_h = []
                dmue_a = []
                amax_h = []
                amax_a = []
                ind = np.arange(len(labels))
                for key in labels:
                    dmue_h.append(full_dmue[key]['H']/full_nvib[key]['H'])
                    dmue_a.append(full_dmue[key]['A']/full_nvib[key]['A'])
                    amax_h.append(full_amax[key]['H'])
                    amax_a.append(full_amax[key]['A'])

            fig = plt.figure(figid, figsize=figsize, dpi=300)
            fig.subplots_adjust(hspace=0.001)

            ax1 = fig.add_subplot(111)
            if compmeths:
                ax2 = ax1.twinx()
                ax3 = ax1.twiny()
                ax1.barh(ind-width/2, dmue_h[1], width, color=c_MAE['H'])
                ax1.barh(ind+width/2, amax_h[1], width, color=c_MAX['H'])
                ax1.barh(ind+3*width/2, amax_a[1], width, color=c_MAX['A'])
                ax1.barh(ind+5*width/2, dmue_a[1], width, color=c_MAE['A'])
                ax2.barh(ind-width/2, dmue_h[0], width, color=c_MAE['H'])
                ax2.barh(ind+width/2, amax_h[0], width, color=c_MAX['H'])
                ax2.barh(ind+3*width/2, amax_a[0], width, color=c_MAX['A'])
                ax2.barh(ind+5*width/2, dmue_a[0], width, color=c_MAE['A'])
                ax3.barh(ind-width/2, ddmue_h, width, edgecolor='black',
                         lw=.4, color='None')
                ax3.barh(ind+width/2, damax_h, width, edgecolor='black',
                         lw=.4, color='None')
                ax3.barh(ind+3*width/2, damax_a, width, edgecolor='black',
                         lw=.4, color='None')
                ax3.barh(ind+5*width/2, ddmue_a, width, edgecolor='black',
                         lw=.4, color='None')
                ax1.set_xlim(-err_max, err_max)
                ax3.set_xlim(-derr_max, derr_max)
                ax3.set_xlabel(r'Difference between errors / cm$^{-1}$')
                ax1.set_xlabel(r'Error wrt Experiment / cm$^{-1}$')
                ax1.grid(axis='y', lw=.2, ls='--', c='gray')
                rlab = [item.replace(r'_', r'\_') for item in rlabels]
                ax2.set_yticks(ind+width)
                ax2.set_yticklabels(tuple(rlab), fontsize=8)
                llab = [item.replace(r'_', r'\_') for item in llabels]
                ax1.set_yticks(ind+width)
                ax1.set_yticklabels(tuple(llab), fontsize=8)
                ax2.axvline(0.0, lw=.2, color='black')
            else:
                hmax = ax1.bar(ind-width/2, amax_h, width, color=c_MAX['H'])
                hmue = ax1.bar(ind+width/2, dmue_h, width, color=c_MAE['H'])
                amax = ax1.bar(ind+3*width/2, amax_a, width, color=c_MAX['A'])
                amue = ax1.bar(ind+5*width/2, dmue_a, width, color=c_MAE['A'])
                autolabel(ax1, hmue, c_MAE['H'])
                autolabel(ax1, hmax, c_MAX['H'])
                autolabel(ax1, amue, c_MAE['A'])
                autolabel(ax1, amax, c_MAX['A'])
                # ymax = max(*dmue_h, *amax_h, *amax_a, *dmue_a)
                # ax1.set_ylim(top=ymax*1.1)
                # ymax = ax1.get_ylim()[1]
                # ax1.set_ylim(top=ymax+10)
                ax1.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
                ax1.grid(axis='x', lw=.2, ls='--', c='gray')
                ax1.set_xticks(ind+width)
                lab_meth = [item.replace(r'_', r'\_') for item in labels]
                ax1.set_xticklabels(tuple(lab_meth), fontsize=8, rotation=25,
                                    ha="right")
            ax1.legend((r'$|$MAX$|$ ($\omega$)', r'MAE ($\omega$)',
                        r'$|$MAX$|$ ($\nu$)', r'MAE ($\nu$)'),
                       loc='best', ncol=2)

            if imgfile is not None:
                plt.savefig(imgfile, bbox_inches='tight')
        else:
            lvl = levels[0]
            if lvl == 'H':
                label = r'$\omega'
            else:
                label = r'$\nu'
            width = 0.4
            dmue = []
            amax = []
            for key in labels:
                dmue.append(full_dmue[key][lvl]/full_nvib[key][lvl])
                amax.append(full_amax[key][lvl])
            figsize = (11, 7)
            fig = plt.figure(figid, figsize=figsize, dpi=300)
            fig.subplots_adjust(hspace=0.001)

            ax1 = fig.add_subplot(111)
            ax1.bar(ind-width, dmue, width, color=c_MAE[lvl])
            ax1.bar(ind+width, amax, width, color=c_MAX[lvl])
            ax1.legend((f'MAE ({label})', f'$|$MAX$|$ ({label})'))
            ax1.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
            ax1.grid(axis='x', lw=.2, ls='--', c='gray')
            ax1.set_xticks(ind)

            lab_meth = [item.replace(r'_', r'\_') for item in labels]
            ax1.set_xticklabels(tuple(lab_meth), fontsize=5, rotation=25,
                                ha="right")
            if imgfile is not None:
                plt.savefig(imgfile, bbox_inches='tight')
        plt.show()
        plt.close(figid)


def build_dist_data(op_mol: str,
                    op_xax: str,
                    absval: bool,
                    labels: Sequence[str],
                    colors: Sequence[str],
                    levels: tp.Optional[Sequence[str]],
                    dcalc: TypeMols,
                    dref: TypeRDat,
                    do_plot: bool) -> None:
    """Build (and display) distribution of data.

    Builds statistical data and displays if requested.
    Also handles call to the function in charge of printing the data.

    Parameters
    ----------
    op_mol
        Option for the molecule
        - all: treat all molecules together (NYI)
        - each: treat each molecule separately
        - name of the molecule of interest
    op_xax
        Option for the scale of the X axis
        - linear: print data based on the value of the quantity
        - range: group quantities by range of values
        - index: use an index (normal mode normally)
    absval
        Use absolute values for the differences.
    labels
        List of all labels, in proper order.
    colors
        List of colors for each method/label.
    levels
        Sequence with the levels of theory of interest
    dcalc
        Computed data as `[mol][label][level][mode] = data`
    dref
        Reference data as `[mol][mode] = data`
    do_plot
        Plot results

    Raises
    ------
    NotImplementedError
        Treated altogether each molecule is not yet possible
    KeyError
        Mismatch between molecules in `dref` and `dcalc`
    """
    def get_diff(d: dict[int, float], mode: int, do_abs: bool):
        """Get and transform differences if needed."""
        val = d.get(mode)
        if val is not None and do_abs:
            val = abs(val)
        return val

    if op_mol.lower() == 'all':
        raise NotImplementedError('Inclusive treatment of molecules NYI')

    lstyle = {'H': '--', 'A': '-'}
    figx = 16
    figy = 8
    if levels is None:
        num_HA = 1
    else:
        num_HA = len(levels)
    nmols = len(set(dcalc) & set(dref))
    if nmols == 0:
        raise KeyError('No overlap between `dref` and `dcalc`')
    if do_plot:
        figid = 1
        figsize = (figx*num_HA, figy*nmols)
        fig = plt.figure(figid, figsize=figsize, dpi=300, tight_layout=True)
        panels = {}
        pan_id = 0
        if levels is not None:
            for mol in dcalc:
                panels[mol] = {}
                for level in levels:
                    pan_id += 1
                    panels[mol][level] = fig.add_subplot(nmols, num_HA, pan_id)
        else:
            for mol in dcalc:
                pan_id += 1
                panels[mol] = fig.add_subplot(nmols, num_HA, pan_id)
        # plt.subplots_adjust(hspace=0.001)

    for mol in dcalc:
        if mol not in dref:
            print(f'Molecule "{mol}" not available in experiment')
        diffs = print_data(mol, labels, dcalc[mol], dref[mol])
        if do_plot:
            if op_xax[:1] == 'l':
                xdat = sorted([(key, dref[mol][key]) for key in dref[mol]],
                              key=lambda x: x[1])
                xaxis = [item[1] for item in xdat]
                xlabel = r'Wavenumbers / cm$^{-1}$'
            elif op_xax[:1] == 'i':
                xdat = sorted([(key, None) for key in dref[mol]],
                              key=lambda x: x[0])
                xaxis = [item[0] for item in xdat]
                xlabel = r'Index'
            else:
                raise NotImplementedError('Range NYI')
            nkey = 0
            for label, color in zip(labels, colors):
                if label not in diffs:
                    continue
                nkey += 1
                pars = {'marker': 'd', 'label': label.replace(r'_', r'\_')}
                if color is not None:
                    pars['color'] = color
                for lvl in diffs[label]:
                    if levels is None:
                        panel = panels[mol]
                    else:
                        panel = panels[mol][lvl]
                    yaxis = []
                    for i, _ in xdat:
                        yaxis.append(get_diff(diffs[label][lvl], i, absval))
                    panel.plot(xaxis, yaxis, ls=lstyle[lvl], **pars)
            if levels is not None:
                for ilv, lvl in enumerate(levels):
                    panel = panels[mol][lvl]
                    fmt = '{} (level: {})'
                    panel.set_title(fmt.format(mol, lvl))
                    if ilv == 0:
                        panel.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
                    panel.set_xlabel(xlabel)
                    panel.set_xticks([int(item) for item in xaxis])
                    if op_xax[:1] == 'l':
                        panel.set_xticklabels(tuple([int(item)
                                                     for item in xaxis]),
                                              rotation='vertical')
                    if ilv == len(levels)-1:
                        panel.legend(bbox_to_anchor=(1.02, 1),
                                     borderaxespad=0.)
            else:
                panel = panels[mol]
                fmt = '{}'
                panel.set_title(fmt.format(mol))
                panel.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
                panel.set_xlabel(xlabel)
                panel.set_xticks([int(item) for item in xaxis])
                if op_xax[:1] == 'l':
                    panel.set_xticklabels(tuple([int(item) for item in xaxis]),
                                          rotation='vertical')
                panel.legend(bbox_to_anchor=(1.02, 1), borderaxespad=0.)
    # plt.savefig(infile[0:-4]+".png",num=figid,figsize=figsize)
    if do_plot:
        plt.show()
        plt.close(figid)


def main():
    """Run the main program."""
    opts = parse_args(sys.argv[1:])

    # Basic parameters
    if MPL_LOADED:
        if opts.fontsize:
            plt.rcParams.update({'font.size': opts.fontsize})
        # Mode-specific MPL parameters
        if opts.mode == 'dist':
            if opts.linewidth:
                plt.rcParams.update({'lines.linewidth': opts.linewidth})

    dmols, dmeth, dbset, dvlvl, dfnam, dtags, dcols = \
        parse_optfile(opts.optfile)
    if opts.mol.lower() not in ('each', 'all'):
        i = 0
        while i < len(dmols):
            if dmols[i] == opts.mol:
                i += 1
            else:
                del dmols[i]
                del dmeth[i]
                del dbset[i]
                del dvlvl[i]
                del dfnam[i]
                del dtags[i]
                del dcols[i]
    if len(dmols) == 0:
        print('ERROR: Molecule not present in option file. Exiting',
              file=sys.stderr)
        sys.exit(1)

    # Extracts computed data
    # First, check the levels
    if opts.level in ('b', 'both'):
        levels = ('H', 'A')
    elif opts.level in ('m', 'mixed'):
        levels = None
    else:
        levels = (opts.level[0].upper(), )
    labels = []
    colors = []
    dcalc = {}
    for i, infile in enumerate(dfnam):
        if not os.path.exists(infile):
            print(f'ERROR: File {infile} not found', file=sys.stderr)
            sys.exit()
        mol = dmols[i]
        if mol not in dcalc:
            dcalc[mol] = {}
        label = dtags[i] or dmeth[i] + '/' + dbset[i]
        try:
            if levels is None:
                if dvlvl[i] in ('A', 'H'):
                    levels_ = (dvlvl[i], )
                elif dvlvl[i] == 'B':
                    levels_ = ('H', 'A')
                else:
                    levels_ = ('?H', '?A')
            else:
                levels_ = levels
            data = get_energies(infile, levels_)
            if levels is None and len(data.keys()) > 1:
                for level, val in data.items():
                    tag = label + f' ({level})'
                    dcalc[mol][tag] = {level: val}
                    if tag not in labels:
                        labels.append(tag)
                        colors.append(dcols[i])
            else:
                dcalc[mol][label] = data
                if label not in labels:
                    labels.append(label)
                    colors.append(dcols[i])
        except IndexError as e:
            print(f'Error while parsing data in file {infile}.')
            print(e)
            sys.exit()
        except QuantityError as e:
            print(f'Error: Unable to process data from file {infile}.')
            print(e)
            sys.exit()

    # Check if reorder necessary
    if opts.sort == 'alpha':
        labels.sort()

    # Extracts reference data
    try:
        dref = parse_reffile(opts.reffile, opts.mol)
    except ValueError as e:
        print(e)
        sys.exit(1)
    except IndexError:
        print('No reference data for the molecule(s) of interest. Exiting.')
        sys.exit(0)

    # Check that there is some overlap between the reference and computed sets
    if not (set(dcalc) & set(dref)):
        msg = 'WARNING: Reference and computed data do not overlap.\n' +\
              '         Nothing to do.'
        print(msg)
        sys.exit(0)

    # Build data
    if opts.mode == 'mean':
        if levels is None:
            print('ERROR: Levels mixing not yet available for averages.')
        if opts.compare is not None:
            if not os.path.exists(opts.compare):
                print('ERROR: Methods comparison file does not exist.')
                sys.exit()
        build_mean_data(opts.mol, labels, levels, dcalc, dref, MPL_LOADED,
                        opts.compare, 'image.pdf')
    elif opts.mode == 'dist':
        build_dist_data(opts.mol, opts.xscale, opts.absval, labels, colors,
                        levels, dcalc, dref, MPL_LOADED)


if __name__ == '__main__':
    main()
