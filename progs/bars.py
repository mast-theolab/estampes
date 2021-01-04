"""
    BARS: Benchmark on Anharmonicity - Results Survey

Simple tool to merge, compare and plot results from a set of anharmonic
  calculations.
"""

import sys
import os
import argparse
import typing as tp

import numpy as np

from estampes.base import ArgumentError, ParseKeyError, QuantityError
from estampes.parser import DataFile, build_qlabel

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
    mpl_loaded = True
    # mpl.rcParams['backend'] = 'SVG'
    # mpl.rcParams['backend.qt4'] = 'PyQt4'
    plt.rcParams['font.family'] = 'sans-serif'
except ImportError:
    mpl_loaded = False


PROGNAME = 'bars'
HELP_OPTF = '''Option file with information on input files for the benchmark.
Each line contains the following columns (separated by ";"):
1. molecule name
2. label for the functional/electronic structure method
3. label for the basis set
4. name of the file with data
5. (optional) associated color (Matplotlib-compatible)'''
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

TypeMol = tp.List[str]
TypeESM = tp.List[str]
TypeGBS = tp.List[str]
TypeInF = tp.List[str]
TypeCol = tp.List[tp.Optional[str]]
TypeOptDat = tp.Tuple[TypeMol, TypeESM, TypeGBS, TypeInF, TypeCol]
# Type for calculated values (for a given method)
TypeCalc = tp.Dict[int, tp.Dict[str, float]]
# Type for methods: [method][mode][level] = value
TypeMeth = tp.Dict[str, TypeCalc]
# Type for statistics: [method][level][mode] = value
TypeDiff = tp.Dict[str, tp.Dict[int, tp.Dict[str, float]]]
# Type for reference data: [mol][mode] = value
TypeRDat = tp.Dict[str, tp.Dict[int, float]]


def parse_args(args: tp.Sequence[str]) -> argparse.Namespace:
    """Defines and checks options given in arguments.

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
    parser.add_argument('--level',
                        type=str.lower, default='both',
                        choices=('a', 'anh', 'anharm', 'h', 'harm', 'b',
                                 'both'),
                        help='Level of theory of interest (harm/anharm/both)')
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
    """Parses option file.

    Parses Option file and returns the data as series of list:
    1. List of molecule names
    2. List of electronic structure methods/functionals
    3. List of Gaussian basis sets
    4. List of input filenames
    5. List of color codes

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
    lnam = []
    lcol = []
    for line in fobj:
        if not line.lstrip().startswith('#'):
            cols = line.split(csep)
            lcols = len(cols)
            if lcols == 5:
                lcol.append(cols[4].strip())
            elif lcols == 4:
                lcol.append(None)
            else:
                raise IndexError('At least 4 columns expected in row')
            lnam.append(cols[3].strip())
            lbas.append(cols[2].strip())
            lmet.append(cols[1].strip())
            lmol.append(cols[0].strip())

    return lmol, lmet, lbas, lnam, lcol


def parse_reffile(fobj: tp.IO,
                  mol: str) -> TypeRDat:
    """Parses reference data file.

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


def get_energies(dfname: str, levels: tp.Sequence[str]) -> TypeCalc:
    """Gets vibrational energies from data file.

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
        Harmonic and anharmonic quantities, sorted by mode index
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
    if set(levels) - {'H', 'A'}:
        raise ArgumentError('Unsupported level of theory.')
    for level in levels:
        dkeys[level] = {
            'assign': build_qlabel('vtrans', level),
            'freq': build_qlabel('vlevel', level),
        }
    keys = [item for key in dkeys for item in dkeys[key].values()]
    # Extract relevant data
    try:
        indata = dfile.get_data(*keys)
    except ParseKeyError as e:
        msg = 'Vibrational data missing' + str(e)
        raise IndexError(msg) from None

    # Check if we have the states are variational or harmonic/DVPT2
    if 'A' in levels:
        state1_nq = indata[dkeys['A']['assign']][1][1][0][1]
        if state1_nq == 0:
            return IndexError('Variational state notation found.')

    data = {}
    for level in levels:
        for i in indata[dkeys[level]['assign']]:
            fsta = indata[dkeys[level]['assign']][i][1]
            if len(fsta) == 1 and fsta[0][1] == 1:
                mode = fsta[0][0]
                if mode not in data:
                    data[mode] = {}
                data[mode][level] = indata[dkeys[level]['freq']][i]

    return data


def print_data(mol: str,
               qcmeths: tp.List[str],
               levels: tp.Sequence[str],
               dcalc: TypeMeth,
               dref: tp.Dict[int, float],
               skip_on_ref: bool = True) -> TypeDiff:
    """Prints data

    Prints data in CSV-compliant files.
    Returns the signed differences, in the form of dictionary
    [method][level][mode] = value

    Parameters
    ----------
    mol
        Molecule name (label)
    qcmeths
        "Dictionary" of all chemical available methods, in proper order
    levels
        Sequence with the levels of theory of interest
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
    fmt_fp = '{{:{}.{}f}}'.format(prec+7, prec)
    blank = (prec+7)*' '
    fmt_func = '{{:{}s}}'.format(2*(prec+7)+3)  # 3 for ' ; '
    fmt_base = ' {:4d}'
    fmt_ref = '{{:{}s}}'.format(prec+7)
    txt_head = ' Mode ; ' + fmt_ref.format('Ref.')
    txt_head2 = '      ; ' + fmt_ref.format(' ')
    fmt_delta = ' {:4s} ; ' + fmt_ref.format(' ')
    diffs = {}
    nmax = 0
    with open('res_{}.csv'.format(mol), 'w') as fobj:
        txt0 = txt_head
        txt1 = txt_head2
        for key in qcmeths:
            if key not in dcalc:  # pass if method not supported
                continue
            txt0 += ' ; ' + fmt_func.format(key)
            if max(dcalc[key]) > nmax:
                nmax = max(dcalc[key])
            diffs[key] = {}
            for lvl in levels:
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
            for key in qcmeths:
                if key not in dcalc:  # pass if method not supported
                    continue
                if i not in dcalc[key]:
                    res.extend([blank, blank])
                    continue
                for lvl in levels:
                    res.append(fmt_fp.format(dcalc[key][i][lvl]))
                    if dodiff:
                        diffs[key][lvl][i] = dcalc[key][i][lvl] - dref[i]
            if (''.join(res)).strip():
                fobj.write(fmt_base.format(i) + ' ; '.join(res) + '\n')

        # Print statistics
        res1 = []
        res2 = []
        for key in qcmeths:
            if key not in dcalc:  # pass if method not supported
                continue
            for lvl in levels:
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
                    qcmeths: tp.List[str],
                    levels: tp.Sequence[str],
                    dcalc: tp.Dict[str, TypeMeth],
                    dref: TypeRDat,
                    do_plot: bool,
                    imgfile: tp.Optional[str] = None) -> tp.NoReturn:
    """Builds (and displays) average statistical data.

    Builds average statistical data and displays if requested.
    Also handles call to the function in charge of printing the data.

    Parameters
    ----------
    op_mol
        Option for the molecule
        - all: treat all molecules together
        - each: treat each molecule separately (NYI)
        - name of the molecule of interest
    qcmeths
        "Dictionary" of all chemical available methods, in proper order
    levels
        Sequence with the levels of theory of interest
    dcalc
        Computed data as `[mol][method][mode][level] = data`
    dref
        Reference data as `[mol][mode] = data`
    do_plot
        Plot results
    imgfile
        Image file (only saved if not None)

    Raises
    ------
    NotImplementedError
        Treated separately each molecule is not yet possible
    KeyError
        Mismatch between molecules in `dref` and `dcalc`
    """
    if op_mol.lower() == 'each':
        raise NotImplementedError('Separate treatment of molecules NYI')

    nmols = len(set(dcalc) & set(dref))
    if nmols == 0:
        raise KeyError('No overlap between `dref` and `dcalc`')

    full_amax = {}
    full_dmue = {}
    full_nvib = {}
    for key in qcmeths:
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
        diffs = print_data(mol, qcmeths, levels, dcalc[mol], dref[mol])
        for key in diffs:
            for lvl in levels:
                for mode in diffs[key][lvl]:
                    val = abs(diffs[key][lvl][mode])
                    full_nvib[key][lvl] += 1
                    full_dmue[key][lvl] += val
                    if val > full_amax[key][lvl]:
                        full_amax[key][lvl] = val

    if do_plot:
        c_MAE = {'H': '#1317ff', 'A': '#bd2626'}
        c_MAX = {'H': '#6baff0', 'A': '#ff86db'}

        figid = 1
        figsize = (11, 7)
        # width2 = 0.4

        ind = np.arange(len(qcmeths))
        if len(levels) == 2:
            width = 0.18
            dmue_h = []
            dmue_a = []
            amax_h = []
            amax_a = []
            for key in qcmeths:
                dmue_h.append(full_dmue[key]['H']/full_nvib[key]['H'])
                dmue_a.append(full_dmue[key]['A']/full_nvib[key]['A'])
                amax_h.append(full_amax[key]['H'])
                amax_a.append(full_amax[key]['A'])

            fig = plt.figure(figid, figsize=figsize, dpi=300)
            fig.subplots_adjust(hspace=0.001)

            ax1 = fig.add_subplot(111)
            ax1.bar(ind-width/2, dmue_h, width, color=c_MAE['H'])
            ax1.bar(ind+width/2, amax_h, width, color=c_MAX['H'])
            ax1.bar(ind+3*width/2, dmue_a, width, color=c_MAE['A'])
            ax1.bar(ind+5*width/2, amax_a, width, color=c_MAX['A'])
            ax1.legend((r'MAE ($\omega$)', r'$|$MAX$|$ ($\omega$)',
                        r'MAE ($\nu$)',
                        r'$|$MAX$|$ ($\nu$)'))
            ax1.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
            ax1.grid(axis='x', lw=.2, ls='--', c='gray')
            ax1.set_xticks(ind+width)

            lab_meth = [item.replace(r'_', r'\_') for item in qcmeths]
            ax1.set_xticklabels(tuple(lab_meth), fontsize=5, rotation=25,
                                ha="right")
            if imgfile is not None:
                plt.savefig(imgfile, num=figid, figsize=figsize)
        else:
            lvl = levels[0]
            if lvl == 'H':
                label = r'$\omega'
            else:
                label = r'$\nu'
            width = 0.4
            dmue = []
            amax = []
            for key in qcmeths:
                dmue.append(full_dmue[key][lvl]/full_nvib[key][lvl])
                amax.append(full_amax[key][lvl])

            fig = plt.figure(figid, figsize=figsize, dpi=300)
            fig.subplots_adjust(hspace=0.001)

            ax1 = fig.add_subplot(111)
            ax1.bar(ind-width, dmue, width, color=c_MAE[lvl])
            ax1.bar(ind+width, amax, width, color=c_MAX[lvl])
            ax1.legend((r'MAE ({})'.format(label),
                        r'$|$MAX$|$ ({})'.format(label)))
            ax1.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
            ax1.grid(axis='x', lw=.2, ls='--', c='gray')
            ax1.set_xticks(ind)

            lab_meth = [item.replace(r'_', r'\_') for item in qcmeths]
            ax1.set_xticklabels(tuple(lab_meth), fontsize=5, rotation=25,
                                ha="right")
            if imgfile is not None:
                plt.savefig(imgfile, num=figid, figsize=figsize)
        plt.show()
        plt.close(figid)


def build_dist_data(op_mol: str,
                    op_xax: str,
                    absval: bool,
                    qcmeths: tp.List[str],
                    qccols: tp.Dict[str, str],
                    levels: tp.Sequence[str],
                    dcalc: tp.Dict[str, TypeMeth],
                    dref: TypeRDat,
                    do_plot: bool) -> None:
    """Builds (and displays) distribution of data.

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
        Use absolute values for the differences
    qcmeths
        "Dictionary" of all chemical available methods, in proper order
    qccols
        Dictionary of colors for each method
    levels
        Sequence with the levels of theory of interest
    dcalc
        Computed data as `[mol][method][mode][level] = data`
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
    def get_diff(d: tp.Dict[int, float], mode: int, do_abs: bool):
        """Gets and transforms differences if needed."""
        val = d.get(mode)
        if val is not None and do_abs:
            val = abs(val)
        return val

    if op_mol.lower() == 'all':
        raise NotImplementedError('Inclusive treatment of molecules NYI')

    lstyle = {'H': '--', 'A': '-'}
    figX = 16
    figY = 8
    NHA = len(levels)
    nmols = len(set(dcalc) & set(dref))
    if nmols == 0:
        raise KeyError('No overlap between `dref` and `dcalc`')

    if do_plot:
        figid = 1
        figsize = (figX*NHA, figY*nmols)
        fig = plt.figure(figid, figsize=figsize, dpi=300, tight_layout=True)
        # plt.subplots_adjust(hspace=0.001)
        pan_id = 0

    for mol in dcalc:
        if mol not in dref:
            print(f'Molecule "{mol}" not available in experiment')
        diffs = print_data(mol, qcmeths, levels, dcalc[mol], dref[mol])
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
            for ilv, lvl in enumerate(levels):
                pan_id += 1
                panel = fig.add_subplot(nmols, NHA, pan_id)
                for key in qcmeths:
                    label = key.replace(r'_', r'\_')
                    if key not in diffs:
                        continue
                    yaxis = []
                    for i, _ in xdat:
                        yaxis.append(get_diff(diffs[key][lvl], i, absval))
                    panel.plot(xaxis, yaxis, ls=lstyle[lvl], marker='d',
                               label=label)
                fmt = '{} (level: {})'
                panel.set_title(fmt.format(mol, lvl))
                if ilv == 0:
                    panel.set_ylabel(r'Error wrt Experiment (cm$^{-1}$)')
                panel.set_xlabel(xlabel)
                panel.set_xticks([int(item) for item in xaxis])
                if op_xax[:1] == 'l':
                    panel.set_xticklabels(tuple([int(item) for item in xaxis]),
                                          rotation='vertical')
                if ilv == len(levels)-1:
                    panel.legend(bbox_to_anchor=(1.02, 1), borderaxespad=0.)
    # plt.savefig(infile[0:-4]+".png",num=figid,figsize=figsize)
    if do_plot:
        plt.show()
        plt.close(figid)


def main():
    opts = parse_args(sys.argv[1:])

    # Basic parameters
    if mpl_loaded:
        if opts.fontsize:
            plt.rcParams.update({'font.size': opts.fontsize})
        # Mode-specific MPL parameters
        if opts.mode == 'dist':
            if opts.linewidth:
                plt.rcParams.update({'lines.linewidth': opts.linewidth})

    dmols, dmeth, dbset, dfnam, dcols = parse_optfile(opts.optfile)
    if opts.mol.lower() not in ('each', 'all'):
        i = 0
        while i < len(dmols):
            if dmols[i] == opts.mol:
                i += 1
            else:
                del(dmols[i])
                del(dmeth[i])
                del(dbset[i])
                del(dfnam[i])
                del(dcols[i])
    if len(dmols) == 0:
        print('ERROR: Molecule not present in option file. Exiting',
              file=sys.stderr)
        sys.exit(1)

    # Extracts computed data
    # First, check the levels
    if opts.level in ('b', 'both'):
        levels = ('H', 'A')
    else:
        levels = (opts.level[0].upper(), )
    qcmeths = []
    qccols = []
    dcalc = {}
    for i, infile in enumerate(dfnam):
        if not os.path.exists(infile):
            print(f'ERROR: File {infile} not found', file=sys.stderr)
            sys.exit()
        mol = dmols[i]
        if mol not in dcalc:
            dcalc[mol] = {}
        meth = dmeth[i] + '/' + dbset[i]
        if meth not in qcmeths:
            qcmeths.append(meth)
            qccols.append(dcols[i])
        try:
            dcalc[mol][meth] = get_energies(infile, levels)
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
        qcmeths.sort()

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
        build_mean_data(opts.mol, qcmeths, levels, dcalc, dref, mpl_loaded)
    elif opts.mode == 'dist':
        build_dist_data(opts.mol, opts.xscale, opts.absval, qcmeths, qccols,
                        levels, dcalc, dref, mpl_loaded)


if __name__ == '__main__':
    main()
