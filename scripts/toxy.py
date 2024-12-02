"""Build and write a XY (CSV) file of a broadened spectrum.

A simple script to extract spectral data from a file and optionally
broaden it.  The new data are stored in a new file or written in the
terminal.
The script is essentially a wrapper to the Spectrum class, providing
basic parsing.  It is intended to do simple tasks, fast.
"""

import sys
import os
import typing as tp

from estampes.base.errors import ArgumentError
from estampes.parser import DataFile
from estampes.base.spectrum import Spectrum, SPEC2DATA
from estampes.base.aliases import broadening_functions, level_theory, \
    spectroscopy

PROGNAME = os.path.basename(sys.argv[0])

LAYOUT_LENGTHS = [4, 9]


def parse_opt(option: str, is_file: bool = False
              ) -> tp.Tuple[str, tp.Optional[tp.Dict[str, str]]]:
    """Split an argument-option and process optional information.

    Takes an argument from the command-line and checks if it has
    optional arguments of the form:

    option(key=value).

    It splits the information and returns it, as tuple, with the
    actual option and a dictionary of sub-options.

    For filenames, since they can contain parentheses, the function
    needs to be informed beforehand, and checks first if file exists.
    If not, it tries to split.
    The function does not formally do any check.  The test on file
    existence is just a simple measure to see if there may be options.

    Positional argument are indicated with a number (their position,
    starting from 0).

    Parameters
    ----------
    option
        Command-line argument to split.
    is_file
        If True, the option is assumed to be file-like (as a path).

    Returns
    -------
    str
        Option name
    dict
        Sub-options, as `option: value`.
    """
    if is_file:
        if os.path.exists(option):
            return option, None

    parsing = option.rsplit('(', maxsplit=1)
    if len(parsing) == 1:
        return option, None
    else:
        # We check that the last non-empty character is a closing parenthesis
        substr = parsing[1].strip()
        if not substr.endswith(')'):
            msg = f'Option "{option}" incorrectly form. \
Expected format: option(key=val)'
            raise ArgumentError(msg)
        substr = substr[:-1]
        true_opt = parsing[0]
        subopts = {}
        for i, item in enumerate(substr.split(',')):
            if '=' not in item:
                subopts[i] = item
            else:
                key, val = item.split('=', maxsplit=1)
                subopts[key] = val
        return true_opt, subopts


def main():
    '''Main function.'''
    if len(sys.argv) not in (LAYOUT_LENGTHS
                             + [i+1 for i in LAYOUT_LENGTHS]):
        print(f'''\
Usage: {PROGNAME} spec level [func hwhm lower upper grain] [outfile].

with
spec:  type of spectroscopy: IR, VCD, RS...
level: level of theory: a(nharm), h(harm), e(lectronic)...
func:  broadening function: l(orentzian), g(aussian)
hwhm:  half-width at half-maximum/half-width at half-height (in cm-1).
lower: lower bound of the X axis (in cm-1)
upper: upper bound of the X axis (in cm-1)
grain: resolution on the X axis (distance between 2 points, in cm-1)
outfile: output file

Alternatively, the broadening and axis parameters can be inverted:
{PROGNAME} spec level [lower upper grain func hwhm] [outfile].
''')
        sys.exit()

    # First argument: input file
    try:
        fname, subopts = parse_opt(sys.argv[1], is_file=True)
    except ArgumentError as err:
        print(f'Error found while parsing filename\n{err}')
        sys.exit(1)
    if not os.path.exists(fname):
        print(f'ERROR: File {fname} does not exist')
        sys.exit(1)
    if subopts is not None:
        print('Error: sub-options for filenames not yet supported.')
        sys.exit(10)
    dfile = DataFile(fname)

    # Second argument: spectroscopic technique
    # - First check if format spec(opt):
    try:
        key, subopts = parse_opt(sys.argv[2])
    except ArgumentError as err:
        print(f'Error found while parsing spectroscopy\n{err}')
        sys.exit(1)
    if subopts is not None:
        print('Error: Spectroscopy-related options NYI.')
        sys.exit(10)
    try:
        spec = spectroscopy[key.casefold()]
    except KeyError:
        print(f'ERROR: Unrecognized spectroscopy, {key}')
        print(f'Available keys: {", ".join(spectroscopy)}')
        sys.exit(1)
    if spec not in SPEC2DATA:
        print('WARNING: unsupported spectroscopy. Spectrum will be normalized')
        yunit = 'normalized'
    else:
        yunit = SPEC2DATA[spec]['unit']

    # Third argument: level of theory
    try:
        key, subopts = arg = sys.argv[3]
    except ArgumentError as err:
        print(f'Error found while parsing level of theory\n{err}')
        sys.exit(1)
    if subopts is not None:
        print('Error: Options for levels of theory NYI.')
        sys.exit(10)
    try:
        level = level_theory[key.casefold()]
    except KeyError:
        print(f'ERROR: Unrecognized level of theory, {arg}')
        print(f'Available keys: {", ".join(level_theory)}')
        sys.exit(1)

    # Check layout after 3rd argument.
    # We can check by testing if the argument is a number of not
    try:
        arg = int(sys.argv[4])
        layout = 2
    except ValueError:
        layout = 1

    # Layout 1: func, hw, lower, upper, grain
    if layout == 1:
        arg_func = sys.argv[4]
        arg_hwhm = sys.argv[5]
        arg_xmin = sys.argv[6]
        arg_xmax = sys.argv[7]
        arg_xres = sys.argv[8]
        if len(sys.argv) == 10:
            arg_fout = sys.argv[9]
        else:
            arg_fout = None
    # Layout 2: lower, upper, grain, func, hw
    elif layout == 2:
        arg_xmin = sys.argv[4]
        arg_xmax = sys.argv[5]
        arg_xres = sys.argv[6]
        arg_func = sys.argv[7]
        arg_hwhm = sys.argv[8]
        if len(sys.argv) == 10:
            arg_fout = sys.argv[9]
        else:
            arg_fout = None
    else:
        print('ERROR: Unrecognized layout.')
        sys.exit(1)

    # -- Parsing: broadening
    # Broadening function
    try:
        func = broadening_functions[arg_func.casefold()]
    except KeyError:
        print(f'ERROR: Unrecognized broadening function, {arg_func}')
        print(f'Available keys: {", ".join(broadening_functions)}')
        sys.exit(1)

    # Half-width at half-maximum
    try:
        hwhm = float(arg_hwhm)
        if hwhm <= 0.0:
            print('ERROR: The HWHM must be greater than 0.')
    except ValueError:
        print('ERROR: HWHM is expected to be a number.')
        sys.exit(1)

    # -- Parsing: axis parameters
    # Lower bound
    try:
        xmin = float(arg_xmin)
    except ValueError:
        print('ERROR: Lower bound is expected to be a number.')
        sys.exit(1)

    # Upper bound
    try:
        xmax = float(arg_xmax)
    except ValueError:
        print('ERROR: Upper bound is expected to be a number.')
        sys.exit(1)

    if xmin > xmax:
        print('NOTE: Lower and upper bounds are inverted, correcting.')
        xmin, xmax = xmax, xmin

    # Resolution/graining
    try:
        xres = float(arg_xres)
    except ValueError:
        print('ERROR: Graining/resolution is expected to be a number.')
        sys.exit(1)

    # -- Output
    if arg_fout is None:
        outfile = sys.stdout
    else:
        if os.path.exists(arg_fout):
            print(f'NOTE: File {arg_fout} already exists.  Overwritten.')
        outfile = open(arg_fout, 'w', encoding='utf-8')

    # Now process data
    spec = Spectrum(dfile, spec, level)
    spec.set_broadening(hwhm, func, yunit, xmin=xmin, xmax=xmax, xres=xres)

    fmt = '{:12.5f} {:15.6e}\n'
    for x, y in zip(spec.xaxis, spec.yaxis):
        outfile.write(fmt.format(x, y))


if __name__ == '__main__':
    main()