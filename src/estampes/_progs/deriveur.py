"""Program DERIVEUR.

DERIVEUR: A Simple Deriving Program

This simple program builds the necessary steps for the numerical
differentiation of chosen quantities.
"""

import os
import sys
import argparse
from collections.abc import Sequence

from estampes.data.physics import PHYSFACT as phfact

from estampes.extras.derive_core import (
    HELP_FMT, HELP_QTY, HELP_TMPL_FNAME, HELP_TMPL_MAIN, HELP_TMPL_OPT,
    KEYS_COORDS, LIM_STEP,
    check_fname_pattern)
from estampes.extras.derive_diff import main_diff
from estampes.extras.derive_shift import main_shift

PROGNAME = "deriveur"


########################################################################
#                                                                      #
#                            USER OPTIONS                              #
#                                                                      #
########################################################################
def parse_cmdargs(options: Sequence[str]) -> argparse.Namespace:
    """Define and check options given in arguments.

    Defines available options in the program and parse a list of
    arguments given in input, and returns the results as a populated
    namespace.

    Parameters
    ----------
    options
        List of options/arguments to parse

    Returns
    -------
    :obj:`argparse.Namespace`
        Returned results from the parsing
    """
    parser = argparse.ArgumentParser(
        prog=PROGNAME,
        formatter_class=argparse.RawTextHelpFormatter)
    msg = f'''Differentiation coordinates.  Possible choices are:
- mass-weighted normal coords: {', '.join(KEYS_COORDS['Q'])}
- reduced normal coords: {', '.join(KEYS_COORDS['q'])}
- Cartesian coords: {', '.join(KEYS_COORDS['X'])}'''
    parser.add_argument('-c', '--coord',
                        type=lambda x: x.lower() if len(x) > 2 else x,
                        choices=[item for enum in KEYS_COORDS.values()
                                 for item in enum],
                        default='mode',
                        metavar='coords',
                        help=msg)
    parser.add_argument('-H', '--fullhelp', action='store_true', default=False,
                        help='Print the full help, with the sub-parsers.')
    parser.add_argument('-g', '--geomfile',
                        help='File with reference geometry. The structure from'
                             ' the reference file will be superposed to it.')
    parser.add_argument('-i', '--indexes',
                        help='Index or list of indexes, separated by commas '
                        + 'or dashes to specify ranges.')
    msg = f'''Filenames pattern for displaced geometries. It can contain:
{HELP_TMPL_FNAME}'''
    parser.add_argument('-p', '--pattern', help=msg)
    msg = 'Value of the differentiation step.  ' \
        + 'The unit should be consistent "--unit".'
    parser.add_argument('-s', '--step', help=msg)
    parser.add_argument('--symm', action='store_true',
                        help='Use symmetry information if available.')
    parser.add_argument('-u', '--unit', type=str.lower, default='ang',
                        choices=('ang', 'aa', 'au', 'a0', 'bohr'),
                        help='Unit for the coordinates.')

    psub = parser.add_subparsers(help='Operating modes of deriveur',
                                 required=True)

    # ------------------------------------------------------------------
    pbuild = psub.add_parser('shift',
                             help='Build derivation points from shifted '
                             + 'geometries',
                             formatter_class=argparse.RawTextHelpFormatter)
    pbuild.set_defaults(mode='shift')
    msg = 'Carry out displacement over several points along derivation ' \
        + 'coordinate.  The steps are read "symmetrically", so ' \
        + '2*multi-steps  will be built.'
    pbuild.add_argument('--multi-steps', type=int, default=1, help=msg)
    msg = f'''\
Template file used as basis to generate the files with displaced geometries.
The file should at least contain 2 keywords:
{HELP_TMPL_MAIN}
Other supported keywords are:
{HELP_TMPL_OPT}
In addition, the file can include the keywords supported for patterns:
{HELP_TMPL_FNAME}
'''
    pbuild.add_argument('-t', '--template', help=msg)
    pbuild.add_argument('reffile',
                        help='Reference file to build derivation points')

    # ------------------------------------------------------------------
    pcalc = psub.add_parser('diff', help='Carry out numerical differentiation',
                            formatter_class=argparse.RawTextHelpFormatter)
    pcalc.set_defaults(mode='diff')
    pcalc.add_argument('reffile',
                       help='File with data for reference geometry')
    # pcalc.add_argument('posfile',
    #                    help='File with data corresponding to positive step')
    # pcalc.add_argument('negfile',
    #                    help='File with data corresponding to negative step')
    pcalc.add_argument('-f', '--format', type=str.lower,
                       choices=('indatanm', 'indatax'), default='indatanm',
                       help=HELP_FMT)
    msg = '''Force orientation of normal modes by multiplying by +/-1.
This should be used primarily for debugging.
Numbers are expected as -1/(+)1, separated by commas, the value is ignored.'''
    pcalc.add_argument('--force-nmorient', help=msg)
    msg = 'Deactivate the correction of the normal modes orientation (for ' \
        + 'backward compatibility if necessary)'
    pcalc.add_argument('--no-nmorient', action='store_true',
                       help=msg)
    pcalc.add_argument('-L', '--Lmatfile',
                       help='File to use for the definition of L, '
                            'transformation matrix from Cartesian to normal '
                            'coordinates')
    pcalc.add_argument('-o', '--output',
                       help='Save output in given file instead of terminal')
    pcalc.add_argument('-q', '--quantity', help=HELP_QTY)
    pcalc.add_argument('--refheader', action='store_true',
                       help='Print reference header')
    # Add keyword to support spectroscopy instead to generate automatically the
    # quantities to parse.

    if {'-H', '--fullhelp'} & set(sys.argv):
        print('''
 Main Options
 ------------''')
        parser.print_help()
        print('''
 Submodule - Shift (build displaced geometries)
 ----------------------------------------------''')
        pbuild.print_help()
        print('''
 Submodule - Derivation (build numerical differentiations)
 ---------------------------------------------------------''')
        pcalc.print_help()
        sys.exit()

    opts = parser.parse_args(options)

    return opts


########################################################################
#                                                                      #
#                                 MAIN                                 #
#                                                                      #
########################################################################

def main():
    """Run the main program."""
    args = parse_cmdargs(sys.argv[1:])
    # Definition of displacement coordinates
    if args.coord in KEYS_COORDS['Q']:
        dcoord = 'Q'
    elif args.coord in KEYS_COORDS['X']:
        dcoord = 'X'
    elif args.coord in KEYS_COORDS['q']:
        dcoord = 'q'
        print('ERROR: Reduced normal coordinates not yet implemented')
        sys.exit()
    else:
        print('ERROR: Coordinate types unsupported')
        sys.exit()
    # Check unit for distances
    if args.unit in ('au', 'a0', 'bohr'):
        au2unit = 1.0
        msg = '''\
NOTE: "Units=AU" must be provided in the Route section of Gaussian
      for the geometries to be correctly read.
'''
        print(msg)
    else:
        au2unit = phfact.bohr2ang
    # Definition of step size
    if dcoord == 'Q':
        dstep = (0.01/phfact.bohr2ang if args.step is None
                 else args.step/au2unit)
    elif dcoord == 'X':
        dstep = (0.001/phfact.bohr2ang if args.step is None
                 else args.step/au2unit)
    elif dcoord == 'q':
        dstep = (0.01/phfact.bohr2ang if args.step is None
                 else args.step/au2unit)
    else:
        print('Internal error: displacement step not defined')
        sys.exit(9)
    if dstep <= LIM_STEP:
        print('ERROR: step too small')
        sys.exit(1)
    if args.geomfile is not None:
        if not os.path.exists(args.geomfile):
            print(f'ERROR: File "{args.geomfile}" not found')
            sys.exit(1)
        # fobj = DataFile(args.geomfile)
        # qlabs = {'atcrd': QLabel('atcrd')}
        # ctarget = np.array(
        #     fobj.get_data(**qlabs)['atcrd'].data).reshape(-1, 3)
        print('Not yet implemented!')
        sys.exit()

    # Check reference file (true for both builder and differentiator)
    if not os.path.exists(args.reffile):
        print(f'ERROR: File "{args.reffile}" not found')
        sys.exit(1)
    elif os.path.splitext(args.reffile)[1].lower() not in ('.fchk', 'fch'):
        print('Only formatted checkpoint files are supported for the '
              + 'reference file for now')
        sys.exit(9)

    # Modes
    if args.mode == 'shift':
        if args.multi_steps <= 0:
            print('ERROR: Number of steps must be strictly positive.')
            sys.exit(1)
        if args.multi_steps >= 10:
            print('WARNING: A large number of steps will be built.')

        if args.template is not None:
            if not os.path.exists(args.template):
                print(f'ERROR: Template file {args.template} does not exist.')
                sys.exit(1)

        main_shift(args.reffile, args.template, dcoord, args.multi_steps,
                   dstep, au2unit, args)

    elif args.mode == 'diff':
        if args.pattern is None:
            print('ERROR: A file pattern is needed to build the differences')
            sys.exit(2)
        tmpl_file = check_fname_pattern(args.pattern, dcoord == 'X')
        main_diff(args.reffile, tmpl_file, args.quantity, dcoord, dstep, args)
    else:
        print('ERROR: Unsupported mode')
        sys.exit()


if __name__ == '__main__':
    main()
