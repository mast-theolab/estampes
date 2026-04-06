"""Construct Overall Raman Scatterings Activity and Invariants from
Resonance Simulations.

Small program to build Raman and ROA tensors and invariants by explicit
    summation of contributions from multiple intermediate states.
"""

import sys
import os
import argparse
import configparser as cfg
import typing as tp

import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

from estampes.base import QLabel
from estampes.parser.base import DataFile
from estampes.base.errors import DataError, QuantityError
from estampes.base.spectro import RamanInvariants, raman_intensities, \
    RAMAN_SETUPS
from estampes.tools.spec import broaden, convert_y
from estampes.visual.plotspec import plot_spec_2D

DEBUG = False

TEMPLATE_INI = '''\
[state:X]
# Data for contributions from state X (should be an integer > 0)
# Resonance Raman data file
rrfile =
# Ground-to-excited electronic transition data file
# Derivatives should be included as well.
excfile =
# Include or not the state in calculations
include = yes/no

[output]
# CSV files were data should be printed
csvfile =
# Display final spectrum
display = yes/no
# File where image is stored
imgfile =

[spectrum]
# Lower bound (in cm-1)
lower =
# Upper bound (in cm-1)
upper =
# Grain (delta x, in cm-1)
grain =
# Broadening function: gaussian, lorentzian
function =
# Half-width at half-maximum (in cm-1)
hwhm =

[raman]
# Do ROA instead of Raman
roa = true/false
# Raman setup
setup =
# Incident frequency, either as index or as a value in cm-1
omega =
# Damping constant, only relevant if FFR results contain imaginary term
gamma =

[transitions]
# Section on vibrational transitions to include
# Behavior to adopt in presence of red-dim calculations
# Possible values:
# - pure: only proceed if the reddim is fully consistent between RR jobs
# - cautious: only proceed with states included in all RR jobs
# - inclusive: include as many states as possible, even if partially done
reddim =
'''

OK_BROADEN = ('gaussian', 'lorentzian')
OK_SPECTRO = ('RR', 'RROA')

DEF_PARAMS = {
    'output': {
        'csvfile': None,
        'display': False,
        'imgfile': None,
    },
    'spec': {
        'lower': 0.0,
        'upper': 4000.0,
        'grain': 2.0,
        'broad': 'lorentzian',
        'hwhm': 4.0,
    },
    'raman': {
        'roa': False,
        'setup': 'SCP(180)',
        'omega': None,
        'gamma': None,
    },
    'states': {
        'reddim': None
    }
}

DEF_RRINFO = {
    'rrfile': None,
    'excfile': None,
    'state': None
}


def parse_cmdline() -> argparse.Namespace:
    """Parse command line and return options.

    Parses options from the command line and sets default values
      otherwise.
    """
    parser = argparse.ArgumentParser()

    msg = 'Gaussian output file.  Each file is appended to the others.'
    parser.add_argument(
        '-f', '--rrfile', action='append', dest='rrfiles',
        help=msg
    )
    msg = '''Incident frequency. By default, automatically chosen.
    Possible values are:
    - auto: choose the first one available in each file (default)
    - positive integer: choose incident frequency num. X
    - real value: must correspond exactly to an existing one.
                  Use "--print-omega" to read them'''
    parser.add_argument(
        '-W', '--omega', default=DEF_PARAMS['raman']['omega'],
        help=msg
    )
    msg = 'Print all available incident frequencies, in the exact way ' \
        + 'they are given in the output file'
    parser.add_argument(
        '--print-omega', action='store_true', default=False,
        help=msg
    )
    parser.add_argument(
        '-t', '--type', type=str.upper, choices=OK_SPECTRO,
        default='RROA' if DEF_PARAMS['raman']['roa'] else 'RR',
        help='Spectroscopy: RR, RROA'
    )
    parser.add_argument(
        '-s', '--setup', choices=RAMAN_SETUPS['full'],
        default=DEF_PARAMS['raman']['setup'],
        help=f'Raman setup (default: {DEF_PARAMS["raman"]["setup"]})'
    )
    parser.add_argument(
        '-o', '--output', default=DEF_PARAMS['output']['csvfile'],
        help='Output file to save the generated CSV file'
    )
    parser.add_argument(
        '--lower', type=float, default=DEF_PARAMS['spec']['lower'],
        help='Lower bound'
    )
    parser.add_argument(
        '--upper', type=float, default=DEF_PARAMS['spec']['upper'],
        help='Upper bound'
    )
    parser.add_argument(
        '--grain', type=float, default=DEF_PARAMS['spec']['grain'],
        help='Discretization parameter of the X axis (δx)'
    )
    parser.add_argument(
        '--broaden', type=str.lower, choices=OK_BROADEN,
        default=DEF_PARAMS['spec']['broad'],
        help='Broadening function'
    )
    parser.add_argument(
        '--hwhm', type=float, default=DEF_PARAMS['spec']['hwhm'],
        help='Half-width at half-maximum'
    )
    parser.add_argument(
        '-D', '--display', action='store_true',
        default=DEF_PARAMS['output']['display'],
        help='Display the spectrum using Matplotlib'
    )
    msg = '''How to treat reduced-dimensionality cases with different number \
of vibrational states:
- cautious: only include states present in all files
- inclusive: include all states, even if only partially treated.'''
    parser.add_argument(
        '--reddim', type=str.lower, choices=('cautious', 'inclusive'),
        default=DEF_PARAMS['states']['reddim'],
        help=msg
    )
    parser.add_argument(
        '-i', '--inifile',
        help='Option file. If present, all other options are ignored.'
    )
    parser.add_argument(
        '--gen-ini', action='store_true', default=False,
        help='Generate a basic option file.'
    )

    opts = parser.parse_args()

    return opts


def parse_cmdargs(opts: argparse.Namespace
                  ) -> tp.Tuple[tp.Dict[str, tp.Dict[str, tp.Any]],
                                tp.List[tp.Dict[str, tp.Any]]]:
    """Parse command line arguments.

    Parses command line arguments and returns the processed options in
    two data structures.

    - params: simulation parameters
    - rr_info: data on RR-related files.

    Parameters
    ----------
    opts
        Result of the parsing by `argparse`.

    Returns
    -------
    dict
        Dictionary with main simulation parameters.
    list
        List containing information on RR-specific computations.
    """
    # We can do a simple reference to DEF_PARAMS since there is no risk of
    #   conflict with other definitions
    params = DEF_PARAMS
    params['output']['csvfile'] = opts.output
    params['output']['display'] = opts.display
    params['spec']['lower'] = opts.lower
    params['spec']['upper'] = opts.upper
    params['spec']['grain'] = opts.grain
    params['spec']['broad'] = opts.broaden
    params['spec']['hwhm'] = opts.hwhm
    params['raman']['roa'] = opts.type == 'RROA'
    params['raman']['setup'] = opts.setup
    params['raman']['omega'] = opts.omega
    params['states']['reddim'] = opts.reddim

    rr_info = []
    for file in opts.rrfiles:
        rr_info.append(DEF_RRINFO.copy())
        rr_info[-1]['rrfile'] = file

    return params, rr_info


def parse_inifile(inifile: str
                  ) -> tp.Tuple[tp.Dict[str, tp.Dict[str, tp.Any]],
                                tp.List[tp.Dict[str, tp.Any]]]:
    """Parse option file with INI structure.

    Parses options stored in option file with INI-style structure and
    returns the processed options in two data structures.

    - params: simulation parameters
    - rr_info: data on RR-related files.

    Parameters
    ----------
    inifile
        File containing the options.

    Returns
    -------
    dict
        Dictionary with main simulation parameters.
    list
        List containing information on RR-specific computations.
    """
    if not os.path.exists(inifile):
        raise FileNotFoundError(inifile)
    opts = cfg.ConfigParser()
    opts.read(inifile)

    # Extract sections
    secs = {key.strip().lower(): key for key in opts.sections()}

    # We can do a simple reference to DEF_PARAMS since there is no risk of
    #   conflict with other definitions
    params = DEF_PARAMS

    if 'output' in secs:
        optsec = opts[secs['output']]
        eqvs = {'csvfile': 'csvfile', 'display': 'display',
                'imgfile': 'imgfile'}
        for key_ini, key_dict in eqvs.items():
            if (val := optsec.get(key_ini)) is not None:
                params['output'][key_dict] = val
    if 'spectrum' in secs:
        optsec = opts[secs['spectrum']]
        eqvs = {'lower': 'lower', 'upper': 'upper', 'grain': 'grain',
                'hwhm': 'hwhm'}
        for key_ini, key_dict in eqvs.items():
            if (val := optsec.getfloat(key_ini)) is not None:
                params['spec'][key_dict] = val
        # For broadening function, we need to do some check.
        if (val := optsec.get('function')) is not None:
            if val.lower() not in OK_BROADEN:
                raise KeyError('Unsupported broadening function.')
            else:
                params['spec']['broad'] = val.lower()
    if 'raman' in secs:
        optsec = opts[secs['raman']]
        eqvs = {'omega': 'omega', 'gamma': 'gamma',
                'ffrfile': 'ffrfile'}
        for key_ini, key_dict in eqvs.items():
            if (val := optsec.get(key_ini)) is not None:
                params['raman'][key_dict] = val
        if (val := optsec.getboolean('roa')) is not None:
            params['raman']['roa'] = val
        if (val := optsec.get('setup')) is not None:
            if val not in RAMAN_SETUPS['full']:
                raise KeyError('Unsupported Raman setup.')
            else:
                params['raman']['setup'] = val
    if 'transitions' in secs:
        optsec = opts[secs['transitions']]
        eqvs = {'reddim': 'reddim'}
        for key_ini, key_dict in eqvs.items():
            if (val := optsec.get(key_ini)) is not None:
                params['states'][key_dict] = val

    rr_info = []
    for sec in secs:
        if sec.startswith('state:'):
            optsec = opts[secs[sec]]
            if not optsec.getboolean('include', fallback=True):
                continue
            rr_info.append(DEF_RRINFO.copy())
            val = optsec.get('rrfile')
            if val is None:
                raise KeyError(f'Missing RR data file for {secs[sec]}')
            rr_info[-1]['rrfile'] = val
            rr_info[-1]['excfile'] = optsec.get('excfile')
            try:
                num = int(sec.split(':', maxsplit=1)[1])
            except ValueError as err:
                raise KeyError('Wrong state definition. Expected: state:num') \
                    from err
            rr_info[-1]['state'] = num

    return params, rr_info


def get_RR_transinfo(ref_file: str,
                     rr_data: tp.Dict[str, tp.Dict[str, DataFile]],
                     rd_data: tp.Dict[str, DataFile],
                     rd_behavior: tp.Optional[str] = None
                     ) -> tp.Dict[str, tp.Dict[str, tp.Any]]:
    """Get transition information from RR data.

    Extracts the information on the transitions stored in the data files
    from the resonance Raman calculations.
    The function also checks that the data are sufficiently consistent,
    otherwise it raises an error.

    The behavior to adopt in presence of reduced-dimension calculations
    can be tuned.

    Parameters
    ----------
    ref_file
        Reference file for comparison.
    rr_data
        Data from RR files.
    rd_data
        reduced-dimensionality data for each RR file.
    rd_behavior
        Behavior to adopt in presence of red-dim calcs.
        Supported values are:

        cautious: only consider transitions available in all files.
        inclusive: include as many transitions as possible.

    Returns
    -------
    str
        Nature of the error.
    """
    trans_info = {}

    # First check if there is some inconsistency in the data files.
    states_differ = False
    n_files = len(rr_data)
    if n_files > 1:
        for file, data in (item for item in rr_data.items()
                           if item[0] != ref_file):
            if (data['states'].data != rr_data[ref_file]['states'].data
                    or data['energies'].data
                    != rr_data[ref_file]['energies'].data):
                states_differ = True
                if rd_behavior is None:
                    msg = '''\
ERROR: Inconsistency in state definition between files.
       Check that the molecules are the same.'''
                    raise DataError(msg=msg)
                else:
                    msg = '''\
Warning: Inconsistency found in number and nature of states between files.
         Assuming that this is caused by reduced-dimensionality schemes.'''
                    print(msg)
                    if rd_behavior == 'cautious':
                        msg = '''\
States which are only partially treated will be ignored.'''
                        print(msg)
                    else:
                        msg = '''\
All states will be included, even if only partially treated.'''
                        print(msg)
                    break

    # Now process the transition information
    # /!\: the format could be adapted to each case with some extra processing
    fmt_nq = '{:4d}'
    fmt_st = fmt_nq+'^{:1d}'
    if states_differ:
        for i_file, (file, data) in enumerate(rr_data.items()):
            for key, val in data['states'].data.items():
                st1, st2 = val
                if st1 == (0, 0):
                    ss1 = fmt_nq.format(0)
                elif rd_data[file] is None:
                    ss1 = ','.join([fmt_st.format(i, n) for i, n in st1])
                else:
                    rd = rd_data[file]
                    ss1 = ','.join([fmt_st.format(rd.get('state1')[i], n)
                                    for i, n in st1])
                if st2 == (0, 0):
                    ss2 = fmt_nq.format(0)
                elif rd_data[file] is None:
                    ss2 = ','.join([fmt_st.format(i, n) for i, n in st2])
                else:
                    rd = rd_data[file]
                    ss2 = ','.join([fmt_st.format(rd.get('state2')[i], n)
                                    for i, n in st2])
                key_state = f'{ss1} -> {ss2}'
                if key_state not in trans_info:
                    trans_info[key_state] = {
                        'ener': data['energies'].data[key],
                        'info': []
                    }
                trans_info[key_state]['info'].append((i_file, key))
        if rd_behavior == 'cautious':
            # We need to copy the keys to then be able to remove from dict
            all_keys = list(trans_info.keys())
            for key in all_keys:
                if len(trans_info[key]['info']) < n_files:
                    del trans_info[key]
    else:
        for key, val in rr_data[ref_file]['states'].data.items():
            st1, st2 = val
            if st1 == (0, 0):
                ss1 = fmt_nq.format(0)
            elif rd_data[ref_file] is None:
                ss1 = ','.join([fmt_st.format(i, n) for i, n in st1])
            else:
                rd = rd_data[ref_file]
                ss1 = ','.join([fmt_st.format(rd.get('state1')[i], n)
                                for i, n in st1])
            if st2 == (0, 0):
                ss2 = fmt_nq.format(0)
            elif rd_data[ref_file] is None:
                ss2 = ','.join([fmt_st.format(i, n) for i, n in st2])
            else:
                rd = rd_data[ref_file]
                ss2 = ','.join([fmt_st.format(rd.get('state2')[i], n)
                                for i, n in st2])
            key_state = f'{ss1} -> {ss2}'
            if key_state not in trans_info:
                trans_info[key_state] = {
                    'ener': rr_data[ref_file]['energies'].data[key],
                    'info': [(ifile, key) for ifile in range(n_files)]
                }

    return trans_info


def set_RR_incfreq(ref_file: str,
                   rr_data: tp.Dict[str, tp.Dict[str, DataFile]],
                   list_omegas: tp.Sequence[tp.Sequence[str]],
                   user_omega: tp.Optional[tp.Union[str, int]] = None
                   ) -> str:
    """Set the incident frequency based on resonance Raman data.

    Sets the reference incident frequency depending on the user choice.
    If missing, the first incident frequency is used.
    Otherwise, either a reference value or the index within a list of
    incident frequencies stored in the data files can be provided.

    Parameters
    ----------
    ref_file
        Reference file for comparison.
    rr_data
        Data from RR files.
    list_omegas
        List of incident frequencies stored in each RR data file.
    user_omega
        Choice of incident frequency.

    Returns
    -------
    str
        Incident frequency to extract data.
    """

    # Check that there is some overlap
    valid_omegas = set(list_omegas[0])
    for i in range(1, len(list_omegas)):
        valid_omegas &= set(list_omegas[i])
    valid_omegas = sorted(list(valid_omegas))

    if not valid_omegas:
        msg = 'Inconsistency in incident frequencies between the files.'
        raise DataError(msg=msg)

    # Now choose the incident frequencies.
    if user_omega is None:
        omega = valid_omegas[0]
    elif '.' in user_omega:
        if user_omega in valid_omegas:
            omega = user_omega
        else:
            raise QuantityError('omega', 'Incident frequency not found.')
    else:
        # Index case: check that it is a correct integer
        # We assume that the indexes are based on the whole file, using the
        #   first as reference.
        try:
            i = int(user_omega) - 1
            try:
                omega = rr_data[ref_file]['incfrq'].get('keys')[i]
            except IndexError as err:
                raise QuantityError(
                    'omega',
                    f'Cannot found incident frequency of index {i+1}') from \
                    err
            if omega not in valid_omegas:
                raise QuantityError(
                    'omega',
                    f'Incident frequency found at index {i+1} ({omega}) '
                    + ' is not available in all files')
        except ValueError as err:
            raise KeyError('Unrecognized value for the incident frequency') \
                from err

    return omega


def main():
    """Main program.

    Runs the main program.
    """
    # Basic Checks for Options Inconsistency
    # --------------------------------------
    args = parse_cmdline()
    # Check first if user has requested read/generate an option file.
    if args.gen_ini:
        print(TEMPLATE_INI)
        sys.exit()
    elif args.inifile is not None:
        try:
            params, rr_info = parse_inifile(args.inifile)
        except FileNotFoundError:
            print(f'Option file "{args.inifile}" not found.')
            sys.exit(1)
        except KeyError as err:
            print(f'Error while parsing option file: {err}')
            sys.exit(1)
    else:
        params, rr_info = parse_cmdargs(args)

    if not rr_info:
        print('Missing Gaussian output files.  Nothing to do')
        sys.exit(1)

    for item in rr_info:
        if not os.path.exists(item['rrfile']):
            print(f'ERROR: {item["rrfile"]} does not exist')
            sys.exit(1)

    # Data Extraction for Resonance Raman
    # -----------------------------------
    # Define list of quantities based on spectroscopy:
    qlabs = {
        'incfrq': QLabel(quantity=1300, level='VE'),
        'alpha': QLabel(quantity=1301, level='VE'),
        'states': QLabel(quantity='vtrans', descriptor='RR'),
        'energies': QLabel(quantity='vlevel', descriptor='RR'),
        'E00': QLabel(quantity='FCDat', descriptor='E(0-0)'),
        'estate': QLabel(quantity='FCDat', descriptor='ExcState'),
        'geomIS': QLabel(quantity='FCDat', descriptor='GeomIS')
    }
    qlab_rd = QLabel(quantity='FCDat', descriptor='RedDim')

    if params['raman']['roa']:
        qlabs['Grom'] = QLabel(quantity=1302, level='VE')
        qlabs['Gcal'] = QLabel(quantity=1303, level='VE')
        qlabs['Arom'] = QLabel(quantity=1304, level='VE')
        qlabs['Acal'] = QLabel(quantity=1305, level='VE')

    # extract info
    rr_data = {}
    rd_data = {}
    for item in rr_info:
        dfile = DataFile(item['rrfile'])
        rr_data[item['rrfile']] = dfile.get_data(**qlabs)
        rd_data[item['rrfile']] = dfile.get_data(rd=qlab_rd,
                                                 error_noqty=False)['rd']

    # Set reference file for comparison
    reffile = rr_info[0]['rrfile']

    # List incident frequencies
    incfrqs = []
    for data in rr_data.values():
        incfrqs.append(data['incfrq'].get('keys'))

    # Print list of incident frequencies and exit
    # Only available through commandline so we directly use args.
    if args.print_omega:
        for info, incfrq in zip(rr_info, incfrqs):
            print(f'File: {os.path.basename(info['rrfile'])}')
            print(f'-> Incident frequencies: {", ".join(incfrq)}.')
        sys.exit()

    # Extraction of transition information
    try:
        trans_info = get_RR_transinfo(reffile, rr_data, rd_data,
                                      params['states']['reddim'])
    except DataError as err:
        print(err)
        sys.exit(3)

    # Setup of reference incident frequency
    try:
        omega = set_RR_incfreq(reffile, rr_data, incfrqs,
                               params['raman']['omega'])
    except (DataError, QuantityError) as err:
        print(f'ERROR: {err}')
        sys.exit(2)
    except KeyError as err:
        print(f'ERROR: {err}')
        print(' Wrong value for -W/--omega')
        sys.exit(2)

    # Processing RR Data
    # ------------------
    e_trans = {}
    rr_invs = {}
    rr_int = {}

    # for key in sorted(trans_info):
    #     print(key, trans_info[key])

    for state in sorted(trans_info):
        e_trans[state] = trans_info[state]['ener']
        f_info = trans_info[state]['info']
        if DEBUG:
            print(f'Transition {state}, energy: {e_trans[state]}')
        alpha = None
        indedip_mdip = None
        indedip_equa = None
        edip_indmdip = None
        edip_indequa = None
        for ifile, i in f_info:
            rrfile = rr_info[ifile]['rrfile']
            if alpha is None:
                alpha = np.array(rr_data[rrfile]['alpha'].data[omega][i])
                if params['raman']['roa']:
                    indedip_mdip = np.array(
                        rr_data[rrfile]['Grom'].data[omega][i])
                    indedip_equa = np.array(
                        rr_data[rrfile]['Arom'].data[omega][i])
                    edip_indmdip = np.array(
                        rr_data[rrfile]['Gcal'].data[omega][i])
                    edip_indequa = np.array(
                        rr_data[rrfile]['Acal'].data[omega][i])
            else:
                alpha += np.array(rr_data[rrfile]['alpha'].data[omega][i])
                if params['raman']['roa']:
                    indedip_mdip += np.array(
                        rr_data[rrfile]['Grom'].data[omega][i])
                    indedip_equa += np.array(
                        rr_data[rrfile]['Arom'].data[omega][i])
                    edip_indmdip += np.array(
                        rr_data[rrfile]['Gcal'].data[omega][i])
                    edip_indequa += np.array(
                        rr_data[rrfile]['Acal'].data[omega][i])

        rr_invs[state] = RamanInvariants(
            alpha=alpha,
            indedip_mdip=indedip_mdip,
            indedip_equa=indedip_equa,
            edip_indmdip=edip_indmdip,
            edip_indequa=edip_indequa,
            E_inc=float(omega),
            E_scat=float(omega)-e_trans[state])
        if DEBUG:
            print(f' alpha^2 = {rr_invs[state].alpha2}')
            print(f' beta_s(alpha)^2 = {rr_invs[state].b_s_alph}')
            print(f' beta_a(alpha)^2 = {rr_invs[state].b_a_alph}')
            if params['raman']['roa']:
                print(f' alpha.G{{roman}} = {rr_invs[state].alp_romG}')
                print(f' beta_s(G{{roman}}) = {rr_invs[state].b_s_romG}')
                print(f' beta_a(G{{roman}}) = {rr_invs[state].b_a_romG}')
                print(f' alpha.G{{script}} = {rr_invs[state].alp_calG}')
                print(f' beta_s(G{{script}}) = {rr_invs[state].b_s_calG}')
                print(f' beta_a(G{{script}}) = {rr_invs[state].b_a_calG}')
                print(f' beta_s(A{{roman}}) = {rr_invs[state].b_s_romA}')
                print(f' beta_a(A{{roman}}) = {rr_invs[state].b_a_romA}')
                print(f' beta_s(A{{script}}) = {rr_invs[state].b_s_calA}')
                print(f' beta_a(A{{script}}) = {rr_invs[state].b_a_calA}')

        rr_int[state] = raman_intensities(
            rr_invs[state], setup=params['raman']['setup'],
            do_ROA=params['raman']['roa'], do_FFR=False
        )

    yfact, xfunc, ycorr = convert_y('ROA', 'I:cm^3/mol/sr', 'ROA:bohr^6',
                                    incfreq=float(omega))
    if params['spec']['lower'] >= params['spec']['upper']:
        print('ERROR: Incorrect energy range.')
    if params['spec']['grain'] <= 0:
        print('ERROR: delta x must be positive.')
    if params['spec']['hwhm'] <= 0:
        print('ERROR: hwhm must be strictly positive.')
    axes = {}
    axes['x'] = np.arange(params['spec']['lower'], params['spec']['upper'],
                          params['spec']['grain'])
    xvals = list(e_trans.values())
    yvals = list(rr_int.values())
    axes['y'] = broaden(xvals, yvals, axes['x'], params['spec']['broad'],
                        params['spec']['hwhm'], yfact, xfunc, ycorr)

    if params['output']['csvfile']:
        with open(params['output']['csvfile'], 'w', encoding='utf-8') as fobj:
            for x, y in zip(axes['x'], axes['y']):
                fobj.write(f'{x:12.4f}  {y:16.8e}\n')
    else:
        print('      x               y')
        for x, y in zip(axes['x'], axes['y']):
            print(f'{x:12.4f}  {y:16.8e}')

    if params['output']['display']:
        _, ax = plt.subplots(1, 1)
        plot_spec_2D(axes, ax, 'Wavenumbers / cm$^{-1}$',
                     'Intensity / cm$^3$ mol$^{-1}$ sr$^{-1}$')
        plt.show()


if __name__ == '__main__':
    main()
