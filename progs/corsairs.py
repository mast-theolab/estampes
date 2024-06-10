"""Construct Overall Raman Scatterings Activity and Invariants from
Resonance Simulations.

Small program to build Raman and ROA tensors and invariants by explicit
    summation of contributions from multiple intermediate states.
"""

import sys
import os
import argparse

import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

from estampes.base import QLabel
from estampes.parser.base import DataFile
from estampes.base.spectro import RamanInvariants, raman_intensities, \
    RAMAN_SETUPS
from estampes.tools.spec import broaden, convert_y
from estampes.visual.plotspec import plot_spec_2D

DEBUG = False


def get_options() -> argparse.Namespace:
    """Parse and return options

    Parses options from the command line and sets default values
      otherwise.
    """
    parser = argparse.ArgumentParser()

    msg = 'Gaussian output file.  Each file is appended to the others.'
    parser.add_argument(
        '-f', '--file', action='append', dest='logfiles',
        help=msg
    )
    msg = '''Incident frequency. By default, automatically chosen.
    Possible values are:
    - auto: choose the first one available in each file (default)
    - positive integer: choose incident frequency num. X
    - real value: must correspond exactly to an existing one.
                  Use "--print-omega" to read them'''
    parser.add_argument(
        '-W', '--omega',
        help=msg
    )
    msg = 'Print all available incident frequencies, in the exact way ' \
        + 'they are given in the output file'
    parser.add_argument(
        '--print-omega', action='store_true', default=False,
        help=msg
    )
    parser.add_argument(
        '-t', '--type', type=str.upper, choices=('RR', 'RROA'),
        help='Raman setup (default:SCP(180)u)'
    )
    parser.add_argument(
        '-s', '--setup', choices=RAMAN_SETUPS['full'], default='SCP(180)',
        help='Raman setup (default:SCP(180))'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file to save the generated CSV file'
    )
    parser.add_argument(
        '--lower', type=float, default=0.0,
        help='Lower bound'
    )
    parser.add_argument(
        '--upper', type=float, default=4000.0,
        help='Upper bound'
    )
    parser.add_argument(
        '--grain', type=float, default=2.0,
        help='Discretization parameter of the X axis (δx)'
    )
    parser.add_argument(
        '--broaden', type=str.lower, choices=('gaussian', 'lorentzian'),
        default='lorentzian',
        help='Broadening function'
    )
    parser.add_argument(
        '--hwhm', type=float, default=4.0,
        help='Half-width at half-maximum'
    )
    parser.add_argument(
        '-D', '--display', action='store_true', default=False,
        help='Display the spectrum using Matplotlib'
    )
    msg = '''How to treat reduced-dimensionality cases with different number \
of vibrational states:
- cautious: only include states present in all files
- inclusive: include all states, even if only partially treated.'''
    parser.add_argument(
        '--reddim', type=str.lower, choices=('cautious', 'inclusive'),
        help=msg
    )

    opts = parser.parse_args()

    return opts


def main():
    """Main program.

    Runs the main program.
    """
    opts = get_options()
    if not opts.logfiles:
        print('Missing Gaussian output files.  Nothing to do')
        sys.exit(1)

    for logfile in opts.logfiles:
        if not os.path.exists(logfile):
            print(f'ERROR: {logfile} does not exist')
            sys.exit()

    # Define list of quantities based on spectroscopy:
    qlabs = {
        'incfrq': QLabel(quantity=1300, level='VE'),
        'alpha': QLabel(quantity=1301, level='VE'),
        'states': QLabel(quantity='vtrans', descriptor='RR'),
        'energies': QLabel(quantity='vlevel', descriptor='RR'),
        'E00': QLabel(quantity='FCDat', descriptor='E(0-0)')
    }
    qlab_rd = QLabel(quantity='FCDat', descriptor='RedDim')

    if opts.type == 'RROA':
        qlabs['Grom'] = QLabel(quantity=1302, level='VE')
        qlabs['Gcal'] = QLabel(quantity=1303, level='VE')
        qlabs['Arom'] = QLabel(quantity=1304, level='VE')
        qlabs['Acal'] = QLabel(quantity=1305, level='VE')

    # extract info
    qdata = {}
    rd_info = {}
    for logfile in opts.logfiles:
        dfile = DataFile(logfile)
        qdata[logfile] = dfile.get_data(**qlabs)
        rd_info[logfile] = dfile.get_data(rd=qlab_rd, error_noqty=False)['rd']

    # List incident frequencies
    incfrqs = []
    for logfile in opts.logfiles:
        incfrqs.append(qdata[logfile]['incfrq'].get('keys'))

    # Print list of incident frequencies and exit
    if opts.print_omega:
        for i, logfile in enumerate(opts.logfiles):
            print(f'File: {os.path.basename(logfile)}')
            print(f'-> Incident frequencies: {", ".join(incfrqs[i])}.')
        sys.exit()

    # Check that there is some overlap
    valid_omegas = set(incfrqs[0])
    for i in range(1, len(incfrqs)):
        valid_omegas &= set(incfrqs[i])
    valid_omegas = sorted(list(valid_omegas))

    if not valid_omegas:
        print(
            'ERROR: Inconsistency in incident frequencies between the files.')
        sys.exit(2)

    # Now choose the incident frequencies.
    if opts.omega is None:
        omega = valid_omegas[0]
    elif '.' in opts.omega:
        if opts.omega in valid_omegas:
            omega = opts.omega
        else:
            print('ERROR: Incident frequency not found.')
            sys.exit(2)
    else:
        # Index case: check that it is a correct integer
        # We assume that the indexes are based on the whole file, using the
        #   first as reference.
        try:
            i = int(opts.omega) - 1
            try:
                omega = qdata[opts.logfiles[0]]['incfrq'].get('keys')[i]
            except IndexError:
                print(f'ERROR: Cannot found incident frequency of index {i+1}')
            if omega not in valid_omegas:
                fmt = 'ERROR: Incident frequency found at index {} ({}) ' \
                     + ' is not available in all files'
                print(fmt.format(i+1, omega))
        except ValueError:
            print('ERROR: Unrecognized value for -W/--omega')

    do_roa = opts.type == 'RROA'

    # Extra check: check that the states are consistent or do some post-process
    # We choose the first file as reference
    reffile = opts.logfiles[0]
    states_differ = False
    for logfile in opts.logfiles[1:]:
        if (qdata[logfile]['states'].data != qdata[reffile]['states'].data
                or qdata[logfile]['energies'].data
                != qdata[reffile]['energies'].data):
            states_differ = True
            if opts.reddim is None:
                msg = '''\
ERROR: Inconsistency in state definition between files.
       Check that the molecules are the same.'''
                print(msg)
                sys.exit(3)
            else:
                msg = '''\
Warning: Inconsistency found in number and nature of states between files.
         Assuming that this is caused by reduced-dimensionality schemes.'''
                print(msg)
                if opts.reddim == 'cautious':
                    msg = '''\
States which are only partially treated will be ignored.'''
                    print(msg)
                else:
                    msg = '''\
All states will be included, even if only partially treated.'''
                    print(msg)
                break

    st_info = {}
    e_trans = {}
    nfiles = len(opts.logfiles)
    # TODO: the format could be adapted to each case with some extra processing
    fmt_nq = '{:4d}'
    fmt_st = fmt_nq+'^{:1d}'
    if states_differ:
        for ifile, logfile in enumerate(opts.logfiles):
            for key, val in qdata[logfile]['states'].data.items():
                st1, st2 = val
                if st1 == (0, 0):
                    ss1 = fmt_nq.format(0)
                elif rd_info[logfile] is None:
                    ss1 = ','.join([fmt_st.format(i, n) for i, n in st1])
                else:
                    rd = rd_info[logfile]
                    ss1 = ','.join([fmt_st.format(rd.get('state1')[i], n)
                                    for i, n in st1])
                if st2 == (0, 0):
                    ss2 = fmt_nq.format(0)
                elif rd_info[logfile] is None:
                    ss2 = ','.join([fmt_st.format(i, n) for i, n in st2])
                else:
                    rd = rd_info[logfile]
                    ss2 = ','.join([fmt_st.format(rd.get('state2')[i], n)
                                    for i, n in st2])
                key_state = f'{ss1} -> {ss2}'
                if key_state not in st_info:
                    st_info[key_state] = {
                        'ener': qdata[logfile]['energies'].data[key],
                        'info': []
                    }
                st_info[key_state]['info'].append((ifile, key))
        if opts.reddim == 'cautious':
            # We need to copy the keys to then be able to remove from dict
            all_keys = list(st_info.keys())
            for key in all_keys:
                if len(st_info[key]['info']) < nfiles:
                    del st_info[key]
    else:
        for key, val in qdata[reffile]['states'].data.items():
            st1, st2 = val
            if st1 == (0, 0):
                ss1 = fmt_nq.format(0)
            elif rd_info[reffile] is None:
                ss1 = ','.join([fmt_st.format(i, n) for i, n in st1])
            else:
                rd = rd_info[reffile]
                ss1 = ','.join([fmt_st.format(rd.get('state1')[i], n)
                                for i, n in st1])
            if st2 == (0, 0):
                ss2 = fmt_nq.format(0)
            elif rd_info[reffile] is None:
                ss2 = ','.join([fmt_st.format(i, n) for i, n in st2])
            else:
                rd = rd_info[reffile]
                ss2 = ','.join([fmt_st.format(rd.get('state2')[i], n)
                                for i, n in st2])
            key_state = f'{ss1} -> {ss2}'
            if key_state not in st_info:
                st_info[key_state] = {
                    'ener': qdata[logfile]['energies'].data[key],
                    'info': [(ifile, key) for ifile in range(nfiles)]
                }
    rr_invs = {}
    rr_int = {}

    # for key in sorted(st_info):
    #     print(key, st_info[key])

    for state in sorted(st_info):
        e_trans[state] = st_info[state]['ener']
        f_info = st_info[state]['info']
        if DEBUG:
            print(f'Transition {state}, energy: {e_trans[state]}')
        alpha = None
        indedip_mdip = None
        indedip_equa = None
        edip_indmdip = None
        edip_indequa = None
        for ifile, i in f_info:
            logfile = opts.logfiles[ifile]
            if alpha is None:
                alpha = np.array(qdata[logfile]['alpha'].data[omega][i])
                if do_roa:
                    indedip_mdip = np.array(
                        qdata[logfile]['Grom'].data[omega][i])
                    indedip_equa = np.array(
                        qdata[logfile]['Arom'].data[omega][i])
                    edip_indmdip = np.array(
                        qdata[logfile]['Gcal'].data[omega][i])
                    edip_indequa = np.array(
                        qdata[logfile]['Acal'].data[omega][i])
            else:
                alpha += np.array(qdata[logfile]['alpha'].data[omega][i])
                if do_roa:
                    indedip_mdip += np.array(
                        qdata[logfile]['Grom'].data[omega][i])
                    indedip_equa += np.array(
                        qdata[logfile]['Arom'].data[omega][i])
                    edip_indmdip += np.array(
                        qdata[logfile]['Gcal'].data[omega][i])
                    edip_indequa += np.array(
                        qdata[logfile]['Acal'].data[omega][i])

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
            if do_roa:
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
            rr_invs[state], setup=opts.setup, do_ROA=do_roa, do_FFR=False
        )

    yfact, xfunc = convert_y('ROA', 'I:cm^3/mol/sr', 'ROA:bohr^6',
                             incfreq=float(omega))
    if opts.lower >= opts.upper:
        print('ERROR: Incorrect energy range.')
    if opts.grain <= 0:
        print('ERROR: delta x must be positive.')
    if opts.hwhm <= 0:
        print('ERROR: hwhm must be strictly positive.')
    axes = {}
    axes['x'] = np.arange(opts.lower, opts.upper, opts.grain)
    xvals = list(e_trans.values())
    yvals = list(rr_int.values())
    axes['y'] = broaden(xvals, yvals, axes['x'], opts.broaden, opts.hwhm,
                        yfact, xfunc)

    if opts.output:
        with open(opts.output, 'w', encoding='utf-8') as fobj:
            for x, y in zip(axes['x'], axes['y']):
                fobj.write(f'{x:12.4f}  {y:16.8e}\n')
    else:
        print('      x               y')
        for x, y in zip(axes['x'], axes['y']):
            print(f'{x:12.4f}  {y:16.8e}')

    if opts.display:
        _, ax = plt.subplots(1, 1)
        plot_spec_2D(axes, ax, 'Wavenumbers / cm$^{-1}$',
                     'Intensity / cm$^3$ mol$^{-1}$ sr$^{-1}$')
        plt.show()


if __name__ == '__main__':
    main()
