"""Program: SOAR.

SOAR: Structural Overlap Analysis and Representation.

Simple program to analyze differences in structures.
"""
import argparse
import typing as tp

import numpy as np

import matplotlib.pyplot as plt

from estampes.base import QLabel
from estampes.parser import DataFile
from estampes.tools.vib import build_dusch_J, build_dusch_K
from estampes.visual.plotmat import plot_jmat, plot_kvec

PROGNAME = 'soar.py'


def parse_args(args: tp.Optional[tp.Sequence[str]] = None
               ) -> argparse.Namespace:
    """Define and check options given in arguments.

    Defines available options in the program, parses a list of arguments
    given in input, and returns the results as a populated namespace.

    Parameters
    ----------
    args
        List of options/arguments to parse.

    Returns
    -------
    obj:`argparse.Namespace`
        Returned results from the parsing.
    """
    parser = argparse.ArgumentParser(
        prog=PROGNAME,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('reffile',
                        help='File with data from the reference structure.')
    parser.add_argument('newfile',
                        help='File with data from the other structure.')

    if args is not None:
        opts = parser.parse_args(args)
    else:
        opts = parser.parse_args()

    return opts


def main():
    """Run the main program."""
    dopts = parse_args()

    dfile_ref = DataFile(dopts.reffile)
    dfile_new = DataFile(dopts.newfile)

    keys = {
        # 'E': QLabel(quantity=1),
        'mass': QLabel(quantity='AtMas'),
        'crd': QLabel(quantity='AtCrd')
    }
    Lmat_ref, freq_ref = dfile_ref.get_hess_data()
    Lmat_new, freq_new = dfile_new.get_hess_data()

    data_ref = dfile_ref.get_data(**keys)
    data_new = dfile_new.get_data(**keys)

    jmat = build_dusch_J(Lmat_new, Lmat_ref)
    kvec = build_dusch_K(Lmat_ref, np.array(data_ref['mass'].data),
                         np.array(data_new['crd'].data)
                         - np.array(data_ref['crd'].data))

    figsize = (10, 8)
    fig, subp = plt.subplots(2, 1, tight_layout=True)
    fig.set_size_inches(figsize)
    _ = plot_jmat(jmat, subp[0])
    _ = plot_kvec(kvec, subp[1])
    # fig.colorbar(plot)

    print('''\
 num. | Ref. freq. |  New freq.
 -----+------------+-----------''')
    i = 0
    for x, y in zip(freq_ref, freq_new):
        i += 1
        print(' {:4d} | {:10.4f} | {:10.4f}'.format(i, x, y))

    plt.show()


if __name__ == '__main__':
    main()
