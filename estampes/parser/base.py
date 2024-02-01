"""General module for the computational chemistry programs.

Provides general routines and wrappers for software-independent parsing
operations.

Notes
-----
* The submodules of parsers should provide a function get_data, which
  will be used in the main class.
"""

import os
import typing as tp

from estampes.base import TypeQData
from estampes.parser import csv, xyz
from estampes.parser.gaussian import glog, fchk


# ==============
# Module Classes
# ==============


class DataFile(object):
    """Create a data file object.

    Main class to build data files to facilitate data management.
    Acts as a wrapper to the format-/program-specific wrappers.
    """

    def __init__(self, filename: str,
                 filetype: tp.Optional[str] = None):
        """Build the DataFile instance.

        Builds the DataFile instance.
        The filetype detection is automatic based on the file extension
        but can be overridden.

        Parameters
        ----------
        filename : str
            Filename.
        filetype : str, optional
            Filetype.
            Supported: 'fchk', 'glog', 'xyz'
        """
        if filetype is None:
            ftype = os.path.splitext(filename)[1][1:].lower()
        else:
            ftype = filetype.lower()
        if ftype in ('fchk', 'fch'):
            self._dfile = fchk.FChkIO(filename)
            self._parser = fchk
        elif ftype in ('glog', 'log', 'out'):
            self._dfile = glog.GLogIO(filename)
            self._parser = glog
        elif ftype in ('csv', 'txt'):
            self._dfile = csv.FileCSV(filename)
            self._parser = csv
        elif ftype == 'xyz':
            self._dfile = xyz.FileXYZ(filename)
            self._parser = xyz
        else:
            raise NotImplementedError('Unsupported filetype')

    @property
    def filename(self) -> str:
        """Name of the file used to extract data."""
        return self._dfile.filename

    @property
    def version(self) -> tp.Tuple[str, tp.Any]:
        """Version of the program used to generate the file.

        The version is a tuple with:

        # The program name
        # Program/Format-specific version information
        """
        return self._dfile.full_version

    def get_data(self,
                 *qlabels,
                 error_noqty: bool = True,
                 **keys4qlabels) -> TypeQData:
        """Get data from Data File object.

        Wrapper to internal get_data functions.
        """
        return self._parser.get_data(self._dfile, *qlabels,
                                     error_noqty=error_noqty,
                                     **keys4qlabels)

    def get_hess_data(self,
                      get_evec: bool = True,
                      get_eval: bool = True,
                      pre_data: tp.Optional[TypeQData] = None
                      ) -> tp.Tuple[tp.Any]:
        """Get or generate Hessian-related data.

        Extracts or builds the Hessian related data.
        Acts as a wrapper to parser-specific functions and
        post-processing functions.
        """
        return self._parser.get_hess_data(self._dfile, get_evec=get_evec,
                                          get_eval=get_eval, pre_data=pre_data)
