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

from estampes.base import QLabel, QDataType
from estampes.base.errors import (
    ArgumentError, InternalError, ParsingError, QuantityError)
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

    def __init__(self, filename: str, filetype: str | None = None):
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
        self._ftype = None
        if ftype in ('fchk', 'fch'):
            self._dfile = fchk.FChkIO(filename)
            self._ftype = 'fchk'
        elif ftype in ('glog', 'log', 'out'):
            self._dfile = glog.GLogIO(filename)
            self._ftype = 'glog'
        elif ftype in ('csv', 'txt'):
            self._dfile = csv.FileCSV(filename)
            self._ftype = 'csv'
        elif ftype == 'xyz':
            self._dfile = xyz.FileXYZ(filename)
            self._ftype = 'xyz'
        else:
            raise NotImplementedError('Unsupported filetype')

    @property
    def filename(self) -> str:
        """Name of the file used to extract data."""
        return self._dfile.filename

    @property
    def version(self) -> tuple[str, tp.Any]:
        """Version of the program used to generate the file.

        The version is a tuple with:

        # The program name
        # Program/Format-specific version information
        """
        return self._dfile.full_version

    def get_data(self,
                 *qlabels: str | QLabel,
                 error_noqty: bool = True,
                 **keys4qlabels) -> QDataType:
        """Get data from Data File object.

        Wrapper to internal get_data functions.
        """
        return self._dfile.get_data(*qlabels, error_noqty=error_noqty,
                                    **keys4qlabels)

    def get_hess_data(self,
                      get_evec: bool | str = True,
                      get_eval: bool | str = True,
                      get_rmas: bool | str = False,
                      get_lweigh: bool | str = False,
                      pre_data: QDataType | None = None,
                      force_calc: bool | None = None
                      ) -> dict[str, tp.Any]:
        """Get or build Hessian data (eigenvectors and values).

        This function retrieves or builds the eigenvectors and
        eigenvalues.
        Contrary to ``get_data`` which only looks for available data,
        this functions looks for alternative forms to build necessary
        data.
        It also returns a Numpy array instead of Python lists.
        Preloaded data can be provided to avoid duplicating extraction
        queries.  They are expected in the same format as given by
        `DataFile` methods.

        The default behavior depends on the datafile:

        * for Gaussian FChk files, final data have sufficient precision
        * for Gaussian log files, final data have insufficient precision
          and must be regenerated

        For `get_` keywords, if a string is provided, it is used as key
        in the returned dictionary.

        Parameters
        ----------
        get_evec
            Return the eigenvectors (dimensionless).
        get_eval
            Return the eigenvalues.
        get_rmas
            Return the reduced masses if available.
        get_lweigh
            Return the mass-weighted eigenvectors (L) matrix, generally
        pre_data
            Database with quantities already loaded from previous queries.
        force_calc
            Force the computation of the eigenvectors and eigenvalues.

        Returns
        -------
        dict
            Dictionary containing the quantities requested.
            Possible quantities are:

            * eigenvectors (default: `lmat`)
            * eigenvalues (default: `freq`)
            * reduced masses (default: `redmas`, if available).
            * mass-weighted eigenvectors (default: `lmweigh`)

        Raises
        ------
        ValueError
            Inconsistent values given in input.
        IOError
            Error if file object not set but needed.
        IndexError
            Quantity not found.

        Notes
        -----
        * Numpy is needed to run this function
        * Data can be given in argument or will be extracted from the
          datafile object.
        * The eigenvectors are returned as (nvib, 3*natoms).
        """
        import numpy as np
        from estampes.tools.math import square_ltmat
        from estampes.tools.vib import build_vibrations, convert_hess_evec

        if not (get_evec or get_eval or get_lweigh):
            raise ArgumentError('Nothing to do')
        if get_rmas and not (get_evec or get_lweigh):
            raise ArgumentError(
                'Reduced masses only available with evec or lweigh')

        natoms = None
        nvib = None
        atmas = None
        atcrd = None
        hessvec = None
        hessval = None
        fccart = None
        key_ffx = None
        key_evec = None
        eigvec = None
        eigval = None
        redmas = None
        lmweig = None
        if force_calc is None:
            force_calc = False if self._ftype == 'fchk' else True
        # force_calc = True
        do_calc = force_calc
        if pre_data is not None:
            for key in pre_data:
                qlabel = pre_data[key].qlabel
                if qlabel.label == 'natoms':
                    natoms = pre_data['natoms'].data
                elif qlabel.label == 'nvib':
                    nvib = pre_data['nvib'].data
                elif qlabel.label == 'atmas':
                    atmas = np.array(pre_data[key].data)
                elif qlabel.label == 'atcrd':
                    atcrd = np.array(pre_data[key].data)
                elif qlabel.label == 'hessvec':
                    key_evec = key
                    hessvec = np.array(pre_data[key].data)
                elif qlabel.label == 'hessdat' and qlabel.kind == 'freq':
                    hessval = np.array(pre_data[key].data)
                elif (qlabel.label == 1 and qlabel.derord == 2
                      and qlabel.dercrd == 'X'):
                    key_ffx = key
                    fccart = np.array(pre_data[key].data)

        if force_calc and (
                fccart is not None and
                atcrd is not None and
                atmas is not None):
            if key_ffx is None or pre_data is None:
                raise InternalError('pre_data should have been properly set.')
            if isinstance(pre_data[key_ffx].shape, str):
                shape = pre_data[key_ffx].shape.lower()  # type: ignore
            else:
                raise InternalError('Undetermined shape of FFX data.')
            if (len(fccart.shape) == 1 or shape == 'lt'):
                res = build_vibrations(
                    square_ltmat(pre_data[key_ffx].data),
                    atmas, atcrd, True, bool(get_evec), bool(get_eval),
                    bool(get_rmas), bool(get_lweigh), nvib=nvib)
                if get_evec:
                    eigvec = res['evec']
                if get_eval:
                    eigval = res['freq']
                if get_rmas:
                    redmas = res['redmas']
                if get_lweigh:
                    lmweig = res['lmweigh']

        elif not force_calc:
            if get_evec and hessvec is not None:
                if key_evec is None or pre_data is None:
                    raise InternalError('pre_data[key_evec] should be set')
                # Check if eigenvectors need to be corrected
                calc_rmas = (get_rmas and
                             pre_data[key_evec].dtype == 'L.M^{-1/2}')
                eigvec = convert_hess_evec(
                    hessvec, atmas, natoms,
                    pre_data[key_evec].dtype,  # type: ignore
                    do_norm=not calc_rmas)
                if calc_rmas:
                    redmas = (eigvec**2).sum(axis=1)
                    eigvec /= np.sqrt(redmas[:, np.newaxis])
                if get_lweigh:
                    lmweig = np.array(hessvec)
            if get_eval and hessval is not None:
                eigval = hessval

        calc_evec = (get_evec or get_lweigh) and eigvec is None
        calc_eval = get_eval and eigval is None
        if calc_evec or calc_eval:
            # We are missing data, now extracting data and recomputing.
            read_data = {}
            if calc_evec or calc_eval:
                read_data['natoms'] = QLabel(quantity='natoms')
                read_data['d2EdX2'] = QLabel(quantity=1, derorder=2,
                                             dercoord='X')
                read_data['atmas'] = QLabel(quantity='atmas')
                read_data['atcrd'] = QLabel(quantity='atcrd')
                read_data['nvib'] = QLabel(quantity='nvib')
            if calc_evec:
                read_data['hessvec'] = QLabel(quantity='hessvec', level='H')
            if calc_eval:
                read_data['hessval'] = QLabel(quantity='hessdat',
                                              descriptor='freq', level='H')
            tmp_data = self.get_data(error_noqty=False, **read_data)
            if tmp_data is None:
                msg = 'Failed to extract data necessary to build normal ' \
                    + 'modes-related quantities'
                raise ParsingError(msg)
            do_calc = force_calc and tmp_data['d2EdX2'] is not None
            if not do_calc:
                if calc_evec:
                    if (tmp_data['hessvec'] is not None
                            and tmp_data['atmas'] is not None):
                        hessvec = tmp_data['hessvec'].data
                        atmas = np.array(tmp_data['atmas'].data)
                        natoms = tmp_data['natoms'].data
                        calc_rmas = (get_rmas and
                                     tmp_data['hessvec'].dtype == 'L.M^{-1/2}')
                        eigvec = convert_hess_evec(
                            hessvec, atmas, natoms,
                            tmp_data['hessvec'].dtype,  # type: ignore
                            do_norm=not calc_rmas)
                        if calc_rmas:
                            redmas = (eigvec**2).sum(axis=1)
                            eigvec /= np.sqrt(redmas[:, np.newaxis])
                        if get_lweigh:
                            lmweig = np.reshape(hessvec, (-1, 3*natoms))
                    else:
                        do_calc = True
                if calc_eval:
                    if tmp_data['hessval'] is not None:
                        eigval = np.array(tmp_data['hessval'].data)
                    else:
                        do_calc = True

            if do_calc:
                if (tmp_data['d2EdX2'] is not None
                        and tmp_data['atmas'] is not None):
                    nvib = None if tmp_data['nvib'] is None \
                            else tmp_data['nvib'].data
                    res = build_vibrations(
                        square_ltmat(tmp_data['d2EdX2'].data),
                        np.array(tmp_data['atmas'].data),
                        np.array(tmp_data['atcrd'].data),
                        True, bool(get_evec), bool(get_eval), bool(get_rmas),
                        bool(get_lweigh), nvib=nvib)
                    if get_evec:
                        eigvec = res['evec']
                    if get_eval:
                        eigval = res['freq']
                    if get_rmas:
                        redmas = res['redmas']
                    if get_lweigh:
                        lmweig = res['lmweigh']
                else:
                    if calc_evec:
                        msg = 'Unable to retrieve force constants eigenvectors'
                        raise QuantityError(msg)
                    if calc_eval:
                        msg = 'Unable to retrieve normal mode wavenumbers'
                        raise QuantityError(msg)

        result = {}
        if get_evec:
            if eigvec is None:
                raise QuantityError('Missing eigenvectors')
            key = get_evec if isinstance(get_evec, str) else 'evec'
            result[key] = eigvec
        if get_eval:
            if eigval is None:
                raise QuantityError('Missing eigenvalues')
            key = get_eval if isinstance(get_eval, str) else 'eval'
            if key in result:
                raise ArgumentError(
                    'Overlap in keywords for returned values in get_hess_data')
            result[key] = eigval
        if get_rmas:
            if redmas is None:
                raise QuantityError('Missing reduced masses')
            key = get_rmas if isinstance(get_rmas, str) else 'redmas'
            if key in result:
                raise ArgumentError(
                    'Overlap in keywords for returned values in get_hess_data')
            result[key] = redmas
        if get_lweigh:
            if lmweig is None:
                raise QuantityError('Missing mass-weighted eigenvectors')
            key = get_lweigh if isinstance(get_lweigh, str) else 'lmweigh'
            if key in result:
                raise ArgumentError(
                    'Overlap in keywords for returned values in get_hess_data')
            result[key] = lmweig

        return result
