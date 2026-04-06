"""Provide basic functions common to the parsers.

This modules provides some functions redundant between parsers.
"""
import typing as tp

from estampes.base import ArgumentError, QLabel

_tp_QLab = tp.Union[str, QLabel]


def parse_qlabels(qlabs_list: tp.Sequence[_tp_QLab],
                  qlabs_dict: tp.Mapping[str, _tp_QLab]
                  ) -> tp.Union[tp.Dict[str, QLabel], tp.Dict[str, str]]:
    """Parse lists of qlabels and generate unique mapping.

    Parses lists of qlabels, as sequences or mappings, and generates a
    unique mapping "key: qlabel_object".

    Parameters
    ----------
    qlabs_list
        List of qlabels as strings (old) or QLabel objects (new format)
    qlabs_dict
        Dictionary of old and/or new qlabels, with prefedined key.

    Returns
    -------
    dict
        Mapping key: QLabel object.
    dict
        Mapping for overlapping keys.
        User-given keys in `qlab_dict` can actually refer to the same
        `qlabel` as another one.
        The original key is given in this dictionary.
    """
    def qlab2qlabel(qlab: _tp_QLab) -> tp.Tuple[str, QLabel]:
        """Convert a qlabel to the new QLabel format.

        Checks if qlab is in the old string format or the new QLabel
        instance format and returns the correct form in all cases.
        """
        if isinstance(qlab, str):
            try:
                qstr = qlab
                qobj = QLabel(qlab)
            except ArgumentError as err:
                raise ArgumentError(f'Wrong qlabel {qlab}') from err
        elif isinstance(qlabel, QLabel):
            qobj = qlab
            qstr = str(qlab)
        else:
            raise ArgumentError(f'Unrecognized format for qlabel {qlabel}')

        return qstr, qobj

    qty_dict = {}
    dupl_key = {}
    uniq_qlabs = {}
    for qlabel in qlabs_list:
        qstr, qobj = qlab2qlabel(qlabel)
        if qstr not in uniq_qlabs:
            uniq_qlabs[qstr] = qstr
            qty_dict[qstr] = qobj
        else:
            dupl_key[qstr] = uniq_qlabs[qstr]
    for key, qlabel in qlabs_dict.items():
        qstr, qobj = qlab2qlabel(qlabel)
        if qstr not in uniq_qlabs:
            uniq_qlabs[qstr] = key
            qty_dict[key] = qobj
        else:
            dupl_key[key] = uniq_qlabs[qstr]

    return qty_dict, dupl_key


def reshape_dblock(dblock: tp.Sequence[tp.Any],
                   dshape: tp.Sequence[int]
                   ) -> tp.List[tp.Any]:
    """Reshapes a data block.

    Reshapes a data block, normally extracted from fchk,
    assuming a Fortran-like array was stored in memory.
    If the shape is smaller than the length of `dblock`, the function
    tries to add a dimension, checking the consistency of this new
    dimension.

    Parameters
    ----------
    dblock
        Data block to reshape.
    dshape
        Shape of the returned data block.

        `0`/`-1` can be used for automatic definition of a dimension.

    Returns
    -------
    list
        Reshaped data block.

    Raises
    ------
    ValueError
        Raised if shape inconsistent with the size of `dblock`.
    """
    lshape = len(dshape)
    lblock = len(dblock)
    sshape = sum(dshape)
    # Check consistency between shape and data block size
    sshape = 1
    nzero = []
    for i, dim in enumerate(dshape):
        if dim <= 0:
            nzero.append(i)
        else:
            sshape *= dim
    _dshape = list(dshape)
    if nzero:
        if len(nzero) > 1:
            raise ValueError('Shape inconsistent with data block size')
        _dshape[nzero[0]] = lblock // sshape
        if _dshape[nzero[0]]*sshape != lblock:
            raise ValueError('Shape inconsistent with data block size')
    elif sshape != lblock:
        # If shape inconsistent, check if last dim was implicit
        # If so, product of other dimensions must be multiple with the
        #   full size
        new_dim = lblock // sshape
        if sshape*new_dim == lblock:
            lshape += 1
        else:
            raise ValueError('Shape inconsistent with data block size')
    # Reshape data block
    if lshape == 1:
        data = dblock
    elif lshape == 2:
        dim1 = _dshape[0]
        data = []
        for i in range(0, lblock, dim1):
            data.append([dblock[i+j] for j in range(dim1)])
    elif lshape == 3:
        dim1 = _dshape[0]
        dim2 = _dshape[1]
        dim12 = dim1*dim2
        data = []
        for i in range(0, lblock, dim12):
            data.append([])
            for j in range(0, dim12, dim1):
                data[-1].append([dblock[i+j+k] for k in range(dim1)])
    else:
        raise NotImplementedError()
    return data
