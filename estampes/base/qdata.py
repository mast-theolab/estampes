"""Define the Quantity Data (QData) class.

Defines the structure of the QData class, used to return information
on data extracted from files.
"""
import typing as tp
from collections.abc import Sequence

from estampes.base.errors import ArgumentError
from estampes.base.qlabel import QLabel


class QData():
    """Define the Quantity Data class.

    The Quantity Data (QData) class provides the structure to build
    objects with the necessary information to process data returned by
    ESTAMPES parsers.
    """

    def __init__(self, qlabel: tp.Union[QLabel, str]) -> None:
        """Initialize QData instance.

        Initializes the QData instance.

        Parameters
        ----------
        qlabel
            Quantity label, as a QLabel object or string (old style)
        Notes
        -----
        The construction of the QData object is intended to be incremental.
        __init__ does the minimal work, the content is then built by
        calling the general setter with the field(s) of interest.
        """
        if isinstance(qlabel, str):
            try:
                self.__qlabel = QLabel(qlabel)
            except ArgumentError as err:
                raise ArgumentError('Wrong format of the qlabel') from err
        elif isinstance(qlabel, QLabel):
            self.__qlabel = qlabel
        else:
            raise ArgumentError('Unknown format of the qlabel')
        self.__data = None
        self.__dtype = None
        self.__unit = None
        self.__shape = None
        self.__basefields = {
            'data': 'stored data',
            'unit': 'unit of the stored data',
            'dtype': 'data type, like quantity stored if variants exist',
            'shape': 'Shape of the data field (scalar, array)'
        }
        self.__extrafields = {}

    @property
    def data(self) -> tp.Any:
        """Return the data content."""
        return self.__data

    @property
    def dtype(self) -> tp.Optional[str]:
        """Return the data type."""
        return self.__dtype

    @property
    def unit(self) -> tp.Optional[str]:
        """Return the unit of the stored data."""
        return self.__unit

    @property
    def shape(self) -> tp.Optional[str]:
        """Return the shape of the stored data."""
        return self.__shape

    @property
    def qlabel(self) -> str:
        """Return the string representation of QLabel."""
        return str(self.__qlabel)

    def reset(self) -> None:
        """Reset content of QData object.

        Resets the content of the QData object and delete all extra
        fields.

        NOTE: The QLabel field is kept unchanged.
        """
        self.__data = None
        self.__dtype = None
        self.__unit = None
        self.__shape = None
        for field in self.__extrafields:
            delattr(self, field)
        self.__extrafields = {}

    def set(self, *, data: tp.Optional[tp.Any] = None,
            dtype: tp.Optional[str] = None,
            unit: tp.Optional[str] = None,
            shape: tp.Optional[str] = None,
            **fields: tp.Dict[str, tp.Any]) -> None:
        """Set content of field(s).

        Sets the content of one or more fields.

        Parameters
        ----------
        data
            Main data field.
        dtype
            Data type, as string.
        unit
            Unit of stored data.
        fields
            Extra fields
        """
        if data is not None:
            if self.__data is not None:
                raise AttributeError('Data already set. Please reset first.')
            self.__data = data

        if dtype is not None:
            if self.__dtype is not None:
                raise AttributeError(
                    'Data type already set. Please reset first.')
            self.__dtype = dtype

        if unit is not None:
            if self.__unit is not None:
                raise AttributeError('Unit already set. Please reset first.')
            self.__unit = unit

        if shape is not None:
            if self.__shape is not None:
                raise AttributeError('Shape already set. Please reset first.')
            self.__shape = shape

        if fields:
            for field, content in fields.items():
                field_ = f'__{field}'
                if getattr(self, field_) is not None:
                    raise AttributeError(
                        f'{field} already set. Please reset first.')
                setattr(self, field_, content)

    def add_field(self, field: str, *,
                  value: tp.Optional[tp.Any] = None,
                  desc: tp.Optional[str] = None) -> None:
        """Add non-standard field to QData.

        Adds a non-standard field to the QData object.
        Optionally, it is possible to set the field directly with a
        value.
        If the field exists, the operation is simply ignored.

        Parameters
        ----------
        field
            Non-standard field to add to the object.
        value
            Optional value (default: None).
        desc
            Description of the field.
        """
        if not isinstance(field, str):
            raise ArgumentError('Only string-style fields are allowed.')
        field_ = f'__{field}'
        if not hasattr(self, field_):
            self.__extrafields[field] = desc
            setattr(self, field_, value)

    def get(self,
            field: str, *,
            default: tp.Optional[tp.Any] = None) -> tp.Any:
        """Return the content of any non-standard field."""
        if default is None:
            return getattr(self, f'__{field}', None)
        else:
            return getattr(self, f'__{field}', default)

    def list_fields(self) -> tp.List[str]:
        """Return the list of all fields."""
        return self.__basefields | self.__extrafields

    def extra_fields(self) -> tp.List[str]:
        """Return the extra fields."""
        return {field: getattr(self, f'__{field}')
                for field in self.__extrafields}

    def copy(self, *, only: tp.Optional[Sequence[str]] = None,
             exclude: tp.Optional[Sequence[str]] = None) -> 'QData':
        """Return a copy of the QData object.

        Returns a complete or partial copy of the QData object.

        Parameters
        ----------
        only
            List of fields to include (all others are excluded).
        exclude
            List of fields to exclude (all others are included).

        Returns
        -------
        QData
            New partial or full copy of the object.

        Notes
        -----
        `only` and `exclude` are mutually exclusive.
        """
        if only and exclude:
            raise ArgumentError(
                'Copy() does not support both only and excluded fields.')
        if only:
            fields = []
            for field in only:
                if field not in self.list_fields():
                    raise ArgumentError(f'Unknown field to include: {field}')
                fields.append(field)
        elif exclude:
            fields = list(self.list_fields())
            for field in exclude:
                try:
                    fields.remove(field)
                except ValueError as err:
                    raise ArgumentError(f'Unknown field to exclude: {field}') \
                        from err
        else:
            fields = list(self.list_fields())
        new = QData(qlabel=self.__qlabel)
        if 'data' in fields:
            new.set(data=self.__data)
        if 'dtype' in fields:
            new.set(dtype=self.__dtype)
        if 'unit' in fields:
            new.set(unit=self.__unit)
        if 'shape' in fields:
            new.set(shape=self.__shape)
        for key, val in self.extra_fields():
            if key in fields:
                new.add_field(key, value=val)

        return new


TypeQData = tp.Dict[str, QData]
