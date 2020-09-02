"""Module providing basic error classes

A basic module providing exceptions for ESTAMPES.

Classes
-------
ArgumentError
    Generates an error for inconsistency/errors in arguments.
ParsingError
    Basic container for parsing-specific errors.
ParseDataError
    Generates an error if quantity/property not available/supported.
ParseKeyError
    Generates an error if keyword not found.
QuantityError
    Generates an error if quantity is not supported.
"""

import typing as tp


# ==============
# Module Classes
# ==============

class ArgumentError(Exception):
    """Generates an error for inconsistency/errors in arguments.

    Generates an error if inconsistency/errors are found in arguments in
      calls to functions/methods.

    Parameters
    ----------
    name : str
        Name of the argument
    msg : str, optional
        Message to be printed instead of default one
    """
    def __init__(self, name: str, msg: tp.Optional[str] = None) -> None:
        if msg is None:
            msg = f'Error in argument: {name}'
        super(ArgumentError, self).__init__(msg)


class ParsingError(Exception):
    """Basic container for parsing-specific errors."""


class ParseDataError(ParsingError):
    """Generates an error if quantity/property not available/supported.

    Generates an error if a quantity is not available or not yet supported.

    Parameters
    ----------
    name : str
        Name of the quantity to search
    prog : str, optional
        Name of the program
    msg : str, optional
        Message to be printed instead of default one
    """
    def __init__(self, name: str, msg: tp.Optional[str] = None,
                 prog: str = 'Gaussian') -> None:
        if msg is None:
            msg = f'Unsupported quantity in {prog} file: {name}'
        super(ParseDataError, self).__init__(msg)


class ParseKeyError(ParsingError):
    """Generates an error if keyword not found.

    Generates an error if a keyword is not found.

    Parameters
    ----------
    key : str
        Name of the quantity to search
    msg : str, optional
        Message to be printed instead of default one
    prog : str, optional
        Name of the program
    """
    def __init__(self, key: str, msg: tp.Optional[str] = None,
                 prog: str = 'Gaussian') -> None:
        if msg is None:
            msg = f'Keyword not found in {prog} file: {key}'
        super(ParseKeyError, self).__init__(msg)


class QuantityError(Exception):
    """Generates an error if quantity is not supported.

    Generates an error if a quantity is not supported either as a whole
      or in specific conditions.

    Parameters
    ----------
    name : str
        Name of the quantity.
    msg : str, optional
        Message to be printed instead of default one.
    """
    def __init__(self, name: str, msg: tp.Optional[str] = None) -> None:
        if msg is None:
            msg = f'Unsupported quantity: {name}'
        super(QuantityError, self).__init__(msg)
