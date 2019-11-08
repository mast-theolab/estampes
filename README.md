ESTAMPES
========

Description
-----------
A prototypical program for spectral analysis, which stands for _Experimental Support Toolbox for the Analysis, Modelling, Plotting and Elucidation of Spectra_.

The toolbox is separated in 3 components:
- Low-level modules for data parsing, which relies as much as possible on internal Python modules.
- Intermediate-level variables, functions and classes providing basic data of physical-chemical interest: constants, conversion functions, quantity properties, parameters for display...
- higher-level stand-along scripts
- the ESTAMPES gui (in development)

Parser
------
The parser submodule provides a generic interface to quantum chemical programs through a single class, `DataFile`, which handles software-specific operations and returns data in a common format.
The class provides 2 types of methods:
- `get_data`: the basic function to retrieve data, only checks if data are available.
- more evolved functions (e.g., `get_hess_data`), which perform more advanced operations and can reconstruct data
The module provides also 2 functions to build (`build_qlabel`) and parse (`parse_qlabel`) the label structure used by parsers to correctly identify the quantities to get.
The module contains software-specific submodules, which can be used directly but may make any script or program relying on them software-dependent.

More details can be found in `doc/readme.adoc`

Scripts
-------
Deriveur:
    A simple script to build displaced geometry for numerical differentiation