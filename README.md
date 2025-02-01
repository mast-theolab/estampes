ESTAMPES
========

Description
-----------
A prototypical program for spectral analysis, which stands for _Experimental Support Toolkit for Analysis, Modelling, Processing and Elaboration of Spectra_.

The toolkit is separated in 3 components:
* Low-level modules for data parsing, which relies as much as possible on internal Python modules.
* Intermediate-level variables, functions and classes providing basic data of physical-chemical interest: constants, conversion functions, quantity properties, parameters for display...
* higher-level stand-along scripts

Some proof-of-concepts, working scripts are also provided in the package.

Installation
------------
The version is not yet available on *PyPI* and must be installed manually in root directory of the package:

> **Rolling release**
> 1. Clone the repository
> 2. In the local repository, run:
> ```
> pip install -e .
> ```


> **Directly with pip**
>
> Run
> ```
> pip install -e git+https://github.com/jbloino/estampes.git#egg=estampes
> ```

Requirements
------------
* Low-level modules do not require modules outside of the standard Python installation, to facilitate data parsing
* Visualization tools require additional tools
    * 3D visualization (`visual.molview`): `PySide6`, `NumPy`
    * Plot charting (`visual.plotmat`, `visual.plotspec`): `NumPy`
* Scripts:
    * **ESParser**: `NumPy`, `PySide6` (optional)

Parser
------
The parser submodule provides a generic interface to quantum chemical programs through a single class, `DataFile`, which handles software-specific operations and returns data in a common format.
The class provides 2 types of methods:
* `get_data`: the basic function to retrieve data, only checks if data are available.
* more evolved functions (e.g., `get_hess_data`), which perform more advanced operations and can reconstruct data
The module provides also 2 functions to build (`build_qlabel`) and parse (`parse_qlabel`) the label structure used by parsers to correctly identify the quantities to get.
The module contains software-specific submodules, which can be used directly but may make any script or program relying on them software-dependent.

More details can be found in `doc/readme.adoc`

Scripts
-------
**ESParser**:
    A simple script to parse and plot data related to vibrational, vibronic or electronic spectroscopy, or to display basic molecular information.