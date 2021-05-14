# Changelog

This file documents the changes to the ESTAMPES API and the associated tools.

ESTAMPES API is divided in 3 components:
- **API**: core API, in submodules `base`, `parser`.  Can work with a pure Python installation.  `NumPy` may be needed for higher-level capabilities
- **LIB**: library tools and physical databases complementing **API**.  Most capabilities are implemented on top of Python's core features.  Dependency on `NumPy` can be present when coupled with high-level core APIs.
- **VIZ**: visualization API.  It depends on visualization library, like `Qt` or `Matplotlib`

Tools have each one their own trigram
- **APP**: Estampes Suite Program (new program)
- **BLS**: _Ballast_
- **BRS**: _Bars_
- **ESP**: _ESParser_
- **GUI**: Main graphical interface

Other blocks are:
- **ETC**: Miscellaneous features
- **DOC**: Documentation


## Unreleased

### Added
- **API** - `parser/xyz`: New function `parse_xyz` to simply parse a single configuration in an opened XYZ file.

### Fixed
- **API** - Incorrect parsing of a *qlabel* made of a single keyword.
- **API** - Incorrect definition of a tuple in `parse_qlabel`.
- **API** - `parser/xyz`: Fixed a problem of reference to a private variable in a module.
- **API** - `parser/xyz`: Fixed incorrect construction of the data in `get_data`.

## [0.4.0] - 2021-05-07

### Notice
* This version introduces a number of breaking changes to the library and API.  Updates of external programs are likely necessary.  Be warned:
    * `convert_y` returns a function for the contribution of x to y.
    * `broaden` expects a function (or `None`) to include contributions from x into the definition of y.
    * The internal *qlabel* has changed: the old sub-option is now divided into a quantity-specific option and the level of theory (`H`, `A`).

### Added
- **LIB** - `convert_y` supports sub-options (e.g., incident frequencies for Raman).
- **LIB** - `convert_y` (tools/spec) now supports scaling factors in units.
- **LIB** - New function `property_units` in `data/property` to convert quantities between common unit systems (e.g., SI, cgs).
- **LIB** - New function `convert_expr` in `tools/char` to convert simple "human" mathematical expressions to Python-valid ones.
- **LIB** - New submodule `char` in `tools/` for string-related operations.

### Fixed
- **API** - Scaling factor in unit of electronic rotatory strength extracted from Gaussian log transferred to quantity itself to facilitate conversions.
- **API** - `parse_qlabel` correctly parses the sub-option of `qty_tag=2` (same as `atcrd`).
- **PY** - Changed all import aliases inside ESTAMPES to have the form `from xxx import yyy as zzz`.  This should reduce the number of problems with *PIP* installations.

### Changed
- **DOC** - Added more documentation.  The `doc/` folder has also been reorganized to be clearer.
- **LIB** - `convert_y` (`tools/spec`) now returns a function or `None* for the inclusion of x in the definition of y.  This gives more flexibility than the previous power exponent.
- **API** - Change in the definition of *qlabel* to give more flexibility: quantity-specific options and level of theory are now separated.


## [0.3.2] - 2021-04-30

### Added
- **API** - Added support of electronic dipole strengths and rotatory strengths from Gaussian log files.

### Fixed
- **API** - Fixed error in definition of electronic transition in `parse_qlabel`.


## [0.3.1] - 2021-04-23

### Added

- **ESP** - Added support of vibrational spectroscopies.
- **BRS** - Support of hybrid schemes and composite analyses mixing different level of calculation.
- **BLS** - Scaling now supports more operations (see documentation).
- **VIZ** - New function `format_label` in `visual/plotspec.py` to format unit labels for display.
- **LIB** - New `xunit` and `yunit` attributes in class `Spectrum`.
- **API** - Support of electronic energies in Gaussian log.
- **DOC** - Initial documentation for _ESParser_.

### Fixed

- **LIB** - Fixed support of multiple spectra in Gaussian log files.
- **LIB** - Fixed deprecation warning from `Matplotlib` on X/Y scale.
- **API** - Anharmonic vibrational energies correctly extracted in absence of intensities (DCPT2, HDCPT2).
- **API** - Support numbers like `2e6` as floating point numbers.

### Changed

- **BLS** - Scaling is now done *before* shifting.
- **LIB** - Class `Spectrum` can now support `DataFile` objects in the instance creation.
- **LIB** - Dictionaries for quantities related to spectroscopies separated by types (electronic, vibrational...)


## [0.3.0] - 2021-01-04

### Added

- **APP** - Added new program `Bars` (available as `esbars` through **pip**).
- **BLS** - Added possibility to deactivate the display of the figure in a new window.
- **BLS** - The final figure can be saved in BALLAST.
- **BLS** - Generated curves in BALLAST can now be saved to CSV files.
- **BLS** - Automatic inclusion of a line at Y=0 if visible positive and negative bands are detected in a spectrum.
- **LIB** - Added broadening data for VCD and IR from rotational/dipole strengths.
- **API** - Added partial support of VCD and IR spectroscopies in the `Spectrum` class.
- **API** - Added parsing of dipole and rotatory strengths, frequencies and vibrational transition levels in Gaussian log file.

### Fixed

- **ESP** - Vibronic TD spectra were incorrectly processed.
- **API** - Parsing of broadening parameters in Gaussian output files for vibronic calculations with low progressions or any WARNING message.

## [0.2.3] - 2020-10-04

### Fixed

- **API** - Fixed issues in the definition of the bounds of the X axis for the broadening in the `set_broadening` method of the `Spectrum` class.

## [0.2.2] - 2020-10-04

### Added

- **ESP** - Added support of title on figures related to vibronic spectroscopies.
- **BLS** - Added parameters for the X axis for the broadening of spectra.
- **BLS** - Geometry can be defined proportionally to the number of rows/columns.
- **BLS** - Added option to add a panel label.
- **BLS** - Ranges of subplots are supported in the definition of the curves.
- **BLS** - Support for multiple spectra in a figure.
- **VIZ** - Simple mechanism to add a panel label in an angle of a spectrum in `SpecLayout`.

### Fixed

- **BLS** - Fixed an error in the baseline shift for negative or multi-sign spectra.
- **BLS** - Reduced computational cost for curves plotted on multiple subplots.
- **API** - Unable to parse Gaussian version if the day of the release date was a single digit.
- **API** - The number of columns is supported by the constructor of `SpecLayout`, so all supported properties can be given now when creating an object.
- **API** - Fixed an error in the setup of the canvas title in class `SpecLayout`.

### Changed

- **ESP** - Tighter layout in the figure.
- **BLS** - Tighter layout in the figures.

## [0.2.1] - 2020-09-26

### Added

- **BLS** - Added the support of `linewidth` in the INI file.
- **VIZ** - Added possibility to draw a grid on matrices (`plot_jmat`, `plot_cmat`) and vectors (`plot_kvec`) to improve readability.
- **API** - The `Spectrum` class supports a linewidth parameter for the visualization aspect.

### Fixed

- **ESP** - `molview` in _ESParser_ was broken after update of the data structure.
- **VIZ** - Format when following the cursor is improved on vectors or matrices.
- **VIZ** - Fixed normal modes numbering to start at 1 for matrices and vectors.
- **API** - Wrong corrected unit for the integrated intensities of vibronic calculations in `glog.py`.
- **LIB** - Fixed issues in the broadening code.
- **LIB** - Units in `convert_y` are case-sensitive.
- **LIB** - Fixed some inconsistencies in the aliasing of the units in `convert_y`.

### Changed

- **VIZ** - The numbering for the matrices start on the bottom right to be consistent with K vector-like representations.


## [0.2.0] - 2020-09-22

### Added

- **APP** - Added new program `Ballast` (available as `esballast` through **pip**).
- **VIZ** - Added `SpecLayout` class to manipulate the layout information of spectra.
- **LIB** - Added broadening function
- **LIB** - Added Gaussian and Lorentzian distribution functions.
- **API** - Added `Spectrum` class in `estampes.base.spectrum` to handle extraction and processing of spectra-related data.
- **API** - Added simple parser for CSV files containing spectroscopic data ("XY" files).

### Changed

- **API** - Added type for color specification formats.
- **API** - Changed data structure to be always a dictionary.  The new structure is more flexible this way.
- **API** - Added quantity identifier `AnySpc` to store generic spectral data (spectrum, axis information).
- **API** - Added quantity identifier `VTrans` to store vibrational transition information.
- **API** - `FCDat:SpcLeg` and `FCDat:BShape` are merged into `FCDat:SpcPar` to facilitate the construction of correct parameters (ex: integrated intensity vs broadened intensity).


## [0.1.4] - 2020-09-21

### Fixed

- **VIZ** - Fixed errors with the definition of colors in 2D spectra.

## [0.1.3] - 2020-09-16

### Fixed

- **ESP** - Display of full J matrix in vibronic spectra.


## [0.1.2] - 2020-09-09

### Changed

- Changed versioning system from calendar to semantic.


## [0.1.1] - 2020-09-03

### Fixed

- **API** - Wrong reshaping if a leading dimension was unknown.
- **API** - Error in reshaping command in the FChk low-level parser.
- **API** - Formatted checkpoint files extensions not recognized.

## [0.1.0] - 2020-09-02

This is the first public release of the ESTAMPES library, with a proof-of-concept script.

### Added
- Initial module architecture
- **API** - Added basic parsing infrastructure: construction of `qlabels` and `DataFile` class
- **API** - Added preliminary support of Gaussian output and formatted checkpoint files.
- **LIB** - Added basic atomic database for visualization and calculation.
- **LIB** - Added basic physical constants (CODATA values) and conversion functions for common conversions (`phys_fact`).
- **LIB** - Added generator of standard data for common molecular properties (`property_data`).
- **LIB** - Added conversion tools for atomic labels/numbers
- **LIB** - Added basic conversion tools for storage units and storage structures.
- **LIB** - Added some basic mathematical functions.
- **LIB** - Added basic function to orient univocally an array of normal coordinates.
- **VIZ** - Added some low-level visualization routines for Matplotlib or Qt.
- **ESP** - Initial version of ESParser.
- **ESP** - Possibility to display the molecular geometry as a 3D figure. 
- **ESP** - Some vibronic data can be parsed and displayed
