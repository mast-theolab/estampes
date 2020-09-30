# Changelog

This file documents the changes to the ESTAMPES API and the associated tools.

ESTAMPES API is divided in 3 components:
- **API**: core API, in submodules `base`, `parser`.  Can work with a pure Python installation.  `NumPy` may be needed for higher-level capabilities
- **LIB**: library tools and physical databases complementing **API**.  Most capabilities are implemented on top of Python's core features.  Dependency on `NumPy` can be present when coupled with high-level core APIs.
- **VIZ**: visualization API.  It depends on visualization library, like `Qt` or `Matplotlib`

Tools have each one their own trigram
- **APP**: Estampes Suite Program (new program)
- **BLS**: _Ballast_
- **ESP**: _ESParser_
- **GUI**: Main graphical interface

Other blocks are:
- **ETC**: Miscellaneous features
- **DOC**: Documentation

## Unreleased


### Added

- **BLS** - Ranges of subplots are supported in the definition of the curves.
- **BLS** - Support for multiple spectra in a figure.

### Fixed

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
