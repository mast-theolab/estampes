# Changelog

This file documents the changes to the ESTAMPES API and the associated tools.

ESTAMPES API is divided in 3 components:
- **API**: core API, in submodules `base`, `parser`.  Can work with a pure Python installation.  `NumPy` may be needed for higher-level capabilities
- **LIB**: library tools and physical databases complementing **API**.  Most capabilities are implemented on top of Python's core features.  Dependency on `NumPy` can be present when coupled with high-level core APIs.
- **VIZ**: visualization API.  It depends on visualization library, like `Qt` or `Matplotlib`

Tools have each one their own trigram
- **APP**: Estampes Suite Program (new program)\
- **2XY**: _ToXY_
- **BLS**: _Ballast_
- **BRS**: _Bars_
- **CRS**: _Corsairs_
- **ESP**: _ESParser_
- **MRG**: _Mirage_
- **OAR**: _Soar_
- **GUI**: Main graphical interface

Other blocks are:
- **ETC**: Miscellaneous features
- **DOC**: Documentation
- **DEV**: Development conventions/cleanup


## Unreleased

### Added
- **MRG**: New option `--vibs` provides a list of normal modes, providing a way to look at different vibrations in a single run.
- **VIZ**: `MolUI` can now export XYZ structures, considering orientation, translation and/or vibrational displacements.
- **VIZ**: New shortcut 'CTRL+G' (CMD+G on Mac) to export the XYZ structures.
- **VIZ**: New getter for atomic coordinates can return the atomic coordinates or displaced atomic coordinates with `get_at_crd(True)`.
- **VIZ**: `MolUI` can now update molecular vibrations during animation.
- **VIZ**: `Molecule.redraw()` can redraw the molecule from the stored geometries.
- **VIZ**: `Molecule.animation_status` provides a getter/setter on the animation status (running, paused, stopped...).

### Fixed
- **VIZ**: Fixed warning message about unset label when moving the cursor out of the bar plot in `plot_kvec`.
- **VIZ**: Fixed warning message "invalid value encountered in divide" when the vibration to display in `visual.VibView` or `visual.povrender` had a null component on one atom.
- **LIB**: `DataFile.get_hess_data` failed to return the eigenvectors/eigenvalues from a Gaussian log file if the Cartesian force constants were missing.
- **LIB**: `DataFile.get_hess_data` failed to read vibrational energies when generated from the Cartesian force constants through `tools.vib.build_vibrations`.
- **API**: Fixed broken extraction of non-transition quantities in excited states.

### Changed
- **MRG**: The list of frequencies cannot be edited anymore and remains more readable when the window is resized.

## [0.6.0] - 2025-02-16

### Important
1. The logic in the POV-Ray description scene has been changed in this version to better match the interactive visualizer.  This means that the Z axis is now the vertical and Y points toward the camera, while the old previous, standard form was to have Y upwards and Z pointing in the direction opposite to the camera.  Note that to keep operations more intuitive, the left-handed representation is fully maintained.

### Added
- **2XY**: New simple script to generate convoluted band-shapes or extract spectral data from files supported by ESTAMPES.
- **BLS**: Added option `RamanSetup` in the curve specifications to choose the Raman/ROA setup.
- **BLS**: Added option `RamanLaser` in the curve specifications to choose the Raman/ROA incident frequencies.  Alternatively, `RamanWInc` can be used.  The incident frequencies is expected in cm-1.
- **BLS**: A filetype can now be provided to the file parser if the file extensions is ambiguous.
- **BLS**: Ballast can now support multiple files in a single curve, for instance to produce the spectrum of mixtures.  The format is "file1 @ weight1 & file2 @ weight2" (spaces are irrelevant between symbols).
- **MRG**: Mirage was updated to support the new representations of normal modes.
- **MRG**: Mirage now recognizes colors specifications as sub-options of `--vid-model`.
- **VIZ**: A summary of the commands available in the molecular viewer with `MolWin` is now printed by default (can be deactivated with parameter `skip_guide=True` given to constructor).  It can be recalled at any moment with 'CTRL+SHIFT+H' (CMD+SHIFT+H on Mac).
- **VIZ**: `MolWin` can now animate vibrations with 'CTRL+A' (CMD+A on Mac). The animation can be paused and resumed.
- **VIZ**: New method of `Molecule` to control animated atomic displacements, "animate".  While intended for vibrations, it can be used for any kind of displacements.
- **VIZ**: `MolWin` can now export screenshots with 'CTRL+P' (CMD+P on Mac). The export is still experimental and can give offset colors.
- **VIZ**: `MolWin` now supports a new shortcut 'CTRL+E' (COMMAND+E on Mac) to export a POV-Ray file describing the represented in the same orientation (currently translations are not considered).
- **VIZ**: `Molecule`, `VibMode` and `POVBuilder` now support explicit description of the representation model and material.  Old keywords like `rad_atom_as_bond` are deprecated and have been removed.
- **VIZ**: `Molecule` now supports other "materials" (WIP).
- **VIZ**: The description scene in POV-Ray has been changed to match the modified convention.  This way, the interactive viewer and POV-Ray describe by default the same scene.
- **VIZ**: `Molecule` and `VibView` have a new method to retrieve visualization parameters, `viz_options`.
- **VIZ**: Added two new representations of normal modes: as arrows centered on the atoms (midarrows), as dual arrows showing the "positive" and "negative" displacements (dualarrows).
- **LIB**: New module `base.aliases` to provide common aliases for keywords (e.g., spectroscopy, level of theory) regularly used within ESTAMPES.
- **LIB**: The `Spectrum` class can now support multiple datasets and relative weights (for now purely numerical).
- **LIB**: New function `tools.vib.build_vibrations` constructs the normal coordinates vector and frequencies Cartesian force constants.  The function is inspired by the white paper of Gaussian (https://gaussian.com/vib/), transcribed in Python by M. Fusè.
- **LIB**: New function `tools.vib.convert_hess_evec` can convert mass-weighted eigenvectors to dimensionless eigenvectors and fix the shape of the matrix.
- **LIB**: New function `tools.vib.norm_evec` normalizes each eigenvector.
- **LIB**: New dictionary with X11 RGB colors names available in `data.colors`.
- **LIB**: New module to store color names: `data.colors`.  It currently supports X11, XKCD and the default color palette of Matplotlib.
- **LIB**: New function to parse an option value (for instance in `--arg=value`) of the form "value(option,key=option)".
- **LIB**: New alias `COLOR_NAMES` containing all RGB color dictionaries dictionaries in `data.colors`.
- **LIB**: New conversion function from a color specification to RGB tuple, `data.colors.to_rgb_list`.
- **LIB**: Default colors of vibrational modes are now stored in `data.visuals.VIBCOLS`.
- **LIB**: `visual.povrender.write_pov_vib` now fully supports color specifications.
- **LIB**: `visual.vibview.VidMode` now fully supports color specifications.
- **API**: Multi-jobs (explicit or internal linked jobs) in Gaussian log files are now supported.  ESTAMPES assumes that the jobs are all related to the same task (for instance opt+freq) and does not support selective extractions from a specific job (e.g., linked freq jobs with different levels of theory) for now.
- **API**: Added support of anharmonic X matrix from Gaussian log files.

### Fixed
- **BLS**: Fixed the construction of the grid of labels with multiple plots, which led to an IndexError.
- **BLS**: `gen-longini` now displays the comments below and not anymore inline so they can be properly parsed by the INI parser.
- **LIB**: Fixed conversion of "hybrid" Raman/ROA activities in `convert_y`.
- **LIB**: `Spectrum.get_xunit()` and `Spectrum.get_yunit()` could fail if no unit was defined for a given spectrum (for instance with experimental CSV files) and a conversion to LaTeX format was requested.
- **LIB**: Fixed handedness for chiral molecules in generated POV-Ray file.
- **LIB**: The `molecule` object was not displayed in the generated POV-Ray scene if no vibration was requested.
- **API**: Fixed parsing of harmonic Raman and ROA activities with multiple incident frequencies.
- **API**: Fixed parsing of energy second derivatives from Gaussian log files.
- **API**: ESTAMPES could fail to read the last diagonal element of the Cartesian force constants from a Gaussian log file if the last element was in a separate block in the output (case where 3*natoms-1 is a muliple of 5).

### Changed
- **LIB**: `spec.convert_y` now returns 3 objects, distinguishing conversions to be applied on values from the X axis and instead on values from the original X values (e.g., Raman/ROA hybrid units).
- **LIB**: `spec.broaden` now supports an additional argument, `xconv`, for conversions specific to X values, that cannot be done using the X axis.
- **LIB**: `Spectrum` now does an approximate search on the incident frequencies if the value in input is truncated to the integer part.
- **API**: Method `dfile.get_hess_data` has been completely rewritten to provide more accurate results by properly projecting out rotations and translations, and avoid code duplication.

### Removed
- **API**: Functions `get_hess_data` in parsers' modules are now deprecated and have been removed.  The `dfile.get_hess_data` method should be used instead or `tools.vib.build_vibrations` for a more general form.


## [0.5.2] - 2024-09-17

### Added
- **CRS**: Added `Corsairs` program to build Raman tensors and invariants from multiple resonance simulations.
- **BLS**: Ballast can now display automatically X and Y labels if not provided by the user.
- **MRG**: Added possibility to choose the output file to store the POVRay description content.
- **MRG**: Added option to silence help messages when building POV-Ray files.
- **LIB**: The `get_xunit` and `get_yunit` methods of `base.spectrum.Spectrum` now support additional arguments to request a LaTeX-compatible format for the unit.
- **LIB**: Added method to check if the x and y labels have been set in a `visual.plotspec.SpecLayout` instance.  To force "set" empty labels, it is necessary to use the dedicated setters.
- **LIB**: Added function `tools.char.unit_to_tex` to convert internal compact unit notation to the TeX format.  The function is quite basic and only supports properly built strings.
- **LIB**: The `write_pov` method of `visual.povrender.POVBuilder` now prints the file(s) currently written and help messages for users unfamiliar with POVRay.  The messages can be deactivated by setting the parameter `verbose` to `False`.
- **LIB**: `plot_kvec` can display two shift vectors as images with respect to the origin for comparison.  The secondary vector is printed with negative values.
- **LIB**: `orient_modes()` now support eigenvectors arrays given in 2  or 3 dimensions, in the latter case as [NVib, NAtoms, 3].
- **LIB**: The `POVBuilder` class can now support data provided directly to the constructor without having to extract them from input files.
- **API**: Added support of resonance Raman activity as `QLabel(quantity="ramact", descriptor="RR")`.
- **API**: Support of electronic transition property derivatives from the formatted checkpoint file.
- **API**: The index of the excited electronic state for vibronic calculations can be extracted from Gaussian log files.
- **API**: `DataError` class for general errors related to operations on extracted data.
- **API**: Extraction of vibrational energy levels and states from Gaussian formatted checkpoint files.

### Fixed
- **ESP**: Failure to plot J matrix on recent versions of Python.
- **BLS**: Ballast failed to parse the `Layout` block if present.
- **BLS**: Ballast failed to read all parameters in `Layout` block, which could result in the curves not being displayed.
- **BLS**: `--gen-ini` and `--gen-longini` do not print anymore an extra `\` at the start of the file.
- **MRG**: Failure to run if no vibrational mode was provided.
- **MRG**: No bond was displayed in the 3D window.
- **LIB**: Fixed test to check if an atom has been clicked in `visual.molview`.
- **LIB**: Fixed routine to compute the inertia moment.
- **LIB**: Fixed parsing of CSV files in `Spectrum` class.
- **LIB**: Fixed missing square root to the mass when display the shift vector K.
- **API**: Added support of imaginary frequencies in Gaussian log file when reading energy levels or Hessian frequencies.
- **API**: Fixed parsing of anharmonic static Raman from Gaussian log files.
- **API**: Fixed `FChkIO.write_data`.
- **API**: Fixed parsing of energy levels (`vlevel`) and transition information (`vtrans`) for resonance Raman calculations.
- **API**: Updated parsing of resonance-Raman properties (13xx) with new Gaussian output that reports Gamma together with Omega.
- **API**: Fixed support of incident frequencies for frequency-dependent properties in Gaussian output file.
- **API**: Improved support of frequency-dependent properties in Gaussian formatted checkpoint file.
- **API**: Fixed parsing of dipole and rotatory strengths from the Gaussian formatted checkpoint file.
- **API**: Electronic transition moments of electric quadrupole from Gaussian formatted checkpoint files are now correctly supported.

### Changed
- **MRG**: Options to choose the graphical output (rendering) are gather in a single argument, `--render`.  `--render=display` is equivalent to the old `-D/--display` and `--render=povray` to the old `--render`.
- **LIB**: The default coordinates systems for the derivatives is Cartesian (`X`) in `QLabel`.
- **LIB**: `visual.povrender` now build a POVRay object containing the whole description of a vibration.
- **LIB**: `visual.povrender` can display together the vibration representation and the molecule to facilitate the transformation operations and make them consistent.
- **API**: `QData.get()` does not raise an AttributeError if the the quantity is the field is not found, but simply `None`.
- **API**: It is possible to specify a default value to `QData.get()` in a way similar to the dictionary method.  The `default` keyword must be specified to set it.
- **API**: Improved error messages for missing/unsupported quantities.
- **API**: The rotation matrix for _QLabel_ `92` is now given as a list of list to respect the matrix structure.


## [0.5.1] - 2024-05-15

### Notice
The `visual` module has been heavily refactored, which may break scripts using the molecular visualization primitives.

### Added
- **MRG**: Added representation models: balls-and-sticks, sticks, spheres.
- **MRG**: Added possibility to tune the threshold factor to identify bonds.
- **MRG**: Added possibility to scale atom/bond radii directly from command-line.
- **MRG**: Added direct visualization, which becomes the default.
- **MRG**: Added possibility to visualize normal modes if available in input file.
- **MRG**: Added keyword to request the production of a POVRay input file.
- **MRG**: Added possibility to choose the material (only for rendering).
- **API**: Support of first to third derivatives of frequency-dependent properties with respect to normal coordinates from Gaussian log files.
- **API**: Fixed function `get_hess_data` function for formatted checkpoint files.
- **API**: Harmonic ROA/Raman setup data were not properly set if anharmonic data were present in Gaussian log files.
- **LIB**: Axes labels can now be changed in `plot_jmat`.
- **LIB**: Added functions in `tools.math` to compute bond lengths, angles and dihedral angles between atomic coordinates.
- **LIB**: New module `data.visual` contains the constants previously stored in `visual.molview`.
- **LIB**: New module `visual.vibview` handling the interactive representation of vibrations.
- **LIB**: New class `POVBuilder` in `visual.povrender` to facilitate the creation of POV-Ray files.
- **LIB**: New function `visual.povrender.write_pov_vib` to build the representation of vibrations.

### Fixed
- **BLS**: The program incorrectly probed for aliases to keywords while analyzing the INI file.
- **ESP**: ESParser was not properly updated to match the new behavior of `DataFile.get_data`.
- **LIB**: A typo led the Spectrum class to incorrectly check for the incident frequency in compatible spectroscopies.
- **LIB**: Superposition failed when a subset of modes was chosen, with a shift between the two structures.
- **API**: FChk failed to extract Hessian data, erroneously signaling missing data.
- **API**: GLog parser could fail to extract multiple quantities with overlapping information.
- **API**: Fixed reading of transition energies from Gaussian fchk file.
- **API**: Fixed support of electronic transition moments of properties.
- **API**: Support of energy from Gaussian fchk files.
- **API**: `parser.gaussian.fchk` did not return the right data structure.
- **API**: *QLabel* now properly supports a numerical reference state.
- **API**: Inverted logic caused the default derivative order for properties to be `None` instead of `0` in QLabel instances.
- **API**: Fixed unit definition when parsing IR intensity from Gaussian log file.
- **API**: *QLabel* now properly supports **tuples** as reference state specification.
- **API**: The default reference state for a _qlabel_ is back to 'c' to avoid issues with some routines in parsers that required a state to be set.

### Changed
- **API**: The GLog parser has been split into several sub-modules to facilitate maintenance and development.
- **LIB**: `MolWin` has been moved from `visual.molview` to `visual.molui`.
- **LIB**: All primitives to construct an input for POV-Ray have been moved from `visual.molview` to `visual.povrender`.
- **LIB**: Refactoring of the code in `visual.povrender` to make it more manual.


## [0.5.0] - 2024-02-01

### Compatibility issues
- This version introduces breaking changes at the API level!

### Added
- **API** - New *QLabel* class provided by `base.qlabel`, to handle _qlabel_ specifications.
- **API** - New *QData* class provided by `base.qdata` to structure and handle extracted data (_qdata_).
- **API** - `parser.base.build_qlabel` and `parser.base.parse_qlabel` are obsolete and have been removed.
- **API** - New submodule `parser.functions` to contain basic functions for the parser.

### Fixed
- **MRG** - The Mirage program is now properly installed as a `esmirage` executable.

### Changed
- **API** - `reshape_dblock` has been moved from `parser.base` to `parser.functions`

## [0.4.5] - 2024-01-29

### Added
- **OAR** - New script/program Soar to analyse structural difference through different metrics (for now the normal-modes overlaps) and display them graphically.
- **MRG** - New script/program Mirage will provide some facilities for molecular visualization: display, rendering...
- **LIB** - PovRay-centric routines in `visual.molview` now full supports atomic numbers as atomic labels.
- **API** - Support of vibrational reduced masses as a `qlabel` and from available parsers.

### Fixed
- **LIB** - Error when trying to truncate the broadening in `tools.spec.broaden`.
- **API** - Wrong datablock parsed for excited-state quantities read from Gaussian log file.
- **API** - Broken support of non-active modes from Gaussian anharmonic log files.
- **API** - Support of energy gradient and second derivatives in fchk file.

### Changed
- **BLS** - The filename is now displayed when a problem is found with available data.
- **LIB** - Function `data.atom.atomic_data` now respects the types of data in `atoms` and will preserve them as keys for the generated dictionary.  This simplifies working with atomic numbers or labels.
- **API** - `qlabel` "hessval" is renamed "hessdat" to underline the fact that it can be used to extract different types of data connected to the Hessian matrix.


## [0.4.4] - 2023-08-21

### Added
- **VIZ** - Added support of click events in molecule object for the 3D visualization.
- **LIB** - New function `tools.vib.build_dusch_J` to compute the Duschinsky matrix J.
- **LIB** - New function `tools.vib.build_dusch_K` to compute the shift vector K associated to the Duschinsky transformation.
- **LIB** - New function `tools.mol.inertia_mom` to compute the inertia moment from atomic masses and coordinates.
- **LIB** - Conversion factors from dipole strengths in atomic units to intensity for OPA.
- **LIB** - The `Spectrum` class now supports ECD and OPA spectra from rotatory strengths and dipole strengths, respectively.
- **LIB** - Conversion factors from rotatory strengths to intensity for ECD.
- **API** - New type for atomic masses.
- **API** - Support extraction of non-adiabatic couplings from Gaussian fchk file.
- **API** - New key in returned dictionary from `DataFile.get_data`: `qlabel`, which recalls the actual *qlabel* used.
- **API** - The `get_data` method of `DataFile` supports aliases of the type (`alias=qlabel`) to be used in the returned data dictionary.
- **API** - New `get_hess_data` method in `DataFile` to get the eigenvectors (L) to transform from mass-weighted Cartesian to normal coordinates and frequencies.  This method acts as a wrapper and tries different ways to build the data.
- **API** - Added qlabel *FCData:RedDim* and support in `parser.gaussian.glog` to extract information on normal-modes numbering in reduced-dimensionality schemes.
- **API** - Lengths and velocity gauges can be specified in the *qlabel* for dipole and rotatory strengths.
- **API** - The parser for the Gaussian log now indicates the *qlabel* for the missing quantity in the raised error.
- **API** - Added Make file to build the API documentation.

### Fixed
- **LIB** - Corrected the invariants related to the electric dipole-induced electric quadrupole in `spectro.RamanInvariants`.
- **API** - Improved parsing of version from Gaussian fchk file to support cases of files generated through utilities instead of main suite.
- **API** - Parsing of Hessian eigenvectors and eigenvalues from Gaussian fchk file.
- **API** - Reading the electronic energy from Gaussian fchk files returned an array instead of a scalar.
- **API** - `parser.gaussian.fchk` incorrectly built list of auxiliary keywords needed to process transition/excited-state properties.
- **API** - Parsing of ground-to-excited rotatory strengths in Gaussian log files.
- **API** - Fixed broken support of *AtCrd:all* in Gaussian log files.
- **API** - Fixed cases where *NAtoms* was not properly read from Gaussian log files.
- **API** - Improved internal code documentation to be compatible with Sphinx.

## [0.4.3] - 2023-02-24

### Added
- **BLS** - Added the generation of template INI files with `--gen-ini` and `--gen-longini` to facilitate the creation of input.
- **LIB** - Added possibility to overwrite original X and Y axes in class `Spectrum` through method `overwrite_axis`.
- **LIB** - New `base.spectro` module for spectroscopy-related classes and methods.  It can provide tools to build quantities or observables from more basic quantities.
- **LIB** - New class `base.spectro.RamanInvariants` to build Raman/ROA invariants from tensors.
- **LIB** - New function `base.spectro.raman_intensities` to compute the intensity based on the Raman setup.
- **API** - Added support of Raman/ROA activities given in atomic units in `tool.spec.convert_y`.
- **API** - Parsing of vibronic Resonance Raman/RROA from Gaussian log file.

### Fixed
- **ESP** - ESParser can now work without PySide installed (molecular visualization is then unavailable.)
- **ESP** - ESParser incorrectly stated that the Duschinsky matrix was identity if only the full matrix was printed and not the one after reduction of the dimension of the system (default case if the printing of J is not requested).
- **BLS** - Added correct error message if an input file for a curve is not found, instead of a useless Python traceback.
- **BLS** - Fixed legend being printed even if user explicitly deactivated it.
- **BLS** - Fixed error with `legend = auto`.
- **DEV** - Molecule-centric types are not anymore provided by `visual.molview` but are now more logically in `base.types`.
- **DEV** - Incorrect use of `typing.NoReturn` could confuse syntax analyzers (ex: Pylance).
- **API** - ROA intensity scaled by 45 to match what is done with Raman.

### Changed
- **BLS** - Removed comma as delimiter in CSV files generated by Ballast, to facilitate parsing by other programs.
- **VIZ** - Molecular visualization now requires Qt for Python 6 (PySide 6), as the support of Qt5/PySide2 is becoming scarce (no support of new hardware architectures).
- **API** - Commas get precedence over tabs as delimiters in CSV files (`parser/csv`).


## [0.4.2] - 2022-06-14

### Added
- **API** - **EXPERIMENTAL**: `parser/gaussian/glog` now supports extracting Raman and ROA activity.  Note that because of some peculiarities of Raman and how Gaussian handles it, the parser may still fail in some cases.
- **API** - `parser/gaussian/glog` now fully supports the band assignments of FC calculations.
- **API** - `parser/gaussian/glog` now supports GVPT2 variational coefficients.
- **API** - `parser/gaussian/glog` now supports the electric dipole.
- **API** - `parser/gaussian/glog` supports the qlabel `intens:IR`.
- **API** - Types for atomic labels and bonds.
- **API** - Added qlabel `Intens:` to extract intensity-related quantities for spectroscopies.
- **API** - Added support of reduced-dimensionality VPT2 in Gaussian log files.  Note that there are limitations on the use: `vtrans` and `vlevels` must be provided together for it to work.
- **VIZ** - Added class `MolWin` (`estampes.visual.molview`) to build a 3D representation of one or more molecules.
- **VIZ** - Added possibility to reverse vertically the normal mode listing for `plot_jmat`, `plot_cmat` and `plot_kvec` with `top_down`.
- **VIZ** - Added possibility to reverse the normal mode ordering for `plot_jmat`, `plot_cmat` and `plot_kvec` with `top_down`, so that normal mode 1 is the highest energy.
- **LIB** - The `superpose` function (`estampes.tools.math`) now supports a mask to apply the superposition procedure on a subset of atoms.
- **LIB** - New tool library: `mol`, for molecule-centric operations.
- **LIB** - Added support of Raman (static and dynamic) and ROA spectroscopies in class `Spectrum`.  The spectroscopy keyword `RS` now refers to the dynamic one.
- **LIB** - Added collector dictionary `params` in constructor of `Spectrum` to support spectroscopy-related keywords, like the incident frequency and spectroscopic setup  of Raman and ROA.
- **LIB** - Added some conversion factors in `tools.spec.convert_y` for Raman and ROA spectroscopies.  As there is some arbitrariness in the formulas, the conversion may not be final.
- **LIB** - Added conversion factor from integrated intensity to molar absorption coefficient for IR in `tools.spec.convert_y`

### Fixed
- **API** - Fixed parser of sub-option for `vptdat` in `parse_qlabel`.
- **API** - `parser/xyz` now properly supports the _qlabels_.
- **LIB** - The translation vector returned by `tools.math.superpose` was not properly mass-weighted.
- **VIZ** - Fixed function to generate a POV-Ray input file.  The input is now properly generated.
- **VIZ** - Alignment problem in the plotted KVec, with the bars offset by +1 compared to the mode numbering.
- **VIZ** - When moving the cursor, the normal mode index and K values did not match for the shift vector (K) in the interactive mode.
- **DOC** - Header doc and type hints of `tools.math.superpose`.

### Changed
- **LIB** - `list_bonds` is now part of a dedicated module, `tools.mol`, and not in `visual.molview` to facilitate accesses without the need for Qt libs.
- **LIB** - `convert_labsymb` now fixes the case of the atomic symbol if the conversion is from `str` to `str`.
- **LIB** - Class `Spectrum` uses qlabel `intens:IR` instead or `dipstr` to extract IR intensities, since the later is not systematically available.
- **DEV** - Improved some basic type hints, adding support for Numpy is present.

### Removed
- **API** _ `parser/xyz`: `read_data` does not support anymore the geometry index for files with multiple geometries.  This caused a non-standard behavior of the wrapper, and is partially superseded by the proper support of _qlabels_.  More parameters are accessible by using directly the internal class without passing by the `DataFile` wrapper class.


## [0.4.1] - 2021-07-14

### Added
- **API** - `parser/xyz`: Function `parse_xyz` can now ignore empty lines before a new configuration (not permitted by the parser itself, however).
- **API** - `parser/xyz`: New function `parse_xyz` to simply parse a single configuration in an opened XYZ file.

### Fixed
- **EST** - Fixed incorrect edit in the code which prevented using several features of ESTAMPES in 0.4.0 (Issue [#1](https://github.com/jbloino/estampes/issues/1))
- **API** - Fixed parsing of atomic masses in Gaussian log files.
- **API** - Incorrect parsing of a *qlabel* made of a single keyword.
- **API** - Incorrect definition of a tuple in `parse_qlabel`.
- **API** - `parser/xyz`: Fixed error if XYZ file contained empty line in the end
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
