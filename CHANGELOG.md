# Changelog


All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------------------------------------------------------------------------


## v1.2.4

### Changed
  - Cleanup of return values of plotting functions, all return `nothing` now.
  - no implicit `show` call in potting functions; do this explicitly if needed

## v1.2.3

### Changed
  - Density of states and mobilities are by default now set one, while band edge energies are set to zero.
  - If the user does not define `data.bulkRecombination`, by default recombination is set off. Caution: For semiconductor applications with electrons and holes, this method does method may still need to be initialized, e.g., when working with the Schottky barrier lowering boundary model.


## v1.2.2

### Added
  - In `equilibrium_solve!()`, we have now two additional inputs: `verbose` and `yabstol` to control the secant method for finding vacancy energy levels

## v1.2.1

#### Fixed
  - Fixed broken math mode in docs
  - fixed some plotting, by properly including parameters
  - fixed unnecessary photogeneration loop in notebook

## v1.2.0

### Added
  - Modified `equilibrium_solve!()`. Directly computes now the correct energy level such that the average vacancy density matches the user-defined target, if the argument is `vacancyEnergyCalculation = true`. Check the [package documentation](https://wias-pdelib.github.io/ChargeTransport.jl/stable/PSC/) for more information
  - Added method integrated_density, which computes the integrated carrier density for a given species `icc` and region `ireg`
  - Extended Data structure with `data.regionVolumes`, which stores the volume (measure) of each subregion

### Changed
  - loop for increasing photogeneration rate is now also included internally, see Ex103. There is now no need to do this by the user.

## v1.1.0

### Fixed
  - fixed broken v1.0.0 release

## v1.0.0

### Added
  - `Params` can be constructed with `numberOfRegions`, `numberOfBoundaryRegions` and `numberOfCarriers`
  - `Params` can be constructed directly from problem specific parameter structs

### Changed
  - parameter files are replaced by parameter structs with explicit parameter access: `p = parameter_set(); p.foo` to access parameter `foo`;
    you can specify in the examples which `parameter_set` is used.
  - all examples scripts are overhauled with the new parameter set usage
  - global unit factors are removed: we rely on local unit factors from `LessUnitFul.jl`, provided by `@local_unitfactors` and `ufac""`
  - new globally available dimensionless `constants` object, containing the default physical constants
  - new application specific `teSCA_constants`, `pdelib_constants`, `unit_constants` are also available
  - `Data` needs a `constants` object as a key word argument, defaults to standard constants
  - notebook folder name from `pluto-examples` to `notebooks`

### Removed
  - Thermal voltage `UT` is no longer part of the `Params`, since this value may depend on different definitions of the elementary charge `q`
  - methods taking `UT` as an argument take the temperature now instead
  - exported global physical constants
  - `enable_trap_carrier!()` method as the underlying model and discretization were not correctly set up

## v0.6.0 July 23, 2025

#### Fixed
  - documentation via Documenter.jl

## v0.5.0 July 17, 2025

#### Added
  - functions for internal data handling: datadir, examplesdir, parametersdir
  - possibility to add a stimulated recombination term for laser applications
  - struct ParamsOptical to hold fields for laser applications
#### Fixed
  - correct inclusion of parameter files independent of folder, from which they are started
  - corrected definition of numberOfNodes in function ParamsNodal
#### Changed
  - adjusted argument inival for function equilibrium_solve! to be inserted if desired

## v0.4.0 May 26, 2025

#### Changed
  - adjusted global (constants and units) are mutable globals now for compatibility with julia 1.2 (will change with a proper export in upcoming 1.0 release)


## v0.3.0 April 29, 2025

#### Added
  - code quality checks
  - Runic code formatting
  - post-process methods
