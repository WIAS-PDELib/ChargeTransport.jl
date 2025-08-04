# Changelog


All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------------------------------------------------------------------------
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
  - `Data` needs a `constants` object as an argument
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
