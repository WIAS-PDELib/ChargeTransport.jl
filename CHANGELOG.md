# Changelog


All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------------------------------------------------------------------------
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
