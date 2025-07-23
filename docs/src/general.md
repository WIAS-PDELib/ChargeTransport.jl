ChargeTransport.jl -- Simulating charge transport in semiconductors
================================

[![Build status](https://github.com/WIAS-PDELib/ChargeTransport.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/WIAS-PDELib/ChargeTransport.jl/actions/workflows/ci.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ChargeTransport.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ChargeTransport.jl/dev)


`ChargeTransport.jl` simulates charge transport in semiconductors. To this end, it discretizes
the semiconductor drift-diffusion equations via the Voronoi finite volume method as implemented in [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl).

### Special features

- heterostructures
- 1D, 2D and 3D simulations
- stationary and transient simulations
- IV curves and scan protocols
- an arbitrary amount of charge carriers (electrons, heavy holes, light holes, ions, ...)
- thermodynamically consistent, physics preserving numerical methods
- different charge carrier statistics per species (Boltzmann, Blakemore, Fermi-Dirac)
- Auger, radiative, Shockley-Read-Hall recombination including transient traps
- uniform and Beer-Lambert generation

Installation and first steps
================================
The installation can easily be done via the Julia REPL with the following commands

```julia
julia> using Pkg
julia> Pkg.add("ChargeTransport")
```

We recommend have a look at the example files:

```@contents
Pages = [
    "GaAs.md",
    "PSC.md",
    ]
Depth = 2
```

You can load an example as follows

```julia
julia> include("Ex103_PSC.jl")
julia> Ex103_PSC.main()
julia> Ex103_PSC.main(plotting = true) # show plots
```
Since the examples are encapsulated into modules, you can load as many examples as you wish. If you would like to modify one of the examples, consider using [Revise.jl](https://github.com/timholy/Revise.jl) and `includet`.

Literature
===========

The simulations in the following papers are based on ChargeTransport.jl:

[1.] D. Abdel, P. Farrell and J. Fuhrmann. [Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation.](https://link.springer.com/article/10.1007/s11082-021-02803-4) Optical and Quantum Electronics **53**, 163 (2021).

[2.] D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).

[3.] D. Abdel, C. Chainais-Hillairet, P. Farrell and M. Herda. [Numerical analysis of a finite volume scheme for charge transport in perovskite solar cells.](https://doi.org/10.1093/imanum/drad034) IMA Journal of Numerical Analysis (2023).

[4.] D. Abdel, N. E. Courtier and P. Farrell. [Volume exclusion effects in perovskite charge transport modeling.](https://doi.org/10.1007/s11082-023-05125-9) Optical and Quantum Electronics **55**, 884 (2023).

[5.] B. Spetzler, D. Abdel, F. Schwierz, M. Ziegler and P. Farrell. [The Role of Vacancy Dynamics in Two-Dimensional Memristive Devices.](https://doi.org/10.1002/aelm.202300635) Advanced Electronic Materials (2023).

[6.] D. Abdel, A. Glitzky and M. Liero. [Analysis of a drift-diffusion model for perovskite solar cells.](https://doi.org/10.3934/dcdsb.2024081) Discrete and Continuous Dynamical Systems - Series B (2024).

[7.] D. Abdel, M. Herda, M. Ziegler, C. Chainais-Hillairet, B. Spetzler, P. Farrell. [Numerical analysis and simulation of lateral memristive devices: Schottky, ohmic, and multi-dimensional electrode models.](https://doi.org/10.48550/arXiv.2412.15065) submitted (2024).

[8.] B. Spetzler, E. Spetzler, S. Zamankhani, D. Abdel, P. Farrell, K.-U. Sattler, M. Ziegler. [Physics-Guided Sequence Modeling for Fast Simulation and Design Exploration of 2D Memristive Devices.](https://doi.org/10.48550/arXiv.2505.13882
) (2025).

[9.] D. Abdel, J. Relle, T. Kirchartz, P. Jaap, J. Fuhrmann, S. Burger, C. Becker, K. Jäger, P. Farrell. [Unravelling the mystery of enhanced open-circuit voltages in nanotextured perovskite solar cells.](https://doi.org/10.48550/arXiv.2506.10691) submitted (2025).
