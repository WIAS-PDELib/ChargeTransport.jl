ChargeTransport.jl -- Simulating charge transport in semiconductors
================================

[![Build status](https://github.com/PatricioFarrell/ChargeTransport.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/PatricioFarrell/ChargeTransport.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ChargeTransport.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ChargeTransport.jl/dev)


`ChargeTransport.jl` simulates charge transport in semiconductors. To this end, it discretizes
the semiconductor drift-diffusion equations via the Voronoi finite volume method as implemented in [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl).

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

We recommend to have a look at the example files:

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

[1] D. Abdel, P. Farrell and J. Fuhrmann. [Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation.](https://link.springer.com/article/10.1007/s11082-021-02803-4) Optical and Quantum Electronics 53 (163) (2021).

[2] D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).

