module ChargeTransport

# for documentation
using DocStringExtensions: DocStringExtensions, SIGNATURES, TYPEDEF,
    TYPEDFIELDS, TYPEDSIGNATURES
# grid initializer
using ExtendableGrids: ExtendableGrids, BFaceNodes, BFaceRegions, CellNodes,
    CellRegions, Coordinates, ExtendableGrid,
    NumBFaceRegions, NumCellRegions, append!, dim_space,
    num_cellregions, num_nodes, subgrid
using ForwardDiff: ForwardDiff
# visualizer wrapper
using GridVisualize: GridVisualize, GridVisualizer, reveal, scalarplot!
# for interpolation of data
using Interpolations: Interpolations, Gridded, Linear
# printing
using Printf: @printf
# for interpolation of data
using Roots: Roots, find_zero
# for generating sparse arrays
using SparseArrays: SparseArrays, spzeros
# PDE solver with a FVM spatial discretization
using VoronoiFVM: VoronoiFVM, ContinuousQuantity, DiscontinuousQuantity,
    TestFunctionFactory, boundary_dirichlet!, fbernoulli_pm, physics!,
    unknown_indices, value

# for internal data handling (naming "inspired" by DrWatson)
# we do not export these functions
datadir(args...) = joinpath(pkgdir(ChargeTransport), "data", args...)
examplesdir(args...) = joinpath(pkgdir(ChargeTransport), "examples", args...)
parametersdir(args...) = joinpath(pkgdir(ChargeTransport), "parameter_files", args...)

include("ct_constants.jl")

export kB, Planck_constant, mₑ, q, ε0
##################################################################

include("ct_units.jl")

export K, J, A, V, m, s, C, kg, Hz, kHz, W, kW
export cm, mm, μm, nm, ms, μs, ns, ps, eV
export tiny_penalty_value
##################################################################

# parameter files (add new parameter sets to the list below)
for parameter_set in [
        :Params_Laser_simple,
        :Params_PSC_PCBM_MAPI_Pedot,
        :Params_PSC_TiO2_MAPI_spiro,
    ]
    include(parametersdir("$(parameter_set).jl"))
    @eval export $parameter_set
end

include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk
export FermiDiracOneHalfTeSCA, FermiDiracZero
##################################################################

include("ct_datatypes.jl")

export StandardFuncSet
export QType
export QFModelType, DiscontQF, ContQF

export OuterBoundaryModelType, OuterBoundaryModelType, InterfaceModelType
export OhmicContact, SchottkyContact, SchottkyBarrierLowering, MixedOhmicSchottkyContact
export InterfaceNone, InterfaceRecombination

export OhmicContactModelType, OhmicContactDirichlet, OhmicContactRobin

export ModelType, Transient, Stationary

export FluxApproximationType
export ScharfetterGummel, ExcessChemicalPotential, DiffusionEnhanced, DiffusionEnhancedModifiedDrift, GeneralizedSG
export ScharfetterGummelGraded, ExcessChemicalPotentialGraded
export ExcessChemicalPotentialDiffusive, ConcentrationGradient, DensityProduct

export InEquilibrium, OutOfEquilibrium

export SRHModelType, SRHWithoutTrapsType, SRHWithTrapsType
export SRHOff, SRHWithoutTrapsStationary
export SRHTrapsTransient, SRHTrapsStationary
export AuxModelSRHType, SRH2SpeciesPresentTrapDens

export GenerationModelType
export GenerationNone, GenerationBeerLambert, GenerationUniform, GenerationUserDefined
export LaserModelType, LaserModelOff, LaserModelOn
export BarrierLoweringType
export BarrierLoweringOn, BarrierLoweringOff
##################################################################

include("ct_physics.jl")

export get_BEE, get_DOS, etaFunction, get_density
export breaction!, bstorage!, reaction!, storage!, flux!
export zeroVoltage
export BeerLambert
export storage!
export reaction!, SRHRecombination!, RadiativeRecombination!, SRRecombination!, Photogeneration!
##################################################################

include("ct_system.jl")

export Params, ParamsNodal, ParamsOptical, Data, System
export BulkRecombination, set_bulk_recombination

export enable_ionic_carrier!
export enable_trap_carrier!, add_trap_density_Poisson!

export equilibrium_solve!
export enable_species!, enable_boundary_species!
export solve, solve!
export unknowns, NewtonControl, SolverControl
export TestFunctionFactory, integrate, testfunction

export gridplot

export set_contact!
export compute_open_circuit_voltage
export electroNeutralSolution
export show_params, show_paramsoptical, trap_density!
export get_current_val, charge_density

##################################################################

include("ct_plotting.jl")

export set_plotting_labels
export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann
export plot_solution, plot_IV

end # module
