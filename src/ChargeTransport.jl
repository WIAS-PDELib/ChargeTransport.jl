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
    NewtonControl, SolverControl, TestFunctionFactory,
    boundary_dirichlet!, fbernoulli_pm, physics!,
    unknown_indices, value


include("ct_constants.jl")

export kB, Planck_constant, mₑ, q, ε0
##################################################################

include("ct_units.jl")

export K, J, A, V, m, s, C, kg, Hz, kHz, W, kW
export cm, mm, μm, nm, ms, μs, ns, ps, eV
export tiny_penalty_value
##################################################################

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

export Params, ParamsNodal, Data, System
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
export electroNeutralSolution, print_jacobi
export show_params, trap_density!
export get_current_val, charge_density

##################################################################

include("ct_plotting.jl")

export set_plotting_labels
export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann
export plot_solution, plot_IV

end # module
