module ChargeTransport

using VoronoiFVM
using ExtendableGrids
using Printf
using DocStringExtensions
using SparseArrays
using Roots


include("ct_constants.jl")

export kB, Planck_constant, mₑ, q, ε0


include("ct_units.jl")

export K, J, A, V, m, s, C, kg
export cm, mm, μm, nm, ms, μs, ns, ps, eV
export tiny_penalty_value

##################################################################

include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk
export FermiDiracOneHalfTeSCA, FermiDiracZero

##################################################################

include("ct_datatypes.jl")

export StandardFuncSet, QType

export OuterBoundaryModel, OuterBoundaryModel, InterfaceModel
export OhmicContact, SchottkyContact
export InterfaceModelNone, InterfaceModelSurfaceReco, InterfaceModelTangentialFlux, InterfaceModelDiscontqF
export InterfaceModelIonCharge, InterfaceModelSurfaceRecoAndTangentialFlux

export ModelType, Transient, Stationary

export FluxApproximations
export ScharfetterGummel, ExcessChemicalPotential, DiffusionEnhanced, GeneralizedSG
export ScharfetterGummelGraded, ExcessChemicalPotentialGraded

export InEquilibrium, OutOfEquilibrium

export SRHWithoutTraps, SRHWithTraps
export SRHOff, SRHWithoutTrapsStationary, SRHTrapsTransient

export GenerationModel, GenerationNone, GenerationBeerLambert, GenerationUniform

export ScanProtocolType, LinearScanProtocol
##################################################################

include("ct_physics.jl")

export breaction!, bstorage!, reaction!, storage!, flux!

##################################################################

include("ct_system.jl")

export Params, ParamsNodal, Data, System
export BulkRecombination, set_bulk_recombination, IonicChargeCarriers, enable_ionic_carriers
export equilibrium_solve!
export set_contact!
export compute_densities!, compute_energies!, electroNeutralSolution!, print_jacobi
export show_params, trap_density!
export set_time_mesh, get_current_val, charge_density
export enable_traps!

##################################################################

include("ct_plotting.jl")

export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann
export plot_solution, plot_IV


end # module