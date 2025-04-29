##########################################################
##########################################################
"""
$(TYPEDEF)

A struct holding all necessary information for building bulk recombination.
With help of this constructor we can read out the indices the user chooses for
electron and hole quasi Fermi potentials.

$(TYPEDFIELDS)

"""
mutable struct BulkRecombination

    """
    Index for FVM construction of electron quasi Fermi potential.
    """
    iphin::Int64

    """
    Index for FVM construction of hole quasi Fermi potential.
    """
    iphip::Int64

    """
    Boolean for present Auger recombination in bulk.
    """
    bulk_recomb_Auger::Bool

    """
    Boolean for present radiative recombination in bulk.
    """
    bulk_recomb_radiative::Bool

    """
    DataType for present SRH recombination in bulk. This needs to be a Type due to cases
    with or without mobile traps.
    """
    bulk_recomb_SRH::SRHModelType

    """
    Data type with which you can include a stationary trap density
    to the right-hand side of the Poisson equation. This stationary trap
    density corresponds to the number of unoccupied trap states.
    """
    SRH_2species_trap::AuxModelSRHType

    BulkRecombination() = new()

end


"""
$(SIGNATURES)

Corresponding constructor for the bulk recombination model.
"""
function set_bulk_recombination(;
        iphin = 1, iphip = 2,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    bulkRecombination = BulkRecombination()

    bulkRecombination.iphin = iphin
    bulkRecombination.iphip = iphip

    bulkRecombination.bulk_recomb_Auger = bulk_recomb_Auger
    bulkRecombination.bulk_recomb_radiative = bulk_recomb_radiative

    if bulk_recomb_SRH == true
        bulkRecombination.bulk_recomb_SRH = SRHStationary
        bulkRecombination.SRH_2species_trap = SRHStationary
    else
        bulkRecombination.bulk_recomb_SRH = SRHOff
        bulkRecombination.SRH_2species_trap = SRHOff
    end

    return bulkRecombination

end

###########################################################
###########################################################

"""
$(TYPEDEF)

A struct holding all information necessary for enabling traps in the SRH recombination.
With help of this constructor we can read out the index the user chooses for trap quasi Fermi
potentials and the respective regions in which they are defined.

$(TYPEDFIELDS)

"""
mutable struct TrapCarrier

    """
    Index of trap carrier user defines.
    """
    trapCarrier::Int64

    """
    Corresponding regions where trap carrier is assumed to be present.
    """
    regions::Array{Int64, 1}

    TrapCarrier() = new()

end

"""
$(SIGNATURES)
This method takes the user information concerning present trap charge carriers,
builds a struct of Type TrapCarrier and add this struct to the trapCarrierList.
"""
function enable_trap_carrier!(; data = data::Data, trapCarrier::Int64, regions::Array{Int64, 1})

    enableTraps = TrapCarrier()

    enableTraps.trapCarrier = trapCarrier
    enableTraps.regions = regions

    if data.modelType == Transient
        data.bulkRecombination.bulk_recomb_SRH = SRHTrapsTransient
    else
        data.bulkRecombination.bulk_recomb_SRH = SRHTrapsStationary
    end

    return push!(data.trapCarrierList, enableTraps)

end

"""
$(TYPEDEF)
Auxiliary struct with the charge number zt and the effective density of trap states Nt
for a stationary trap density which corresponds to the number of unoccupied states.
$(TYPEDFIELDS)
"""

mutable struct AuxiliaryStationaryTrapValues

    """
    Index of traps.
    """
    zt::Int64

    """
    Array with the corresponding effective density of trap states.
    """
    Nt::Array{Float64, 1}

    AuxiliaryStationaryTrapValues() = new()

end

"""
$(SIGNATURES)
This method includes traps for a simplified model, where the trap carriers are not
considered as additional carrier with an own continuity equation. In this case the trap
density is additionally added to the right-hand side of Poisson equation.
"""
function add_trap_density_Poisson!(; data = data::Data, zt = 1::Int64, Nt = 5.0e20 * ones(Float64, data.params.numberOfRegions)::Array{Float64, 1})

    data.bulkRecombination.SRH_2species_trap = SRH2SpeciesPresentTrapDens
    aux_trap_values = AuxiliaryStationaryTrapValues()
    aux_trap_values.zt = zt
    aux_trap_values.Nt = Nt
    return data.AuxTrapValues = aux_trap_values

end
###########################################################
###########################################################

"""
$(TYPEDEF)

A struct holding all information necessary on the ionic charge carriers which are
the index of the charge carrier and the respective region in which they are defined.
This struct along with all information necessary will be stored in an Array ionicCarrierList.
Note that it is possible to use ions as well as ion vacancies.

$(TYPEDFIELDS)

"""

mutable struct IonicCarrier

    """
    Index for data construction of ionic charge carrier
    """
    ionicCarrier::Int64

    """
    Corresponding regions where the ionic charge carrier is assumed to be present.
    """
    regions::Array{Int64, 1}

    IonicCarrier() = new()

end


"""
$(SIGNATURES)

This method takes the user information concerning present ionic charge carriers,
builds a struct of Type IonicCarrier and add this struct to the ionicCarrierList.
"""
function enable_ionic_carrier!(data; ionicCarrier::Int64, regions::Array{Int64, 1})

    enableIons = IonicCarrier()

    enableIons.ionicCarrier = ionicCarrier
    enableIons.regions = regions

    return push!(data.ionicCarrierList, enableIons)

end

###########################################################
###########################################################

###########################################################
###########################################################

"""
$(TYPEDEF)

A struct holding all information necessary for Schottky barrier lowering boundary
conditions. The implementation of this type of boundary condition needs two additional
species, see the explanation in breaction!(args ..., ::Type{SchottkyBarrierLowering}) for
further information.
$(TYPEDFIELDS)

"""

mutable struct BarrierLoweringSpecies

    """
    Datatype which gives information whether barrier lowering is turned on or off.
    """
    BarrierLoweringOn::BarrierLoweringType

    """
    Index of additional electric potential for the case with standard Schottky contacts.
    """
    ipsiStandard::QType

    """
    Additional species, where the projected gradient of the electric potential without
    Schottky barrier lowering is stored.
    """
    ipsiGrad::QType

    """
    Boundary region numbers, where Schottky barrier lowering boundary conditions are defined.
    """
    breg::Array{Int64, 1}

    """
    This quantity is needed to define the generic operator.
    """

    idx::Union{VoronoiFVM.SparseSolutionIndices, LinearIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}}

    BarrierLoweringSpecies() = new()

end


###########################################################
###########################################################
"""
$(TYPEDEF)

A struct holding the physical region dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct Params

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    """
    Number of nodes used for the disretization of the domain ``\\mathbf{\\Omega}``.
    """
    numberOfNodes::Int64

    """
    Number of subregions ``\\mathbf{\\Omega}_k`` within the domain ``\\mathbf{\\Omega}``.
    """
    numberOfRegions::Int64

    """
    Number of boundary regions ``(\\partial \\mathbf{\\Omega})_k`` such that
    `` \\partial \\mathbf{\\Omega} = \\cup_k (\\partial \\mathbf{\\Omega})_k``.
    Note that here are inner and outer boundaries calculated.
    """
    numberOfBoundaryRegions::Int64

    """
    Number of moving charge carriers.
    """
    numberOfCarriers::Int64

    """
    Parameter for the direction of illumination. If illumination is coming from the left,
    then set this value to 1. Otherwise, if the illumination comes from the right,
    set this value to -1.
    """
    invertedIllumination::Int64

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    """
    A given constant temperature.
    """
    temperature::Float64

    """
    The thermal voltage, which reads  ``U_T = k_B T / q``.
    """
    UT::Float64

    """
    The parameter of the Blakemore statistics (needed for the generalizedSG flux).
    """
    γ::Float64

    """
    Prefactor of electro-chemical reaction of internal boundary conditions.
    """
    r0::Float64

    """
    Prefactor for stationary SRH recombination.
    """
    prefactor_SRH::Float64

    """
    Parameter for the shift of generation peak of the Beer-Lambert generation profile.
    """
    generationPeak::Float64

    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################

    """
    An array for the given Schottky barriers at present Schotkky contacts.
    """
    SchottkyBarrier::Array{Float64, 1}

    """
    An array containing a constant value for the applied voltage.
    """
    contactVoltage::Array{Float64, 1}


    """
    An array containing a constant value for the electric potential
    in case of Dirichlet boundary conditions.
    """
    bψEQ::Array{Float64, 1}

    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    """
    An array with the corresponding charge numbers
    ``z_\\alpha`` for all carriers ``\\alpha``.
    """
    chargeNumbers::Array{Float64, 1}


    ###############################################################
    ####    number of boundary regions x number of carriers    ####
    ###############################################################
    """
    An array with the corresponding boundary band-edge energy values
    ``E_\\alpha`` in each region for each carrier ``\\alpha``.
    """
    bBandEdgeEnergy::Array{Float64, 2}

    """
    An array with the corresponding boundary effective density of states values
    ``N_\\alpha`` for each carrier ``\\alpha``.
    """
    bDensityOfStates::Array{Float64, 2}


    """
    A 2D array with the corresponding boundary mobility values `` \\mu_\\alpha``
    in each boundary region for each carrier ``\\alpha``.
    """
    bMobility::Array{Float64, 2}

    """
    A 2D array with the corresponding boundary doping values for each carrier ``\\alpha``.
    """
    bDoping::Array{Float64, 2}

    """
    A 2D array with the corresponding boundary velocity values for each carrier ``\\alpha``,
    when assuming Schottky contacts.
    """
    bVelocity::Array{Float64, 2}

    """
    An array to define the reaction coefficient at internal boundaries.

    """
    bReactionCoefficient::Array{Float64, 2}


    ###############################################################
    ####   number of bregions x 2 (for electrons and holes!)   ####
    ###############################################################
    """
    A 2D array with the corresponding recombination surface boundary velocity values
    for electrons and holes.
    """
    recombinationSRHvelocity::Array{Float64, 2}


    """
    A 2D array with the corresponding recombination surface boundary density values
    for electrons and holes.
    """
    bRecombinationSRHTrapDensity::Array{Float64, 2}


    """
    A 2D array with the corresponding recombination surface recombination velocities.
    """
    bRecombinationSRHLifetime::Array{Float64, 2}

    """
    A 2D array containing the equilibrium density of electric charge carriers at the boundary.
    """
    bDensityEQ::Array{Float64, 2}


    ###############################################################
    ####        number of regions x number of carriers         ####
    ###############################################################
    """
    A 2D array with the corresponding doping values for each carrier ``\\alpha`` on each region.
    """
    doping::Array{Float64, 2}

    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    densityOfStates::Array{Float64, 2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    bandEdgeEnergy::Array{Float64, 2}

    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    mobility::Array{Float64, 2}


    ###############################################################
    #### number of regions x 2 (for electrons and holes only!) ####
    ###############################################################
    """
    A 2D array with the corresponding SRH lifetimes ``\\tau_n, \\tau_p``
    for electrons and holes.
    """
    recombinationSRHLifetime::Array{Float64, 2}

    """
    A 2D array with the corresponding time-independent SRH trap densities
    ``n_{\\tau}, p_{\\tau}`` for electrons and holes.
    """
    recombinationSRHTrapDensity::Array{Float64, 2}

    """
    A 2D array with the corresponding Auger coefficients for electrons and holes.
    """
    recombinationAuger::Array{Float64, 2}


    ###############################################################
    ####                   number of regions                   ####
    ###############################################################
    """
    A region dependent dielectric constant.
    """
    dielectricConstant::Array{Float64, 1}

    """
    A region dependent image force dielectric constant.
    """
    dielectricConstantImageForce::Array{Float64, 1}

    """
    A region dependent array for the prefactor in the generation process which is the
    incident photon flux.
    """
    generationIncidentPhotonFlux::Array{Float64, 1}
    """
    A region dependent array for an uniform generation rate.
    """
    generationUniform::Array{Float64, 1}

    """
    A region dependent array for the absorption coefficient in the generation process.
    """
    generationAbsorption::Array{Float64, 1}

    """
    A region dependent array for the radiative recombination rate.
    """
    recombinationRadiative::Array{Float64, 1}

    ###############################################################
    Params() = new() # standard constructor

end


"""
$(TYPEDEF)

A struct holding the physical nodal, i.e. space-dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct ParamsNodal

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    """
    A node dependent dielectric constant.
    """
    dielectricConstant::Array{Float64, 1}

    """
    A 1D array with the corresponding doping values on each node.
    """
    doping::Array{Float64, 1}
    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha`` for each carrier
    ``\\alpha`` on each node.
    """
    mobility::Array{Float64, 2}

    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha`` for
    each carrier ``\\alpha`` on each node.
    """
    densityOfStates::Array{Float64, 2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha`` for each carrier
    ``\\alpha`` on each node.
    """
    bandEdgeEnergy::Array{Float64, 2}

    ###############################################################
    ParamsNodal() = new()

end

"""
$(TYPEDEF)

A struct holding all data information including model and numerics information,
but also all physical parameters for a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct Data{TFuncs <: Function, TVoltageFunc <: Function, TGenerationData <: Union{Array{Float64, 1}, Array{Float64, 2}, Array{Float64, 3}, Function}}

    ###############################################################
    ####                   model information                   ####
    ###############################################################
    """
    An array with the corresponding distribution function ``\\mathcal{F}_\\alpha`` for all
    carriers ``\\alpha``.
    """
    F::Array{TFuncs, 1}

    """
    An datatype containing the information, whether at least on quasi Fermi potential is
    assumed to be continuous or discontinuous.
    """
    qFModel::QFModelType

    """
    An array of DataTypes with the type of boundary model for each boundary
    (interior and exterior).
    """
    boundaryType::Array{BoundaryModelType, 1}

    """
    An array containing predefined functions for the applied bias in dependence of time
    at each outer boundary.
    """
    contactVoltageFunction::Array{TVoltageFunc, 1}

    """
    A struct containing information concerning the bulk recombination model.
    """
    bulkRecombination::BulkRecombination

    """
    A function/Array containing the user-specific photogeneration rate. 
    It can be a function which is specified in the user example 
    or an array which is read in and calculated with,
    e.g., an external software.
    """
    generationData::TGenerationData
    ###############################################################
    ####        Information on present charge carriers         ####
    ###############################################################

    """
    An array containing information on whether charge carriers are continuous or
    discontinuous. This is needed for building the AbstractQuantities which handle the
    indices of charge carriers on different regions.
    """
    isContinuous::Array{Bool, 1}

    """
    This list stores all charge carriers with the correct type needed for VoronoiFVM.
    """
    chargeCarrierList::Array{QType, 1}


    """
    This list stores all electric carrier indices, i.e. the one of electrons and holes.
    """
    electricCarrierList::Array{Int64, 1}

    """
    This list contains all defined ionic carriers as a struct of Type IonicCarrier with
    all needed information on the ionic carriers (can be either ions or ion vacancies).
    """
    ionicCarrierList::Array{IonicCarrier, 1}

    """
    This list contains all defined trap carriers for the SRH recombination
    as a struct of Type TrapCarrier with all needed information on the trap carriers.
    """
    trapCarrierList::Array{TrapCarrier, 1}

    """
    A struct which contains auxiliary trap values for the stationary setting.
    """
    AuxTrapValues::AuxiliaryStationaryTrapValues

    """
    This variable stores the index of the electric potential. Based on the user choice we have
    with this new type the opportunity to simulate discontinuous unknowns.
    """
    index_psi::QType

    """
    This is a struct containing all information necessary to simulate Schottky Barrier Lowering.
    """
    barrierLoweringInfo::BarrierLoweringSpecies

    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    """
    A DataType for the flux discretization method.
    """
    fluxApproximation::Array{FluxApproximationType, 1}

    """
    A DataType for equilibrium or out of equilibrium calculations.
    """
    calculationType::CalculationType

    """
    A DataType for transient or stationary calculations.
    """
    modelType::ModelType

    """
    A DataType for for generation model.
    """
    generationModel::GenerationModelType

    """
    An embedding parameter used to solve the nonlinear Poisson problem, where for
    λ1 = 0 the right hand-side is set to zero whereas for
    for λ1 = 1 we have a full space charge density.
    """
    λ1::Float64

    """
    An embedding parameter for the generation rate.
    """
    λ2::Float64

    """
    An embedding parameter for an electrochemical reaction.
    """
    λ3::Float64

    """
    Possibility to change the implementation of the ohmic contact boundary model
    for the electric potential (Dirichlet or Robin)
    """
    ohmicContactModel::OhmicContactModelType

    ###############################################################
    ####             Templates for DOS and BEE                 ####
    ###############################################################

    """
    Within this template, information concerning the band-edge energy
    of each carrier is stored locally which saves allocations.
    We have two of such templates due to the two point flux approximation schemes.
    """
    tempBEE1::Array{Float64, 1}

    """
    See the description of tempBEE1.
    """
    tempBEE2::Array{Float64, 1}

    """
    Within this template, information concerning the effective DOS
    of each carrier is stored locally which saves allocations.
    We have two of such templates due to the two point flux approximation schemes.
    """
    tempDOS1::Array{Float64, 1}

    """
    See the description of tempDOS2.
    """
    tempDOS2::Array{Float64, 1}

    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    """
    A struct holding all region dependent parameter information. For more information see
    struct Params.
    """
    params::Params

    """
    A struct holding all space dependent parameter information. For more information see
    struct ParamsNodal.
    """
    paramsnodal::ParamsNodal

    ###############################################################
    Data{TFuncs, TVoltageFunc, TGenerationData}() where {TFuncs, TVoltageFunc, TGenerationData} = new()

end


"""
$(TYPEDEF)

A struct holding all information necessary for a drift-diffusion type system.

$(TYPEDFIELDS)

"""
mutable struct System

    """
    A struct holding all data information, see Data
    """
    data::Data

    """
    A struct holding system information for the finite volume system.
    """
    fvmsys::VoronoiFVM.AbstractSystem

    ###############################################################
    System() = new()

end


##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Simplified constructor for Params which only takes the grid and the numberOfCarriers as argument.

"""
function Params(grid, numberOfCarriers)

    numberOfNodes = num_nodes(grid)
    numberOfRegions = grid[NumCellRegions]
    numberOfBoundaryRegions = grid[NumBFaceRegions]
    ###############################################################

    params = Params()

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    params.numberOfNodes = numberOfNodes
    params.numberOfRegions = numberOfRegions
    params.numberOfBoundaryRegions = numberOfBoundaryRegions
    params.numberOfCarriers = numberOfCarriers
    params.invertedIllumination = 1                       # we assume that light enters from the left.

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    params.temperature = 300 * K
    params.UT = (kB * 300 * K) / q # thermal voltage
    params.γ = 0.27                # parameter for Blakemore statistics
    params.r0 = 0.0                 # r0 prefactor electro-chemical reaction
    params.prefactor_SRH = 1.0
    params.generationPeak = 0.0                 # parameter which shifts Beer-Lambert generation peak

    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################
    params.SchottkyBarrier = spzeros(Float64, numberOfBoundaryRegions)
    params.contactVoltage = spzeros(Float64, numberOfBoundaryRegions)
    params.bψEQ = spzeros(Float64, numberOfBoundaryRegions)

    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    params.chargeNumbers = spzeros(Float64, numberOfCarriers)

    ###############################################################
    ####     number of carriers x number of boundary regions   ####
    ###############################################################
    params.bBandEdgeEnergy = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDensityOfStates = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bMobility = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDoping = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bVelocity = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bReactionCoefficient = 1.0e15 / s * ones(numberOfCarriers, numberOfBoundaryRegions)

    ###############################################################
    ####   2 x number of bregions (for electrons and holes!)   ####
    ###############################################################
    params.recombinationSRHvelocity = spzeros(Float64, 2, numberOfBoundaryRegions)
    params.bRecombinationSRHTrapDensity = spzeros(Float64, 2, numberOfBoundaryRegions)
    params.bRecombinationSRHLifetime = spzeros(Float64, 2, numberOfBoundaryRegions)
    params.bDensityEQ = spzeros(Float64, 2, numberOfBoundaryRegions)

    ###############################################################
    ####        number of carriers x number of regions         ####
    ###############################################################
    params.doping = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.densityOfStates = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.bandEdgeEnergy = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.mobility = spzeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    #### 2 x number of regions (for electrons and holes only!) ####
    ###############################################################
    params.recombinationSRHLifetime = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationSRHTrapDensity = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationAuger = spzeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    ####                   number of regions                   ####
    ###############################################################
    params.dielectricConstant = spzeros(Float64, numberOfRegions)
    params.dielectricConstantImageForce = spzeros(Float64, numberOfRegions)
    params.generationUniform = spzeros(Float64, numberOfRegions)
    params.generationIncidentPhotonFlux = spzeros(Float64, numberOfRegions)
    params.generationAbsorption = spzeros(Float64, numberOfRegions)
    params.recombinationRadiative = spzeros(Float64, numberOfRegions)

    ###############################################################
    return params

end


"""
$(TYPEDSIGNATURES)

Simplified constructor for ParamsNodal which only takes the grid
and the numberOfCarriers as argument.

"""
function ParamsNodal(grid, numberOfCarriers)

    numberOfNodes = length(grid[Coordinates])

    ###############################################################

    paramsnodal = ParamsNodal()

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    paramsnodal.dielectricConstant = spzeros(Float64, numberOfNodes)
    paramsnodal.doping = spzeros(Float64, numberOfNodes)

    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    paramsnodal.mobility = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.densityOfStates = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.bandEdgeEnergy = spzeros(Float64, numberOfCarriers, numberOfNodes)

    ###############################################################
    return paramsnodal

end


"""
$(TYPEDSIGNATURES)

Simplified constructor for Data which only takes the grid
and the numberOfCarriers as argument. Here, all necessary information
including the physical parameters, but also some numerical information
are located.

"""
function Data(grid, numberOfCarriers; contactVoltageFunction = [zeroVoltage for i in 1:grid[NumBFaceRegions]], generationData = [0.0], statfunctions::Type{TFuncs} = StandardFuncSet) where {TFuncs}

    numberOfBoundaryRegions = grid[NumBFaceRegions]

    ###############################################################
    # save the type of the inserted contact voltage function
    TypeVoltageFunc = Union{}

    for ii in eachindex(contactVoltageFunction)
        TypeVoltageFunc = Union{TypeVoltageFunc, typeof(contactVoltageFunction[ii])}
    end

    # save the type of generation data
    TypeGenerationData = typeof(generationData)

    # construct a data struct
    data = Data{TFuncs, TypeVoltageFunc, TypeGenerationData}()

    ###############################################################
    ####                   model information                   ####
    ###############################################################

    data.F = TFuncs[ Boltzmann for i in 1:numberOfCarriers]
    data.qFModel = ContQF
    data.boundaryType = BoundaryModelType[InterfaceNone for i in 1:numberOfBoundaryRegions]
    data.contactVoltageFunction = contactVoltageFunction
    data.generationData = generationData

    # bulkRecombination is a struct holding the input information
    data.bulkRecombination = set_bulk_recombination(
        iphin = 1, iphip = 2,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ###############################################################
    ####        Information on present charge carriers         ####
    ###############################################################
    # default values for most simple case
    data.isContinuous = Bool[true for ii in 1:numberOfCarriers]
    data.chargeCarrierList = QType[ii  for ii in 1:numberOfCarriers]
    data.electricCarrierList = Int64[ii for ii in 1:2]                       # electrons and holes
    data.ionicCarrierList = IonicCarrier[]
    data.trapCarrierList = TrapCarrier[]
    data.AuxTrapValues = AuxiliaryStationaryTrapValues()
    data.index_psi = numberOfCarriers + 1
    data.barrierLoweringInfo = BarrierLoweringSpecies()
    data.barrierLoweringInfo.BarrierLoweringOn = BarrierLoweringOff # set in general case barrier lowering off

    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    data.fluxApproximation = FluxApproximationType[ScharfetterGummel for i in 1:numberOfCarriers]
    data.calculationType = OutOfEquilibrium      # do performances InEquilibrium or OutOfEquilibrium
    data.modelType = Stationary            # indicates if we need additional time dependent part
    data.generationModel = GenerationNone        # generation model
    data.λ1 = 1.0                   # λ1: embedding parameter for NLP
    data.λ2 = 1.0                   # λ2: embedding parameter for G
    data.λ3 = 1.0                   # λ3: embedding parameter for electro chemical reaction
    data.ohmicContactModel = OhmicContactDirichlet # OhmicContactRobin also possible

    ###############################################################
    ####             Templates for DOS and BEE                 ####
    ###############################################################

    data.tempBEE1 = spzeros(Float64, numberOfCarriers)
    data.tempBEE2 = spzeros(Float64, numberOfCarriers)
    data.tempDOS1 = spzeros(Float64, numberOfCarriers)
    data.tempDOS2 = spzeros(Float64, numberOfCarriers)

    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    data.params = Params(grid, numberOfCarriers)
    data.paramsnodal = ParamsNodal(grid, numberOfCarriers)

    ###############################################################

    return data

end

###########################################################
###########################################################

"""
$(SIGNATURES)

System constructor which builds all necessary information needed based on the input parameters
with special regard to the quasi Fermi potential model. This is the main struct in which all
information on the input data, but also on the solving system, are stored.

"""
function System(grid, data; kwargs...)

    # We have currently two cases, where we use the discontinuous qF framework:
    # 1. interface charge carriers are defined
    # 2. the user chooses by themselves at least one discontinuous qF

    if all(data.isContinuous) == false
        data.qFModel = DiscontQF
    end

    # At this point, we choose a system based on usual integer indexing or quantity indexing.
    ctsys = build_system(grid, data, data.qFModel; kwargs...)

    return ctsys

end


"""
$(TYPEDSIGNATURES)

The core of the system constructor. Here, the system for continuous quasi Fermi potentials is build.

"""
function build_system(grid, data, ::Type{ContQF}; kwargs...)

    #################################################################################
    ##### Set the recombinations parameters correctly based on user information #####

    # put Auger, radiative and SRH recombination on or off (based on user information)
    if data.bulkRecombination.bulk_recomb_Auger == false
        data.params.recombinationAuger .= 0.0
    end

    if data.bulkRecombination.bulk_recomb_radiative == false
        data.params.recombinationRadiative .= 0.0
    end

    if data.bulkRecombination.bulk_recomb_SRH == SRHOff
        data.params.prefactor_SRH = 0.0
        # need to define at least one entry within each region to be non-zero. Otherwise get a NaN expression in reaction.
        for ireg in 1:grid[NumCellRegions]
            data.params.recombinationSRHTrapDensity[1, ireg] = 1.0
            data.params.recombinationSRHLifetime[1, ireg] = 1.0
        end
    end

    #################################################################################
    #####    Check, if Schottky barrier lowering conditions applicable or not   #####

    boundaryReg = Int64[]
    for ibreg in eachindex(data.boundaryType)
        if data.boundaryType[ibreg] == SchottkyBarrierLowering
            push!(boundaryReg, ibreg)
        end
    end

    if dim_space(grid) > 1 && !isempty(boundaryReg)
        error("Schottky Barrier Lowering so far only implemented in 1D.")
    elseif dim_space(grid) == 1 && length(boundaryReg) == 1
        error("Schottky Barrier Lowering only working for two contacts.")
    elseif dim_space(grid) == 1 && !isempty(boundaryReg)
        data.barrierLoweringInfo.BarrierLoweringOn = BarrierLoweringOn
    end

    #################################################################################
    #####        Set carrier lists correctly based on user information          #####
    #####    Build system for VoronoiFVM and enable carriers accordingly        #####
    ctsys = System()
    ctsys.data = data

    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOff
        physics = VoronoiFVM.Physics(
            data = data,
            flux = flux!,
            reaction = reaction!,
            storage = storage!,
            breaction = breaction!,
            bstorage = bstorage!,
            bflux = bflux!
        )
    else # in this case we add the generic operator
        physics = VoronoiFVM.Physics(
            data = data,
            flux = flux!,
            reaction = reaction!,
            storage = storage!,
            breaction = breaction!,
            bstorage = bstorage!,
            bflux = bflux!,
            generic = generic_operator!
        )
    end

    ctsys.fvmsys = VoronoiFVM.System(grid, physics; kwargs...)

    data = ctsys.fvmsys.physics.data

    ######################################
    # continuous case = integer indexing
    data.chargeCarrierList = collect(1:data.params.numberOfCarriers)
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p
    data.electricCarrierList = [iphin, iphip]
    num_species_sys = data.params.numberOfCarriers + 1
    data.index_psi = num_species_sys

    # electrons and holes
    for icc in data.electricCarrierList
        enable_species!(ctsys, icc, 1:data.params.numberOfRegions)
    end

    # if ionic carriers are present
    for icc in data.ionicCarrierList
        enable_species!(ctsys, icc.ionicCarrier, icc.regions)
    end

    # if trap carriers are present
    for icc in data.trapCarrierList
        enable_species!(ctsys, icc.trapCarrier, icc.regions)
    end

    # we need no loop for interface carriers, since in this case there are not present.

    # enable lastly the electric potential on whole domain
    enable_species!(ctsys, data.index_psi, 1:data.params.numberOfRegions)

    ######################################
    # Fill in boundary parameters. Note the convention that left boundary = 1, right boundary = 2
    # and that first region = 1, second region = 2
    bregionLeft = 1
    bregionRight = 2
    regionLeft = 1
    regionRight = data.params.numberOfRegions

    for icc in data.chargeCarrierList
        if iszero(data.paramsnodal.densityOfStates[icc, :])
            data.params.bDensityOfStates[icc, bregionLeft] = data.params.densityOfStates[icc, regionLeft]
            data.params.bDensityOfStates[icc, bregionRight] = data.params.densityOfStates[icc, regionRight]
        end

        if iszero(data.paramsnodal.bandEdgeEnergy[icc, :])
            data.params.bBandEdgeEnergy[icc, bregionLeft] = data.params.bandEdgeEnergy[icc, regionLeft]
            data.params.bBandEdgeEnergy[icc, bregionRight] = data.params.bandEdgeEnergy[icc, regionRight]
        end

        if iszero(data.paramsnodal.doping)
            data.params.bDoping[icc, bregionLeft] = data.params.doping[icc, regionLeft]
            data.params.bDoping[icc, bregionRight] = data.params.doping[icc, regionRight]
        end

    end

    ######################################
    # add here additional electric potential and boundary species in case of Schottky
    # barrier lowering conditions
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn

        data.barrierLoweringInfo.ipsiStandard = data.index_psi + 1
        data.barrierLoweringInfo.ipsiGrad = data.index_psi + 2
        data.barrierLoweringInfo.breg = boundaryReg

        enable_species!(ctsys, data.barrierLoweringInfo.ipsiStandard, 1:data.params.numberOfRegions)
        enable_boundary_species!(ctsys, data.barrierLoweringInfo.ipsiGrad, boundaryReg)

        # for detection of number of species
        VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species_sys)

        data.barrierLoweringInfo.idx = unknown_indices(unknowns(ctsys))

    end

    # for detection of number of species
    VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species_sys)

    return ctsys

end

"""
$(TYPEDSIGNATURES)

The core of the system constructor. Here, the system for discontinuous quasi Fermi potentials is build.

"""
function build_system(grid, data, ::Type{DiscontQF}; kwargs...)

    #################################################################################
    ##### Set the recombinations parameters correctly based on user information #####

    # put Auger, radiative and SRH recombination on or off (based on user information)
    if data.bulkRecombination.bulk_recomb_Auger == false
        data.params.recombinationAuger .= 0.0
    end

    if data.bulkRecombination.bulk_recomb_radiative == false
        data.params.recombinationRadiative .= 0.0
    end

    if data.bulkRecombination.bulk_recomb_SRH == SRHOff
        data.params.prefactor_SRH = 0.0
        # need to define at least one entry within each region to be non-zero. Otherwise get a NaN expression in reaction.
        for ireg in 1:grid[NumCellRegions]
            data.params.recombinationSRHTrapDensity[1, ireg] = 1.0
            data.params.recombinationSRHLifetime[1, ireg] = 1.0
        end
    end

    #################################################################################
    ##### Set carrier lists correctly based on user information #####

    fvmsys = VoronoiFVM.System(grid; kwargs...)

    #########################################
    # electrons and holes
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    data.chargeCarrierList[iphin] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = iphin)
    data.chargeCarrierList[iphip] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = iphip)
    data.electricCarrierList = [iphin, iphip]

    #########################################
    # if ionic carriers are present
    for icc in data.ionicCarrierList
        enable_species!(ctsys, icc.ionicCarrier, icc.regions)
    end

    #########################################
    # if trap carriers are present
    for icc in data.trapCarrierList
        enable_species!(ctsys, icc.trapCarrier, icc.regions)
    end

    #########################################

    data.index_psi = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions)


    #########################################
    # Fill in boundary parameters. Note the convention that left boundary = 1, right boundary = 2
    # and that first region = 1, second region = 2
    bregionLeft = 1
    bregionRight = 2
    regionLeft = 1
    regionRight = data.params.numberOfRegions
    for icc in data.chargeCarrierList
        data.params.bDensityOfStates[icc, bregionLeft] = data.params.densityOfStates[icc, regionLeft]
        data.params.bBandEdgeEnergy[icc, bregionLeft] = data.params.bandEdgeEnergy[icc, regionLeft]

        data.params.bDensityOfStates[icc, bregionRight] = data.params.densityOfStates[icc, regionRight]
        data.params.bBandEdgeEnergy[icc, bregionRight] = data.params.bandEdgeEnergy[icc, regionRight]

    end

    #########################################
    # DA: Note that Schottky barrier lowering is for the discontinuous case not implemented yet.

    #################################################################################
    #####                 Build system for VoronoiFVM                           #####

    physics = VoronoiFVM.Physics(
        data = data,
        flux = flux!,
        reaction = reaction!,
        breaction = breaction!,
        storage = storage!,
        bstorage = bstorage!,
        bflux = bflux!
    )

    # add the defined physics to system
    physics!(fvmsys, physics)

    ctsys = System()
    ctsys.fvmsys = fvmsys
    ctsys.data = data

    return ctsys

end

###########################################################
###########################################################

function show_params(ctsys::System)

    params = ctsys.data.params
    for name in fieldnames(typeof(params))[1:end]
        @printf("%30s = ", name)
        println(getfield(params, name))
    end

    return
end

function Base.show(io::IO, this::ParamsNodal)
    for name in fieldnames(typeof(this))[1:end]
        @printf("%30s = ", name)
        println(io, getfield(this, name))
    end
    return
end

###########################################################
###########################################################

"""
$(TYPEDSIGNATURES)

Master function which applies the voltage ``\\Delta u``at the
boundary ibreg for the chosen contact model.

"""

set_contact!(ctsys, ibreg, ; Δu) = __set_contact!(ctsys, ibreg, Δu, ctsys.data.boundaryType[ibreg])

# For schottky contacts
function __set_contact!(ctsys, ibreg, Δu, ::Type{SchottkyContact})

    ctsys.fvmsys.physics.data.params.contactVoltage[ibreg] = Δu
    return ctsys.data.params.contactVoltage[ibreg] = Δu

end

# For internal boundaries, do nothing
function __set_contact!(ctsys, ibreg, Δu, ::InterfaceModelType)
    return
end

# For schottky contacts with barrier lowering
function __set_contact!(ctsys, ibreg, Δu, ::Type{SchottkyBarrierLowering})

    # set Schottky barrier and applied voltage
    return ctsys.data.params.contactVoltage[ibreg] = Δu

end


function __set_contact!(ctsys, ibreg, Δu, ::Type{OhmicContact})

    ctsys.fvmsys.physics.data.params.contactVoltage[ibreg] = Δu
    return ctsys.data.params.contactVoltage[ibreg] = Δu

end


function __set_contact!(ctsys, ibreg, Δu, ::Type{MixedOhmicSchottkyContact})

    ctsys.fvmsys.physics.data.params.contactVoltage[ibreg] = Δu
    return ctsys.data.params.contactVoltage[ibreg] = Δu

end

###########################################################
###########################################################
# Wrappers for methods of VoronoiFVM

enable_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_species!(ctsys.fvmsys, ispecies, regions)
enable_boundary_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_boundary_species!(ctsys.fvmsys, ispecies, regions)

unknowns(ctsys::System) = VoronoiFVM.unknowns(ctsys.fvmsys)

solve(ctsys::System; kwargs...) = VoronoiFVM.solve(ctsys.fvmsys; kwargs...)
## DA: This one will be deleted soon (in VoronoiFVM):
solve!(solution, initialGuess, ctsys, ; control = control, tstep = tstep) = VoronoiFVM.solve!(solution, initialGuess, ctsys.fvmsys, control = control, tstep = tstep)

TestFunctionFactory(ctsys::System) = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
integrate(ctsys::System, tf, solution, inival, Δt) = VoronoiFVM.integrate(ctsys.fvmsys, tf, solution, inival, Δt)
integrate(ctsys::System, tf, solution) = VoronoiFVM.integrate(ctsys.fvmsys, tf, solution)
integrate(ctsys::System, F::Function, U; kwarges...) = VoronoiFVM.integrate(ctsys.fvmsys, F, U; kwarges...)
testfunction(factory::VoronoiFVM.TestFunctionFactory, bc0, bc1) = VoronoiFVM.testfunction(factory::VoronoiFVM.TestFunctionFactory, bc0, bc1)


# Solver Control and Newton Control are the same
function NewtonControl()

    control = VoronoiFVM.SolverControl()
    control.handle_exceptions = true # put by default handle exceptions to true

    return control
end

function SolverControl()

    control = VoronoiFVM.SolverControl()
    control.handle_exceptions = true # put by default handle exceptions to true

    return control
end
###########################################################
###########################################################

# Wrappers for GridVisualize

gridplot(grid::ExtendableGrid; Plotter, kwargs...) = GridVisualize.gridplot(grid::ExtendableGrid; Plotter, kwargs...)

###########################################################
###########################################################

"""
$(TYPEDSIGNATURES)

Functions which calculates the equilibrium solution in case of non-present fluxes and zero bias.

"""

function equilibrium_solve!(ctsys::System; control = VoronoiFVM.NewtonControl(), nonlinear_steps = 20.0, inival = nothing)

    ctsys.fvmsys.physics.data.calculationType = InEquilibrium
    grid = ctsys.fvmsys.grid

    data = ctsys.fvmsys.physics.data
    params = ctsys.fvmsys.physics.data.params
    paramsnodal = ctsys.fvmsys.physics.data.paramsnodal
    bnode = grid[BFaceNodes]
    ipsi = data.index_psi

    # We set zero voltage for each charge carrier at all outer boundaries for equilibrium calculations.
    for ibreg in grid[BFaceRegions]
        set_contact!(ctsys, ibreg, Δu = 0.0)
    end

    # initialize solution and starting vectors
    if inival === nothing
        inival = unknowns(ctsys)
        inival .= 0.0
    end

    sol = unknowns(ctsys)

    # we slightly turn a linear Poisson problem to a nonlinear one with these variables.
    I = collect(nonlinear_steps:-1:0.0)
    LAMBDA = 10 .^ (-I)
    if ctsys.fvmsys.physics.data.boundaryType[1] != SchottkyBarrierLowering
        prepend!(LAMBDA, 0.0)
    end

    for i in eachindex(LAMBDA)

        if control.verbose == "n"
            println("λ1 = $(LAMBDA[i])")
        end
        ctsys.fvmsys.physics.data.λ1 = LAMBDA[i]
        try
            sol = VoronoiFVM.solve(ctsys.fvmsys, inival = inival, control = control)
        catch
            error("try to adjust nonlinear_steps, currently set to $(nonlinear_steps) or adjust Newton control parameters.")
        end

        inival = sol

    end

    for ibreg in grid[BFaceRegions]
        # here we assume that in multidimensions, we receive a constant value of the electric potential at the boundary
        # check for applications, where this is not the case
        bψVal = view(sol[ipsi, :], subgrid(grid, [ibreg], boundary = true))[1]
        params.bψEQ[ibreg] = bψVal
    end

    # calculate equilibrium densities (especially needed for Schottky boundary conditions)
    for icc in data.electricCarrierList
        for ibreg in grid[BFaceRegions]
            Ncc = params.bDensityOfStates[icc, ibreg] + paramsnodal.densityOfStates[icc, bnode[ibreg]]
            Ecc = params.bBandEdgeEnergy[icc, ibreg] + paramsnodal.bandEdgeEnergy[icc, bnode[ibreg]]

            eta = params.chargeNumbers[icc] / params.UT * ((sol[icc, bnode[ibreg]] - sol[ipsi, bnode[ibreg]]) + Ecc / q)
            params.bDensityEQ[icc, ibreg] = Ncc * data.F[icc](eta)
        end
    end

    # set now calculationType to outOfEquilibrium for further calculations
    data.calculationType = OutOfEquilibrium

    # save changes on fvmsys of VoronoiFVM likewise in ctsys.data
    ctsys.data = ctsys.fvmsys.physics.data

    return sol

end

###########################################################
###########################################################

"""
Calculates current for time dependent problem.
"""
function get_current_val(ctsys, U, Uold, Δt) # DA: But caution, still need some small modification!

    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

    # left outer boundary = 1; right outer boundary = 2 (caution with order)
    tf = testfunction(factory, [1], [2])
    I = integrate(ctsys, tf, U, Uold, Δt)

    current = 0.0
    for ii in eachindex(I)
        current = current + I[ii]
    end

    # DA: caution I[ipsi] not completely correct. In our examples, this does not effect something,
    # but we need derivative here.
    return current
end
###########################################################
###########################################################

"""
Calculates current for stationary problem.
"""
function get_current_val(ctsys, U)

    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

    # left outer boundary = 1; right outer boundary = 2 (caution with order)
    tf = testfunction(factory, [1], [2])
    I = VoronoiFVM.integrate(ctsys.fvmsys, tf, U)

    current = 0.0
    for ii in eachindex(I)
        current = current + I[ii]
    end

    return current

end

###########################################################
###########################################################

"""

$(SIGNATURES)

For given bias vector and given IV vector this method calculates the open circuit voltage
for solar cells under illumination.
"""

function compute_open_circuit_voltage(bias::Array{Float64, 1}, IV::Array{Float64, 1})

    # http://juliamath.github.io/Interpolations.jl/latest/control/#Gridded-interpolation-1
    interpolated_IV = Interpolations.interpolate((bias,), IV, Gridded(Linear()))

    return find_zero(interpolated_IV, (bias[1], bias[end]))
end


"""

$(TYPEDSIGNATURES)

Compute the electro-neutral solution for the Boltzmann approximation.
It is obtained by setting the left-hand side in
the Poisson equation equal to zero and solving for ``\\psi``.
The charge carriers may obey different statistics functions.
Currently, this one is not well tested for the case of charge carriers beyond electrons and holes.
"""
function electroNeutralSolution(ctsys)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    params = data.params
    paramsnodal = data.paramsnodal

    if params.numberOfCarriers > 2
        error("this method is currently only working for electrons and holes")
    end

    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    psi0Vector = zeros(num_nodes(grid))
    psi0Values = zeros(num_cellregions(grid))
    cellnodes = grid[CellNodes]
    cellregions = grid[CellRegions]

    for ireg in 1:num_cellregions(grid)

        Ec = params.bandEdgeEnergy[iphin, ireg]
        Ev = params.bandEdgeEnergy[iphip, ireg]
        T = params.temperature
        Nc = params.densityOfStates[iphin, ireg]
        Nv = params.densityOfStates[iphip, ireg]
        C = params.doping[iphin, ireg] - params.doping[iphip, ireg]       # N_D - N_A
        Nintr = sqrt(Nc * Nv * exp((Ec - Ev) / (-kB * T)))

        psi0Values[ireg] = (Ec + Ev) / (2 * q) - 0.5 * (kB * T / q) * log(Nc / Nv) + (kB * T / q) * asinh(C / (2 * Nintr))

    end

    for icell in 1:size(cellnodes, 2)
        for inode in 1:size(cellnodes, 1)
            psi0Vector[cellnodes[inode, icell]] = psi0Values[cellregions[icell]]
        end
    end

    return psi0Vector

end

"""

$(TYPEDSIGNATURES)

Compute the charge density for each region separately.
"""
function charge_density(ctsys, sol)
    return VoronoiFVM.integrate(ctsys.fvmsys, reaction!, sol)[ctsys.data.index_psi, :]
end


"""

$(TYPEDSIGNATURES)

Compute the charge density, i.e. the right-hand side of Poisson's equation.

"""
function charge_density(psi0, phi, UT, EVector, chargeNumbers, dopingVector, dosVector, FVector)
    # https://stackoverflow.com/questions/45667291/how-to-apply-one-argument-to-arrayfunction-1-element-wise-smartly-in-julia
    return sum(-chargeNumbers .* dopingVector) + sum(chargeNumbers .* dosVector .* (etaFunction(psi0, phi, UT, EVector, chargeNumbers) .|> FVector))
end

"""

$(TYPEDSIGNATURES)

First try of debugger. Print the Jacobi matrix for a given node, i.e. the number of node in
the grid and not the exact coordinate. This is only done for the one dimensional case so far.
"""
function printJacobi(node, sys)
    ctdata = data(sys)
    numberOfNodes = ctdata.numberOfNodes
    return if node == 1
        println(sys.matrix[1:3, 1:9])
    elseif node == numberOfNodes
        println(sys.matrix[(3 * numberOfNodes - 2):(3 * numberOfNodes), (3 * numberOfNodes - 8):(3 * numberOfNodes)])
    else
        println(sys.matrix[(3 * node - 2):(3 * node), (3 * node - 5):(3 * node + 3)])
    end
end
