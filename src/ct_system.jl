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
    else
        bulkRecombination.bulk_recomb_SRH = SRHOff
    end

    return bulkRecombination

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

    push!(data.ionicCarrierList, enableIons)

    return

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
    An array for the given Schottky barriers at present Schottky contacts.
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
$(TYPEDSIGNATURES)

Simplified constructor for Params which only takes the numberOfRegions, numberOfBoundaryRegions and numberOfCarriers as argument.

"""
function Params(numberOfRegions, numberOfBoundaryRegions, numberOfCarriers)

    @local_unitfactors K s

    params = Params()

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    params.numberOfRegions = numberOfRegions
    params.numberOfBoundaryRegions = numberOfBoundaryRegions
    params.numberOfCarriers = numberOfCarriers
    params.invertedIllumination = 1                       # we assume that light enters from the left.

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    params.temperature = 300 * K
    params.γ = 0.27                # parameter for Blakemore statistics
    params.r0 = 0.0                 # r0 prefactor electro-chemical reaction
    params.prefactor_SRH = 1.0
    params.generationPeak = 0.0                 # parameter which shifts Beer-Lambert generation peak

    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################
    params.SchottkyBarrier = zeros(Float64, numberOfBoundaryRegions)
    params.contactVoltage = zeros(Float64, numberOfBoundaryRegions)
    params.bψEQ = zeros(Float64, numberOfBoundaryRegions)

    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    params.chargeNumbers = zeros(Float64, numberOfCarriers)

    ###############################################################
    ####     number of carriers x number of boundary regions   ####
    ###############################################################
    params.bBandEdgeEnergy = zeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDensityOfStates = zeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bMobility = zeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDoping = zeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bVelocity = zeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bReactionCoefficient = 1.0e15 / s * ones(numberOfCarriers, numberOfBoundaryRegions)

    ###############################################################
    ####   2 x number of bregions (for electrons and holes!)   ####
    ###############################################################
    params.recombinationSRHvelocity = zeros(Float64, 2, numberOfBoundaryRegions)
    params.bRecombinationSRHTrapDensity = zeros(Float64, 2, numberOfBoundaryRegions)
    params.bRecombinationSRHLifetime = zeros(Float64, 2, numberOfBoundaryRegions)
    params.bDensityEQ = zeros(Float64, 2, numberOfBoundaryRegions)

    ###############################################################
    ####        number of carriers x number of regions         ####
    ###############################################################
    params.doping = zeros(Float64, numberOfCarriers, numberOfRegions)
    params.densityOfStates = zeros(Float64, numberOfCarriers, numberOfRegions)
    params.bandEdgeEnergy = zeros(Float64, numberOfCarriers, numberOfRegions)
    params.mobility = zeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    #### 2 x number of regions (for electrons and holes only!) ####
    ###############################################################
    params.recombinationSRHLifetime = zeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationSRHTrapDensity = zeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationAuger = zeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    ####                   number of regions                   ####
    ###############################################################
    params.dielectricConstant = zeros(Float64, numberOfRegions)
    params.dielectricConstantImageForce = zeros(Float64, numberOfRegions)
    params.generationUniform = zeros(Float64, numberOfRegions)
    params.generationIncidentPhotonFlux = zeros(Float64, numberOfRegions)
    params.generationAbsorption = zeros(Float64, numberOfRegions)
    params.recombinationRadiative = zeros(Float64, numberOfRegions)

    ###############################################################
    return params

end

"""
$(TYPEDSIGNATURES)

Deprecated!

Simplified constructor for Params which only takes the grid and the numberOfCarriers as argument.

"""
function Params(grid::ExtendableGrid, numberOfCarriers)
    @warn "Creating Params with a grid is deprecated and will be removed in future versions of ChargeTransport. Please call `Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)`"
    return Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)
end

###########################################################
###########################################################

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
$(TYPEDSIGNATURES)

Simplified constructor for ParamsNodal which only takes the grid
and the numberOfCarriers as argument.

"""
function ParamsNodal(grid, numberOfCarriers)

    numberOfNodes = num_nodes(grid) # = length(grid[Coordinates][1,:])

    ###############################################################

    paramsnodal = ParamsNodal()

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    paramsnodal.dielectricConstant = zeros(Float64, numberOfNodes)
    paramsnodal.doping = zeros(Float64, numberOfNodes)

    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    paramsnodal.mobility = zeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.densityOfStates = zeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.bandEdgeEnergy = zeros(Float64, numberOfCarriers, numberOfNodes)

    ###############################################################
    return paramsnodal

end

###########################################################
###########################################################


"""
$(TYPEDEF)

A struct holding the physical parameters for the Helmholtz equation simulation in a laser.
$(TYPEDFIELDS)
"""
mutable struct ParamsOptical

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    """
    The wavelength for the laser on hand.
    """
    laserWavelength::Float64

    """
    The laser power.
    """
    power::Float64

    ###############################################################
    ####                   number of regions                   ####
    ###############################################################
    """
    A region dependent array for the absorption coefficient in the
    absorption function in the medium.
    """
    absorption_0::Array{Float64, 1}

    """
    A region dependent array for the gain model coefficient.
    """
    gain_0::Array{Float64, 1}

    """
    A region dependent array for the refractive index coefficient.
    """
    refractiveIndex_0::Array{Float64, 1}

    """
    A region dependent array for the second refractive index coefficient.
    """
    refractiveIndex_d::Array{Float64, 1}

    """
    A region dependent array for the refractive index exponent.
    """
    refractiveIndex_γ::Array{Float64, 1}

    ###############################################################
    ####                 number of eigenvalues                 ####
    ###############################################################
    """
    An array of the eigenvalues.
    """
    eigenvalues::Array{Complex{Float64}, 1}

    ###############################################################
    ####        number of carriers x number of regions         ####
    ###############################################################
    """
    A 2D array with the corresponding free carrier absorption values.
    """
    absorptionFreeCarriers::Array{Float64, 2}

    ###############################################################
    ####        number of nodes x number of eigenvalues        ####
    ###############################################################
    """
    A 2D array with the corresponding eigenvector for eah eigenvalue.
    """
    eigenvectors::Array{Complex{Float64}, 2}

    ###############################################################
    ####        number of carriers + 1 x number of nodes       ####
    ###############################################################
    """
    A 2D array with the calculated solutions ``\\varphi_n``,
    ``\\varphi_p`` and``\\psi`` in all the nodes.
    """
    oldSolution::Array{Float64, 2}

    ###############################################################
    ParamsOptical() = new()

end


"""
$(TYPEDSIGNATURES)
Simplified constructor for ParamsOptical which only takes the grid,
numberOfCarriers and numberOfEigenvalues as argument.
"""
function ParamsOptical(grid, numberOfCarriers, numberOfEigenvalues)

    numberOfNodes = num_nodes(grid)
    numberOfRegions = grid[NumCellRegions]
    ###############################################################

    paramsoptical = ParamsOptical()

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    paramsoptical.laserWavelength = 0.0
    paramsoptical.power = 0.0

    ###############################################################
    ####                   number of regions                   ####
    ###############################################################
    paramsoptical.absorption_0 = zeros(Float64, numberOfRegions)
    paramsoptical.gain_0 = zeros(Float64, numberOfRegions)
    paramsoptical.refractiveIndex_0 = zeros(Float64, numberOfRegions)
    paramsoptical.refractiveIndex_d = zeros(Float64, numberOfRegions)
    paramsoptical.refractiveIndex_γ = zeros(Float64, numberOfRegions)

    ###############################################################
    ####                 number of eigenvalues                 ####
    ###############################################################
    paramsoptical.eigenvalues = zeros(Complex, numberOfEigenvalues)

    ###############################################################
    ####        number of carriers x number of regions         ####
    ###############################################################
    paramsoptical.absorptionFreeCarriers = zeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    ####        number of nodes x number of eigenvalues        ####
    ###############################################################
    paramsoptical.eigenvectors = zeros(Complex, numberOfNodes, numberOfEigenvalues)

    ###############################################################
    ####        number of carriers + 1 x number of nodes       ####
    ###############################################################
    paramsoptical.oldSolution = zeros(Float64, numberOfCarriers + 1, numberOfNodes)

    ###############################################################
    return paramsoptical
end

###########################################################
###########################################################


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
    A datatype containing the information, whether at least on quasi Fermi potential is
    assumed to be continuous or discontinuous.
    """
    qFModel::QFModelType

    """
    An array with the measure of each region of the domain.
    """
    regionVolumes::Array{Float64, 1}

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
    A function/array containing the user-specific photogeneration rate.
    It can be a function which is specified in the user example
    or an array which is read in and calculated with,
    e.g., an external software.
    """
    generationData::TGenerationData

    """
    A datatype defining whether the user wants to use the laser model or not.
    """
    laserModel::LaserModelType

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

    """
    A struct holding the physical parameters for the Helmholtz equation simulation in a laser.
    """
    paramsoptical::ParamsOptical

    """
    A struct holding the dimensionless physical constants used for the simulations.
    """
    constants::Constants

    ###############################################################
    Data{TFuncs, TVoltageFunc, TGenerationData}() where {TFuncs, TVoltageFunc, TGenerationData} = new()

end


"""
$(TYPEDSIGNATURES)

Simplified constructor for Data which only takes the grid
and the numberOfCarriers as argument. Here, all necessary information
including the physical parameters, but also some numerical information
are located.

"""
function Data(grid, numberOfCarriers; constants = ChargeTransport.constants, contactVoltageFunction = [zeroVoltage for i in 1:grid[NumBFaceRegions]], generationData = [0.0], statfunctions::Type{TFuncs} = StandardFuncSet, numberOfEigenvalues = 0) where {TFuncs}

    numberOfBoundaryRegions = grid[NumBFaceRegions]
    numberOfRegions = grid[NumCellRegions]

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

    data.regionVolumes = zeros(numberOfRegions)
    for ireg in 1:numberOfRegions
        subg = subgrid(grid, [ireg])

        # calculate measure of region
        mOmega = 0.0
        for icellVol in subg[ExtendableGrids.CellVolumes]
            mOmega = mOmega + icellVol
        end

        data.regionVolumes[ireg] = mOmega

    end

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

    if numberOfEigenvalues == 0
        data.laserModel = LaserModelOff        # by default, no laser model is used
    else
        data.laserModel = LaserModelOn         # this is needed to define the stimulated recombination in ct_physics.jl
    end

    ###############################################################
    ####        Information on present charge carriers         ####
    ###############################################################
    # default values for most simple case
    data.isContinuous = Bool[true for ii in 1:numberOfCarriers]
    data.chargeCarrierList = QType[ii  for ii in 1:numberOfCarriers]
    data.electricCarrierList = Int64[ii for ii in 1:2]                       # electrons and holes
    data.ionicCarrierList = IonicCarrier[]
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

    data.tempBEE1 = zeros(Float64, numberOfCarriers)
    data.tempBEE2 = zeros(Float64, numberOfCarriers)
    data.tempDOS1 = zeros(Float64, numberOfCarriers)
    data.tempDOS2 = zeros(Float64, numberOfCarriers)

    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    data.params = Params(grid[NumCellRegions], numberOfBoundaryRegions, numberOfCarriers)
    data.paramsnodal = ParamsNodal(grid, numberOfCarriers)
    data.paramsoptical = ParamsOptical(grid, numberOfCarriers, numberOfEigenvalues)

    ###############################################################

    data.constants = constants

    return data

end

###########################################################
###########################################################


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

    q = data.constants.q
    k_B = data.constants.k_B
    T = data.params.temperature

    # if ionic carriers are present
    for iicc in data.ionicCarrierList
        enable_species!(ctsys, iicc.ionicCarrier, iicc.regions)

        for ireg in iicc.regions

            icc = iicc.ionicCarrier # species number chosen by user

            if data.params.bandEdgeEnergy[icc, ireg] == 0.0 && ireg > 1

                ## give some initial value of ionic energy level
                En1 = data.params.bandEdgeEnergy[iphin, ireg - 1]
                Ep1 = data.params.bandEdgeEnergy[iphip, ireg - 1]
                Nn1 = data.params.densityOfStates[iphin, ireg - 1]
                Np1 = data.params.densityOfStates[iphip, ireg - 1]
                C1 = data.params.doping[iphin, ireg - 1] - data.params.doping[iphip, ireg - 1]
                Nintr1 = sqrt(Nn1 * Np1 * exp((En1 - Ep1) / (-k_B * T)))

                psi1 = (En1 + Ep1) / (2 * q) - 0.5 * (k_B * T / q) * log(Nn1 / Np1) + (k_B * T / q) * asinh(C1 / (2 * Nintr1))
                ###
                En2 = data.params.bandEdgeEnergy[iphin, ireg + 1]
                Ep2 = data.params.bandEdgeEnergy[iphip, ireg + 1]
                Nn2 = data.params.densityOfStates[iphin, ireg + 1]
                Np2 = data.params.densityOfStates[iphip, ireg + 1]
                C2 = data.params.doping[iphin, ireg + 1] - data.params.doping[iphip, ireg + 1]
                Nintr2 = sqrt(Nn2 * Np2 * exp((En2 - Ep2) / (-k_B * T)))

                psi2 = (En2 + Ep2) / (2 * q) - 0.5 * (k_B * T / q) * log(Nn2 / Np2) + (k_B * T / q) * asinh(C2 / (2 * Nintr2))

                Na = data.params.densityOfStates[icc, ireg]
                za = data.params.chargeNumbers[icc]
                Ca = data.params.doping[icc, ireg]

                Ea = trunc((k_B * T * log((Ca / Na) / (1 - Ca / Na)) + za * q * (psi1 + psi2) / 2) / q, digits = 3) * q

                data.params.bandEdgeEnergy[icc, ireg] = Ea

            end

        end

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

function show_paramsoptical(ctsys::System)              # ZA: find command to shorten very long vectors in terminal output

    paramsoptical = ctsys.data.paramsoptical
    for name in fieldnames(typeof(paramsoptical))[1:end]
        @printf("%30s = ", name)
        println(display(getfield(paramsoptical, name)))
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
    ctsys.data.params.contactVoltage[ibreg] = Δu
    return

end

# For internal boundaries, do nothing
function __set_contact!(ctsys, ibreg, Δu, ::InterfaceModelType)
    return
end

# For schottky contacts with barrier lowering
function __set_contact!(ctsys, ibreg, Δu, ::Type{SchottkyBarrierLowering})

    # set Schottky barrier and applied voltage
    ctsys.data.params.contactVoltage[ibreg] = Δu
    return

end


function __set_contact!(ctsys, ibreg, Δu, ::Type{OhmicContact})

    ctsys.fvmsys.physics.data.params.contactVoltage[ibreg] = Δu
    ctsys.data.params.contactVoltage[ibreg] = Δu
    return

end


function __set_contact!(ctsys, ibreg, Δu, ::Type{MixedOhmicSchottkyContact})

    ctsys.fvmsys.physics.data.params.contactVoltage[ibreg] = Δu
    ctsys.data.params.contactVoltage[ibreg] = Δu
    return

end

###########################################################
###########################################################
# Wrappers for methods of VoronoiFVM

enable_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_species!(ctsys.fvmsys, ispecies, regions)
enable_boundary_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_boundary_species!(ctsys.fvmsys, ispecies, regions)

unknowns(ctsys::System) = VoronoiFVM.unknowns(ctsys.fvmsys)

solve(ctsys::System; kwargs...) = VoronoiFVM.solve(ctsys.fvmsys; kwargs...)

VoronoiFVM.TestFunctionFactory(ctsys::System) = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
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

Function which calculates the equilibrium solution in case of non-present fluxes and zero bias.

"""

function equilibrium_solve!(ctsys::System; inival = VoronoiFVM.unknowns(ctsys.fvmsys, inival = 0.0), control = VoronoiFVM.NewtonControl(), nonlinear_steps = 20.0)

    ctsys.fvmsys.physics.data.calculationType = InEquilibrium
    grid = ctsys.fvmsys.grid

    data = ctsys.fvmsys.physics.data
    params = ctsys.fvmsys.physics.data.params
    paramsnodal = ctsys.fvmsys.physics.data.paramsnodal
    bnode = grid[BFaceNodes]
    ipsi = data.index_psi
    (; k_B, q) = data.constants


    # We set zero voltage for each charge carrier at all outer boundaries for equilibrium calculations.
    for ibreg in grid[BFaceRegions]
        set_contact!(ctsys, ibreg, Δu = 0.0)
    end


    sol = inival

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

            eta = params.chargeNumbers[icc] / (k_B * params.temperature / q) * ((sol[icc, bnode[ibreg]] - sol[ipsi, bnode[ibreg]]) + Ecc / q)
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

$(TYPEDSIGNATURES)

Calculates the energy value for the vacancies via the Brent method.
This method is a root-finding method, making use of the biscetion method, secant method, as well as the inverse quadratic interpolation, see
R. P. Brent. “An algorithm with guaranteed convergence for finding a zero of a function”. In: Computer J. 14. (1971), p. 422-425.

We will use this method to calculate suitable values for vacancy energy levels and internally modify the corresponding parameter.
"""
function calculate_Ea!(ctsys::System, ytol::Float64 = 1.0e-4, xtol::Float64 = 1.0e-5, maxiter::Int64 = 50; inival = VoronoiFVM.unknowns(ctsys.fvmsys, inival = 0.0), control = VoronoiFVM.NewtonControl())

    data = ctsys.fvmsys.physics.data
    params = data.params

    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p
    T = params.temperature
    k_B = data.constants.k_B
    q = data.constants.q

    x0, y0 = 0.0, 0.0
    x1, y1 = 0.0, 0.0

    for iicc in data.ionicCarrierList

        for ireg in iicc.regions

            icc = iicc.ionicCarrier # species number chosen by user

            # --- define function to be minimized ---
            mOmega = data.regionVolumes[ireg]
            Avgncc(sol) = integrated_density(ctsys, sol = sol, icc = icc, ireg = ireg) / mOmega
            Ca = params.doping[icc, ireg]

            # difference between integral and doping
            F(sol) = (Avgncc(sol) - Ca) / Ca

            # --- for initial value of Ea ---
            Ec = params.bandEdgeEnergy[iphin, ireg - 1]
            Ev = params.bandEdgeEnergy[iphip, ireg - 1]
            Nc = params.densityOfStates[iphin, ireg - 1]
            Nv = params.densityOfStates[iphip, ireg - 1]
            C = params.doping[iphin, ireg - 1] - params.doping[iphip, ireg - 1]
            Nintr = sqrt(Nc * Nv * exp((Ec - Ev) / (-k_B * T)))

            psiL = (Ec + Ev) / (2 * q) - 0.5 * (k_B * T / q) * log(Nc / Nv) + (k_B * T / q) * asinh(C / (2 * Nintr))
            ####################
            Ec = params.bandEdgeEnergy[iphin, ireg + 1]
            Ev = params.bandEdgeEnergy[iphip, ireg + 1]
            Nc = params.densityOfStates[iphin, ireg + 1]
            Nv = params.densityOfStates[iphip, ireg + 1]
            C = params.doping[iphin, ireg + 1] - params.doping[iphip, ireg + 1]
            Nintr = sqrt(Nc * Nv * exp((Ec - Ev) / (-k_B * T)))

            psiR = (Ec + Ev) / (2 * q) - 0.5 * (k_B * T / q) * log(Nc / Nv) + (k_B * T / q) * asinh(C / (2 * Nintr))

            Na = params.densityOfStates[icc, ireg]
            za = params.chargeNumbers[icc]
            x0 = trunc((k_B * T * log((Ca / Na) / (1 - Ca / Na)) + za * q * 0.475 * (psiL + psiR)) / q, digits = 4) * q
            x1 = trunc((k_B * T * log((Ca / Na) / (1 - Ca / Na)) + za * q * 0.51 * (psiL + psiR)) / q, digits = 4) * q

            params.bandEdgeEnergy[icc, ireg] = x0

            ii = 0
            Eafix = false
            while !Eafix && ii <= maxiter
                try
                    ii = ii + 1

                    sol0 = solve(ctsys, inival = inival, control = control)

                    Eafix = true
                    x0 = params.bandEdgeEnergy[icc, ireg]
                    y0 = F(sol0)
                catch
                    params.bandEdgeEnergy[icc, ireg] = trunc(1.005 * (params.bandEdgeEnergy[icc, ireg] / q), digits = 4) * q
                end
            end

            params.bandEdgeEnergy[icc, ireg] = x1
            ii = 0
            Eafix = false
            while !Eafix && ii <= maxiter
                try
                    ii = ii + 1

                    sol1 = solve(ctsys, inival = inival, control = control)

                    Eafix = true
                    x1 = params.bandEdgeEnergy[icc, ireg]
                    y1 = F(sol1)
                catch
                    params.bandEdgeEnergy[icc, ireg] = trunc(1.005 * (params.bandEdgeEnergy[icc, ireg] / q), digits = 4) * q
                end
            end

            if y0 * y1 >= 0
                error("Root not bracketed: y1 and y2 must have opposite signs.")
            end

            if abs(y0) < abs(y1)
                # swap lower and upper bounds
                x0, x1 = x1, x0
                y0, y1 = y1, y0
            end

            # used to track history of earlier approximations
            x3 = x0

            # previous approximation (initially set to x0)
            x2 = x0
            y2 = y0

            bisection = true # flag to save bisec info in last step

            for _ in 1:maxiter

                # --- Stopping criterion in x-space (bracket width) ---
                if abs(x1 - x0) / q < xtol
                    return x1
                end

                # --- Interpolation step ---
                # If three distinct y values, inverse quadratic interpolation.
                # Otherwise, secant method
                if abs(y0 - y2) > ytol && abs(y1 - y2) > ytol
                    x = x0 * y1 * y2 / ((y0 - y1) * (y0 - y2)) +
                        x1 * y0 * y2 / ((y1 - y0) * (y1 - y2)) +
                        x2 * y0 * y1 / ((y2 - y0) * (y2 - y1))
                    x = trunc(x / q, digits = 4) * q # only allow for four digits
                else
                    x = x1 - y1 * (x1 - x0) / (y1 - y0)
                    x = trunc(x / q, digits = 4) * q # only allow for four digits
                end

                # --- Acceptance check for interpolation ---
                # If the interpolation step not well-behaved, then fall back to a bisection step
                delta = ytol
                min1 = abs(x - x1)
                min2 = abs(x1 - x2)
                min3 = abs(x2 - x3)
                if (x < (3x0 + x1) / 4 && x > x1) ||  # outside bracket
                        (bisection && min1 >= min2 / 2) ||  # no improvement
                        (!bisection && min1 >= min3 / 2) || # no improvement compared to previous
                        (bisection && min2 < delta) ||    # step size too small
                        (!bisection && min3 < delta)      # step size too small

                    xnew = trunc((x0 + x1) / (2q), digits = 4) * q
                    if xnew == x0 || xnew == x1
                        xnew = (x0 + x1) / 2   # use the exact midpoint
                    end
                    x = xnew

                    bisection = true
                else
                    bisection = false
                end

                # --- Evaluate at the Ea candidate point ---
                params.bandEdgeEnergy[icc, ireg] = x

                ii = 0
                Eafix = false
                y = 0.0
                while !Eafix && ii <= maxiter
                    try
                        ii = ii + 1

                        sol = solve(ctsys, inival = inival, control = control)

                        Eafix = true
                        x = params.bandEdgeEnergy[icc, ireg]
                        y = F(sol)
                    catch
                        println("       error got catched, so Ea is updated manually")
                        newx = round(0.999 * (params.bandEdgeEnergy[icc, ireg] / q), digits = 4) * q
                        if newx == params.bandEdgeEnergy[icc, ireg]
                            newx -= 1.0e-4 * q   # or another small step
                        end
                        println("       use new energy: ", newx/q)
                        println(" ..")
                        params.bandEdgeEnergy[icc, ireg] = newx
                    end
                end

                println(" ")
                @show x / q
                @show y

                # --- Stopping criterion in y-space ---
                if abs(y) < ytol
                    return x
                end

                # --- Update history ---
                x3 = x2
                x2 = x1

                # --- Update the bracketing interval ---
                if sign(y0) != sign(y)
                    # Root lies between x0 and x
                    x1 = x
                    y1 = y
                else
                    # Root lies between x and x1
                    x0 = x
                    y0 = y
                end

                # Ensure |f(x0)| ≥ |f(x1)| again (x1 is best approximation)
                if abs(y0) < abs(y1)
                    x0, x1 = x1, x0
                    y0, y1 = y1, y0
                end

            end # Brent method

            error("Max iteration exceeded")

        end # each present region for carrier

    end # ionic carrier list

    return
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

function integrated_density(ctsys; sol, icc, ireg)

    saveType = deepcopy(ctsys.data.modelType)
    ctsys.data.modelType = Transient

    integral = ctsys.data.params.chargeNumbers[icc] * ChargeTransport.integrate(ctsys, storage!, sol)[icc, ireg] / ctsys.data.constants.q

    ctsys.data.modelType = saveType

    return integral
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
    (; k_B, q) = data.constants


    params = data.params

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
        Nintr = sqrt(Nc * Nv * exp((Ec - Ev) / (-k_B * T)))

        psi0Values[ireg] = (Ec + Ev) / (2 * q) - 0.5 * (k_B * T / q) * log(Nc / Nv) + (k_B * T / q) * asinh(C / (2 * Nintr))

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
function charge_density(psi0, phi, temperature, EVector, chargeNumbers, dopingVector, dosVector, FVector)
    # https://stackoverflow.com/questions/45667291/how-to-apply-one-argument-to-arrayfunction-1-element-wise-smartly-in-julia
    return sum(-chargeNumbers .* dopingVector) + sum(chargeNumbers .* dosVector .* (etaFunction(psi0, phi, temperature, EVector, chargeNumbers) .|> FVector))
end
