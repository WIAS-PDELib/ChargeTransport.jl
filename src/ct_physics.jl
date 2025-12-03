##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the effective DOS for interior nodes (analogously for boundary nodes and edges).
"""
function get_DOS!(icc::QType, node::VoronoiFVM.Node, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, node.region] + data.paramsnodal.densityOfStates[icc, node.index]

    return

end

# Defining locally the effective DOS for boundary nodes.
function get_DOS!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempDOS1[icc] = data.params.bDensityOfStates[icc, bnode.region] + data.paramsnodal.densityOfStates[icc, bnode.index]

    return

end

# Defining locally the effective DOS for edges.
function get_DOS!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[1]]
    data.tempDOS2[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[2]]

    return
end

# Calculate the DOS on a given interior region.
function get_DOS(icc::QType, ireg::Int, ctsys)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    return data.params.densityOfStates[icc, ireg] .+ view(data.paramsnodal.densityOfStates[icc, :], subgrid(grid, [ireg])) # view nodal dependent DOS on respective grid
end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the band-edge energy for interior nodes (analogously for boundary nodes and edges).
"""
function get_BEE!(icc::QType, node::VoronoiFVM.Node, data)

    return data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, node.region] + data.paramsnodal.bandEdgeEnergy[icc, node.index]

end

# Defining locally the band-edge energy for boundary nodes.
function get_BEE!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempBEE1[icc] = data.params.bBandEdgeEnergy[icc, bnode.region] + data.paramsnodal.bandEdgeEnergy[icc, bnode.index]

    return

end

# Defining locally the band-edge energy for edges.
function get_BEE!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[1]]
    data.tempBEE2[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[2]]

    return

end

# Calculate the band-edge energy on a given interior region.
function get_BEE(icc::QType, ireg::Int, ctsys)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    return data.params.bandEdgeEnergy[icc, ireg] .+ view(data.paramsnodal.bandEdgeEnergy[icc, :], subgrid(grid, [ireg]))

end
##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the statistics function for interior nodes.
"""
function etaFunction!(u, node::VoronoiFVM.Node, data, icc)

    get_BEE!(icc, node::VoronoiFVM.Node, data)

    E = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / (data.constants.k_B * data.params.temperature) * ((u[icc] - u[data.index_psi]) * data.constants.q + E)

end

"""
$(TYPEDSIGNATURES)

The argument of the statistics function for boundary nodes.
"""
function etaFunction!(u, bnode::VoronoiFVM.BNode, data, icc) # bnode.index refers to index in overall mesh


    get_BEE!(icc, bnode::VoronoiFVM.BNode, data)
    E = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / (data.constants.k_B * data.params.temperature) * ((u[icc] - u[data.index_psi]) * data.constants.q + E)
end


"""
$(TYPEDSIGNATURES)

The argument of the statistics function for edges.
"""

function etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    get_BEE!(icc, edge::VoronoiFVM.Edge, data)

    E1 = data.tempBEE1[icc];  E2 = data.tempBEE2[icc]

    return etaFunction(u[data.index_psi, 1], u[icc, 1], data.params.temperature, E1, data.params.chargeNumbers[icc], data.constants),
        etaFunction(u[data.index_psi, 2], u[icc, 2], data.params.temperature, E2, data.params.chargeNumbers[icc], data.constants)
end

"""
$(TYPEDSIGNATURES)

The argument of the statistics function for a given solution on a given interior region.
"""
function etaFunction(sol, ireg::Int, ctsys, icc::QType)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    Ecc = get_BEE(icc, ireg, ctsys)
    # view solution on respective grid
    solcc = view(sol[icc, :], subgrid(grid, [ireg]))
    solpsi = view(sol[data.index_psi, :], subgrid(grid, [ireg]))

    return @. data.params.chargeNumbers[icc] / (data.constants.k_B * data.params.temperature) * ((solcc - solpsi) * data.constants.q + Ecc)
end


"""
$(TYPEDSIGNATURES)

The argument of the statistics function for given ``\\varphi_\\alpha``
and ``\\psi``

``z_\\alpha / (k_B  T)   ( (\\varphi_\\alpha - \\psi) * q + E_\\alpha ).``

The parameters ``E_\\alpha`` and ``z_\\alpha`` are given as vectors.
This function may be used to compute the charge density, i.e. the
right-hand side of the Poisson equation.
"""
function etaFunction(psi, phi, temperature, E, z, constants)
    return @. z / (constants.k_B * temperature) * ((phi - psi) * constants.q + E)
end

##########################################################
##########################################################

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for interior nodes.

"""
function get_density!(u, node::VoronoiFVM.Node, data, icc)

    get_DOS!(icc, node, data)

    Ncc = data.tempDOS1[icc]
    eta = etaFunction!(u, node, data, icc) # calls etaFunction!(u,node::VoronoiFVM.Node,data,icc)

    return Ncc * data.F[icc](eta)

end

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for interior nodes.

"""
function get_density!(u, bnode::VoronoiFVM.BNode, data, icc)

    get_DOS!(icc, bnode, data)

    Ncc = data.tempDOS1[icc]
    eta = etaFunction!(u, bnode, data, icc) # calls etaFunction!(u,node::VoronoiFVM.BNode,data,icc)

    return Ncc * data.F[icc](eta)

end


"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for edges.

"""
function get_density!(u, edge::VoronoiFVM.Edge, data, icc)

    get_DOS!(icc, edge, data)

    Ncc1 = data.tempDOS1[icc]
    Ncc2 = data.tempDOS2[icc]

    eta1, eta2 = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    return Ncc1 * data.F[icc](eta1), Ncc2 * data.F[icc](eta2)

end


"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for given interior region corresponding
to a homogeneous set of parameters.

"""
function get_density(sol, ireg::Int, ctsys, icc::QType)

    data = ctsys.fvmsys.physics.data

    Ncc = get_DOS(icc, ireg, ctsys)
    eta = etaFunction(sol, ireg, ctsys, icc)

    return Ncc .* data.F[icc].(eta)

end


"""
$(TYPEDSIGNATURES)

The densities for given potentials  ``\\varphi_\\alpha``
and ``\\psi``

"""
function get_density(sol, data, icc, ireg, ; inode)

    N = data.params.densityOfStates[icc, ireg]
    E = data.params.bandEdgeEnergy[icc, ireg]
    z = data.params.chargeNumbers[icc]

    eta = etaFunction(sol[data.index_psi, inode], sol[icc, inode], data.params.temperature, E, z, data.constants)

    return N .* data.F[icc].(eta)
end
###########################################################
###########################################################

function emptyFunction()
end

"""
Function in case of an applied voltage equal to zero at one boundary.
"""
zeroVoltage(t) = 0.0

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master breaction! function. This is the function which enters VoronoiFVM and hands over
for each boundary the chosen boundary model.

"""
breaction!(f, u, bnode, data) = breaction!(f, u, bnode, data, data.boundaryType[bnode.region])

#################################################################################################

breaction!(f, u, bnode, data, ::Type{OhmicContact}) = breaction!(f, u, bnode, data, data.calculationType)

# in case of equilibrium conditions, we choose the initial ohmic contact boundary model
breaction!(f, u, bnode, data, ::Type{InEquilibrium}) = breaction!(f, u, bnode, data, OhmicContactRobin)

breaction!(f, u, bnode, data, ::Type{OutOfEquilibrium}) = breaction!(f, u, bnode, data, data.ohmicContactModel)

"""
$(TYPEDSIGNATURES)

Creates ohmic boundary conditions via a penalty approach with penalty parameter ``\\delta``.
For example, the right-hand side for the electrostatic potential ``\\psi`` is implemented as

``f[\\psi]  = -q/\\delta   ( (p - N_a) - (n - N_d) )``,

assuming a bipolar semiconductor. In general, we have for some given charge number ``z_\\alpha``

``f[\\psi] =  -q/\\delta  \\sum_\\alpha{ z_\\alpha  (n_\\alpha - C_\\alpha) },``

where ``C_\\alpha`` corresponds to some doping w.r.t. the species ``\\alpha``.

The boundary conditions for electrons and holes are dirichlet conditions, where

`` \\varphi_{\\alpha} = U```

with ``U`` as an applied voltage.
"""
function breaction!(f, u, bnode, data, ::Type{OhmicContactRobin})

    params = data.params
    paramsnodal = data.paramsnodal

    ipsi = data.index_psi

    # electrons and holes entering right hand-side for BC of ipsi
    for icc in data.electricCarrierList         # Array{Int64, 1}

        icc = data.chargeCarrierList[icc]  # Array{QType, 1}
        ncc = get_density!(u, bnode, data, icc)

        # subtract doping
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * (params.bDoping[icc, bnode.region])
        # add charge carrier
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

    end

    # if ionic carriers are present
    for iicc in data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if bnode.cellregions[1] ∈ iicc.regions    # bnode.cellregions = [bnode.region, 0] for outer boundary.
            icc = iicc.ionicCarrier           # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            ncc = get_density!(u, bnode, data, icc)

            # subtract doping
            f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * (params.bDoping[icc, bnode.region])
            # add charge carrier
            f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

        end

    end
    # if trap carriers are present
    for iicc in data.trapCarrierList
        # add trap carriers only in defined regions (otherwise get NaN error)
        if bnode.cellregions[1] ∈ iicc.regions    # bnode.cellregions = [bnode.region, 0] for outer boundary.
            icc = iicc.trapCarrier           # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            ncc = get_density!(u, bnode, data, icc)

            # subtract doping
            f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * (params.bDoping[icc, bnode.region])
            # add charge carrier
            f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

        end
    end

    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]
    f[ipsi] = - data.λ1 * 1 / tiny_penalty_value * data.constants.q * f[ipsi]

    # electrons and holes boundary condition
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    Δu = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    boundary_dirichlet!(f, u, bnode, species = iphin, region = bnode.region, value = Δu)
    boundary_dirichlet!(f, u, bnode, species = iphip, region = bnode.region, value = Δu)
    return

end


"""
$(TYPEDSIGNATURES)

Creates ohmic boundary conditions via Dirichlet BC for the electrostatic potential ``\\psi``

``\\psi  = \\psi_0 + U``,

where ``\\psi_0`` contains some given value and ``U`` is an applied voltage.

``f[\\psi] =  -q/\\delta  \\sum_\\alpha{ z_\\alpha  (n_\\alpha - C_\\alpha) },``

where ``C_\\alpha`` corresponds to some doping w.r.t. the species ``\\alpha``.

The boundary conditions for electrons and holes are dirichlet conditions, where

`` \\varphi_{\\alpha} = U.```

"""
function breaction!(f, u, bnode, data, ::Type{OhmicContactDirichlet})

    params = data.params

    # DA: we get here an issue with the allocation, if we pass into boundary_dirichlet! something which is not of type Int,
    # as e.g. an AbstractQuantity
    ipsi = params.numberOfCarriers + 1 # data.index_psi
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    Δu = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)
    ψ0 = params.bψEQ[bnode.region]

    boundary_dirichlet!(f, u, bnode, species = iphin, region = bnode.region, value = Δu)
    boundary_dirichlet!(f, u, bnode, species = iphip, region = bnode.region, value = Δu)
    boundary_dirichlet!(f, u, bnode, species = ipsi, region = bnode.region, value = ψ0 + Δu)
    return

end

"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions. For the electrostatic potential we assume

``\\psi = - \\phi_S/q + U, ``

where  ``\\phi_S`` corresponds to a given value (non-negative Schottky barrier) and ``U`` to the applied voltage.
The quantity ``\\phi_S`` needs to be specified in the main file.
For electrons and holes we assume the following

``f[n_\\alpha]  =  z_\\alpha q v_\\alpha (n_\\alpha - n_{\\alpha, 0})``,

where ``v_{\\alpha}`` can be treated as a surface recombination mechanism and is given. The parameter
``n_{\\alpha, 0}`` is the equilibrium density of the charge carrier ``\\alpha`` and can be
calculated via

``n_{\\alpha, 0}= N_\\alpha \\mathcal{F}_\\alpha \\Bigl( - z_\\alpha/ U_T (E_c - E_\\alpha) - \\phi_S) / q \\Bigr). ``

"""

function breaction!(f, u, bnode, data, ::Type{SchottkyContact})

    params = data.params
    paramsnodal = data.paramsnodal
    ipsi = data.index_psi
    iphin = data.bulkRecombination.iphin
    Ec = params.bBandEdgeEnergy[iphin, bnode.region]

    for icc in data.electricCarrierList       # Array{Int64, 1}

        icc = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used

        get_DOS!(icc, bnode, data);  get_BEE!(icc, bnode, data)
        Ni = data.tempDOS1[icc]
        Ei = data.tempBEE1[icc]
        etaFix = - params.chargeNumbers[icc] / (data.constants.k_B * params.temperature) * (((Ec - Ei) - params.SchottkyBarrier[bnode.region]))

        ncc = get_density!(u, bnode, data, icc)

        f[icc] = params.chargeNumbers[icc] * data.constants.q * params.bVelocity[icc, bnode.region] * (ncc - Ni * data.F[icc](etaFix))

    end

    # function evaluation causes allocation!!!
    Δu = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    ipsiIndex = length(data.chargeCarrierList) + 1 # This is necessary, since passing something other than an Integer in boundary_dirichlet!() causes allocations
    boundary_dirichlet!(f, u, bnode, species = ipsiIndex, region = bnode.region, value = (- (params.SchottkyBarrier[bnode.region] - Ec) / data.constants.q) + Δu)
    return

end

###########################################################################
###########################################################################
"""
$(TYPEDSIGNATURES)
A mixed Schottky-Ohmic boundary type condition, where we impose on the electric potential (Schottky)

``\\psi = - \\phi_S/q + U, ``

with  ``\\phi_S`` as given value (non-negative Schottky barrier) and ``U`` to the applied voltage.
The quantity ``\\phi_S`` needs to be specified in the main file.
For electrons and holes we assume the following (Ohmic)

`` \\varphi_{\\alpha} = U``.
"""

function breaction!(f, u, bnode, data, ::Type{MixedOhmicSchottkyContact})

    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p
    ipsiIndex = length(data.chargeCarrierList) + 1 # This is necessary, since passing something other than an Integer in boundary_dirichlet!() causes allocations

    params = data.params
    Ec = params.bBandEdgeEnergy[iphin, bnode.region]
    Δu = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)


    # electric potential BC
    boundary_dirichlet!(f, u, bnode, species = ipsiIndex, region = bnode.region, value = (- (params.SchottkyBarrier[bnode.region] - Ec) / data.constants.q) + Δu)

    # electrons and holes boundary condition
    boundary_dirichlet!(f, u, bnode, species = iphin, region = bnode.region, value = Δu)
    boundary_dirichlet!(f, u, bnode, species = iphip, region = bnode.region, value = Δu)
    return

end

###########################################################################
###########################################################################

"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions with additional lowering which are modelled as

`` \\psi = - \\phi_S/q  + \\sqrt{ -\\frac{ q  \\nabla_{\\boldsymbol{\\nu}} \\psi_\\mathrm{R}}{4\\pi \\varepsilon_\\mathrm{i}}} + U``,

where `` \\psi_\\mathrm{R}`` denotes the electric potential with standard Schottky contacts and the same space charge density as `` \\psi`` and where ``\\varepsilon_\\mathrm{i}}}`` corresponds to the image force dielectric constant.

To solve for this additional boundary conditions the projected gradient ``\\nabla_{\\boldsymbol{\\nu}} \\psi_\\mathrm{R} `` is stored within a boundary species and calculated in the method generic_operator!().
"""
function breaction!(f, u, bnode, data, ::Type{SchottkyBarrierLowering})

    params = data.params
    ipsi = data.index_psi
    iphin = data.bulkRecombination.iphin
    Ec = params.bBandEdgeEnergy[iphin, bnode.region]
    ipsiStandard = data.barrierLoweringInfo.ipsiStandard
    ipsiGrad = data.barrierLoweringInfo.ipsiGrad

    q = data.constants.q

    if data.calculationType == OutOfEquilibrium
        for icc in data.electricCarrierList       # Array{Int64, 1}

            icc = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used
            ncc = get_density!(u, bnode, data, icc)

            f[icc] = params.chargeNumbers[icc] * q * params.bVelocity[icc, bnode.region] * (ncc - params.bDensityEQ[icc, bnode.region])

        end
    end

    # function evaluation causes allocation!!!
    Δu = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    if u[ipsiGrad] < 0
        PsiS = sqrt(- q / (4 * pi * params.dielectricConstantImageForce[bnode.cellregions[1]]) * u[ipsiGrad]) #bnode.cellregions[1] ∈ iicc.regions    # bnode.cellregions = [bnode.region, 0] for outer boundary.
    else
        PsiS = 0.0
    end

    f[ipsiStandard] = 1 / tiny_penalty_value * (u[ipsi] + (params.SchottkyBarrier[bnode.region] - Ec) / q - PsiS - Δu)
    f[ipsi] = 1 / tiny_penalty_value * (u[ipsiStandard] + (params.SchottkyBarrier[bnode.region] - Ec) / q - Δu)
    return

end


# This breaction! function is chosen when no interface model is chosen.
breaction!(f, u, bnode, data, ::Type{InterfaceNone}) = emptyFunction()


function breaction!(f, u, bnode, data, ::Type{InterfaceRecombination})

    params = data.params
    (; q, k_B) = data.constants

    if data.calculationType == InEquilibrium
        return
    end

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    n = get_density!(u, bnode, data, iphin)
    p = get_density!(u, bnode, data, iphip)

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (k_B * params.temperature))
    excessDensTerm = n * p * (1.0 - exponentialTerm)

    if params.recombinationSRHvelocity[iphip, bnode.region] ≈ 0.0
        vp = 1.0e30
    else
        vp = 1.0 / params.recombinationSRHvelocity[iphip, bnode.region]
    end

    if params.recombinationSRHvelocity[iphin, bnode.region] ≈ 0.0
        vn = 1.0e30
    else
        vn = 1.0 / params.recombinationSRHvelocity[iphin, bnode.region]
    end

    kernelSRH = 1.0 / (vp * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + vn * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region]))

    for icc in data.electricCarrierList
        icc = data.chargeCarrierList[icc]
        f[icc] = q * params.chargeNumbers[icc] * kernelSRH * excessDensTerm
    end

    return
end

###########################################################################
###########################################################################

# Robin boundary conditions for Gate contact
function breaction!(f, u, bnode, data, ::Type{GateContact})

    params = data.params
    ipsi = data.index_psi

    # Robin Condition psi
    f[ipsi] = (1 / params.dielectricConstant[bnode.region]) * ((params.oxidePermittivity / params.oxideThickness) * (u[ipsi] - params.contactVoltage[bnode.region] - params.additionalVoltage) - params.surfacechargeDensity)

    # Dirichlet Conditions phin and phip (electrons and holes)
    # Homogeneous Neumann boundary conditions by default
    return
end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Generic operator to save the projected gradient of electric potential
(for system with standard Schottky contacts). Note that this currently
only working in one dimension!

"""
function generic_operator!(f, u, fvmsys, data)

    f .= 0

    coord = fvmsys.grid[Coordinates]
    n = length(coord)
    barrierLoweringInfo = fvmsys.physics.data.barrierLoweringInfo
    ipsiStandard = barrierLoweringInfo.ipsiStandard
    ipsiGrad = barrierLoweringInfo.ipsiGrad

    idx = barrierLoweringInfo.idx

    f[idx[ipsiGrad, 1]] = u[idx[ipsiGrad, 1]] - (-1) * (u[idx[ipsiStandard, 2]] - u[idx[ipsiStandard, 1]]) / (coord[2] - coord[1])
    f[idx[ipsiGrad, n]] = u[idx[ipsiGrad, n]] - (u[idx[ipsiStandard, n]] - u[idx[ipsiStandard, n - 1]]) / (coord[n] - coord[n - 1])
    return


end


##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bstorage! function. This is the function which enters VoronoiFVM and hands over
for each boundary the time-dependent part of the chosen boundary model.

"""
bstorage!(f, u, bnode, data) = bstorage!(f, u, bnode, data, data.modelType)

bstorage!(f, u, bnode, data, ::Type{Stationary}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{Transient}) = bstorage!(f, u, bnode, data, data.boundaryType[bnode.region])


bstorage!(f, u, bnode, data, ::Type{InterfaceNone}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{InterfaceRecombination}) = emptyFunction()

# No bstorage! is used, if an ohmic and schottky contact model is chosen.
bstorage!(f, u, bnode, data, ::OuterBoundaryModelType) = emptyFunction()

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bflux! function. This is the function which enters VoronoiFVM and hands over
for each boundary the flux within the boundary.

"""
bflux!(f, u, bedge, data) = emptyFunction()


##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master reaction! function. This is the function which enters VoronoiFVM and hands over
reaction terms for concrete calculation type and bulk recombination model.

"""
reaction!(f, u, node, data) = reaction!(f, u, node, data, data.calculationType)

"""
$(TYPEDSIGNATURES)
Reaction in case of equilibrium, i.e. no generation and recombination is considered.
"""
function reaction!(f, u, node, data, ::Type{InEquilibrium})

    ipsi = data.index_psi
    # RHS of Poisson
    RHSPoisson!(f, u, node, data, ipsi)
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        RHSPoisson!(f, u, node, data, ipsiStandard)
    end

    # zero reaction term for all icc (stability purpose)
    for icc in data.electricCarrierList # Array{Int64, 1}
        icc = data.chargeCarrierList[icc] # Array{QType 1}
        f[icc] = u[icc]
    end

    for iicc in data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions
            icc = iicc.ionicCarrier           # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            f[icc] = u[icc]
        end
    end
    for iicc in data.trapCarrierList
        # add trap carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions
            icc = iicc.trapCarrier            # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            f[icc] = u[icc]
        end
    end

    return
end

function StimulatedRecombination(u, node, data)

    params = data.params
    paramsoptical = data.paramsoptical
    ireg = node.region
    (; k_B, q, Planck_constant) = data.constants

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]
    ipsi = data.index_psi

    n = get_density!(u, node, data, iphin)
    p = get_density!(u, node, data, iphip)

    hbar = Planck_constant / (2 * pi)
    c0 = 299_792_458
    k0 = 2 * pi / paramsoptical.laserWavelength
    ω0 = k0 * c0
    kBT = k_B * params.temperature

    Ec = get_BEE!(iphin, node, data)
    Ev = get_BEE!(iphip, node, data)

    n0 = paramsoptical.refractiveIndex_0[ireg]
    nd = paramsoptical.refractiveIndex_d[ireg]
    γn = paramsoptical.refractiveIndex_γ[ireg]

    g0 = paramsoptical.gain_0[ireg]

    eValue = paramsoptical.eigenvalues[1]
    beta = sqrt(-eValue)
    eVector = paramsoptical.eigenvectors[node.index, 1]

    power = paramsoptical.power
    expTerm1 = exp((-q * u[iphin] - Ec + q * u[ipsi]) / kBT)
    expTerm2 = exp((Ev + q * u[iphip] - q * u[ipsi]) / kBT)
    expTerm3 = exp(((-q * (u[iphin] - u[iphip])) - (hbar * ω0)) / kBT) - 1
    gainDenominator = (1 + expTerm1) * (1 + expTerm2)
    gain = (g0 / gainDenominator) * expTerm3

    refractive = n0 - (nd * ((n + p) / 2))^γn
    RstimValue = ((refractive * gain) / (hbar * ω0)) * power * (((abs.(eVector)) .^ 2) / (real(beta) / k0))

    return RstimValue

end


function addRecombination!(f, u, node, data)

    params = data.params
    ireg = node.region
    (; q, k_B) = data.constants

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]

    n = get_density!(u, node, data, iphin)
    p = get_density!(u, node, data, iphip)

    taun = params.recombinationSRHLifetime[iphin, ireg]
    n0 = params.recombinationSRHTrapDensity[iphin, ireg]
    taup = params.recombinationSRHLifetime[iphip, ireg]
    p0 = params.recombinationSRHTrapDensity[iphip, ireg]

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (k_B * data.params.temperature))
    excessDensTerm = n * p * (1.0 - exponentialTerm)

    # calculate recombination kernel. If user adjusted Auger, radiative or SRH recombination,
    # they are set to 0. Hence, adding them here, has no influence since we simply add by 0.0.
    kernelRad = params.recombinationRadiative[ireg]
    kernelAuger = (params.recombinationAuger[iphin, ireg] * n + params.recombinationAuger[iphip, ireg] * p)
    kernelSRH = params.prefactor_SRH / (taup * (n + n0) + taun * (p + p0))
    kernel = kernelRad + kernelAuger + kernelSRH

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################
    f[iphin] = q * params.chargeNumbers[iphin] * kernel * excessDensTerm
    f[iphip] = q * params.chargeNumbers[iphip] * kernel * excessDensTerm
    return nothing
end

function addStimulatedRecombination!(f, u, node, data, ::Type{LaserModelOff})
    return nothing
end

function addStimulatedRecombination!(f, u, node, data, ::Type{LaserModelOn})

    params = data.params
    q = data.constants.q

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # calculate stimulatedRecombination
    stimulatedRecombination = StimulatedRecombination(u, node, data)

    f[iphin] = f[iphin] + q * params.chargeNumbers[iphin] * stimulatedRecombination
    f[iphip] = f[iphip] + q * params.chargeNumbers[iphip] * stimulatedRecombination

    return nothing
end

function addGeneration!(f, u, node, data)

    generationTerm = generation(data, node, data.generationModel)

    for icc in data.electricCarrierList
        icc = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used and depicted here
        f[icc] = f[icc] - data.constants.q * data.params.chargeNumbers[icc] * generationTerm
    end

    return
end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of Poisson equation, i.e. which builds
the space charge density.
"""
function RHSPoisson!(f, u, node, data, ipsi)

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################

    # electrons and holes entering right hand-side of Poisson in each layer
    for icc in data.electricCarrierList          # Array{Int64, 1}

        icc = data.chargeCarrierList[icc]   # Array{QType, 1}
        ncc = get_density!(u, node, data, icc)

        f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * (data.params.doping[icc, node.region])  # subtract doping
        f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier

    end

    for iicc in data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions

            icc = iicc.ionicCarrier           # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})
            ncc = get_density!(u, node, data, icc)

            f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * (data.params.doping[icc, node.region])  # subtract doping
            f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier
        end
    end
    for iicc in data.trapCarrierList
        # add trap carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions

            icc = iicc.trapCarrier           # species number chosen by user
            icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            ncc = get_density!(u, node, data, icc)

            f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * (data.params.doping[icc, node.region])  # subtract doping
            f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier
        end
    end

    f[ipsi] = f[ipsi] - data.paramsnodal.doping[node.index]
    f[ipsi] = - data.constants.q * data.λ1 * f[ipsi]

    return

end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of electric charge carriers.
"""
function RHSContinuityEquations!(f, u, node, data)

    # dependent on user information concerning recombination
    if data.bulkRecombination.bulk_recomb
        addRecombination!(f, u, node, data)
    end
    # dependent on user information concerning trap dynamics
    addTrapCaptureEscape!(f, u, node, data)
    # dependent on user information concerning laser model
    addStimulatedRecombination!(f, u, node, data, data.laserModel)
    # dependent on user information concerncing generation
    addGeneration!(f, u, node, data)
    return nothing

end


"""
$(TYPEDSIGNATURES)
Sets up the right-hand sides. Assuming a bipolar semiconductor
the right-hand side for the electrostatic potential becomes
  ``f[ψ]  = - q ((p - N_a) - (n - N_d) ) = - q  \\sum  n_\\alpha  (n_\\alpha - C_\\alpha) ``
for some doping ``C_\\alpha`` w.r.t. to the species ``\\alpha``.
The right-hand sides for the charge carriers read as
``f[n_\\alpha] =  - z_\\alpha  q (G -  R) ``
for all charge carriers ``n_\\alpha``.
The recombination includes radiative, Auger and Shockley-Read-Hall
recombination. For latter recombination process the stationary simplification is implemented.
The recombination is only implemented for electron and holes and assumes
that the electron index is 1 and the hole index is 2.
"""
function reaction!(f, u, node, data, ::Type{OutOfEquilibrium})

    ipsi = data.index_psi
    RHSPoisson!(f, u, node, data, ipsi)               # RHS of Poisson
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn # additional Poisson in case of Barrier Lowering
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        RHSPoisson!(f, u, node, data, ipsiStandard)
    end

    # First, set RHS to zero for all icc
    for icc in eachindex(data.chargeCarrierList) # Array{Int61, 1}
        icc = data.chargeCarrierList[icc]   # Array{QType, 1}
        f[icc] = 0.0
    end

    # Then, add RHS of continuity equations based on user information
    RHSContinuityEquations!(f, u, node, data) # RHS of Charge Carriers with special treatment of recombination

    return

end


"""
$(SIGNATURES)
Compute trap densities for a given trap energy.
[Currently, only done for the Boltzmann statistics and for region dependent parameters.]
"""
function trap_density(icc, ireg, params, Et, constants)

    return params.densityOfStates[icc, ireg] * exp(params.chargeNumbers[icc] * (params.bandEdgeEnergy[icc, ireg] - Et) / (constants.k_B * params.temperature))
end

# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations with a uniform generation rate.
function generation(data, node, ::Type{GenerationUniform})

    return data.λ2 * data.params.generationUniform[node.region]
end


# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations obeying the Beer-Lambert law.
# only works in 1D till now; adjust node, when multidimensions
function generation(data, node, ::Type{GenerationBeerLambert})

    params = data.params
    ireg = node.region
    node = node.coord[node.index]

    return data.λ2 .* params.generationIncidentPhotonFlux[ireg] .* params.generationAbsorption[ireg] .* exp.(- params.invertedIllumination .* params.generationAbsorption[ireg] .* (node .- params.generationPeak))

end


# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations with a user defined generation rate.
# only works in 1D till now; adjust node, when multidimensions
function generation(data, node, ::Type{GenerationUserDefined})

    return data.λ2 .* data.generationData[node.index]

end

generation(data, node, ::Type{GenerationNone}) = 0.0

"""
$(SIGNATURES)
Beer-Lambert function for the visualization of this type of photogeneration profile.
"""

function BeerLambert(ctsys, ireg, node)

    data = ctsys.fvmsys.physics.data
    params = data.params

    return params.generationIncidentPhotonFlux[ireg] .* params.generationAbsorption[ireg] .* exp.(- params.invertedIllumination .* params.generationAbsorption[ireg] .* (node .- params.generationPeak))

end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master storage! function. This is the function which enters VoronoiFVM and hands over
a storage term, if we consider transient problem.
"""
storage!(f, u, node, data) = storage!(f, u, node, data, data.modelType)

storage!(f, u, node, data, ::Type{Stationary}) = emptyFunction()


storage!(f, u, node, data, ::Type{Transient}) = storage!(f, u, node, data, data.calculationType)

storage!(f, u, node, data, ::Type{InEquilibrium}) = emptyFunction()
"""
$(TYPEDSIGNATURES)
The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have
``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha``
and for the electrostatic potential
``f[ψ] = 0``.
"""
function storage!(f, u, node, data, ::Type{OutOfEquilibrium})

    params = data.params
    ipsi = data.index_psi
    q = data.constants.q

    for icc in data.electricCarrierList       # Array{Int64, 1}

        icc = data.chargeCarrierList[icc] # get correct index in chargeCarrierList
        ncc = get_density!(u, node, data, icc)
        f[icc] = q * params.chargeNumbers[icc] * ncc

    end

    for iicc in data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # Here we do not need to check, if carrier is present in a specific region.
        # This is directly handled by VoronoiFVM.
        icc = iicc.ionicCarrier           # species number chosen by user
        icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

        ncc = get_density!(u, node, data, icc)
        f[icc] = q * params.chargeNumbers[icc] * ncc
    end
    for iicc in data.trapCarrierList
        icc = iicc.trapCarrier            # species number chosen by user
        icc = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

        ncc = get_density!(u, node, data, icc)
        f[icc] = q * params.chargeNumbers[icc] * ncc
    end

    f[ipsi] = 0.0

    return

end


function DensityProduct(f, u, node, data)

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]

    return f[1] = get_density!(u, node, data, iphin) * get_density!(u, node, data, iphip)

end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master flux functions which enters VoronoiFVM. Flux discretization scheme is chosen in two steps. First, we need
to see, if we are in or out of equilibrium. If, InEquilibrium, then
no flux is passed. If outOfEquilibrium, we choose the flux approximation
which the user chose for each charge carrier. For the displacement flux we use a finite difference approach.
"""
flux!(f, u, edge, data) = flux!(f, u, edge, data, data.calculationType)


# Finite difference discretization of the displacement flux.
function displacementFlux!(f, u, edge, data)

    params = data.params
    paramsnodal = data.paramsnodal

    ipsi = data.index_psi
    nodel = edge.node[2]   # left node
    nodek = edge.node[1]   # right node
    ireg = edge.region

    dpsi = u[ipsi, 2] - u[ipsi, 1]

    dielConst = params.dielectricConstant[ireg] + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek]) / 2
    f[ipsi] = - dielConst * dpsi

    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn # additional Poisson in case of Barrier Lowering
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        f[ipsiStandard] = - dielConst * (u[ipsiStandard, 2] - u[ipsiStandard, 1])
    end

    return

end


function flux!(f, u, edge, data, ::Type{InEquilibrium})
    ## discretization of the displacement flux (LHS of Poisson equation)
    displacementFlux!(f, u, edge, data)
    return
end

function flux!(f, u, edge, data, ::Type{OutOfEquilibrium})

    ## discretization of the displacement flux (LHS of Poisson equation)
    displacementFlux!(f, u, edge, data)

    for icc in data.electricCarrierList   # correct index of electric carriers of Type Int64
        chargeCarrierFlux!(f, u, edge, data, icc, data.fluxApproximation[icc])
    end

    for icc in data.ionicCarrierList
        icc = icc.ionicCarrier       # correct index number chosen by user of Type Int64
        chargeCarrierFlux!(f, u, edge, data, icc, data.fluxApproximation[icc])
    end

    return
end


# The classical Scharfetter-Gummel flux scheme. This also works for space-dependent
# band-edge energy, but not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ScharfetterGummel})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; q, k_B) = data.constants

    j0 = k_B * params.temperature / q * params.mobility[icc, ireg]

    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    bp, bm = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi * q - bandEdgeDiff) / (k_B * params.temperature))
    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end

# The classical Scharfetter-Gummel flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the
# discretization scheme is modified.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ScharfetterGummelGraded})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants

    mobility = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek]) / 2
    j0 = (k_B * params.temperature / q) * mobility

    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
        bp, bm = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi * q - bandEdgeDiff) / (k_B * params.temperature))
    else
        bp, bm = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi * q - bandEdgeDiff) / (k_B * params.temperature) - (log(paramsnodal.densityOfStates[icc, nodel]) - log(paramsnodal.densityOfStates[icc, nodek])))
    end

    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end

# The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ExcessChemicalPotential})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants

    j0 = (k_B * params.temperature / q) * params.mobility[icc, ireg]

    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal = etaFunction!(u, edge, data, icc)

    Q = params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff) / (k_B * params.temperature)) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak))
    bp, bm = fbernoulli_pm(Q)

    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end

# Reconstructing the concentration gradients
function ConcentrationGradient(f, u, edge, data)

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]

    nnk, nnl = get_density!(u, edge, data, iphin)
    npk, npl = get_density!(u, edge, data, iphip)

    f[iphin] = - data.params.chargeNumbers[iphin] * (nnl - nnk)
    f[iphip] = - data.params.chargeNumbers[iphip] * (npl - npk)
    return

end

# The excess chemical potential flux discretization scheme without force term for postprocessing
function ExcessChemicalPotentialDiffusive(f, u, edge, data)

    params = data.params
    paramsnodal = data.paramsnodal

    for icc in data.chargeCarrierList

        nodek = edge.node[1]   # left node
        nodel = edge.node[2]   # right node
        ireg = edge.region
        (; k_B, q) = data.constants

        j0 = (k_B * params.temperature / q) * params.mobility[icc, ireg]

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        etak, etal = etaFunction!(u, edge, data, icc)

        Q = params.chargeNumbers[icc] * ((- bandEdgeDiff) / (k_B * params.temperature)) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak))
        bp, bm = fbernoulli_pm(Q)

        ncck, nccl = get_density!(u, edge, data, icc)

        f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    end

    return
end

# The excess chemical potential flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the discretization scheme is modified.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ExcessChemicalPotentialGraded})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants


    mobility = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek]) / 2
    j0 = (k_B * params.temperature / q) * mobility

    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal = etaFunction!(u, edge, data, icc)

    if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
        Q = params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff) / (k_B * params.temperature)) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak))
    else
        Q = params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff) / (k_B * data.temperature)) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(paramsnodal.densityOfStates[icc, nodel]) - log(paramsnodal.densityOfStates[icc, nodek])))
    end

    bp, bm = fbernoulli_pm(Q)
    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end

# The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme
# is used for the regularization of the removable singularity. This also works for
# space-dependent band-edge energy, but not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{DiffusionEnhanced})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants


    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    if (log(data.F[icc](etal)) - log(data.F[icc](etak))) ≈ 0.0 # regularization idea coming from Pietra-Jüngel scheme
        gk = exp(etak) / data.F[icc](etak)
        gl = exp(etal) / data.F[icc](etal)
        g = 0.5 * (gk + gl)
    else
        g = (etal - etak) / (log(data.F[icc](etal)) - log(data.F[icc](etak)))
    end

    j0 = (k_B * params.temperature / q) * params.mobility[icc, ireg] * g

    bp, bm = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi * q - bandEdgeDiff) / (k_B * params.temperature * g))
    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end

# The diffusion enhanced scheme by Bessemoulin-Chatard for fluxes not based on a nonlinear diffusion
# but on a modified drift.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{DiffusionEnhancedModifiedDrift})

    params = data.params
    paramsnodal = data.paramsnodal

    icc = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants


    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    if (log(data.F[icc](etal)) - log(data.F[icc](etak))) ≈ 0.0 # regularization idea coming from Pietra-Jüngel scheme
        gk = exp(etak) / data.F[icc](etak)
        gl = exp(etal) / data.F[icc](etal)
        g = 0.5 * (gk + gl)
    else
        g = (etal - etak) / (log(data.F[icc](etal)) - log(data.F[icc](etak)))
    end

    j0 = (k_B * params.temperature / q) * params.mobility[icc, ireg]

    bp, bm = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi * q - bandEdgeDiff) / (k_B * params.temperature * g))
    ncck, nccl = get_density!(u, edge, data, icc)

    f[icc] = - params.chargeNumbers[icc] * q * j0 * (bm * nccl - bp * ncck)

    return

end


# # The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation
# # which arise when considering the generalized Scharfetter-Gummel scheme in case of Blakemore
# # statistics. Hence, it should be exclusively worked with, when considering the Blakemore
# # statistics. This also works for space-dependent band-edge energy, but not for
# # space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{GeneralizedSG})

    max_iter = 300          # for Newton solver
    it = 0            # number of iterations (newton)
    damp = 0.1          # damping factor

    params = data.params
    paramsnodal = data.paramsnodal

    # DA: we get issues with allocations, when allowing non Integer icc.
    #icc          = data.chargeCarrierList[icc]
    ipsi = data.index_psi
    nodek = edge.node[1]   # left node
    nodel = edge.node[2]   # right node
    ireg = edge.region
    (; k_B, q) = data.constants

    j0 = (k_B * params.temperature / q) * params.mobility[icc, ireg]

    dpsi = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
    etak, etal = etaFunction!(u, edge, data, icc)

    # use Sedan flux as starting guess
    Q = params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff) / (k_B * params.temperature)) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak))
    bp, bm = fbernoulli_pm(Q)
    ncck, nccl = get_density!(u, edge, data, icc)

    jInitial = (bm * nccl - bp * ncck)

    implicitEq(j::Real) = (fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff)) / (k_B * params.temperature) + params.γ * j)[2] * exp(etal) - fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi * q - bandEdgeDiff) / (k_B * params.temperature)) - params.γ * j)[1] * exp(etak)) - j

    delta = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
    oldup = 1.0
    while (it < max_iter)
        Fval = implicitEq(jInitial)
        dFval = ForwardDiff.derivative(implicitEq, jInitial)

        if isnan(value(dFval)) || value(abs(dFval)) < delta
            @show value(jInitial), value(Fval), value(dFval)
            error("singular derivative in exact SG scheme")
        end

        update = Fval / dFval
        jInitial = jInitial - damp * update

        if abs(update) < delta
            break
        end
        oldup = value(update)

        it = it + 1
        damp = min(damp * 1.2, 1.0)
    end

    f[icc] = - params.chargeNumbers[icc] * q * j0 * jInitial

    return

end

# ##########################################################
# recombination kernels for calculating recombination currents
# ##########################################################


function SRRecombination!(f, u, bnode, data)

    params = data.params

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    n = get_density!(u, bnode, data, iphin)
    p = get_density!(u, bnode, data, iphip)

    (; k_B, q) = data.constants


    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (k_B * params.temperature))
    excessDensTerm = n * p * (1.0 - exponentialTerm)

    if params.recombinationSRHvelocity[iphip, bnode.region] ≈ 0.0
        vp = 1.0e30
    else
        vp = 1.0 / params.recombinationSRHvelocity[iphip, bnode.region]
    end

    if params.recombinationSRHvelocity[iphin, bnode.region] ≈ 0.0
        vn = 1.0e30
    else
        vn = 1.0 / params.recombinationSRHvelocity[iphin, bnode.region]
    end

    kernelSRH = 1.0 / (vp * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + vn * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region]))

    for icc in data.electricCarrierList
        icc = data.chargeCarrierList[icc]
        f[icc] = q * params.chargeNumbers[icc] * kernelSRH * excessDensTerm
    end

    return
end

function SRHRecombination!(f, u, node, data)

    params = data.params
    ireg = node.region
    (; k_B, q) = data.constants


    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]

    n = get_density!(u, node, data, iphin)
    p = get_density!(u, node, data, iphip)

    taun = params.recombinationSRHLifetime[iphin, ireg]
    n0 = params.recombinationSRHTrapDensity[iphin, ireg]
    taup = params.recombinationSRHLifetime[iphip, ireg]
    p0 = params.recombinationSRHTrapDensity[iphip, ireg]

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (k_B * data.params.temperature))
    excessDensTerm = n * p * (1.0 - exponentialTerm)

    kernelSRH = params.prefactor_SRH / (taup * (n + n0) + taun * (p + p0))
    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################
    f[iphin] = q * params.chargeNumbers[iphin] * kernelSRH * excessDensTerm
    f[iphip] = q * params.chargeNumbers[iphip] * kernelSRH * excessDensTerm

    return

end


function RadiativeRecombination!(f, u, node, data)

    params = data.params
    ireg = node.region
    (; k_B, q) = data.constants


    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin = data.chargeCarrierList[iphin]
    iphip = data.chargeCarrierList[iphip]

    n = get_density!(u, node, data, iphin)
    p = get_density!(u, node, data, iphip)

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (k_B * data.params.temperature))
    excessDensTerm = n * p * (1.0 - exponentialTerm)

    # calculate recombination kernel. If user adjusted Auger, radiative or SRH recombination,
    # they are set to 0. Hence, adding them here, has no influence since we simply add by 0.0.
    kernelRad = params.recombinationRadiative[ireg]
    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################
    f[iphin] = q * params.chargeNumbers[iphin] * kernelRad * excessDensTerm
    f[iphip] = q * params.chargeNumbers[iphip] * kernelRad * excessDensTerm

    return

end

function Photogeneration!(f, u, node, data)

    generationTerm = generation(data, node, data.generationModel)

    for icc in data.electricCarrierList
        icc = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used and depicted here
        f[icc] = data.constants.q * data.params.chargeNumbers[icc] * generationTerm
    end

    return
end

"""
$(TYPEDSIGNATURES)
No trap: Do nothing
"""
function addTrapCaptureEscape!(f, u, node, data, ::Type{NoTrap})
    return nothing
end
"""
$(TYPEDSIGNATURES)
A simple trap with one state that can either be filled or empty
"""
function addTrapCaptureEscape!(f, u, node, data, ::Type{SingleStateTrap})

    (; k_B, q) = data.constants

    capture = data.params.recombinationTrapCaptureRates
    T = data.params.temperature

    for icc in data.chargeCarrierList
        ncc = get_density!(u, node, data, icc)
        Nc = data.params.densityOfStates[icc]

        # Account for non-Boltzmann statistics in detailed balance to compute escape rate.
        # It is assumed that the trap is described using FermiDiracMinusOne. The correction
        # is of F(η)/exp(η) where F is the function used in the carrier state equation.
        if (data.F[icc] == Boltzmann)
            nonBoltzmannReductionFactor = 1.0
        else
            nonBoltzmannReductionFactor = ncc / (Nc * exp(etaFunction!(u, node, data, icc)))
        end

        Ec = data.params.bandEdgeEnergy[icc]
        zc = data.params.chargeNumbers[icc]

        for iitc in data.trapCarrierList
            # add trap carriers only in defined regions (otherwise get NaN error)
            if node.region ∈ iitc.regions
                itc = iitc.trapCarrier            # species number chosen by user
                itc = data.chargeCarrierList[itc] # find correct index within chargeCarrierList
                s = capture[itc, icc, node.region]

                if (s > 0) # Only compute where there is capture
                    zt = data.params.chargeNumbers[itc]

                    ntc = get_density!(u, node, data, itc)
                    Nt = data.params.densityOfStates[itc]
                    Et = data.params.bandEdgeEnergy[itc]

                    # Escape computed from detailed balance assuming traps described using FD-minus one
                    e = s * Nc * exp(zc * (Ec - Et) / (k_B * T)) * nonBoltzmannReductionFactor

                    # Allow for both acceptor and donor trap in one line
                    # e.g.  If acceptor traps (trap charge = -1) then reaction with
                    # conduction band is  r = Nt*( s*n*(1-f) - e*f ). Reaction with
                    # the valence band is r = Nt*( s*p*f - e*(1-f) ).
                    # For donor traps the (1-f) and f swaps, which is done
                    # using sign(zc*zt).
                    captureFactor = (sign(zc * zt) + 1) / 2 - sign(zc * zt) * ntc / Nt
                    escapeFactor = 1 - captureFactor
                    r = Nt * (s * ncc * captureFactor - e * escapeFactor)

                    # For the reaction expression we use the charge of the band as (e.g.) holes can enter
                    # an electron trap from the valence band, and using the trap charge would not capture this.
                    f[icc] = f[icc] + q * zc * r    #
                    f[itc] = f[itc] - q * zc * r    #
                end
            end
        end
    end

    return
end

"""
$(TYPEDSIGNATURES)
Include recombination between bands and traps
"""
function addTrapCaptureEscape!(f, u, node, data)
    return addTrapCaptureEscape!(f, u, node, data, data.bulkRecombination.bulk_recomb_trap)
end
