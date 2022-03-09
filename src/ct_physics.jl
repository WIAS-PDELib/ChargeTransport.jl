##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the effective DOS for interior nodes (analogously for boundary nodes and edges).
"""
function get_DOS!(icc::QType, node::VoronoiFVM.Node, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, node.region] + data.paramsnodal.densityOfStates[icc, node.index]

end

# Defining locally the effective DOS for boundary nodes.
function get_DOS!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempDOS1[icc] = data.params.bDensityOfStates[icc, bnode.region] + data.paramsnodal.densityOfStates[icc, bnode.index]

end

# Defining locally the effective DOS for edges.
function get_DOS!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[1]]
    data.tempDOS2[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[2]]

    return data.tempDOS1, data.tempDOS2
end

"""
$(TYPEDSIGNATURES)

Defining locally the band-edge energy for interior nodes (analougesly for boundary nodes and edges).
"""
function get_BEE!(icc::QType, node::VoronoiFVM.Node, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, node.region] + data.paramsnodal.bandEdgeEnergy[icc, node.index]

end

# Defining locally the band-edge energy for boundary nodes.
function get_BEE!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempBEE1[icc] = data.params.bBandEdgeEnergy[icc, bnode.region] + data.paramsnodal.bandEdgeEnergy[icc, bnode.index]

end

# Defining locally the band-edge energy for edges.
function get_BEE!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[1]]
    data.tempBEE2[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[2]]

    return data.tempBEE1, data.tempBEE2

end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the distribution function for interior nodes.
"""
function etaFunction(u, node::VoronoiFVM.Node, data, icc)

    get_BEE!(icc, node::VoronoiFVM.Node, data)

    E  = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[data.index_psi]) + E / q )

end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for floats.
"""
function etaFunction(u, data, node, region, icc, in_region::Bool)

    if in_region == true
        E  = data.params.bandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    elseif in_region == false
        E  = data.params.bBandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    end

    return etaFunction(u[data.index_psi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc])
end

"""
$(TYPEDSIGNATURES)

The argument of the distribution function for boundary nodes.
"""
function etaFunction(u, bnode::VoronoiFVM.BNode, data, icc) # bnode.index refers to index in overall mesh

    get_BEE!(icc, bnode::VoronoiFVM.BNode, data)
    E  = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[data.index_psi]) + E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for edges.
"""

function etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

    get_BEE!(icc, edge::VoronoiFVM.Edge, data)

    E1 = data.tempBEE1[icc];  E2 = data.tempBEE2[icc]

    return etaFunction(u[data.index_psi, 1], u[icc, 1], data.params.UT, E1, data.params.chargeNumbers[icc]), etaFunction(u[data.index_psi, 2], u[icc, 2], data.params.UT, E2, data.params.chargeNumbers[icc])
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for given ``\\varphi_\\alpha``
and ``\\psi``

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q ).``

The parameters ``E_\\alpha`` and ``z_\\alpha`` are given as vectors.
This function may be used to compute the charge density, i.e. the
right-hand side of the Poisson equation.
"""
function etaFunction(psi, phi, UT, E, z)
    z ./ UT .* ( (phi - psi) .+ E / q )
end

##########################################################
##########################################################

function emptyFunction()
end


"""
$(TYPEDSIGNATURES)
Master breaction! function. This is the function which enters VoronoiFVM and hands over
for each boundary the chosen boundary model.

"""
breaction!(f, u, bnode, data) =  breaction!(f, u, bnode, data, data.boundaryType[bnode.region])


"""
$(TYPEDSIGNATURES)

Creates ohmic boundary conditions via a penalty approach with penalty parameter ``\\delta``.
For example, the right-hand side for the electrostatic potential ``\\psi`` is implemented as

``f[\\psi]  = -q/\\delta   ( (p - N_a) - (n - N_d) )``,

assuming a bipolar semiconductor. In general, we have for some given charge number ``z_\\alpha``

``f[\\psi] =  -q/\\delta  \\sum_\\alpha{ z_\\alpha  (n_\\alpha - C_\\alpha) },``

where ``C_\\alpha`` corresponds to some doping w.r.t. the species ``\\alpha``.

The boundary conditions for the charge carriers are set in the main file. Hence,

``f[n_\\alpha] = 0```

for all charge carriers ``\\alpha``.
"""
function breaction!(f, u, bnode, data, ::Type{OhmicContact})

    params      = data.params
    paramsnodal = data.paramsnodal

    # parameters
    ipsi        = data.index_psi  # final index for electrostatic potential

    for icc ∈ data.chargeCarrierList[1:2] # quantities or integer indices

        get_DOS!(icc, bnode, data)
        Ni      = data.tempDOS1[icc]
        eta     = etaFunction(u, bnode, data, icc) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc)

        # subtract doping
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )
        # add charge carrier
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * Ni * data.F[icc](eta)

        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0

    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]

    f[ipsi] = - data.λ1 * 1 / tiny_penalty_value *  q * f[ipsi]

end

"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions in a first attempt. For the electrostatic potential we assume

``\\psi = \\phi_S/q + U, ``

where  ``\\phi_S`` corresponds to a given value (Schottky barrier) and ``U`` to the applied voltage. For now,
the quantitity ``\\phi_S`` needs to be specified in the main file.
For the charge carriers we assume the following

``f[n_\\alpha]  =  z_\\alpha q v_\\alpha (n_\\alpha - n_{\\alpha, 0})``,

where ``v_{\\alpha}`` can be treated as a surface recombination mechanism and is given. The parameter
``n_{\\alpha, 0}`` is a given value, calculated by the statistical relation, when assuming
no electrical field and a quasi Fermi level equal to the Schottky barrier ``\\phi_S``, i.e.

``n_{\\alpha, 0}= z_\\alpha/ U_T (E_\\alpha - \\phi_S) / q. ``

"""
function breaction!(f, u, bnode, data,  ::Type{SchottkyContact})

    params        = data.params
    paramsnodal   = data.paramsnodal

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin
    iphip       = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    looplist    = (iphin, iphip)
    ipsi        = data.index_psi

    for icc in 1:length(looplist)

        get_DOS!(icc, bnode, data);  get_BEE!(icc, bnode, data)
        Ni     = data.tempDOS1[icc]
        Ei     = data.tempBEE1[icc]
        etaFix = params.chargeNumbers[icc] / params.UT * (  (- params.SchottkyBarrier[bnode.region] + Ei ) / q  )
        eta    = params.chargeNumbers[icc] / params.UT * (  (u[icc]  - u[ipsi]) + Ei / q )

        f[icc] =  data.λ1 * params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * (  Ni  * (data.F[icc](eta) - data.F[icc](etaFix)  ))

    end

end


# This breaction! function is chosen when no interface model is chosen.
breaction!(f, u, bnode, data, ::Type{InterfaceModelNone}) = emptyFunction()


# breaction term for surface recombination.
breaction!(f, u, bnode, data, ::Type{InterfaceModelSurfaceRecoAndTangentialFlux}) = breaction!(f, u, bnode, data, InterfaceModelSurfaceReco)

# breaction term for tangential flux.
breaction!(f, u, bnode, data, ::Type{InterfaceModelTangentialFlux}) = emptyFunction()


function breaction!(f, u, bnode, data, ::Type{InterfaceModelSurfaceReco})
    if data.calculationType == InEquilibrium
        return
    end

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin # integer index of φ_n
    iphip       = data.bulkRecombination.iphip # integer index of φ_p
    looplist    = (iphin, iphip)

    ipsi        = data.index_psi

    params      = data.params
    paramsnodal = data.paramsnodal

    get_DOS!(iphin, bnode, data); get_DOS!(iphip, bnode, data)
    Nc   = data.tempDOS1[iphin]
    Nv   = data.tempDOS1[iphip]
    etan = etaFunction(u, bnode, data, iphin) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc)
    etap = etaFunction(u, bnode, data, iphip) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc)

    n    = Nc * data.F[iphin](etan)
    p    = Nv * data.F[iphip](etap)

    exponentialTerm = exp((q * u[iphin] - q  * u[iphip] ) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    kernelSRH = 1.0 / (  1.0/params.recombinationSRHvelocity[iphip, bnode.region] * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + 1.0/params.recombinationSRHvelocity[iphin, bnode.region] * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region] ) )

    for icc in 1:length(looplist)
        f[icc] =  q * params.chargeNumbers[icc] * kernelSRH *  excessDensTerm
    end


end


function breaction!(f, u, bnode, data, ::Type{InterfaceModelDiscontqFNoReaction})

    if data.calculationType == InEquilibrium
        return emptyFunction()
    end

    ipsi = data.index_psi

    #indices (∈ N) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin # integer index of φ_n
    iphip       = data.bulkRecombination.iphip # integer index of φ_p

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin] # = Quantity or integer
    iphip       = data.chargeCarrierList[iphip] # = Quantity or integer

    params      = data.params
    paramsnodal = data.paramsnodal

    for icc in (iphin, iphip)

        xx         =  [1.0e14 1.0e14; 1.0e14 1.0e14]

        etan1 = params.chargeNumbers[icc] / params.UT * ( (u[icc, 1] - u[ipsi]) + params.bandEdgeEnergy[icc, bnode.cellregions[1]] / q ) # left
        etan2 = params.chargeNumbers[icc] / params.UT * ( (u[icc, 2] - u[ipsi]) + params.bandEdgeEnergy[icc, bnode.cellregions[1]] / q ) # right

        n1 = params.densityOfStates[icc, bnode.cellregions[1]] * data.F[icc](etan1)
        n2 = params.densityOfStates[icc, bnode.cellregions[1]] * data.F[icc](etan2)

        react     = q * params.chargeNumbers[icc] * xx[icc, bnode.region-2]  * (n1 - n2)

        f[icc, 1] =   react
        f[icc, 2] = - react
    end


end


# breaction term for case where qF are discontinuous.
function breaction!(f, u, bnode, data, ::Type{InterfaceModelDiscontqF})

    if data.calculationType == InEquilibrium
        return emptyFunction()
    end

    ipsi = data.index_psi

    #indices (∈ N) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin # integer index of φ_n
    iphip       = data.bulkRecombination.iphip # integer index of φ_p

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin] # = Quantity or integer
    iphip       = data.chargeCarrierList[iphip] # = Quantity or integer

    params      = data.params
    paramsnodal = data.paramsnodal

    iphin_b1 = data.interfaceCarriers.interfaceIndex[iphin]
    iphip_b1 = data.interfaceCarriers.interfaceIndex[iphip]

    iphin_b1 = data.chargeCarrierList[iphin_b1]
    iphip_b1 = data.chargeCarrierList[iphip_b1]

    zn  = params.chargeNumbers[iphin]
    zp  = params.chargeNumbers[iphip]
    UT  = params.UT

    # take values from intrinsic layer
    Nc_l  = params.densityOfStates[iphin, 1]
    Nv_l  = params.densityOfStates[iphip, 1]
    Ec_l  = params.bandEdgeEnergy[iphin, 1]
    Ev_l  = params.bandEdgeEnergy[iphip, 1]

    Nc_r  = params.densityOfStates[iphin, 2]
    Nv_r  = params.densityOfStates[iphip, 2]
    Ec_r  = params.bandEdgeEnergy[iphin, 2]
    Ev_r  = params.bandEdgeEnergy[iphip, 2]

    mun = params.mobility[iphin, 2]
    mup = params.mobility[iphip, 2]

    Nc_b = params.bDensityOfStates[iphin_b1, bnode.region]
    Nv_b = params.bDensityOfStates[iphip_b1, bnode.region]
    Ec_b = params.bBandEdgeEnergy[iphin_b1, bnode.region]
    Ev_b = params.bBandEdgeEnergy[iphip_b1, bnode.region]

    # since k0 is scaled with q*z_i it's consistent with the storage and the bstorage dimensions
    k0n = q * zn * UT * mun * Nc_r/data.d
    k0p = q * zp * UT * mup * Nv_r/data.d

    # left values
    etan1 = zn / UT * ( (u[iphin, 1] - u[ipsi]) + Ec_l / q ) # left
    etap1 = zp / UT * ( (u[iphip, 1] - u[ipsi]) + Ev_l / q ) # left

    n1    = Nc_l * data.F[iphin](etan1)
    p1    = Nv_l * data.F[iphip](etap1)

    # interface values
    etan_b1 = zn / UT * ( (u[iphin_b1] - u[ipsi]) + Ec_b / q ) # interface
    etap_b1 = zp / UT * ( (u[iphip_b1] - u[ipsi]) + Ev_b / q ) # interface

    n_b1    = Nc_b * data.F[iphin](etan_b1)
    p_b1    = Nv_b * data.F[iphip](etap_b1)

    # right values
    etan2 = zn / UT * ( (u[iphin, 2] - u[ipsi]) + Ec_r / q ) # right
    etap2 = zp / UT * ( (u[iphip, 2] - u[ipsi]) + Ev_r / q ) # right

    n2    = Nc_r * data.F[iphin](etan2)
    p2    = Nv_r * data.F[iphip](etap2)

    Knleft  = exp(- zn/ (kB * data.params.temperature) * (Ec_l - Ec_b))
    Knright = exp(- zn/ (kB * data.params.temperature) * (Ec_r - Ec_b))
    Kpleft  = exp(- zp/ (kB * data.params.temperature) * (Ev_l - Ev_b))
    Kpright = exp(- zp/ (kB * data.params.temperature) * (Ev_r - Ev_b))

    ##############################################################
    reactn1  = k0n * (Knleft^(1/2) * n1/Nc_l - Knleft^(- 1/2) * n_b1/Nc_b)
    reactn2  = k0n * (Knright^(1/2) * n2/Nc_r - Knright^(- 1/2) * n_b1/Nc_b)

    f[iphin, 1] =   reactn1
    f[iphin, 2] =   reactn2
    f[iphin_b1] = - reactn1 - reactn2

    reactp1     = k0p * (Kpleft^(1/2) * p1/Nv_l - Kpleft^(- 1/2) * p_b1/Nv_b)
    reactp2     = k0p * (Kpright^(1/2) * p2/Nv_r - Kpright^(- 1/2) * p_b1/Nv_b)

    f[iphip, 1] =   reactp1
    f[iphip, 2] =   reactp2
    f[iphip_b1] = - reactp1 - reactp2

end


"""
$(TYPEDSIGNATURES)

Electrochemical reaction between interface and bulk ionic species.
This function enters in the internal boundary reaction
in case of an ion charge interface model.

"""
function electrochemicalReaction(data, u, iphia, ipsi, iphiaJunction, ipsiJunction, β, κ, DOS, E) # (1.4.9)

    params             = data.params

    etaExp             = params.chargeNumbers[iphia] / params.UT * ( (u[iphia] - u[iphiaJunction]) + E / q )
    expTerm            =  exp( β * etaExp ) - exp( (β - 1) * etaExp)

    etaInterfaceAnion  = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaJunction] - u[ipsiJunction]) + E / q )
    etaAnion           = params.chargeNumbers[iphia] / params.UT * ( (u[iphia] - u[ipsi]) + E / q )

    densFactor         = ( (DOS^(2/3) * data.F[iphia](etaInterfaceAnion) )^(1/2) * (DOS * data.F[iphia](etaAnion) )^(- 1/2) )^κ

    return densFactor * expTerm

end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bstorage! function. This is the function which enters VoronoiFVM and hands over
for each boundary the time-dependent part of the chosen boundary model.

"""
bstorage!(f, u, bnode, data) = bstorage!(f, u, bnode, data, data.modelType)

bstorage!(f, u, bnode, data, ::Type{Stationary})  = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{Transient}) = bstorage!(f, u, bnode, data, data.boundaryType[bnode.region])


bstorage!(f, u, bnode, data, ::Type{InterfaceModelNone}) = emptyFunction()


function bstorage!(f, u, bnode, data, ::Type{InterfaceModelDiscontqF})
    params = data.params
    ipsi   = data.index_psi

    for icc ∈ data.chargeCarrierList

        get_DOS!(icc, bnode, data)

        Ni     = data.tempDOS1[icc]
        eta    = etaFunction(u, bnode, data, icc) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc)
        f[icc] = q * params.chargeNumbers[icc] * Ni * Boltzmann(eta)

    end

    f[ipsi] = 0.0

end


bstorage!(f, u, bnode, data, ::Type{InterfaceModelDiscontqFNoReaction}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{InterfaceModelSurfaceReco}) = emptyFunction()

# No bstorage! is used, if an ohmic and schottky contact model is chosen.
bstorage!(f, u, bnode, data, ::OuterBoundaryModelType) =  emptyFunction()

bstorage!(f, u, bnode, data, ::Type{InterfaceModelSurfaceRecoAndTangentialFlux}) = bstorage!(f, u, bnode, data, InterfaceModelTangentialFlux)

function bstorage!(f, u, bnode, data, ::Type{InterfaceModelTangentialFlux})

    params      = data.params
    paramsnodal = data.paramsnodal

    #indices (∈ IN) of electron and hole quasi Fermi potentials specified by user (they pass it through recombination)
    iphin       = data.bulkRecombination.iphin # integer index of φ_n
    iphip       = data.bulkRecombination.iphip # integer index of φ_p

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin] # = Quantity or integer
    iphip       = data.chargeCarrierList[iphip] # = Quantity or integer
    looplist    = (iphin, iphip)

    ipsi        = data.index_psi

    for icc = 1:length(looplist)

        get_DOS!(icc, bnode, data)
        Ni     = data.tempDOS[icc]
        eta    = etaFunction(u, bnode, data, icc) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc)
        f[icc] = q * params.chargeNumbers[icc] * Ni * data.F[icc](eta)

    end

    f[ipsi] =  0.0

end
##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bflux! function. This is the function which enters VoronoiFVM and hands over
for each boundary the flux within the boundary.

"""
bflux!(f, u, bedge, data) = bflux!(f, u, bedge, data, data.calculationType)


# In case of equilibrium, the bflux shall not enter.
bflux!(f, u, bedge, data, ::Type{InEquilibrium})             = emptyFunction()


# Out of equilibrium, we need to additionally check for boundary type.
bflux!(f, u, bedge, data, ::Type{OutOfEquilibrium})          = bflux!(f, u, bedge, data, data.boundaryType[bedge.region])

bflux!(f, u, bedge, data, ::Type{InterfaceModelNone})        = emptyFunction()

bflux!(f, u, bedge, data, ::Type{OhmicContact})              = emptyFunction()
bflux!(f, u, bedge, data, ::Type{SchottkyContact})           = emptyFunction()
bflux!(f, u, bedge, data, ::Type{InterfaceModelSurfaceReco}) = emptyFunction()



# Cases, where we have a tangential flux.
bflux!(f, u, bedge, data, ::Type{InterfaceModelSurfaceRecoAndTangentialFlux}) = bflux!(f, u, bedge, data, data.fluxApproximation)

bflux!(f, u, bedge, data, ::Type{InterfaceModelTangentialFlux}) = bflux!(f, u, bedge, data, data.fluxApproximation)


# excess chemical potential flux discretization scheme for inner boundaries.
function bflux!(f, u, bedge, data, ::Type{ExcessChemicalPotential})

    params      =   data.params
    paramsnodal =   data.paramsnodal


    nodel       =   bedge.node[2]
    nodek       =   bedge.node[1]
    ireg        =   bedge.region

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin
    iphip       = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.index_psi                  # final index for electrostatic potential
    looplist    = (iphin, iphip)

    # ############################################################
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc = 1:length(looplist)

        j0           = params.chargeNumbers[icc] * q * params.bMobility[icc, ireg] * params.UT * params.bDensityOfStates[icc, ireg]

        # need to add this to the other etaFunctions
        Ek           = params.bBandEdgeEnergy[icc, bedge.region] + paramsnodal.bandEdgeEnergy[icc, nodek]
        El           = params.bBandEdgeEnergy[icc, bedge.region] + paramsnodal.bandEdgeEnergy[icc, nodel]
        etak         = etaFunction(u[ipsi, 1], u[icc, 1], data.params.UT, Ek, params.chargeNumbers[icc])
        etal         = etaFunction(u[ipsi, 2], u[icc, 2], data.params.UT, El, params.chargeNumbers[icc])

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q            = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm       = fbernoulli_pm(Q)

        f[icc]       = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )

    end

end

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

    # RHS of Poisson
    RHSPoisson!(f, u, node, data)

    # zero reaction term for all icc (stability purpose)
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        f[icc] = u[icc]
    end

end


function addRecombination!(f, u, node, data, ::SRHWithoutTrapsType)

    ireg  = node.region

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin
    iphip       = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]

    get_DOS!(iphin, node, data); get_DOS!(iphip, node, data)

    Nc    = data.tempDOS1[iphin];      Nv = data.tempDOS1[iphip]

    etan  = etaFunction(u, node, data, iphin)
    etap  = etaFunction(u, node, data, iphip)
    n     = Nc * data.F[iphin](etan)
    p     = Nv * data.F[iphip](etap)

    taun  = data.params.recombinationSRHLifetime[iphin, ireg]
    n0    = data.params.recombinationSRHTrapDensity[iphin, ireg]
    taup  = data.params.recombinationSRHLifetime[iphip, ireg]
    p0    = data.params.recombinationSRHTrapDensity[iphip, ireg]

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (kB * data.params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    # calculate recombination kernel. If user adjusted Auger, radiative or SRH recombination,
    # they are set to 0. Hence, adding them here, has no influence since we simply add by 0.0.
    kernelRad   = data.params.recombinationRadiative[ireg]
    kernelAuger = (data.params.recombinationAuger[iphin, ireg] * n + data.params.recombinationAuger[iphip, ireg] * p)
    kernelSRH   = data.params.prefactor_SRH / ( taup * (n + n0) + taun * (p + p0) )
    kernel      = kernelRad + kernelAuger + kernelSRH
    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################
    f[iphin] = q * data.params.chargeNumbers[iphin] * kernel * excessDensTerm
    f[iphip] = q * data.params.chargeNumbers[iphip] * kernel * excessDensTerm

end


function addRecombination!(f, u, node, data, ::SRHWithTrapsType)

    params      = data.params
    ireg        = node.region

    # indices (∈ IN) used by user
    iphin       = data.bulkRecombination.iphin
    iphip       = data.bulkRecombination.iphip
    itrap       = data.enableTraps.traps

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    itrap       = data.chargeCarrierList[itrap]

    get_DOS!(iphin, node, data); get_DOS!(iphip, node, data)

    Nc = data.tempDOS1[iphin];      Nv = data.tempDOS1[iphip]
    Nt = params.densityOfStates[itrap, ireg]

    n     = Nc * data.F[iphin](etaFunction(u, node, data, iphin))
    p     = Nv * data.F[iphip](etaFunction(u, node, data, iphip))
    t     = Nt * data.F[itrap](etaFunction(u, node, data, itrap))

    taun  = params.recombinationSRHLifetime[iphin, ireg]
    n0    = params.recombinationSRHTrapDensity[iphin, ireg]
    taup  = params.recombinationSRHLifetime[iphip, ireg]
    p0    = params.recombinationSRHTrapDensity[iphip, ireg]

    # Rn, Rp agree up to sign with *On the Shockley-Read-Hall Model: Generation-Recombination
    # in Semiconductors* in SIAM Journal on Applied Mathematics, Vol. 67, No. 4 (2007), pp. 1183-1201.
    # The sign is chosen according to *Supporting Information: Consistent Device Simulation Model
    # Describing Perovskite Solar Cells in Steady-State, Transient and Frequency Domain* in ACS (2018)
    if params.chargeNumbers[itrap] == -1
        Rn =  1 / taun * (n * (1-t/Nt) - n0 * t/Nt)
        Rp =  1 / taup * (p * t/Nt     - p0 * (1-t/Nt))
    elseif params.chargeNumbers[itrap] == 1
        Rn =  1 / taun * (n * t/Nt     - n0 * (1-t/Nt))
        Rp =  1 / taup * (p * (1-t/Nt) - p0 * t/Nt)
    end

    f[iphin] = q * params.chargeNumbers[iphin] * Rn
    f[iphip] = q * params.chargeNumbers[iphip] * Rp
    f[itrap] = q * params.chargeNumbers[itrap] * (Rp-Rn)

end

function addGeneration!(f, u, node, data)

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin       = data.bulkRecombination.iphin
    iphip       = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    looplist    = (iphin, iphip)

    generationTerm = generation(data, node.region, node.coord[node.index], data.generationModel)

    for icc = 1:length(looplist)
        f[icc] = f[icc] - q * data.params.chargeNumbers[icc] * generationTerm
    end


end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of Poisson equation, i.e. which builds
the space charge density.
"""
function RHSPoisson!(f, u, node, data)

    ipsi = data.index_psi
    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList[1:2] # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}

        get_DOS!(icc, node, data)

        Ni      = data.tempDOS1[icc]
        eta     = etaFunction(u, node, data, icc)

        f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * ( data.params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * Ni * data.F[icc](eta)   # add charge carrier

    end

    f[ipsi]     = f[ipsi] - data.paramsnodal.doping[node.index]

    f[ipsi]     = - q * data.λ1 * f[ipsi]

end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of electric charge carriers.
"""
function RHSContinuityEquations!(f, u, node, data)

    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        f[icc] = u[icc]
    end

    # dependent on user information concerncing recombination
    addRecombination!(f, u, node, data, data.bulkRecombination.bulk_recomb_SRH)
    # dependent on user information concerncing generation
    addGeneration!(f, u, node, data)

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

    # set RHS to zero for all icc (stability purpose)
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        f[icc]  = u[icc]
    end

    RHSPoisson!(f, u, node, data)             # RHS of Poisson
    RHSContinuityEquations!(f, u, node, data) # RHS of Charge Carriers with special treatment of recombination

end


"""
$(SIGNATURES)
Compute trap densities for a given trap energy.
[Currently, only done for the Boltzmann statistics and for region dependent parameters.]
"""
function trap_density!(icc, ireg, data, Et)
    params      = data.params

    params.densityOfStates[icc, ireg] * exp( params.chargeNumbers[icc] * (params.bandEdgeEnergy[icc, ireg] - Et) / (kB * params.temperature))
end

# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations with a uniform generation rate.
function generation(data, ireg, node, ::Type{GenerationUniform})

    return data.λ2 * data.params.generationUniform[ireg]
end


# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations obeying the Beer-Lambert law.
#only works in 1D till now; adjust node, when multidimensions
function generation(data, ireg, node, ::Type{GenerationBeerLambert})

    params = data.params

    return data.λ2 * params.generationIncidentPhotonFlux[ireg] * params.generationAbsorption[ireg] * exp( - params.generationAbsorption[ireg] * node )

end

generation(data, ireg, node, ::Type{GenerationNone}) = 0.0

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master storage! function. This is the function which enters VoronoiFVM and hands over
a storage term, if we consider transient problem.
"""
storage!(f, u, node, data) = storage!(f, u, node, data, data.modelType)

storage!(f, u, node, data, ::Type{Stationary})  = emptyFunction()


"""
$(TYPEDSIGNATURES)
The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have
``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha``
and for the electrostatic potential
``f[ψ] = 0``.
"""
function storage!(f, u, node, data, ::Type{Transient})

    params = data.params
    ipsi   = data.index_psi

    for icc ∈ data.chargeCarrierList[1:2]

        get_DOS!(icc, node, data)

        Ni     = data.tempDOS1[icc]
        eta    = etaFunction(u, node, data, icc) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc)
        f[icc] = q * params.chargeNumbers[icc] * Ni * data.F[icc](eta)

    end

    f[ipsi] = 0.0

end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master flux functions which enters VoronoiFVM. Flux discretization scheme is chosen in two steps. First, we need
to see, if we are in or out of equilibrium. If, InEquilibrium, then
no flux is passed. If outOfEquilibrium, we choose the flux approximation
which the user chose.
"""
flux!(f, u, edge, data) = flux!(f, u, edge, data, data.calculationType)

function flux!(f, u, edge, data, ::Type{InEquilibrium})

    params      = data.params
    paramsnodal = data.paramsnodal

    ipsi        = data.index_psi # integer index or Quantity
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - (params.dielectricConstant[ireg] + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * ε0 * dpsi

end

flux!(f, u, edge, data, ::Type{OutOfEquilibrium}) = flux!(f, u, edge, data, data.fluxApproximation)


# The classical Scharfetter-Gummel flux scheme. This also works for space-dependent
# band-edge energy, but not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{ScharfetterGummel})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList[1:2]

        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]

        j0           =  params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm       = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT)
        f[icc]       = - j0 * ( bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak) )

    end

end

# The classical Scharfetter-Gummel flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the
# discretization scheme is modified.
function flux!(f, u, edge, data, ::Type{ScharfetterGummelGraded})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 =  params.chargeNumbers[icc] * q * params.UT

        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]

        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobility     = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2

        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            bp, bm   = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT )
        else
            bp, bm   = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc, nodek])) )
        end

        f[icc]       = - j0  * mobility * ( bm  * Nil * data.F[icc](etal) - bp *  Nik * data.F[icc](etak) )

    end

end

# The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{ExcessChemicalPotential})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi


    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0           = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]
        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q            = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff /q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm       = fbernoulli_pm(Q)

        f[icc]       = - j0 * ( bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak) )

    end

end

# The excess chemical potential flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the discretization scheme is modified.
function flux!(f, u, edge, data, ::Type{ExcessChemicalPotentialGraded})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps


    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0           = params.chargeNumbers[icc] * q * params.UT

        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]

        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobilityl    = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodel])
        mobilityk    = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodek])

        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            Q        = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )

        else
            Q        = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc,nodek])) )

        end

        bp, bm       = fbernoulli_pm(Q)
        f[icc]       = - j0  * ( bm  * mobilityl * Nil * data.F[icc](etal) - bp * mobilityk * Nik * data.F[icc](etak) )

    end

end

# The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme
# is used for the regularization of the removable singularity. This also works for
# space-dependent band-edge energy, but not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{DiffusionEnhanced})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    tolReg      =   1.0e-13

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0           = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT
        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]

        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        if abs( (etal - etak)/(etak + etal) ) > tolReg
            g  = (etal - etak ) / ( log(data.F[icc](etal)) - log(data.F[icc](etak)) )
        else
            # regularization idea coming from Pietra-Jüngel scheme
            gk = exp(etak) / data.F[icc](etak)
            gl = exp(etal) / data.F[icc](etal)
            g  = 0.5 * ( gk + gl )
        end

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm       = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / (params.UT * g))
        f[icc]       = - j0 * g * (  bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak))
    end

end

# The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation
# which arise when considering the generalized Scharfetter-Gummel scheme in case of Blakemore
# statistics. Hence, it should be exclusively worked with, when considering the Blakemore
# distribution. This also works for space-dependent band-edge energy, but not for
# space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{GeneralizedSG})

    params        =   data.params
    paramsnodal   =   data.paramsnodal

    max_iteration =   200          # for Newton solver
    it            =   0            # number of iterations (newton)
    damp          =   0.1          # damping factor

    ipsi          =   data.index_psi
    nodel         =   edge.node[2]
    nodek         =   edge.node[1]
    ireg          =   edge.region

    dpsi          =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]       =  - params.dielectricConstant[ireg] * ε0 * dpsi


    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0           = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        get_DOS!(icc, edge, data)

        Nik          = data.tempDOS1[icc]
        Nil          = data.tempDOS2[icc]

        etak, etal   = etaFunction(u, edge, data, icc) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc)

        bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        # use Sedan flux as starting guess
        Q            = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm       = fbernoulli_pm(Q)
        jInitial     = ( bm * Nik * data.F[icc](etal)  - bp * Nil * data.F[icc](etak))

        implicitEq(j::Real) = (fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDiff/q)) /params.UT + params.γ * j )[2] * exp(etal) - fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDiff/q) /params.UT) - params.γ * j)[1] * exp(etak)) - j

        delta        = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
        oldup        = 1.0
        while (it < max_iteration)
            Fval     = implicitEq(jInitial)
            dFval    = ForwardDiff.derivative(implicitEq, jInitial)

            if isnan(value(dFval)) || value(abs(dFval)) < delta
                @show value(jInitial), value(Fval), value(dFval)
                error("singular derivative in exact SG scheme")
            end

            update   = Fval / dFval
            jInitial = jInitial - damp * update

            if abs(update) < delta
                break
            end
            #@show abs(value(update)/oldup)
            oldup = value(update)

            it       = it + 1
            damp     = min(damp * 1.2, 1.0)
        end
        #@show  it
        f[icc]       = - j0 * jInitial
    end

end

##########################################################
##########################################################