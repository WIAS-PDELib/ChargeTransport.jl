
"""
$(SIGNATURES)

Method which determines with input parameters which inner interface model 
was chosen by user.

"""
function inner_interface_model(data::ChargeTransportData)

    countSurfaceReco = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:data.params.numberOfBoundaryRegions
        
        if     data.boundary_type[ireg] ==  interface_model_surface_recombination

            countSurfaceReco = countSurfaceReco + 1

        elseif data.boundary_type[ireg] == interface_model_ion_charge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countSurfaceReco > 0 # build surface_recombination based system

        return interface_model_surface_recombination

    elseif countInterfaceCharge > 0 # build ion interface charge system 

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in 
        #     Example107 works perfectly fine since the choice of species number is already set.
        return interface_model_none

    elseif countSurfaceReco + countInterfaceCharge == 0 # build system without further interface conditions

        return interface_model_none

    end

end

"""
$(SIGNATURES)

Method which determines with input parameters which inner interface model 
was chosen by user.

"""
function inner_interface_model(ctsys::ChargeTransportSystem)


    countSurfaceReco = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:ctsys.data.params.numberOfBoundaryRegions
        
        if     ctsys.data.boundary_type[ireg] ==  interface_model_surface_recombination

            countSurfaceReco = countSurfaceReco + 1

        elseif ctsys.data.boundary_type[ireg] == interface_model_ion_charge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countSurfaceReco > 0 # build surface_recombination based system

        return interface_model_surface_recombination

    elseif countInterfaceCharge > 0 # build ion interface charge system 

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in 
        #     Example107 works perfectly fine since the choice of species number is already set.
        return interface_model_none

    elseif countSurfaceReco + countInterfaceCharge == 0 # build system without further interface conditions

        return interface_model_none

    end

end

"""
$(SIGNATURES)

New system constructor which builds all necessary information needed based
on the input parameters concerning additional interface models.
The idea is that this new constructor will replace the current one.

"""
function ChargeTransportSystem2(grid, data ;unknown_storage)

    ctsys                 = ChargeTransportSystem()

    interface_model       = inner_interface_model(data)
    ctsys                 = build_system(grid, data, unknown_storage, interface_model) 
    
    return ctsys

end

##########################################################
##########################################################
"""
$(SIGNATURES)

The core of the new system constructor. Here, the system for no additional interface model is build.

"""
function build_system(grid, data, unknown_storage, ::Type{interface_model_none})

    num_species  = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1

    ctsys        = ChargeTransportSystem()

    ctsys.data   = data
    
    physics      = VoronoiFVM.Physics(data        = data,
                                      flux        = flux!,
                                      reaction    = reaction!,
                                      breaction   = breaction!,
                                      storage     = storage!,
                                      bstorage    = bstorage!
                                      )

    ctsys.fvmsys = VoronoiFVM.System(grid, physics, unknown_storage = unknown_storage)

    # for detection of number of species. In following versions we can simply delete num_species from physics initialization. 
    VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species) 

    return ctsys

end

"""
$(SIGNATURES)

The core of the new system constructor. Here, the system for additional
surface recombination at inner interfaces is build.

"""
function build_system(grid, data, unknown_storage, ::Type{interface_model_surface_recombination})

    #num_species  = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1

    ctsys        = ChargeTransportSystem()

    fvmsys = VoronoiFVM.System(grid, unknown_storage=unknown_storage)

    data.iphin = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions) # iphin
    data.iphip = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions) # iphip

    for ii = 3:data.params.numberOfCarriers
        ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions) # if ion vacancies are present
    end

    data.ipsi  = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions)    # last quantitiy is psi

    physics      = VoronoiFVM.Physics(data        = data,
                                      flux        = fluxDiscont!,
                                      reaction    = reactionDiscont!,
                                      breaction   = breactionDiscont!
                                      ### DA: insert additionally storage and bstorage
                                      )

    # add the defined physics to system
    physics!(fvmsys, physics)

    # numberOfSpecies within system changes since we consider discontinuous quasi Fermi potentials
    data.params.numberOfSpeciesSystem = VoronoiFVM.num_species(fvmsys)

    ctsys.fvmsys = fvmsys
    ctsys.data   = data

    return ctsys

end

##########################################################
##########################################################

breactionDiscont!(f, u, bnode, data) = breactionDiscont!(f, u, bnode, data, data.boundary_type[bnode.region])


function breactionDiscont!(f, u, bnode, data, ::Type{interface_model_surface_recombination})

    if data.calculation_type == inEquilibrium

        return emptyFunction()

    else
    
        params      = data.params
    
        ########
        etan1 = params.chargeNumbers[1] / params.UT * ( (u[data.iphin, 1] - u[data.ipsi]) + params.bandEdgeEnergy[1, bnode.cellregions[1]] / q )
        etan2 = params.chargeNumbers[1] / params.UT * ( (u[data.iphin, 2] - u[data.ipsi]) + params.bandEdgeEnergy[1, bnode.cellregions[2]] / q )
        
        n1    = params.bDensityOfStates[1, bnode.cellregions[1] ] * data.F[1](etan1)
        n2    = params.bDensityOfStates[1, bnode.cellregions[2] ] * data.F[1](etan2)
        ########
        etap1 = params.chargeNumbers[2] / params.UT * ( (u[data.iphip, 1] - u[data.ipsi]) + params.bandEdgeEnergy[2, bnode.cellregions[1]] / q )
        etap2 = params.chargeNumbers[2] / params.UT * ( (u[data.iphip, 2] - u[data.ipsi]) + params.bandEdgeEnergy[2, bnode.cellregions[2]] / q )
    
        p1    = params.densityOfStates[2, bnode.cellregions[1] ] * data.F[2](etap1)
        p2    = params.densityOfStates[2, bnode.cellregions[2] ] * data.F[2](etap2)
    
       
        d       = 5.0e14   # choose: 5.0e14 to see discontinuity
    
        react_n         = (n1 - n2)/ d
        react_p         = (p1 - p2)/ d
        ########
        f[data.iphin, 1] =   react_n
        f[data.iphin, 2] = - react_n
    
        f[data.iphip, 1] =   react_p
        f[data.iphip, 2] = - react_p

        return f
    end

end
##########################################################
##########################################################



# without recombination!!!!
function reactionDiscont!(f, u, node, data)

    params      = data.params
    paramsnodal = data.paramsnodal
    
    # discontinuous quantities
    etan     = params.chargeNumbers[1] / params.UT * ( (u[data.iphin] - u[data.ipsi]) + params.bandEdgeEnergy[1, node.region] / q )
    etap     = params.chargeNumbers[2] / params.UT * ( (u[data.iphip] - u[data.ipsi]) + params.bandEdgeEnergy[2, node.region] / q )

    f[data.ipsi] = f[data.ipsi] - params.chargeNumbers[1] * (params.doping[1, node.region])  - params.chargeNumbers[2] * (params.doping[2, node.region])  # subtract doping
    f[data.ipsi] = f[data.ipsi] + params.chargeNumbers[1] * params.densityOfStates[1, node.region] * data.F[1](etan) + params.chargeNumbers[2] * params.densityOfStates[2, node.region] * data.F[2](etap)  # add charge carrier

    # if defined, continuous quantities
    otherCarriers = data.iphip.regionspec[data.params.numberOfRegions] + 1

    for icc = otherCarriers:data.ipsi.ispec-1

        E            = params.bandEdgeEnergy[icc, node.region] + paramsnodal.bandEdgeEnergy[icc, node.index]
        eta          = params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[data.ipsi]) + E / q )

        f[data.ipsi] = f[data.ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[data.ipsi] = f[data.ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    f[data.ipsi] = - data.λ1 * q * f[data.ipsi]

    if data.calculation_type == inEquilibrium # these are for stability purposes
        f[data.iphin] = u[data.iphin] - 0.0
        f[data.iphip] = u[data.iphip] - 0.0
    else
        f[data.iphin] =  0.0
        f[data.iphip] =  0.0
    end
 

    return f

end

##########################################################
##########################################################


fluxDiscont!(f, u, edge, data) = fluxDiscont!(f, u, edge, data, data.calculation_type)

function fluxDiscont!(f, u, edge, data, ::Type{inEquilibrium})

    params      = data.params

    dpsi        =   u[data.ipsi, 2] - u[data.ipsi, 1]
    f[data.ipsi] = - params.dielectricConstant[edge.region] * ε0 * dpsi

    return f
    
end


fluxDiscont!(f, u, edge, data, ::Type{outOfEquilibrium}) = fluxDiscont!(f, u, edge, data, data.flux_approximation)



"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function fluxDiscont!(f, u, edge, data, ::Type{ScharfetterGummel})

    params      = data.params
    
    dpsi        =   u[data.ipsi, 2] - u[data.ipsi, 1]
    f[data.ipsi] = - params.dielectricConstant[edge.region] * ε0 * dpsi

    # discontinuous quantities
    j0n         =  params.chargeNumbers[1] * q *  params.UT 
    j0p         =  params.chargeNumbers[2] * q *  params.UT

    etank    = params.chargeNumbers[1] / params.UT * ( (u[data.iphin, 1] - u[data.ipsi, 1]) + params.bandEdgeEnergy[1, edge.region]/q )
    etanl    = params.chargeNumbers[1] / params.UT * ( (u[data.iphin, 2] - u[data.ipsi, 2]) + params.bandEdgeEnergy[1, edge.region]/q )
    
    etapk    = params.chargeNumbers[2] / params.UT * ( (u[data.iphip, 1] - u[data.ipsi, 1]) + params.bandEdgeEnergy[2, edge.region]/q )
    etapl    = params.chargeNumbers[2] / params.UT * ( (u[data.iphip, 2] - u[data.ipsi, 2]) + params.bandEdgeEnergy[2, edge.region]/q ) 
    ############
    
    bpn, bmn = fbernoulli_pm(params.chargeNumbers[1] * dpsi / params.UT)
    bpp, bmp = fbernoulli_pm(params.chargeNumbers[2] * dpsi / params.UT)
   
    f[data.iphin] = - j0n * (params.mobility[1, edge.region] * bmn * params.densityOfStates[1, edge.region]  * data.F[1](etanl) - params.mobility[1, edge.region] * bpn * params.densityOfStates[1, edge.region]  * data.F[1](etank)) 

    f[data.iphip] = - j0p * (params.mobility[2, edge.region] * bmp * params.densityOfStates[2, edge.region]  * data.F[2](etapl) - params.mobility[2, edge.region] * bpp * params.densityOfStates[2, edge.region]  * data.F[2](etapk))

    ############
    # continuous quantities, if defined
    otherCarriers = data.iphip.regionspec[data.params.numberOfRegions] + 1
    nodel         = edge.node[2]; nodek       = edge.node[1]

    for icc = otherCarriers:data.ipsi.ispec-1
        j0                 =  params.chargeNumbers[icc] * q * params.mobility[icc, node.region] * params.UT * params.densityOfStates[icc, node.region]

        Ek                 = params.bandEdgeEnergy[icc, edge.region] + paramsnodal.bandEdgeEnergy[icc, nodek]
        El                 = params.bandEdgeEnergy[icc, edge.region] + paramsnodal.bandEdgeEnergy[icc, nodel]

        etak               = params.chargeNumbers[icc] / params.UT * ( (u[icc, 1] - u[ipsi, 1]) + Ek / q )
        etal               = params.chargeNumbers[icc] / params.UT * ( (u[icc, 2] - u[ipsi, 2]) + El / q )

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm             = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT)
        f[icc]             = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
    end

    return f

end

##########################################################
##########################################################
function breactionDiscont!(f, u, bnode, data, ::Type{ohmic_contact})

    params      = data.params
    paramsnodal = data.paramsnodal 

    # parameters
    α          = 1.0e-10                      # tiny penalty value
 
    innerRegion = bnode.cellregions[1] # tuple with [subRegionNumber 0]


    etan = params.chargeNumbers[1] / params.UT * ( (u[data.iphin.regionspec[innerRegion]] - u[data.ipsi]) + params.bBandEdgeEnergy[1, bnode.region] / q )

    etap = params.chargeNumbers[2] / params.UT * ( (u[data.iphip.regionspec[innerRegion]] - u[data.ipsi]) + params.bBandEdgeEnergy[2, bnode.region] / q )
 
    f[data.ipsi] = f[data.ipsi] - params.chargeNumbers[1] * ( params.bDoping[1, bnode.region] ) - params.chargeNumbers[2] * ( params.bDoping[2, bnode.region] )# subtract doping
    f[data.ipsi] = f[data.ipsi] + params.chargeNumbers[1] * (params.bDensityOfStates[1, bnode.region] + paramsnodal.densityOfStates[1, bnode.index]) * data.F[1](etan) + params.chargeNumbers[2] * (params.bDensityOfStates[2, bnode.region] + paramsnodal.densityOfStates[2, bnode.index]) * data.F[2](etap) # add charge carrier
 
    f[data.ipsi] = - data.λ1 * 1 / α *  q * f[data.ipsi]

    
    return f

end

