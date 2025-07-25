# Parameters from Driftfusion:
# (https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv)

# representing PCBM | MAPI | Pedot:PSS

@kwdef struct Params_PSC_PCBM_MAPI_Pedot

    #####################################################################
    ############################ parameters ############################

    ########## charge carriers ##########

    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    iphia = 3
    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ########## device geometry ##########

    # region numbers
    regionDonor = 1
    regionIntrinsic = 2
    regionAcceptor = 3
    regions = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionDonor = 1
    bregionAcceptor = 2
    bregionJ1 = 3
    bregionJ2 = 4
    numberOfBoundaryRegions = 4

    ## length domains
    h_ndoping = 60.0 * nm
    h_intrinsic = 300.0 * nm
    h_pdoping = 50.0 * nm
    h_total = h_ndoping + h_intrinsic + h_pdoping
    heightLayers = [
        h_ndoping,
        h_ndoping + h_intrinsic,
        h_ndoping + h_intrinsic + h_pdoping,
    ]
    ########## physical values ##########

    ## charge numbers
    zn = -1
    zp = 1
    za = 1

    ## temperature
    T = 300.0 * K

    ## band edge energies
    En = [-3.8, -3.8, -3.0] .* eV
    Ep = [-6.2, -5.4, -5.1] .* eV
    Ea = [0.0, -4.66, 0.0] .* eV

    ## effective densities of density of states
    Nn = [1.0e25, 1.0e25, 1.0e26] ./ (m^3)
    Np = [1.0e25, 1.0e25, 1.0e26] ./ (m^3)
    Na = [0.0, 1.0e26, 0.0] ./ (m^3)


    ## mobilities
    μn = [1.0e-7, 2.0e-3, 1.0e-5] .* (m^2) / (V * s)
    μp = [1.0e-7, 2.0e-3, 1.0e-5] .* (m^2) / (V * s)
    μa = [0.0, 1.0e-14, 0.0] .* (m^2) / (V * s)


    ## relative dielectric permittivity
    ε = [3.0, 23.0, 4.0] .* 1.0

    ## radiative recombination
    r0 = [6.8e-17, 3.6e-18, 6.3e-17] .* cm^3 / s

    ## life times and trap densities
    τn = [1.0e-6, 1.0e-7, 1.0e-6] .* s
    τp = [1.0e-6, 1.0e-7, 1.0e-6] .* s

    ## SRH trap energies
    EI = [-5.0, -4.6, -4.05] .* eV

    ## generation
    incidentPhotonFlux = [0.0, 8.0e20, 0.0] ./ (m^2 * s)
    absorption = [0.0, 1.3e7, 0.0] ./ m
    generationPeak = h_ndoping

    generation_uniform = [0.0, 2.64e27, 0.0] ./ (m^3 * s)

    ## doping
    Cn = 2.09e24 / (m^3)
    Cp = 2.09e24 / (m^3)
    Ca = 1.0e24 / (m^3)

    UT = kB * T / q

end


"""
$(SIGNATURES)

Create a `ChargeTransport.Params` object directly from `Params_PSC_PCBM_MAPI_Pedot`
"""
function Params(p::Params_PSC_PCBM_MAPI_Pedot)

    params = Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        p.numberOfCarriers
    )

    params.temperature = p.T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za

    params.dielectricConstant = p.ε * ε0

    ## effective DOS, band edge energy and mobilities
    params.densityOfStates[p.iphin, :] = p.Nn
    params.densityOfStates[p.iphip, :] = p.Np
    params.densityOfStates[p.iphia, :] = p.Na

    params.bandEdgeEnergy[p.iphin, :] = p.En
    params.bandEdgeEnergy[p.iphip, :] = p.Ep
    params.bandEdgeEnergy[p.iphia, :] = p.Ea

    params.mobility[p.iphin, :] = p.μn
    params.mobility[p.iphip, :] = p.μp
    params.mobility[p.iphia, :] = p.μa


    for ireg in 1:p.numberOfRegions ## interior region data
        ## recombination parameters
        params.recombinationRadiative[ireg] = p.r0[ireg]
        params.recombinationSRHLifetime[p.iphin, ireg] = p.τn[ireg]
        params.recombinationSRHLifetime[p.iphip, ireg] = p.τp[ireg]
        params.recombinationSRHTrapDensity[p.iphin, ireg] = trap_density!(p.iphin, ireg, params, p.EI[ireg])
        params.recombinationSRHTrapDensity[p.iphip, ireg] = trap_density!(p.iphip, ireg, params, p.EI[ireg])
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[p.iphin, p.bregionJ1] = p.Nn[p.regionIntrinsic]
    params.bDensityOfStates[p.iphip, p.bregionJ1] = p.Np[p.regionIntrinsic]

    params.bDensityOfStates[p.iphin, p.bregionJ2] = p.Nn[p.regionIntrinsic]
    params.bDensityOfStates[p.iphip, p.bregionJ2] = p.Np[p.regionIntrinsic]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ1] = p.En[p.regionIntrinsic]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ1] = p.Ep[p.regionIntrinsic]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ2] = p.En[p.regionIntrinsic]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ2] = p.Ep[p.regionIntrinsic]

    ## for surface recombination
    params.recombinationSRHvelocity[p.iphin, p.bregionJ1] = 1.0e1 * cm / s
    params.recombinationSRHvelocity[p.iphip, p.bregionJ1] = 1.0e5 * cm / s

    params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphin, p.regionIntrinsic]
    params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphip, p.regionIntrinsic]

    params.recombinationSRHvelocity[p.iphin, p.bregionJ2] = 1.0e7 * cm / s
    params.recombinationSRHvelocity[p.iphip, p.bregionJ2] = 1.0e1 * cm / s

    params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphin, p.regionIntrinsic]
    params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphip, p.regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[p.iphin, p.regionDonor] = p.Cn
    params.doping[p.iphip, p.regionAcceptor] = p.Cp
    params.doping[p.iphia, p.regionIntrinsic] = p.Ca

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak = p.generationPeak

    ## generation parameters
    params.generationIncidentPhotonFlux = p.incidentPhotonFlux
    params.generationAbsorption = p.absorption
    params.generationUniform = p.generation_uniform

    return params
end
