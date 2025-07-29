# Default parameters of Ionmonger (https://github.com/PerovskiteSCModelling/IonMonger)
# representing TiO2 | MAPI | spiro-OMeTAD

@kwdef struct Params_PSC_TiO2_MAPI_spiro

    #####################################################################
    ############################ parameters ############################

    # used unit factors
    nm = ufac"nm"
    K = ufac"K"
    m = ufac"m"
    V = ufac"V"
    s = ufac"s"

    eV = q * V

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
    numberOfBoundaryRegions = 4

    # boundary region numbers
    bregionDonor = 1
    bregionAcceptor = 2
    bregionJ1 = 3
    bregionJ2 = 4

    ## length domains
    h_ndoping = 100.0 * nm
    h_intrinsic = 400.0 * nm
    h_pdoping = 200.0 * nm
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
    T = 298.0 * K

    ## band edge energies
    En = [-4.0, -3.7, -3.4] .* eV
    Ep = [-5.8, -5.4, -5.1] .* eV
    Ea = [0.0, -4.45, 0.0] .* eV
    Ea_i = Ea[regionIntrinsic]

    ## effective densities of density of states
    Nn = [5.0e25, 8.1e24, 5.0e25] ./ (m^3)
    Np = [5.0e25, 5.8e24, 5.0e25] ./ (m^3)
    Na = [0.0, 1.0e27, 0.0] ./ (m^3)
    Na_i = Na[regionIntrinsic]

    ## mobilities
    μn = [3.89e-4, 6.62e-3, 3.89e-5] .* (m^2) / (V * s)
    μp = [3.89e-4, 6.62e-3, 3.89e-5] .* (m^2) / (V * s)
    μa = [0.0, 3.93e-16, 0.0] .* (m^2) / (V * s)

    ## relative dielectric permittivity
    ε = [10.0, 24.1, 3.0] .* 1.0

    ## radiative recombination
    r0 = [6.8e-17, 3.6e-18, 6.3e-17] .* m^3 / s

    ## life times and trap densities
    τn = [1.0e100, 3.0e-9, 1.0e100] .* s
    τp = [1.0e100, 3.0e-7, 1.0e100] .* s

    ## SRH trap energies
    EI = [-5.0, -4.55, -4.1] .* eV

    ## generation
    incidentPhotonFlux = [0.0, 9.4e20, 0.0] ./ (m^2 * s)
    absorption = [0.0, 1.3e7, 0.0] ./ m
    generationPeak = h_ndoping

    generation_uniform = [0.0, 2.5e27, 0.0] ./ (m^3 * s)

    ## doping
    Cn = 1.0e24 / (m^3)
    Cp = 1.0e24 / (m^3)
    Ca = 1.6e25 / (m^3)
end

"""
$(SIGNATURES)

Create a `ChargeTransport.Params` object directly from `Params_PSC_TiO2_MAPI_spiro`
"""
function Params(p::Params_PSC_TiO2_MAPI_spiro)

    params = Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        p.numberOfCarriers
    )

    params.temperature = p.T
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za

    for ireg in 1:p.numberOfRegions # interior region data

        params.dielectricConstant[ireg] = p.ε[ireg] * ε_0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[p.iphin, ireg] = p.Nn[ireg]
        params.densityOfStates[p.iphip, ireg] = p.Np[ireg]
        params.densityOfStates[p.iphia, ireg] = p.Na[ireg]

        params.bandEdgeEnergy[p.iphin, ireg] = p.En[ireg]
        params.bandEdgeEnergy[p.iphip, ireg] = p.Ep[ireg]
        params.bandEdgeEnergy[p.iphia, ireg] = p.Ea[ireg]

        params.mobility[p.iphin, ireg] = p.μn[ireg]
        params.mobility[p.iphip, ireg] = p.μp[ireg]
        params.mobility[p.iphia, ireg] = p.μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = p.r0[ireg]
        params.recombinationSRHLifetime[p.iphin, ireg] = p.τn[ireg]
        params.recombinationSRHLifetime[p.iphip, ireg] = p.τp[ireg]
        params.recombinationSRHTrapDensity[p.iphin, ireg] = trap_density!(p.iphin, ireg, params, p.EI[ireg])
        params.recombinationSRHTrapDensity[p.iphip, ireg] = trap_density!(p.iphip, ireg, params, p.EI[ireg])

        ## generation parameters
        params.generationIncidentPhotonFlux[ireg] = p.incidentPhotonFlux[ireg]
        params.generationAbsorption[ireg] = p.absorption[ireg]
        params.generationUniform[ireg] = p.generation_uniform[ireg]
    end

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak = p.generationPeak

    ## interior doping
    params.doping[p.iphin, p.regionDonor] = p.Cn
    params.doping[p.iphia, p.regionIntrinsic] = p.Ca
    params.doping[p.iphip, p.regionAcceptor] = p.Cp

    return params
end
