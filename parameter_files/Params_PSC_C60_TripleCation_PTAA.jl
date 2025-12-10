# Parameters from Thiesbrummel et al., nature energy, 2024
# with data publication https://doi.org/10.25446/oxford.24359959
# (https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv)

# representing C60 | triple cation | PTAA

@kwdef struct Params_PSC_C60_TripleCation_PTAA

    #####################################################################
    ############################ parameters ############################

    # physical constants
    constants = ChargeTransport.constants
    q = ChargeTransport.constants.q
    k_B = ChargeTransport.constants.k_B

    # used unit factors
    cm = ufac"cm"
    nm = ufac"nm"
    K = ufac"K"
    m = ufac"m"
    s = ufac"s"
    V = ufac"V"

    eV = constants.q * V

    ########## charge carriers ##########

    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    iphia = 3 # anion vacancy quasi Fermi potential
    ipsi = 4

    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ########## device geometry ##########

    # region numbers
    regionETL = 1
    regionPero = 2
    regionHTL = 3
    regions = [regionETL, regionPero, regionHTL]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionLeft = 1
    bregionRight = 2
    bregionJ1 = 3
    bregionJ2 = 4
    numberOfBoundaryRegions = 4

    heightDev = 750.0 * nm

    ## length domains
    h_ETL = 30.0 * nm       # C60
    h_activePL = 400.0 * nm # perovskite
    h_HTL = 10.0 * nm       # PTAA

    heightLayersPL = [h_ETL, h_ETL + h_activePL, h_ETL + h_activePL + h_HTL]
    h_totalPL = heightLayersPL[end]

    ########## physical values ##########

    ## charge numbers
    zn = -1
    zp = 1
    za = 1

    ## temperature
    T = 300.0 * K

    ## band edge energies
    En = [-3.9, -3.9, -2.5] .* eV
    Ep = [-5.9, -5.53, -5.5] .* eV

    Ea1D = [0.0, -5.267, 0.0] .* eV
    Ea2D = [0.0, -5.287, 0.0] .* eV

    ## effective densities of density of states
    Nn = [1.0e26, 2.2e24, 1.0e26] ./ (m^3)
    Np = [1.0e26, 2.2e24, 1.0e26] ./ (m^3)
    Na = [0.0, 1.0e27, 0.0] ./ (m^3)

    Da = 5.0e-14
    ## mobilities
    μn = [1.0e-6, 1.0e-4, 1.0e-8] .* (m^2) / (V * s)
    μp = [1.0e-6, 1.0e-4, 1.0e-8] .* (m^2) / (V * s)
    μa = [0.0, Da / (k_B * T / q), 0.0] .* (m^2) / (V * s) # 1.9e-12

    ## relative dielectric permittivity
    ε = [5.0, 22.0, 3.5] .* 1.0

    ## radiative recombination
    r0 = [0.0, 3.0e-17, 0.0] .* m^3 / s

    ## SRH life times
    τn = [1.0e100, 2.0e-7, 1.0e100] .* s
    τp = [1.0e100, 2.0e-7, 1.0e100] .* s

    ## SRH trap densities
    ni1 = sqrt(Nn[1] * Np[1] * exp(- (En[1] - Ep[1]) / (k_B * T)))
    ni2 = sqrt(Nn[2] * Np[2] * exp(- (En[2] - Ep[2]) / (k_B * T)))
    ni3 = sqrt(Nn[3] * Np[3] * exp(- (En[3] - Ep[3]) / (k_B * T)))

    nτ = [ni1, ni2, ni3]
    pτ = [ni1, ni2, ni3]

    ## generation
    incidentPhotonFlux = [0.0, 1.42e21, 0.0] ./ (m^2 * s)
    absorption = [0.0, 6.34e6, 0.0] ./ m
    generationPeak = h_ETL + h_activePL
    invertedIllumination = -1

    ## doping
    ## doping of transport layers are just for numerical stability, they are so small, they count as undoped.
    Cn = 1.0e20 / (m^3)
    Cp = 1.0e20 / (m^3)
    Ca = 6.0e22 / (m^3)

    ## metal offsets
    ΔEFLeft = 0.05 * eV  # offset between metal and ETL
    ΔEFRight = 0.05 * eV # offset between metal and HTL

end


"""
$(SIGNATURES)

Create a `ChargeTransport.Params` object directly from `Params_PSC_C60_TripleCation_PTAA`
"""
function Params(p::Params_PSC_C60_TripleCation_PTAA)

    params = Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        p.numberOfCarriers
    )

    params.temperature = p.T
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za

    params.dielectricConstant = p.ε * p.constants.ε_0

    ## effective DOS, band edge energy and mobilities
    params.densityOfStates[p.iphin, :] = p.Nn
    params.densityOfStates[p.iphip, :] = p.Np
    params.densityOfStates[p.iphia, :] = p.Na

    params.bandEdgeEnergy[p.iphin, :] = p.En
    params.bandEdgeEnergy[p.iphip, :] = p.Ep

    params.mobility[p.iphin, :] = p.μn
    params.mobility[p.iphip, :] = p.μp
    params.mobility[p.iphia, :] = p.μa


    for ireg in 1:p.numberOfRegions ## interior region data
        ## recombination parameters
        params.recombinationRadiative[ireg] = p.r0[ireg]
        params.recombinationSRHLifetime[p.iphin, ireg] = p.τn[ireg]
        params.recombinationSRHLifetime[p.iphip, ireg] = p.τp[ireg]
        params.recombinationSRHTrapDensity[p.iphin, ireg] = p.nτ[ireg]
        params.recombinationSRHTrapDensity[p.iphip, ireg] = p.pτ[ireg]
    end

    ##############################################################
    ## inner boundary region data
    params.bDensityOfStates[p.iphin, p.bregionJ1] = p.Nn[p.regionETL]
    params.bDensityOfStates[p.iphip, p.bregionJ1] = p.Np[p.regionPero]

    params.bDensityOfStates[p.iphin, p.bregionJ2] = p.Nn[p.regionPero]
    params.bDensityOfStates[p.iphip, p.bregionJ2] = p.Np[p.regionHTL]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ1] = p.En[p.regionETL]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ1] = p.Ep[p.regionPero]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ2] = p.En[p.regionPero]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ2] = p.Ep[p.regionHTL]

    ############# DA:: AB HIER WEITER SCHAUEN!!!!
    ## for surface recombination

    params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphin, p.regionETL]
    params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphip, p.regionPero]

    params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphin, p.regionPero]
    params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphip, p.regionHTL]

    ## surface recombination velocities will be set in the main folder as
    ## variable parameters
    ##############################################################

    # Schottky barriers for electric potential
    params.SchottkyBarrier[p.bregionLeft] = p.ΔEFLeft
    params.SchottkyBarrier[p.bregionRight] = p.En[3] - p.Ep[3] - p.ΔEFRight

    ## interior doping
    params.doping[p.iphin, p.regionETL] = p.Cn
    params.doping[p.iphip, p.regionHTL] = p.Cp
    params.doping[p.iphia, p.regionPero] = p.Ca

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak = p.generationPeak

    ## generation parameters
    params.generationIncidentPhotonFlux = p.incidentPhotonFlux
    params.generationAbsorption = p.absorption
    params.invertedIllumination = -1

    return params
end
