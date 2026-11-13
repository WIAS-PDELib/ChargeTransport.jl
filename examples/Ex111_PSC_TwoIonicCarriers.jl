#=
# PSC device with two ionic carriers and photogeneration (1D).

Three layer PSC device TiO2 | MAPI | spiro with anion and cation vacancies in the
perovskite layer. Both ionic carriers compensate each other, there is no fixed ionic
background charge. Their band-edge energies are calculated in thermodynamic equilibrium
such that the mean densities equal the prescribed averageDensity.
Afterwards, we perform a forward scan under Beer-Lambert photogeneration.
=#

module Ex111_PSC_TwoIonicCarriers

using ChargeTransport
using ExtendableGrids
using GridVisualize
using LaTeXStrings
using VoronoiFVM

# for convenience
parametersdir = ChargeTransport.parametersdir

# supported Plotters are GLMakie and PythonPlot
# you can set verbose also to true to display some solver information
function main(;
        n = 5,
        Plotter = nothing,
        verbose = false, test = false,
        ########################
        parameter_set = Params_PSC_C60_TripleCation_PTAA, # choose the parameter set
        ########################
        vacancyEnergyCalculation = true,            # assume the vacancy energy level is either given or not
    )

    @local_unitfactors μm cm s ns V K ps Hz W m

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # parameters
    p = parameter_set()
    numberOfCarriers2 = 4 # as we also include the cation vacancies, we have 4 carriers in total


    ######### primary data for I-V scan protocol ##########
    endVoltage = 1.3 * V
    tPrecond = 0.001
    scanrate = 1.0e3 * V / s
    tend = endVoltage / scanrate

    ## Define scan protocol function
    function scanProtocolPrecond(t)
        biasVal = 0.0
        if data.calculationType == InEquilibrium
            biasVal = 0.0
        else
            if 0.0 <= t && t <= tPrecond
                biasVal = endVoltage
            elseif tPrecond <= t  && t <= tPrecond + tend
                biasVal = endVoltage .- scanrate * (t - tPrecond)
            elseif tend + tPrecond < t  && t <= tPrecond + 2 * tend
                biasVal = scanrate * (t - endVoltage / scanrate - tPrecond) # tend = endVoltage / scanrate
            else
                biasVal = 0.0
            end
        end

        # we need this such that during the computation of the equilibrium solution, we are safe that no bias is applied.
        if !data.generationComplete
            biasVal = 0.0
        end
        return biasVal
    end

    # Apply zero voltage on left boundary and scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocolPrecond]

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u = collect(range(0.0, p.h_ETL / 2, step = p.h_ETL / (0.4 * δ)))
    coord_n_g = geomspace(
        p.h_ETL / 2,
        p.h_ETL,
        p.h_ETL / (0.4 * δ),
        p.h_ETL / (0.7 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ETL,
        p.h_ETL + p.h_activePL / k,
        p.h_activePL / (2.8 * δ),
        p.h_activePL / (2.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ETL + p.h_activePL / k,
        p.h_ETL + p.h_activePL,
        p.h_activePL / (2.1 * δ),
        p.h_activePL / (2.8 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ETL + p.h_activePL,
        p.h_ETL + p.h_activePL + p.h_HTL / 2,
        p.h_HTL / (0.2 * δ),
        p.h_HTL / (0.2 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ETL + p.h_activePL + p.h_HTL / 2, p.h_ETL + p.h_activePL + p.h_HTL, step = p.h_HTL / (0.2 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [p.heightLayersPL[1]], p.regionETL, tol = 1.0e-18) # n-doped region   = 1
    cellmask!(grid, [p.heightLayersPL[1]], [p.heightLayersPL[2]], p.regionPero, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [p.heightLayersPL[2]], [p.heightLayersPL[3]], p.regionHTL, tol = 1.0e-18) # p-doped region   = 3

    bfacemask!(grid, [p.heightLayersPL[1]], [p.heightLayersPL[1]], p.bregionJ1, tol = 1.0e-18)
    bfacemask!(grid, [p.heightLayersPL[2]], [p.heightLayersPL[2]], p.bregionJ2, tol = 1.0e-18)

    if Plotter !== nothing
        vis = GridVisualizer(; Plotter, layout = (4, 2), size = (1550, 800))
        gridplot!(vis[1, 1], grid; Plotter, legend = :lt, title = "Grid", xlabel = L"\text{space [m]}", show = true)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################
    data = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction)
    data.modelType = Transient

    ## The default for electrons and holes is Boltzmann. Here, we set it to a more general statistics function
    data.F[p.iphin] = FermiDiracOneHalfTeSCA
    data.F[p.iphip] = FermiDiracOneHalfTeSCA

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )
    data.boundaryType[p.bregionLeft] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionRight] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionJ1] = InterfaceRecombination
    data.boundaryType[p.bregionJ2] = InterfaceRecombination

    ## two compensating ionic carriers: no doping, but a prescribed mean density
    ## It also possible to prescibe a doping for the ionic carriers which will be taken as mean density
    ## Caution: Pay special attention to the space charge density, i.e. RHS of Poisson then
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionPero], averageDensity = [6.0e22 / (m^3)])
    enable_ionic_carrier!(data, ionicCarrier = p.iphic, regions = [p.regionPero], averageDensity = [6.0e22 / (m^3)])

    data.generationModel = GenerationBeerLambert

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params = Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        numberOfCarriers2
    )

    params.temperature = p.T
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za
    params.chargeNumbers[p.iphic] = p.zc

    params.dielectricConstant = p.ε * p.constants.ε_0

    ## effective DOS, band edge energy and mobilities
    params.densityOfStates[p.iphin, :] = p.Nn
    params.densityOfStates[p.iphip, :] = p.Np
    params.densityOfStates[p.iphia, :] = p.Na
    params.densityOfStates[p.iphic, :] = p.Nc

    params.bandEdgeEnergy[p.iphin, :] = p.En
    params.bandEdgeEnergy[p.iphip, :] = p.Ep

    params.mobility[p.iphin, :] = p.μn
    params.mobility[p.iphip, :] = p.μp
    params.mobility[p.iphia, :] = p.μa
    params.mobility[p.iphic, :] = p.μc

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

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak = p.generationPeak

    ## generation parameters
    params.generationIncidentPhotonFlux = p.incidentPhotonFlux
    params.generationAbsorption = p.absorption
    params.invertedIllumination = -1

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################
    control = ChargeTransport.SolverControl()
    control.verbose = verbose
    control.maxiters = 300
    control.max_round = 5
    control.damp_initial = 0.5
    control.damp_growth = 1.21
    control.Δt_max = 5.0e-1

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## calculate equilibrium solution and as initial guess
    ## We recommend calculating the energies for ions for 4 digits or more for better accuracy.
    solution = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = vacancyEnergyCalculation, energyDigits = 4, verboseEnergy = true)
    inival = solution
    solutionEQ = inival

    if Plotter !== nothing
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, p.iphia] = "\$E_a-q\\psi\$"; label_energy[2, p.iphia] = "\$ - q \\varphi_a\$"; label_BEE[p.iphia] = "\$E_a\$"
        label_density[p.iphia] = "\$ n_a \$";      label_solution[p.iphia] = "\$ \\varphi_a\$"
        label_energy[1, p.iphic] = "\$E_c-q\\psi\$"; label_energy[2, p.iphic] = "\$ - q \\varphi_c\$"; label_BEE[p.iphic] = "\$E_c\$"
        label_density[p.iphic] = "\$ n_c \$";      label_solution[p.iphic] = "\$ \\varphi_c\$"

        plot_densities!(vis[1, 2], ctsys, solution, "Initial condition", label_density)
        plot_solution!(vis[2, 1], ctsys, solution, "Initial condition", label_solution)
    end

    if test == false
        println("*** done\n")
    end


    ################################################################################
    if test == false
        println("Loop to increase bias (so that following scan protocol starts near VOC)")
    end
    ################################################################################

    biasLoop = collect(reverse(-range(0.0, endVoltage, length = 21)))

    sol0p0 = ChargeTransport.unknowns(ctsys)
    for istep in 1:length(biasLoop)

        ## turn slowly voltage on
        set_contact!(ctsys, p.bregionRight, Δu = biasLoop[istep])

        ## keep some big time stepping here, so we preserve the vacancy density
        sol0p0 = ChargeTransport.solve(ctsys, inival = inival, control = control, tstep = 1.0e3)
        inival = sol0p0

    end # bias loop

    if !test
        println("*** done\n")
    end

    # Now we have the correct initial condition to start the scan protocol
    ################################################################################
    if test == false
        println("Preconditioning")
    end
    ################################################################################

    control.Δt = 1.0e-4
    control.Δt_min = 1.0e-5
    control.Δt_max = 5.0e-4
    control.Δt_grow = 1.1

    solPrecond = ChargeTransport.solve(ctsys, inival = inival, times = (0.0, tPrecond), control = control)

    if !test
        println("*** done\n")
    end

    ################################################################################
    if !test
        println("Reverse IV Measurement loop")
    end
    ################################################################################

    control.Δt = 4.0e-3 / scanrate
    control.Δt_min = 4.0e-3 / scanrate
    control.Δt_max = 6.0e-2 / scanrate
    control.Δt_grow = 1.15

    solRev = ChargeTransport.solve(ctsys, inival = solPrecond.u[end], times = (tPrecond, tPrecond + tend), control = control)

    if !test
        println("*** done\n")
    end

    ################################################################################
    if !test
        println("Forward IV Measurement loop")
    end
    ################################################################################

    solForw = ChargeTransport.solve(ctsys, inival = solRev.u[end], times = (tPrecond + tend, tPrecond + 2 * tend), control = control)

    if Plotter !== nothing

        plot_densities!(vis[2, 2], ctsys, solForw.u[end], "End time", label_density)
        plot_solution!(vis[3, 1], ctsys, solForw.u[end], "End time", label_solution)
    end

    if !test
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Forward IV curve calculation")
    end
    ################################################################################

    IV = zeros(0) # for saving I-V data
    IVn = zeros(0); IVp = zeros(0)
    IVa = zeros(0); IVc = zeros(0)
    IVψ = zeros(0)

    tvalues = solForw.t
    number_tsteps = length(tvalues)
    biasValues = data.contactVoltageFunction[2].(tvalues[2:end])

    factory = ChargeTransport.TestFunctionFactory(ctsys)
    tf = ChargeTransport.testfunction(factory, [p.bregionLeft], [p.bregionRight])

    for istep in 2:number_tsteps

        Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
        inival = solForw.u[istep - 1]
        solution = solForw.u[istep]

        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)
        IEdgeOld = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, inival)

        push!(IVn, IEdge[p.iphin]); push!(IVp, IEdge[p.iphip])
        push!(IVψ, (IEdge[data.index_psi] - IEdgeOld[data.index_psi]) / Δt)

        push!(IVa, IEdge[p.iphia])
        push!(IVc, IEdge[p.iphic])
        push!(IV, IVψ[istep - 1] + IVn[istep - 1] + IVp[istep - 1] + IVa[istep - 1] + IVc[istep - 1])

    end

    if !test
        println("*** done\n")
    end

    ################################################################################

    if Plotter !== nothing
        scalarplot!(
            vis[3, 2],
            solPrecond.t,
            data.contactVoltageFunction[2].(solPrecond.t);
            clear = false,
            markershape = :cross,
            markersize = 8,
            color = "blue",
            label = "precondition",
            legend = :cc,
            xlabel = L"\text{time [s]}",
            ylabel = L"\text{voltage [V]}",
            title = "Applied bias over time"
        )

        scalarplot!(
            vis[3, 2],
            solRev.t,
            data.contactVoltageFunction[2].(solRev.t);
            clear = false,
            markershape = :cross,
            markersize = 8,
            color = "orange",
            label = "reverse",
            legend = :cc
        )

        scalarplot!(
            vis[3, 2],
            solForw.t,
            data.contactVoltageFunction[2].(solForw.t);
            clear = false,
            markershape = :cross,
            markersize = 8,
            color = "green",
            label = "forward",
            legend = :cc
        )

        scalarplot!(
            vis[4, 1],
            biasValues,
            -IV .* (cm^2) .* 1.0e3;
            linewidth = 2,
            label = "forward",
            legend = :cc,
            color = "blue",
            xlabel = L"\text{applied bias [V]}",
            ylabel = L"current density [$\frac{\text{mA}}{\text{cm}^2}$]",
            title = "Total current",
            clear = false
        )

        reveal(vis)
    end

    if test == false
        println("*** done\n")
    end

################################################################################
    if test == false
        println("Compute fill factor and efficiency")
    end
    ################################################################################

    if test == false
        IV = -IV
        bias = biasValues

        powerDensity = bias .* (IV)           # power density function
        MaxPD, indexPD = findmax(powerDensity)

        open_circuit = compute_open_circuit_voltage(bias, IV)

        IncLightPowerDens = 1000.0 * W / m^2

        fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

        efficiency = 100 * bias[indexPD] * (IV[indexPD]) / (IncLightPowerDens)
        JSC = IV[1] .* (cm)^(2) .* 1.0e3

        println(" ")
        println("The JSC                  is $(round(JSC, digits = 3)) mAcm^{-2}.")
        println("The fill factor          is $(round(fillfactor, digits = 2)) %.")
        println("The efficiency           is $(round(efficiency, digits = 2)) %.")
        println("The open circuit voltage is $(round(open_circuit, digits = 4)) V.")
        println(" ")
    end
    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        integral1 = integrated_density(ctsys, sol = solution, icc = p.iphia, ireg = p.regionPero)
        integral2 = integrated_density(ctsys, sol = solution, icc = p.iphic, ireg = p.regionPero)

        println("Calculated average vacancy density (anions): ", integral1 / data.regionVolumes[p.regionPero])
        println("Calculated average vacancy density (cations): ", integral2 / data.regionVolumes[p.regionPero])
        println(" ")
        if vacancyEnergyCalculation
            vacancyEnergy1 = data.params.bandEdgeEnergy[p.iphia, p.regionPero] / data.constants.q
            vacancyEnergy2 = data.params.bandEdgeEnergy[p.iphic, p.regionPero] / data.constants.q
            println("Value for vacancy energy (anions) is: ", vacancyEnergy1, " eV.")
            println("Value for vacancy energy (cations) is: ", vacancyEnergy2, " eV.")
            println("We recommend to calculate it on a fine grid and saving the values for later use.")
            println(" ")
        end
    end

    testval = sum(filter(!isnan, solutionEQ)) / length(solutionEQ) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -1.0528500585493203
    return main(test = true, vacancyEnergyCalculation = true) ≈ testval
end

end # module
