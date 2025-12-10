#=
# 2D textured perovskite solar cell.
([source code](@__SOURCE_URL__))

Simulating a three layer textured PSC device with mobile ions.
The simulations are performed in 2D on an unstructured grid.

=#

module Ex201_PSC_Textured

using ChargeTransport
using ExtendableGrids
using PyPlot
using VoronoiFVM

function generate_grid(; parameter_set)

    # use the destructuring operator to extract all the necessary parameters
    (;
        h_activePL, heightLayersPL, h_totalPL, h_ETL, h_HTL, regionETL, regionPero, regionHTL, bregionLeft,
        bregionRight, bregionJ1, bregionJ2,
    ) = parameter_set()

    t = 1.0e-13
    ######################

    h_active = h_activePL
    heightLayers = heightLayersPL
    h_total = h_totalPL

    h1 = 4.0e-3
    h2 = 3.0e-3
    h3 = 3.0e-3
    h4 = 3.0e-3

    coord_n1_1 = geomspace(0.0, h_ETL / 2, h2 * h_ETL, h3 * h_ETL, tol = t)
    coord_n1_2 = geomspace(h_ETL / 2, h_ETL, h3 * h_ETL, h1 * h_ETL, tol = t)

    coord_Ac_1 = geomspace(
        h_ETL, h_ETL + 0.15 * h_active,
        h4 * h_ETL, h1 * h_active, tol = t
    )
    coord_Ac_2 = geomspace(
        h_ETL + 0.15 * h_active, h_ETL + 0.5 * h_active,
        h1 * h_active, h3 * h_active, tol = t
    )
    coord_Ac_3 = geomspace(
        h_ETL + 0.5 * h_active, h_ETL + 0.85 * h_active,
        h3 * h_active, h1 * h_active, tol = t
    )
    coord_Ac_4 = geomspace(
        h_ETL + 0.85 * h_active, h_ETL + 1.0 * h_active,
        h1 * h_active, h4 * h_HTL, tol = t
    )

    coord_p_1 = geomspace(
        h_ETL + h_active, h_ETL + h_active + h_HTL / 2,
        h4 * h_HTL, h3 * h_HTL, tol = t
    )
    coord_p_2 = geomspace(
        h_ETL + h_active + h_HTL / 2, h_total,
        h3 * h_HTL, h4 * h_HTL, tol = t
    )

    coordn1 = glue(coord_n1_1, coord_n1_2, tol = t)
    coordAc = glue(coord_Ac_1, coord_Ac_2, tol = t)
    coordAc = glue(coordAc, coord_Ac_3, tol = t)
    coordAc = glue(coordAc, coord_Ac_4, tol = t)
    coordp = glue(coord_p_1, coord_p_2, tol = t)

    coord = glue(coordn1, coordAc, tol = t)
    coord = glue(coord, coordp, tol = t)
    grid = simplexgrid(coord)

    ## cellmask! for setting different inner regions
    cellmask!(grid, [0.0], [heightLayers[1]], regionETL, tol = t) # n-doped region = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionPero, tol = t) # pero region    = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionHTL, tol = t) # p-doped region = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionLeft, tol = t) # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionRight, tol = t) # outer right boundary

    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = t) # left pero boundary
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = t) # right pero boundary

    return grid

end

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(;
        Plotter = PyPlot, plotting = false, verbose = false, test = false,
        parameter_set = Params_PSC_C60_TripleCation_PTAA, # choose the parameter set
        vacancyEnergyCalculation = true,             # assume the vacancy energy level is either given or not
        vETL = 2000 * ufac"cm" / ufac"s", # surface reco velocity at ETL
        vHTL = 500 * ufac"cm" / ufac"s", # surface reco velocity at HTL
    )

    PyPlot.rc("font", family = "sans-serif", size = 14)
    PyPlot.rc("mathtext", fontset = "dejavusans")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    @local_unitfactors V cm m s W

    (; q, ε_0) = ChargeTransport.constants
    eV = q * V

    ######### primary data for I-V scan protocol ##########
    endVoltage = 1.2 * V
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

    # parameter
    p = parameter_set()

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    grid = generate_grid(parameter_set = parameter_set)

    if plotting
        gridplot(grid, Plotter = PyPlot, resolution = (600, 400), linewidth = 0.5, legend = :rc)
    end

    if !test
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in data
    data = Data(grid, p.numberOfCarriers, contactVoltageFunction = contactVoltageFunction)

    ## Possible choices: Stationary, Transient
    data.modelType = Transient

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    data.boundaryType[p.bregionLeft] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionRight] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionJ1] = InterfaceRecombination
    data.boundaryType[p.bregionJ2] = InterfaceRecombination

    data.generationModel = GenerationBeerLambert

    ## Present ionic vacancies in perovskite layer
    ## by default the statistics function is set to FermiDiracMinusOne to limit ion depletion
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionPero])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    data.params = Params(p)

    data.params.recombinationSRHvelocity[p.iphin, p.bregionJ1] = vETL
    data.params.recombinationSRHvelocity[p.iphip, p.bregionJ1] = vETL

    data.params.recombinationSRHvelocity[p.iphin, p.bregionJ2] = vHTL
    data.params.recombinationSRHvelocity[p.iphip, p.bregionJ2] = vHTL

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
    control.damp_initial = 0.5
    control.damp_growth = 1.61 # >= 1
    control.maxiters = 300
    control.max_round = 5
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.Δu_opt = Inf

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    #data.params.bandEdgeEnergy[3, 2] = -5.267 * eV
    solEQ = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = vacancyEnergyCalculation)
    inival = solEQ

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Loop for increasing bias (to start at V near VOC)")
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

    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, p.iphia] = "\$E_a-q\\psi\$"; label_energy[2, p.iphia] = "\$ - q \\varphi_a\$"; label_BEE[p.iphia] = "\$E_a\$"
        label_density[p.iphia] = "\$ n_a \$";      label_solution[p.iphia] = "\$ \\varphi_a\$"
    end

    if plotting
        figure(2)
        plot_densities(PyPlot, ctsys, sol0p0, "Initial condition", label_density)
        figure()
        plot_solution(PyPlot, ctsys, sol0p0, "Initial condition", label_solution)
    end

    if !test
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Preconditioning")
    end
    ################################################################################
    control.Δt = 1.0e-4
    control.Δt_min = 1.0e-5
    control.Δt_max = 1.0e-4
    control.Δt_grow = 1.03

    solPrecond = ChargeTransport.solve(ctsys, inival = inival, times = (0.0, tPrecond), control = control)

    if !test
        println("*** done\n")
    end

    ################################################################################
    if !test
        println("Reverse IV Measurement loop")
    end
    ################################################################################

    control.Δt = 6.0e-4 / scanrate
    control.Δt_min = 6.0e-3 / scanrate
    control.Δt_max = 8.0e-3 / scanrate
    control.Δt_grow = 1.03
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

    if plotting
        figure()
        plot_densities(PyPlot, ctsys, solForw[end], "End time", label_density)
        figure()
        plot_solution(PyPlot, ctsys, solForw[end], "End time", label_solution)
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
    IVa = zeros(0); IVψ = zeros(0)

    ######################
    ISRHn = zeros(0); IRadn = zeros(0)
    ISRHp = zeros(0); IRadp = zeros(0)
    ISRnL = zeros(0)
    ISRnR = zeros(0)
    IGen = zeros(0)

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
        push!(IVψ, (IEdge[p.ipsi] - IEdgeOld[p.ipsi]) / Δt)

        push!(IVa, IEdge[p.iphia])
        push!(IV, IVψ[istep - 1] + IVn[istep - 1] + IVp[istep - 1] + IVa[istep - 1])

        IntSRH = ChargeTransport.integrate(ctsys, SRHRecombination!, solution)
        IntRad = ChargeTransport.integrate(ctsys, RadiativeRecombination!, solution)
        IntSR = ChargeTransport.integrate(ctsys, SRRecombination!, solution, boundary = true)
        IntGen = ChargeTransport.integrate(ctsys, Photogeneration!, solution)

        IntSRHnSum = 0.0; IntRadnSum = 0.0
        IntSRHpSum = 0.0; IntRadpSum = 0.0

        for ii in 1:p.numberOfRegions
            IntSRHnSum = IntSRHnSum - IntSRH[p.iphin, ii]
            IntRadnSum = IntRadnSum - IntRad[p.iphin, ii]

            IntSRHpSum = IntSRHpSum + IntSRH[p.iphip, ii]
            IntRadpSum = IntRadpSum + IntRad[p.iphip, ii]
        end

        IntGenp = IntGen[p.iphip, p.regionPero]
        IntSRnL = - IntSR[p.iphin, p.bregionJ1]
        IntSRnR = - IntSR[p.iphin, p.bregionJ2]

        push!(ISRHn, IntSRHnSum); push!(ISRHp, IntSRHpSum)
        push!(IRadn, IntRadnSum); push!(IRadp, IntRadpSum)
        push!(ISRnL, IntSRnL)
        push!(ISRnR, IntSRnR)
        push!(IGen, IntGenp)

    end

    if plotting

        figure()
        tEnd = tPrecond + 2 * tend

        tt = collect(range(0.0, tEnd, length = 201))
        T = data.contactVoltageFunction[2]
        plot(tt, T.(tt), marker = "o")
        xlabel("time [s]")
        ylabel("voltage [V]")
        PyPlot.tight_layout()

        figure()
        plot(biasValues, -IV .* (cm^2) .* 1.0e3, linewidth = 5, color = "blue", label = "forward")

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]", fontsize = 17)
        PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
        PyPlot.tick_params(which = "both", labelsize = 18)
        PyPlot.tight_layout()

        #####################################
        figure()
        semilogy(biasValues, ISRHn .* (cm^2) .* 1.0e3, linewidth = 5, color = "blue", label = "SRH")
        semilogy(biasValues, IRadn .* (cm^2) .* 1.0e3, linewidth = 5, color = "red", label = "rad")
        semilogy(biasValues, IGen .* (cm^2) .* 1.0e3, linewidth = 5, color = "gold", label = "Gen")

        semilogy(biasValues, ISRnL .* (cm^2) .* 1.0e3, linewidth = 5, color = "black", label = "SR, left")

        semilogy(biasValues, ISRnR .* (cm^2) .* 1.0e3, linewidth = 5, color = "darkgreen", label = "SR, right")

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]", fontsize = 17)
        PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
        PyPlot.tick_params(which = "both", labelsize = 18)
        PyPlot.tight_layout()
    end

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
        integral = integrated_density(ctsys, sol = solEQ, icc = p.iphia, ireg = p.regionPero)
        mOmega = data.regionVolumes[p.regionPero]

        println("Calculated average vacancy density is: ", integral / mOmega)
        println(" ")
        vacancyEnergy = data.params.bandEdgeEnergy[p.iphia, p.regionPero] / q
        println("Value for vacancy energy is: ", vacancyEnergy, " eV. Save this value for later use.")
        println(" ")
    end

   testval = sum(filter(!isnan, solForw.u[end])) / length(solForw.u[end]) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -0.5862627265480347; testvalvacancyEnergyCalculation = -0.5871876928952634
    return main(test = true) ≈ testval && main(test = true, vacancyEnergyCalculation = false) ≈ testvalvacancyEnergyCalculation
end

end # module
