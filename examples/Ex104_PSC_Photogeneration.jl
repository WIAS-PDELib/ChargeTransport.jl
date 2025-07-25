#=
# PSC device with photogeneration rate (1D).
([source code](@__SOURCE_URL__))

Simulating a three layer PSC device TiO2 | MAPI | Pedot with mobile ions where the
ion vacancy accumulation is limited by the Fermi-Dirac integral of order -1.

We perform a linear scan protocol and try out different photogeneration rates.
=#

module Ex104_PSC_Photogeneration

using ChargeTransport
using ExtendableGrids
using PyPlot

# for convenience
parametersdir = ChargeTransport.parametersdir

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(;
        n = 5,
        Plotter = PyPlot,
        plotting = false, verbose = "", test = false,
        ########################
        parameter_set = Params_PSC_TiO2_MAPI_spiro, # choose the parameter set
        ########################
        userdefinedGeneration = false
    ) # you can choose between predefined and user-defined generation profiles

    @local_unitfactors μm cm eV s ns V K ps Hz W m

    if plotting
        Plotter.close("all")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # parameters
    p = parameter_set()

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate = 0.04 * V / s
    number_tsteps = 31
    endVoltage = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = endVoltage / scanrate

    ## Define scan protocol function
    function scanProtocol(t)

        if 0.0 <= t  && t <= tend
            biasVal = 0.0 + scanrate * t
        elseif t > tend  && t <= 2 * tend
            biasVal = scanrate * tend .+ scanrate * (tend - t)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    ## Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocol]

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

    coord_n_u = collect(range(0.0, p.h_ndoping / 2, step = p.h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        p.h_ndoping / 2,
        p.h_ndoping,
        p.h_ndoping / (0.7 * δ),
        p.h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ndoping,
        p.h_ndoping + p.h_intrinsic / k,
        p.h_intrinsic / (2.8 * δ),
        p.h_intrinsic / (2.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ndoping + p.h_intrinsic / k,
        p.h_ndoping + p.h_intrinsic,
        p.h_intrinsic / (2.1 * δ),
        p.h_intrinsic / (2.8 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ndoping + p.h_intrinsic,
        p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2,
        p.h_pdoping / (1.6 * δ),
        p.h_pdoping / (1.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, p.h_ndoping + p.h_intrinsic + p.h_pdoping, step = p.h_pdoping / (1.3 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [p.heightLayers[1]], p.regionDonor, tol = 1.0e-18) # n-doped region   = 1
    cellmask!(grid, [p.heightLayers[1]], [p.heightLayers[2]], p.regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [p.heightLayers[2]], [p.heightLayers[3]], p.regionAcceptor, tol = 1.0e-18) # p-doped region   = 3

    bfacemask!(grid, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1, tol = 1.0e-18)
    bfacemask!(grid, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2, tol = 1.0e-18)

    if plotting
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    if userdefinedGeneration

        subg1 = subgrid(grid, [p.regionDonor]); subg2 = subgrid(grid, [p.regionIntrinsic]); subg3 = subgrid(grid, [p.regionAcceptor])

        gen1 = zeros(length(subg1[Coordinates]) - 1); gen3 = zeros(length(subg3[Coordinates]) - 1)
        gen2 = p.incidentPhotonFlux[p.regionIntrinsic] .* p.absorption[p.regionIntrinsic] .* exp.(- p.absorption[p.regionIntrinsic] .* (subg2[Coordinates] .- p.generationPeak))

        ## we want to get agreement with the region-wise defined photogeneration
        X1 = subg1[Coordinates]; X2 = subg2[Coordinates]; X3 = subg3[Coordinates]

        h1end = X1[end] - X1[end - 1]; h2beg = X2[2] - X2[1]
        h2end = X2[end] - X2[end - 1]; h3beg = X3[2] - X3[1]

        # region reaction multiplies with h2beg/2 ( = | ω_k ∩ region2|) it visits the node only from region2
        # node reaction multiplies with (h1end+h2beg)/2 ( = | ω_k|)  as it visits the node from region 1 and region 2
        # therefore, we need the following weights
        # However, note that | ω_k ∩ region2| is not calculate explicitly but via the simplex
        # components (if we have cellwise loops)
        weight1 = h2beg / (h1end + h2beg) # ( = | ω_k ∩ region2| / | ω_k| )
        weight2 = h2end / (h2end + h3beg)

        gen2[1] = weight1 * gen2[1]; gen2[end] = weight2 * gen2[end]

        generationData = [gen1; gen2'; gen3]

        data = Data(
            grid, p.numberOfCarriers,
            contactVoltageFunction = contactVoltageFunction,
            generationData = generationData
        )
    else

        data = Data(
            grid, p.numberOfCarriers,
            contactVoltageFunction = contactVoltageFunction
        )

    end

    data.modelType = Transient
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )
    data.boundaryType[p.bregionAcceptor] = OhmicContact
    data.boundaryType[p.bregionDonor] = OhmicContact
    data.fluxApproximation .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    if userdefinedGeneration
        data.generationModel = GenerationUserDefined
    else
        data.generationModel = GenerationBeerLambert
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    data.params = Params(p)
    ctsys = System(grid, data, unknown_storage = :sparse)

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    if verbose == true
        control.verbose = verbose
    else
        control.verbose = "eda" # still print the time values
    end
    if test == true
        control.verbose = false # do not show time values in testing case
    end
    control.maxiters = 300
    control.max_round = 5
    control.damp_initial = 0.5
    control.damp_growth = 1.21 # >= 1
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
    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Loop for generation")
    end
    ################################################################################

    # these values are needed for putting the generation slightly on
    I = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need
    ## to fix it in the bias loop, since we have no applied bias. Otherwise we get convergence errors
    ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop

    solutionEQ = inival

    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia] = "\$ n_a \$";      label_solution[iphia] = "\$ \\varphi_a\$"

        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Initial condition", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Initial condition", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    ## put here back the homogeneous Neumann boundary conditions.
    ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 0.0
    ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

    sol = solve(ctsys, inival = inival, times = (0.0, tend), control = control)

    if plotting
        tsol = sol(tend)
        Plotter.figure()
        plot_densities(Plotter, ctsys, tsol, "Densities at end time", label_density)
        Plotter.tight_layout()
        Plotter.figure()
        plot_solution(Plotter, ctsys, tsol, "Solution at end time", label_solution)
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Reverse scan protocol")
    end
    ################################################################################
    inivalReverse = sol(tend)
    solReverse = solve(ctsys, inival = inivalReverse, times = (tend, 2 * tend), control = control)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Curve calculation")
    end
    ################################################################################

    factory = TestFunctionFactory(ctsys)
    tf = testfunction(factory, [p.bregionDonor], [p.bregionAcceptor])

    tvalues = sol.t
    number_tsteps = length(tvalues)
    biasValues = scanProtocol.(tvalues)
    IV = zeros(0)

    for istep in 2:number_tsteps
        Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
        inival = sol.u[istep - 1]
        solution = sol.u[istep]

        I = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii in 1:(p.numberOfCarriers + 1)
            current = current + I[ii]
        end

        push!(IV, current)

    end

    tvaluesReverse = solReverse.t
    number_tstepsReverse = length(tvaluesReverse)
    biasValuesReverse = scanProtocol.(tvaluesReverse)
    IVReverse = zeros(0)

    for istep in 2:number_tstepsReverse
        Δt = tvaluesReverse[istep] - tvaluesReverse[istep - 1] # Time step size
        inival = solReverse.u[istep - 1]
        solution = solReverse.u[istep]

        I = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii in 1:(p.numberOfCarriers + 1)
            current = current + I[ii]
        end

        push!(IVReverse, current)

    end

    if plotting
        Plotter.figure()
        Plotter.plot([tvalues tvaluesReverse], [biasValues biasValuesReverse], marker = "x")
        Plotter.xlabel("time [s]")
        Plotter.ylabel("voltage [V]")
        Plotter.grid()

        Plotter.figure()
        Plotter.plot(biasValues[2:end], -IV, linewidth = 5, label = "forward")
        Plotter.plot(biasValuesReverse[2:end], -IVReverse, linewidth = 5, label = "reverse")
        Plotter.grid()
        Plotter.legend()
        Plotter.xlabel("applied bias [V]")
        Plotter.ylabel("total current [A]")

        Plotter.figure()
        if userdefinedGeneration
            Plotter.plot(coord, data.generationData)
        else
            for ireg in 1:p.numberOfRegions
                subg = subgrid(grid, [ireg])
                Plotter.plot(subg[Coordinates]', BeerLambert(ctsys, ireg, subg[Coordinates])', label = "region $ireg")
            end

        end
        Plotter.legend()
        Plotter.grid()
        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("photogeneration [\$\\frac{1}{cm^3s}\$]")
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute fill factor and efficiency")
    end
    ################################################################################

    bias = biasValues[2:end]
    IV = -IV

    powerDensity = bias .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, IV)

    IncidentLightPowerDensity = 1000.0 * W / m^2

    efficiency = bias[indexPD] * IV[indexPD] / IncidentLightPowerDensity
    fillfactor = (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

    if test == false
        println("The fill factor is $fillfactor %.")
        println("The efficiency  is $efficiency %.")
        println("The open circuit voltage  is $open_circuit.")
    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solutionEQ)) / length(solutionEQ) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -1.052813874410313
    return main(test = true, userdefinedGeneration = false) ≈ testval && main(test = true, userdefinedGeneration = true) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
