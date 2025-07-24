#=
# PSC device on 2D domain (Tensor grid).
([source code](@__SOURCE_URL__))

Simulating a three layer PSC device PCBM | MAPI | Pedot with mobile ions.
The simulations are
performed in 2D on a tensor grid, out of equilibrium and with abrupt interfaces.

=#

module Ex201_PSC_tensorGrid

using ChargeTransport
using ExtendableGrids
using PyPlot

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(;
        n = 3, Plotter = PyPlot, plotting = false, verbose = "", test = false,
        parameter_set = Params_PSC_PCBM_MAPI_Pedot, # choose the parameter set
    )

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # parameter
    p = parameter_set()

    bregionNoFlux = 3
    height = 5.0e-6 * cm

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate = 0.4 * V / s
    number_tsteps = 31
    endVoltage = voltageAcceptor # bias goes until the given voltage at acceptor boundary

    ## with fixed timestep sizes we can calculate the times a priori
    tend = endVoltage / scanrate
    tvalues = range(0, stop = tend, length = number_tsteps)

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

    coord_n_u = collect(range(0.0, p.h_ndoping / 2, step = p.h_ndoping / (0.6 * δ)))
    coord_n_g = geomspace(
        p.h_ndoping / 2, p.h_ndoping,
        p.h_ndoping / (0.7 * δ), p.h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ndoping, p.h_ndoping + p.h_intrinsic / k,
        p.h_intrinsic / (5.1 * δ), p.h_intrinsic / (1.0 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ndoping + p.h_intrinsic / k, p.h_ndoping + p.h_intrinsic,
        p.h_intrinsic / (1.0 * δ), p.h_intrinsic / (5.1 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ndoping + p.h_intrinsic, p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2,
        p.h_pdoping / (1.3 * δ), p.h_pdoping / (0.3 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, p.h_ndoping + p.h_intrinsic + p.h_pdoping, step = p.h_pdoping / (0.6 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord_length = glue(coord, coord_p_u, tol = 10 * t)

    height_L = geomspace(0.0, height / 2, height / (0.4 * δ), height / (0.4 * δ))
    height_R = geomspace(height / 2, height, height / (0.4 * δ), height / (0.4 * δ))
    coord_height = glue(height_L, height_R, tol = 10 * t)

    grid = simplexgrid(coord_length, coord_height)

    ## specify inner regions
    cellmask!(grid, [0.0, 0.0], [p.h_ndoping, height], regionDonor, tol = 1.0e-18)
    cellmask!(grid, [p.h_pdoping, 0.0], [p.h_ndoping + p.h_intrinsic, height], regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid, [p.h_ndoping + p.h_intrinsic, 0.0], [h_total, height], regionAcceptor, tol = 1.0e-18)

    ## specify outer regions
    ## metal interfaces
    bfacemask!(grid, [0.0, 0.0], [0.0, height], bregionDonor)            # BregionNumber = 1
    bfacemask!(grid, [h_total, 0.0], [h_total, height], bregionAcceptor) # BregionNumber = 2

    ## no flux interfaces [xmin, ymin], [xmax, ymax]
    bfacemask!(grid, [0.0, 0.0], [h_total, 0.0], bregionNoFlux)          # BregionNumber = 3
    bfacemask!(grid, [0.0, height], [h_total, height], bregionNoFlux)    # BregionNumber = 3

    if plotting
        gridplot(grid, Plotter = Plotter, resolution = (600, 400), linewidth = 0.5, legend = :lt)
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

    ## Initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor] = OhmicContact

    ## Present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp
    params.chargeNumbers[iphia] = za

    for ireg in 1:numberOfRegions # region data

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn[ireg]
        params.densityOfStates[iphip, ireg] = Np[ireg]
        params.densityOfStates[iphia, ireg] = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg] = En[ireg]
        params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg] = Ea[ireg]

        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]
        params.mobility[iphia, ireg] = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, params, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, params, EI[ireg])

    end

    ## interior doping
    params.doping[iphin, regionDonor] = Cn
    params.doping[iphia, regionIntrinsic] = Ca
    params.doping[iphip, regionAcceptor] = Cp

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 300
    control.max_round = 5
    control.damp_initial = 0.1
    control.damp_growth = 1.21 # >= 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    if plotting # currently, plotting the solution was only tested with PyPlot.
        ipsi = data.index_psi
        X = grid[Coordinates][1, :]
        Y = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[ipsi, :])
        Plotter.title("Electrostatic potential \$ \\psi \$ (Equilibrium)")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        Plotter.tight_layout()
        ################
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[iphin, :])
        Plotter.title("quasi Fermi potential \$ \\varphi_n \$ (Equilibrium)")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    ## for saving I-V data
    IV = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values

    for istep in 2:number_tsteps

        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage; set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t) s")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        inival = solution
    end # time loop

    if plotting
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[ipsi, :])
        Plotter.title("Electrostatic potential \$ \\psi \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        ################
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[iphin, :])
        Plotter.title("quasi Fermi potential \$ \\varphi_n \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        ################
        Plotter.figure()
        Plotter.plot(biasValues, IV .* (cm)^2 / height, label = "", linewidth = 3, marker = "o")
        Plotter.grid()
        Plotter.ylabel("total current [A]") #
        Plotter.xlabel("Applied Voltage [V]")
    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -0.5818799192233242
    return main(test = true) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
