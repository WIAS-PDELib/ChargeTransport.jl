#=
# PSC device with surface recombination (1D).
([source code](SOURCE_URL))

Simulating a three layer PSC device PCBM | MAPI | Pedot with mobile ions with a linear scan protocol.

Here, the surface recombination at internal boundaries is tested.

=#

module Ex106_PSC_SurfaceRecombination

using ChargeTransport
using ExtendableGrids
using PyPlot

# for convenience
parametersdir = ChargeTransport.parametersdir

# you can also use other Plotters, if you add them to the example file
function main(;
        n = 6, Plotter = PyPlot, plotting = false,
        verbose = false, test = false,
        parameter_file = parametersdir("Params_PSC_PCBM_MAPI_Pedot.jl"), # choose the parameter file
    )

    if plotting
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file) # include the parameter file we specified

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate = 1.0 * V / s
    ntsteps = 31
    vend = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = vend / scanrate

    ## with fixed timestep sizes we can calculate the times a priori
    tvalues = range(0, stop = tend, length = ntsteps)

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

    coord_n_u = collect(range(0.0, h_ndoping / 2, step = h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        h_ndoping / 2, h_ndoping,
        h_ndoping / (0.7 * δ), h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        h_ndoping, h_ndoping + h_intrinsic / k,
        h_intrinsic / (5.1 * δ), h_intrinsic / (1.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        h_ndoping + h_intrinsic / k, h_ndoping + h_intrinsic,
        h_intrinsic / (1.1 * δ), h_intrinsic / (5.1 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        h_ndoping + h_intrinsic, h_ndoping + h_intrinsic + h_pdoping / 2,
        h_pdoping / (1.3 * δ), h_pdoping / (0.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(h_ndoping + h_intrinsic + h_pdoping / 2, h_ndoping + h_intrinsic + h_pdoping, step = h_pdoping / (0.8 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionDonor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionAcceptor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = 1.0e-18) # second inner interface

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
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionJ2] = InterfaceRecombination
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

    params = Params(grid, numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp
    params.chargeNumbers[iphia] = za

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
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

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[iphin, bregionJ1] = Nn[regionIntrinsic]
    params.bDensityOfStates[iphip, bregionJ1] = Np[regionIntrinsic]

    params.bDensityOfStates[iphin, bregionJ2] = Nn[regionIntrinsic]
    params.bDensityOfStates[iphip, bregionJ2] = Np[regionIntrinsic]

    params.bBandEdgeEnergy[iphin, bregionJ1] = En[regionIntrinsic]
    params.bBandEdgeEnergy[iphip, bregionJ1] = Ep[regionIntrinsic]

    params.bBandEdgeEnergy[iphin, bregionJ2] = En[regionIntrinsic]
    params.bBandEdgeEnergy[iphip, bregionJ2] = Ep[regionIntrinsic]

    ## for surface recombination
    params.recombinationSRHvelocity[iphin, bregionJ1] = 1.0e1 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ1] = 1.0e5 * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    params.recombinationSRHvelocity[iphin, bregionJ2] = 1.0e7 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ2] = 1.0e1 * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[iphin, regionDonor] = Cn
    params.doping[iphip, regionAcceptor] = Cp
    params.doping[iphia, regionIntrinsic] = Ca

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

    control = SolverControl()
    control.verbose = verbose
    control.damp_initial = 0.9
    control.damp_growth = 1.61 # >= 1
    control.max_round = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

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

    for istep in 2:ntsteps

        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: Δt = $(t)")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival = solution

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        if plotting
            label_solution, label_density, label_energy = set_plotting_labels(data)
            label_solution[iphia] = "\$ \\varphi_a\$"

            Plotter.clf()
            plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu)", label_solution)
            Plotter.pause(0.5)
        end

    end # time loop

    ##res = [biasValues, IV]

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
    return testval


end # main

function test()
    testval = -0.5963272869004673
    return main(test = true) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
