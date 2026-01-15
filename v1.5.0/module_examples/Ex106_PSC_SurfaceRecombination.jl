#=
# PSC device with surface recombination (1D).
([source code](@__SOURCE_URL__))

Simulating a three layer PSC device PCBM | MAPI | Pedot with mobile ions with a linear scan protocol.

Here, the surface recombination at internal boundaries is tested.

=#

module Ex106_PSC_SurfaceRecombination

using ChargeTransport
using ExtendableGrids
using GridVisualize
using LaTeXStrings

function main(;
        n = 6, Plotter = nothing,                   # supported Plotters are GLMakie and PythonPlot
        verbose = false, test = false,
        parameter_set = Params_PSC_PCBM_MAPI_Pedot, # choose the parameter set
        vacancyEnergyCalculation = false,           # assume the vacancy energy level is either given or not
    )

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    @local_unitfactors μm cm s ns V K ps Hz W

    # parameter
    p = parameter_set()

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate = 1.0 * V / s
    vend = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = vend / scanrate

    ## Define scan protocol function
    function linearScanProtocol(t)
        return if t == Inf
            0.0
        else
            scanrate * t
        end
    end

    ## Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, linearScanProtocol]

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
        p.h_ndoping / 2, p.h_ndoping,
        p.h_ndoping / (0.7 * δ), p.h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ndoping, p.h_ndoping + p.h_intrinsic / k,
        p.h_intrinsic / (5.1 * δ), p.h_intrinsic / (1.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ndoping + p.h_intrinsic / k, p.h_ndoping + p.h_intrinsic,
        p.h_intrinsic / (1.1 * δ), p.h_intrinsic / (5.1 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ndoping + p.h_intrinsic, p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2,
        p.h_pdoping / (1.3 * δ), p.h_pdoping / (0.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, p.h_ndoping + p.h_intrinsic + p.h_pdoping, step = p.h_pdoping / (0.8 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [p.heightLayers[1]], p.regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [p.heightLayers[1]], [p.heightLayers[2]], p.regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [p.heightLayers[2]], [p.heightLayers[3]], p.regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], p.bregionDonor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [p.h_total], [p.h_total], p.bregionAcceptor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2, tol = 1.0e-18) # second inner interface

    if Plotter !== nothing
        vis = GridVisualizer(; Plotter, layout = (3, 2), size = (1550, 800))
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

    ## Initialize Data instance and fill in data
    data = Data(grid, p.numberOfCarriers, contactVoltageFunction = contactVoltageFunction)

    ## Possible choices: Stationary, Transient
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

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[p.bregionAcceptor] = OhmicContact
    data.boundaryType[p.bregionJ1] = InterfaceRecombination
    data.boundaryType[p.bregionJ2] = InterfaceRecombination
    data.boundaryType[p.bregionDonor] = OhmicContact

    ## Present ionic vacancies in perovskite layer
    ## by default the statistics function is set to FermiDiracMinusOne to limit ion depletion
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    data.params = Params(p)

    if !vacancyEnergyCalculation
        data.params.bandEdgeEnergy[p.iphia, p.regionIntrinsic] = p.Ea[p.regionIntrinsic]
    end

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

    solution = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = vacancyEnergyCalculation)
    inival = solution

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement")
    end
    ################################################################################

    control.Δt = 5.0e-2
    control.Δt_grow = 1.03
    ## calculation of solution
    sol = ChargeTransport.solve(ctsys, inival = inival, times = (0.0, tend), control = control)

    tvalues = sol.t
    number_tsteps = length(tvalues)

    ## for saving I-V data
    IV = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values

    for istep in 2:number_tsteps

        Δt = tvalues[istep] - tvalues[istep - 1]  # Time step size
        inival = sol.u[istep - 1]
        solution = sol.u[istep]

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)
        push!(IV, current)

        inival = solution

    end # time loop

    solution = sol.u[end]

    biasValues = contactVoltageFunction[p.bregionAcceptor].(tvalues)

    # TODO MO: Ergibt das Sinn? Habe ich für die Title gemacht.
    maxBias = biasValues[end]

    if Plotter !== nothing
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
        ## add labels for anion vacancy
        label_energy[1, p.iphia] = "\$E_a-q\\psi\$"; label_energy[2, p.iphia] = "\$ - q \\varphi_a\$"; label_BEE[p.iphia] = "\$E_a\$"
        label_density[p.iphia] = "\$ n_a \$";      label_solution[p.iphia] = "\$ \\varphi_a\$"

        plot_densities!(vis[1, 2], ctsys, solution, "Applied voltage Δu = $maxBias", label_density)
        plot_solution!(vis[2, 1], ctsys, solution, "Applied voltage Δu = $maxBias", label_solution)

        ###########
        #TODO MO: Im Referenceplot mit PyPlot sind mehr Schritte, wieso?
        scalarplot!(
            vis[2, 2],
            tvalues,
            biasValues;
            markershape = :circle,
            markersize = 8,
            color = :blue,
            xlabel = L"\text{time [s]}",
            ylabel = L"\text{bias [V]}",
            title = "Applied bias over time"
        )
        ###########
        scalarplot!(
            vis[3, 1],
            biasValues[2:end],
            IV .* (cm^2) .* 1.0e3,
            linewidth = 2,
            color = "darkblue",
            xlabel = L"\text{bias [V]}",
            ylabel = L"current density [$\text{mAcm}^{-2}$]",
            title = "Total current"
        )

        reveal(vis)
    end

    if test == false
        println("*** done\n")
    end

    if test == false
        integral = integrated_density(ctsys, sol = solution, icc = p.iphia, ireg = p.regionIntrinsic)

        println("Calculated average vacancy density is: ", integral / data.regionVolumes[p.regionIntrinsic])
        println(" ")
        if vacancyEnergyCalculation
            vacancyEnergy = data.params.bandEdgeEnergy[p.iphia, p.regionIntrinsic] / data.constants.q
            println("Value for vacancy energy is: ", vacancyEnergy, " eV. Save this value for later use.")
            println("We recommend to calculate it on a fine grid.")
            println(" ")
        end
    end

    testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end # main

function test()
    testval = -0.5965444263524541; testvalvacancyEnergyCalculation = -0.5966729068541846
    return main(test = true, vacancyEnergyCalculation = true) ≈ testval && main(test = true, vacancyEnergyCalculation = false) ≈ testvalvacancyEnergyCalculation
end


end # module
