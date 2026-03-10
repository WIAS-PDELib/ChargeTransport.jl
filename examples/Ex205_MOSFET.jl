#=
# Metal-Oxide-Semiconductor Field-Effect Transistor (MOSFET).
([source code](@__SOURCE_URL__))

We simulate charge transport in an n-channel MOSFET consisting of a single material, silicon.
The drain-source current is controlled by an additional voltage applied at the gate contact.
Simulations are performed on a two-dimensional grid.
=#

module Ex205_MOSFET

using ChargeTransport
using VoronoiFVM
using GridVisualize
using ExtendableGrids


# Generates the grid with refinement level `grid_refinement`.
# Resulting node counts:
#   0 →   930
#   1 →  3599
#   2 →  8008
#   3 → 14157
# ...
# Input parameter `p` stores all constants and parameters
function generate_2Dgrid(p; grid_refinement = 0)

    # unit factors
    @local_unitfactors μm

    # gridpoints x-direction
    x0 = -1.5 * μm
    x1 = -1.0 * μm
    x2 = -0.9 * μm
    x3 = -0.8 * μm
    x4 = 0.8 * μm
    x5 = 0.9 * μm
    x6 = 1.0 * μm
    x7 = 1.5 * μm

    # gridpoints y-direction
    y0 = -2.4 * μm
    y1 = -0.5 * μm
    y2 = -0.2 * μm
    y3 = 0.0 * μm

    # refinement in x-direction
    X = collect(range(x0, x7, length = 31 + grid_refinement * 30))

    # refinement in y-direction
    Y1 = collect(range(y0, y1, length = 20 + grid_refinement * 19))
    Y2 = collect(range(y1, y3, length = 11 + grid_refinement * 10))
    Y = glue(Y1, Y2)

    grid = simplexgrid(X, Y)

    # cell regions
    cellmask!(grid, [x0, y2], [x3, y3], p.region_n)
    cellmask!(grid, [x4, y2], [x7, y3], p.region_n)
    cellmask!(grid, [x3, y2], [x4, y3], p.region_p)
    cellmask!(grid, [x0, y0], [x7, y2], p.region_p)

    # boundary regions
    bfacemask!(grid, [x2, y3], [x5, y3], p.bregion_gate)
    bfacemask!(grid, [x0, y3], [x1, y3], p.bregion_drain)
    bfacemask!(grid, [x6, y3], [x7, y3], p.bregion_source)
    bfacemask!(grid, [x0, y0], [x7, y0], p.bregion_bulk)
    bfacemask!(grid, [x0, y0], [x0, y3], p.bregion_neumann)
    bfacemask!(grid, [x7, y0], [x7, y3], p.bregion_neumann)
    bfacemask!(grid, [x1, y3], [x2, y3], p.bregion_neumann)
    bfacemask!(grid, [x5, y3], [x6, y3], p.bregion_neumann)

    return grid
end


# THX at JF!!!
# https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
function tridata(grid::ExtendableGrid)
    coord = grid[Coordinates]
    cellnodes = Matrix(grid[CellNodes])
    return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
end


function main(;
        Plotter = nothing,
        test = false,
        grid_refinement = 0,
        parameter_set = Params_MOSFET
    )

    if Plotter !== nothing && nameof(Plotter) !== :PythonPlot
        error("Plotting only possible for Plotter = PythonPlot")
    end

    if Plotter !== nothing
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # unit factors
    @local_unitfactors μm cm s ns V K

    # import constants and parameters from parameter file
    p = parameter_set()

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    grid = generate_2Dgrid(p, grid_refinement = grid_refinement)

    if Plotter !== nothing
        GridVisualize.gridplot(grid, Plotter = Plotter, resolution = (600, 450), title = "Grid", xlabel = "Width [m]", ylabel = "Height [m]")
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # initialize Data instance
    data = Data(grid, p.numberOfCarriers, constants = p.constants)

    # model
    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true,
    )

    # boundary model
    data.boundaryType[p.bregion_gate] = GateContact
    data.boundaryType[p.bregion_drain] = OhmicContact
    data.boundaryType[p.bregion_source] = OhmicContact
    data.boundaryType[p.bregion_bulk] = OhmicContact

    # flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

    # add parameters to data
    data.params = ChargeTransport.Params(p)

    # Definition ChargeTransport System
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

    control = ChargeTransport.SolverControl()
    control.verbose = ""
    control.maxiters = 70
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.damp_initial = 0.5
    control.max_round = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    solution_eq = equilibrium_solve!(ctsys, control = control)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    biasValues_gate = range(0.0, stop = 5.0, length = 4)
    biasValues_drain = range(0.0, stop = 5.0, length = 43)

    IV_drain = zeros(0)

    inival = copy(solution_eq)
    solution = copy(solution_eq)

    # Bias Loop Gate
    for Δu_gate in biasValues_gate

        if test == false
            println("bias value at gate: Δu = ", Δu_gate, " V")
        end

        set_contact!(ctsys, p.bregion_gate, Δu = Δu_gate)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

    end

    # Bias Loop Drain
    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
    tf = VoronoiFVM.testfunction(factory, [p.bregion_drain], [p.bregion_source])

    for Δu_drain in biasValues_drain

        if test == false
            println("bias value at drain: Δu = ", Δu_drain, " V")
        end

        set_contact!(ctsys, p.bregion_drain, Δu = Δu_drain)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)

        current = IEdge[p.iphin] + IEdge[p.iphip]

        push!(IV_drain, abs(p.zaus * current))
    end

    if test == false
        println("*** done\n")
    end

    if Plotter !== nothing
        plotting(p, grid, ctsys, solution, biasValues_drain, IV_drain, Plotter)
    end

    testval = sum(solution) / length(solution)
    return testval
end # main


# function for testing
function test()
    testval = 0.4992399114009203
    return main(test = true) ≈ testval
end


# function for plotting
function plotting(p, grid, ctsys, solution, biasValues_drain, IV_drain, Plotter)
    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    ########## Surface Electrostatic Potential ##########
    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], solution[p.ipsi, :])
    Plotter.title("Electrostatic Potential (Equilibrium)", fontsize = 16)
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.zlabel("potential [V]", fontsize = 12)
    ax = Plotter.gca()
    ax.tick_params(axis = "x", labelsize = 10)
    ax.tick_params(axis = "y", labelsize = 10)
    ax.tick_params(axis = "z", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.zaxis.get_offset_text().set_fontsize(10)
    current_figure = Plotter.gcf()
    display(current_figure)

    ########## Electron Density ##########
    nn = p.Nc .* exp.(p.zn * (p.constants.q * (solution[p.iphin, :] .- solution[p.ipsi, :]) .+ p.Ec) ./ (p.constants.k_B * p.T))
    nn_cm = nn * 1.0e-6

    Plotter.figure()
    Plotter.tripcolor(
        tridata(grid)..., nn_cm,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e6, vmax = 1.0e20),
        shading = "gouraud",
        rasterized = true
    )
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.title("Electron Density", fontsize = 16)
    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label("Density [\$\\mathrm{cm}^{-3}\$]", fontsize = 12)
    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    display(current_figure)

    ########## Hole Density ##########
    np = p.Nv .* exp.(p.zp * (p.constants.q * (solution[p.iphip, :] .- solution[p.ipsi, :]) .+ p.Ev) ./ (p.constants.k_B * p.T))
    np_cm = np * 1.0e-6

    Plotter.figure()
    Plotter.tripcolor(
        tridata(grid)..., np_cm,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e3, vmax = 1.0e16),
        shading = "gouraud",
        rasterized = true
    )
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.title("Hole Density", fontsize = 16)
    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label("Density [\$\\mathrm{cm}^{-3}\$]", fontsize = 12)
    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    display(current_figure)

    ########## Current Density ##########
    nf = -VoronoiFVM.nodeflux(ctsys.fvmsys, solution)
    j = nf[:, p.iphin, :] ./ (1.0e4 * p.zaus)
    jx = j[1, :]; jy = j[2, :]
    jabs = sqrt.(jx .^ 2 .+ jy .^ 2)

    L = 0.03 * (maximum(YY) - minimum(YY))  # length arrows
    eps = 1.0e-30
    ux = jx ./ (jabs .+ eps);  uy = jy ./ (jabs .+ eps)
    idx = 1:1:length(XX)

    Plotter.figure()
    Plotter.tripcolor(
        tridata(grid)...,
        jabs,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e-8, vmax = 1.0e4),
        shading = "gouraud",
        rasterized = true
    )
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")
    Plotter.quiver(
        XX[idx], YY[idx], (L * ux)[idx], (L * uy)[idx],
        color = "black",
        scale = 1.0,
        angles = "xy",
        scale_units = "xy",
        width = 0.002 # thickness arrows
    )
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.title("Electron Current Density", fontsize = 16)
    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label("|j| [A/cm^2]", fontsize = 12)
    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    display(current_figure)

    ########## IV Curve ##########
    Plotter.figure()
    Plotter.plot(biasValues_drain, IV_drain, marker = "o", markersize = 3)
    Plotter.grid()

    Plotter.title("IV Curve", fontsize = 16)
    Plotter.xlabel("voltage drain [V]", fontsize = 12)
    Plotter.ylabel("total current drain [A]", fontsize = 12)

    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(1.0))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)

    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    return display(current_figure)
end # plotting

end # module
