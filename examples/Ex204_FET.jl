#=
# Field effect transistor.
([source code](@__SOURCE_URL__))

We consider an n-channel Metal-Oxide-Semiconductor (MOS) field effect transistor.
The material is silicon.
=#

module Ex204_FET

using ChargeTransport
using ExtendableGrids
using GLMakie
using GridVisualize
using PyPlot

function main(; plotting = true, Plotter = GLMakie, test = false)

    # for Makie plots
    if plotting
        Plotter.closeall()
    end

    # unit factors
    @local_unitfactors μm cm s ns V K

    # constants
    constants = ChargeTransport.teSCA_constants
    (; q, k_B, ε_0) = constants

    eV = q * V

    ########## charge carriers ##########
    iphin = 1                            # quasi Fermi potential for electrons
    iphip = 2                            # quasi Fermi potential for holes
    ipsi = 3 # needed for PyPlot Plotting
    numberOfCarriers = 2

    ########## device geometry ##########
    # region numbers
    region_gate = 1
    region_drain = 2
    region_source = 3
    region_bulk = 4

    # boundary region numbers
    bregion_gate = 1
    bregion_drain = 2
    bregion_source = 3
    bregion_bulk = 4
    bregion_neutral = 5

    zaus = 15.0e-4 * cm                                 # depth of the device  [cm]
    thickness_ox = 0.044e-4 * cm                        # oxide thickness on gate [cm]

    ########## physical values of Si at room temperature ##########

    Ec = 0.562 * eV                                     # conduction band-edge energy
    Ev = -0.562 * eV                                    # valence band-edge energy
    Nc = 2.86e19 / (cm^3)                                # conduction band density of states
    Nv = 3.1e19 / (cm^3)                                # valence band densitiy of states

    mun = 1200.0 * (cm^2) / (V * s)                     # electron mobility
    mup = 350.0 * (cm^2) / (V * s)                      # hole mobility
    εr = 11.68 * 1.0                                    # relative dielectric permittivity of Si
    εr_ox = 3.9 * 1.0                                   # relative dielectric permittivity of SiO2
    T = 300.0 * K                                       # room temperature

    # Recombination parameters
    Auger = 1.0e-29 * cm^6 / s
    SRH_TrapDensity = 1.0e10 / cm^3
    SRH_LifeTime = 1.0 * ns
    Radiative = 1.0e-10 * cm^3 / s

    # Doping (Vereinfacht, nur eine Dotierung pro Region)
    Na_gate = 1.0e16 / cm^3
    Nd_drain = 1.0e19 / cm^3
    Nd_source = 1.0e19 / cm^3
    Na_bulk = 1.0e15 / cm^3

    # Voltage information
    U_add_gate = 0.55 * V                   # contact voltage on gate [V], siehe unten

    # DA: As these are some surface charges, I think, we need to put either here or in the implementation of the Gate BC the elementary charge.
    qss = 6.0e10 / (cm^2)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # Gridpoints x-direction
    x0 = -1.5 * μm              # beginning left contact (drain)
    x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = 0.9 * μm               # end gate contact
    x4 = 1.0 * μm               # beginning right contact (source)
    x5 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y1 = -0.2 * μm              # beginning n-channel (gate region)
    y2 = 0.0 * μm               # top of the device

    # Refinement x-direction
    X1 = geomspace(x0, x1, 2.0e-7, 5.0e-8)
    X2 = collect(range(x1, x2, length = 3))
    X3 = collect(x2:(0.2 * μm):x3)
    X4 = collect(range(x3, x4, length = 3))
    X5 = geomspace(x4, x5, 5.0e-8, 2.0e-7)
    X_temp = glue(X1, X2)
    X_temp = glue(X_temp, X3)
    X_temp = glue(X_temp, X4)
    X = glue(X_temp, X5)

    # Refinement y-direction
    Y1 = collect(-2.4:0.2:-1.2) .* μm
    Y2 = collect(-1.2:0.1:-0.1) .* μm
    Y3 = geomspace(-0.1 * μm, 0.0 * μm, 4.0e-8, 4.0e-8)
    Y12 = glue(Y1, Y2)
    Y = glue(Y12, Y3)

    grid = simplexgrid(X, Y)

    # cell regions
    cellmask!(grid, [x1, y1], [x4, y2], region_gate)
    cellmask!(grid, [x0, y1], [x1, y2], region_drain)
    cellmask!(grid, [x4, y1], [x5, y2], region_source)
    cellmask!(grid, [x0, y0], [x5, y1], region_bulk)

    # boundary regions

    bfacemask!(grid, [x2, y2], [x3, y2], bregion_gate)
    bfacemask!(grid, [x0, y2], [x1, y2], bregion_drain)
    bfacemask!(grid, [x4, y2], [x5, y2], bregion_source)
    bfacemask!(grid, [x0, y0], [x5, y0], bregion_bulk)
    bfacemask!(grid, [x0, y0], [x0, y2], 0)
    bfacemask!(grid, [x5, y0], [x5, y2], 0)
    bfacemask!(grid, [x1, y2], [x2, y2], 0)
    bfacemask!(grid, [x3, y2], [x4, y2], 0)

    if plotting
        vis = GridVisualizer(; Plotter, layout = (2, 2), size = (1200, 600))
        gridplot!(vis[1, 1], grid; Plotter, title = "Grid", show = true)

        # ohne Visualizer
        #fig_grid = ChargeTransport.gridplot(grid, Plotter = Plotter, title = "Grid")
        #display(fig_grid)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

    # bereits bestehende Parameter
    params.temperature = T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ireg in 1:grid[NumCellRegions] # region data
        params.dielectricConstant[ireg] = εr * ε_0
        params.dielectricConstantOxide[ireg] = εr_ox * ε_0
        params.thicknessOxide[ireg] = thickness_ox
        params.surfacechargeDensityGate[ireg] = qss
        params.additionalVoltageGate[ireg] = U_add_gate

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = mun
        params.mobility[iphip, ireg] = mup

        # recombination parameters, wird das immer übergeben? auch wenn Recombination = false?
        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg] = Auger
        params.recombinationAuger[iphip, ireg] = Auger
    end

    # doping                                                        #     n    p
    params.doping[iphin, region_drain] = Nd_drain                  #   [Nd  0.0;
    params.doping[iphin, region_source] = Nd_source                #    Nd  0.0;
    params.doping[iphip, region_gate] = Na_gate                    #    0.0  Na;
    params.doping[iphip, region_bulk] = Na_bulk                    #    0.0  Na]

    # Initialize Data instance
    data = Data(grid, numberOfCarriers, constants = constants)

    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = false
    )

    # boundary model
    data.boundaryType[bregion_gate] = GateContact
    data.boundaryType[bregion_drain] = OhmicContact
    data.boundaryType[bregion_source] = OhmicContact
    data.boundaryType[bregion_bulk] = OhmicContact

    # flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

    data.params = params

    # Definition ChargeTransport System
    ctsys = System(grid, data, unknown_storage = :sparse)

    # Add doping on boundary (if not added 0 on boundary)
    params.bDoping[:, :] .= 0.0
    params.bDoping[iphin, bregion_drain] = Nd_drain                  #   [Nd  0.0;
    params.bDoping[iphin, bregion_source] = Nd_source                #    Nd  0.0;
    params.bDoping[iphip, bregion_gate] = Na_gate                    #    0.0  Na;
    params.bDoping[iphip, bregion_bulk] = Na_bulk                    #    0.0  Na]

    params.bDensityOfStates[iphin, :] .= Nc
    params.bDensityOfStates[iphip, :] .= Nv

    params.bBandEdgeEnergy[iphin, :] .= Ec
    params.bBandEdgeEnergy[iphip, :] .= Ev

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
    control.verbose = true
    control.maxiters = 70
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.damp_initial = 0.5
    control.max_round = 3
    #control.damp_growth = 1.21 # >= 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # Solution in equilibrium
    solution_eq = equilibrium_solve!(ctsys, control = control)

    if test == false
        println("*** done\n")
    end

    if plotting
        ################################################################################
        # Surface plot equlibrium with Makie
        ################################################################################
        psi_eq = solution_eq[3, :]

        # Get grid coordinates
        coords = grid[Coordinates]
        nx = length(unique(coords[1, :]))
        ny = length(unique(coords[2, :]))

        # reshape to 2D arrays
        X = reshape(coords[1, :], nx, ny)
        Y = reshape(coords[2, :], nx, ny)
        Z = reshape(psi_eq, nx, ny)

        fig = GLMakie.Figure()
        ax = Axis3(
            fig[1, 1];
            xlabel = "length [m]",
            ylabel = "height [m]",
            zlabel = "potential [V]",
            title = "Electrostatic potential (Equilibrium)"
        )

        surface!(ax, X, Y, Z; shading = false)
        wireframe!(ax, X, Y, Z; color = :black, linewidth = 0.5)

        display(GLMakie.Screen(), fig) # open new window

        ################################################################################
        # Density plots equilibrium with PyPlot
        ################################################################################
        # MO: Just one figure plottet at a time
        function tridata(grid::ExtendableGrid)
            coord = grid[Coordinates]
            cellnodes = Matrix(grid[CellNodes])
            return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
        end

        ## visualizing with PyPlot
        vmin = 1.0e15
        vmax = 1.0e25

        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution_eq[iphin, :] .- solution_eq[ipsi, :]) .+ Ec) ./ (k_B * T))
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution_eq[iphip, :] .- solution_eq[ipsi, :]) .+ Ev) ./ (k_B * T))

        fig_n = PyPlot.figure()
        PyPlot.tripcolor(tridata(grid)..., vcat(nn...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        PyPlot.xlabel(" \$x\$ [nm]", fontsize = 17)
        PyPlot.ylabel(" \$y\$ [nm]", fontsize = 17)
        PyPlot.title("n density")
        PyPlot.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        PyPlot.tight_layout()

        display(fig_n)

        fig_h = PyPlot.figure()
        PyPlot.tripcolor(tridata(grid)..., vcat(np...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        PyPlot.xlabel(" \$x\$ [nm]", fontsize = 17)
        PyPlot.ylabel(" \$y\$ [nm]", fontsize = 17)
        PyPlot.title("h density")
        PyPlot.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        PyPlot.tight_layout()

        #display(fig_h)

        ################################################################################
        # Density plots equilibrium with Makie
        ################################################################################
        #=
        grid = ctsys.fvmsys.grid
        numberOfRegions = grid[NumCellRegions]

        # Density plot for electrons
        for ireg in 1:numberOfRegions
            subg = subgrid(grid, [ireg])
            ncc_n = get_density(solution_eq, ireg, ctsys, 1)  # icc = 1 for electrons

            scalarplot!(
                vis[2, 1],
                subg,
                1.0e-6 .* ncc_n;
                clear = false,
                title = "Electron density [1/cm^3]",
                xlabel = "length [m]",
                ylabel = "height [m]",
                yscale = :log
            )
        end

        # Density plot for holes
        for ireg in 1:numberOfRegions
            subg = subgrid(grid, [ireg])
            ncc_h = get_density(solution_eq, ireg, ctsys, 2)  # icc = 2 for holes

            scalarplot!(
                vis[2, 2],
                subg,
                1.0e-6 .* ncc_h;
                clear = false,
                title = "Hole density [1/cm^3]",
                xlabel = "length [m]",
                ylabel = "height [m]",
                yscale = :log
            )
        end
        reveal(vis)
        =#
    end

    #=
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    biasValues_gate = range(0.0, stop = 5.0, length = 4) # 4 Steps in Tesca
    biasValues_drain = range(0.0, stop = 5.0, length = 12) # 12 steps in Tesca
    #IV = zeros(0)

    inival = copy(solution_eq)
    solution = copy(solution_eq)

    # Bias Loop
    for Δu_gate in biasValues_gate

        println("bias value at gate: Δu = ", Δu_gate, " V")

        set_contact!(ctsys, bregion_gate, Δu = Δu_gate)
        set_contact!(ctsys, bregion_source, Δu = 0.0)
        set_contact!(ctsys, bregion_drain, Δu = 0.0)
        set_contact!(ctsys, bregion_bulk, Δu = 0.0)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        #current = get_current_val(ctsys, solution)
        #push!(IV, abs.(zaus * current)) # zaus=wide of device
    end

    for Δu_drain in biasValues_drain

        println("bias value at drain: Δu = ", Δu_drain, " V")

        set_contact!(ctsys, bregion_gate, Δu = 5.0)
        set_contact!(ctsys, bregion_source, Δu = 0.0)
        set_contact!(ctsys, bregion_drain, Δu = Δu_drain)
        set_contact!(ctsys, bregion_bulk, Δu = 0.0)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        #current = get_current_val(ctsys, solution)
        #push!(IV, abs.(zaus * current)) # zaus=wide of device

    end

    if test == false
        println("*** done\n")
    end

    =#

    return
end # function main

end # module
