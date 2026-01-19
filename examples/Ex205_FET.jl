#=
# Field effect transistor.
([source code](@__SOURCE_URL__))

We consider an n-channel Metal-Oxide-Semiconductor (MOS) field effect transistor.
The material is silicon.
=#

module Ex205_FET

using ChargeTransport
using VoronoiFVM          # for IV Curve, Current Density Plot
using GridVisualize       # for Current Density Plot
using ExtendableGrids

# This function is used for the plotting lateron
function tridata(grid::ExtendableGrid)
    coord = grid[Coordinates]
    cellnodes = Matrix(grid[CellNodes])
    return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
end

function main(; Plotter = nothing, test = false)

    if Plotter !== nothing && nameof(Plotter) !== :PythonPlot
        error("Plotting in Ex205_FET is only possible for Plotter = PythonPlot")
    end

    if Plotter !== nothing
        Plotter.close("all")
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
    ipsi = 3                             # electrostatic potential
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

    zaus = 15.0e-4 * cm                                 # depth of the device [cm], see mosS.dio file
    thickness_ox = 0.044e-4 * cm                        # oxide thickness on gate [cm], see mosS.dio file

    ########## physical values of Si at room temperature ##########
    Ec = 0.562 * eV                                     # conduction band-edge energy, Tesca Manual S. 83
    Ev = -0.562 * eV                                    # valence band-edge energy, Tesca Manual S. 83
    Nc = 2.86e19 / (cm^3)                               # conduction band density of states, Tesca Manual S.85
    Nv = 3.1e19 / (cm^3)                                # valence band densitiy of states, Tesca Manual S.85

    mun = 1030.0 * (cm^2) / (V * s)                     # electron mobility, ioffe \leq 1400, Tesca Manual S.92: 1030
    mup = 495.0 * (cm^2) / (V * s)                      # hole mobility \ioffe \leq 450, Tesaca Manual S. 92: 495
    εr = 11.67 * 1.0                                    # relative dielectric permittivity of Si, Tesca Manual S.21: 11.67
    εr_ox = 3.8 * 1.0                                   # relative dielectric permittivity of SiO2 at gate, Tesca Manual S. 21. 3.8
    T = 300.0 * K                                       # room temperature

    # Recombination parameters (from Tesca)
    # Radiative = 0 # Tesca Manual S.127
    Auger_n = 2.8e-31                         # Tesca Manual S. 127: 2.8d-31
    Auger_p = 9.9e-32                         # Tesca Manual S. 127: 9.9d-32
    SRH_TrapDensity_n = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10 1/cm³
    SRH_TrapDensity_p = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10 1/cm³
    SRH_LifeTime_n = 2.0e-4 * s               # Tesca Manual S. 128, tau_n0 = 2e-4 s
    SRH_LifeTime_p = 2.0e-6 * s               # Tesca Manual S. 128, trau_p0 = 2e-6 s
    SRH_Velocity_n = 5.0 * (cm / s)           # Tesca Manual S. 129, VREN = 5.d0 cm/s
    SRH_Velocity_p = 5.0 * (cm / s)           # Tesca Manual S. 129, VREP = 5.d0 cm/s

    # Doping (Simplified, one doping per region, compare gridplot from Tesca)
    Na_gate = 1.0e16 / cm^3
    Nd_drain = 1.0e20 / cm^3
    Nd_source = 1.0e20 / cm^3
    Na_bulk = 1.0e16 / cm^3

    # Voltage information
    U_add_gate = 0.55 * V                   # Contact voltage on gate, see mosS.dio file
    qss = 6.0e10 / (cm^2)                   # Surface charge density on gate, see mosS.dio file

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # Gridpoints x-direction
    x0 = -1.5 * μm              # beginning left contact (drain)
    x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = -0.8 * μm              # end n-doped zone drain
    x4 = 0.8 * μm               # beginning n-doped zone source
    x5 = 0.9 * μm               # end gate contact
    x6 = 1.0 * μm               # beginning right contact (source)
    x7 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y1 = -1.0 * μm              # middle of the device (bulk)
    y2 = -0.2 * μm              # beginning n-channel (gate region)
    y3 = 0.0 * μm               # top of the device

    # Simple refinement in x-direction
    X = collect(range(x0, x7, length = 241))

    # Refinement y-direction
    Y1 = collect(range(y0, y1, length = 8))
    Y2 = collect(range(y1, y2, length = 21))
    Y3 = geomspace(y2, y3, 4.0e-8, 1.0e-10)
    Y12 = glue(Y1, Y2)
    Y = glue(Y12, Y3)

    grid = simplexgrid(X, Y)

    # cell regions
    cellmask!(grid, [x3, y2], [x4, y3], region_gate)
    cellmask!(grid, [x0, y2], [x3, y3], region_drain)
    cellmask!(grid, [x4, y2], [x7, y3], region_source)
    cellmask!(grid, [x0, y0], [x7, y2], region_bulk)

    # boundary regions
    bfacemask!(grid, [x2, y3], [x5, y3], bregion_gate)
    bfacemask!(grid, [x0, y3], [x1, y3], bregion_drain)
    bfacemask!(grid, [x6, y3], [x7, y3], bregion_source)
    bfacemask!(grid, [x0, y0], [x7, y0], bregion_bulk)
    bfacemask!(grid, [x0, y0], [x0, y3], 0)
    bfacemask!(grid, [x7, y0], [x7, y3], 0)
    bfacemask!(grid, [x1, y3], [x2, y3], 0)
    bfacemask!(grid, [x5, y3], [x6, y3], 0)

    if Plotter !== nothing
        Plotter.pygui(true) # Plots as Pop-Ups
        GridVisualize.gridplot(grid, Plotter = Plotter, resolution = (600, 450), legend = :none, title = "Grid", xlabel = "Width [m]", ylabel = "Height [m]")
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
    #paramsnodal = ParamsNodal(grid, numberOfCarriers)

    params.temperature = T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ireg in 1:grid[NumCellRegions] # region data
        params.dielectricConstant[ireg] = εr * ε_0

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = mun
        params.mobility[iphip, ireg] = mup

        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime_n
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime_p
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity_n
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity_p
        params.recombinationAuger[iphin, ireg] = Auger_n
        params.recombinationAuger[iphip, ireg] = Auger_p
    end

    params.recombinationSRHvelocity[iphin, bregion_gate] = SRH_Velocity_n # Surface Recombination just at gate
    params.recombinationSRHvelocity[iphip, bregion_gate] = SRH_Velocity_p

    for ibreg in 1:grid[NumBFaceRegions] #boundary region data
        params.dielectricConstantOxide[ibreg] = εr_ox * ε_0
        params.thicknessOxide[ibreg] = thickness_ox
        params.surfacechargeDensityGate[ibreg] = qss
        params.additionalVoltageGate[ibreg] = U_add_gate
    end

    # doping                                                       #     n    p
    params.doping[iphin, region_drain] = Nd_drain                  #   [Nd  0.0;
    params.doping[iphin, region_source] = Nd_source                #    Nd  0.0;
    params.doping[iphip, region_gate] = Na_gate                    #    0.0  Na;
    params.doping[iphip, region_bulk] = Na_bulk                    #    0.0  Na]

    # Initialize Data instance
    data = Data(grid, numberOfCarriers, constants = constants)

    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann # Tesca Manual S.84 IFERMI = 0

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true,
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

    if Plotter !== nothing
        ################################################################################
        # Doping Profile with PythonPlot
        ################################################################################
        cellregions = grid[CellRegions]
        cellValue = zeros(length(cellregions))

        for i in eachindex(cellregions)
            # determine doping value in cell and number of cell nodes
            # ToDo: - or +? or consider chargeNumbers
            cellValue[i] = (params.doping[1, cellregions[i]] - params.doping[2, cellregions[i]]) * 1.0e-6
        end

        Plotter.figure()
        m = Plotter.tripcolor(
            tridata(grid)...,
            cellValue;
            shading = "flat",
            alpha = 0.75,
            rasterized = true
        )
        # statt colorbar:
        vals = sort(unique(cellValue))
        cmap = m.get_cmap()
        norm = m.norm
        Patch = Plotter.matplotlib.patches.Patch

        handles = [Patch(facecolor = cmap(norm(v)), edgecolor = "k") for v in vals]
        labels = ["$(v)" for v in vals]

        #Plotter.colorbar(orientation = "vertical", label = " Doping [\$\\mathrm{cm}^{-3}\$]", extend = "both")
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Doping Profile", fontsize = 16)

        ax = Plotter.gca()
        ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
        ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
        ax.tick_params(axis = "both", which = "major", labelsize = 10)
        ax.xaxis.get_offset_text().set_fontsize(10)
        ax.yaxis.get_offset_text().set_fontsize(10)
        ax.legend(handles, labels, title = "Doping values [cm⁻³]", loc = "lower right", fontsize = 10, title_fontsize = 8)

        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)
    end

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
    control.verbose = false
    control.maxiters = 70 # bei 69 Convergence Error im 5ten Schritt Equilibrium, ohne den damp_growth
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.damp_initial = 0.5
    control.max_round = 3

    VoronoiFVM.check_allocs!(false) # no allocation warnings, ToDo: Check allocations in the end.
    #control.damp_growth = 1.05 # >= 1 # wenn das 1.21 auch Convergence error beim selben Schritt, wenn kleiner dann weniger iterations (bei 1.05: 64 Iterations)
    #control.method_linear = AMGCLWrap.AMGSolverAlgorithm(blocksize = 2)

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

    if Plotter !== nothing
        ################################################################################
        # Surface plot equilibrium with PythonPlot
        ################################################################################

        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        fig = Plotter.figure()
        ax = fig.add_subplot(111, projection = "3d")
        ax.plot_trisurf(XX[:], YY[:], solution_eq[ipsi, :])

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

        ################################################################################
        # Density plots equilibrium with PythonPlot
        ################################################################################

        # electron density
        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution_eq[iphin, :] .- solution_eq[ipsi, :]) .+ Ec) ./ (k_B * T)) # this is in m, to see this print Nc
        nn_cm = nn * 1.0e-6 # Convert to cm, so it is better comparable to Tesca

        # plot electron density
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., nn_cm,
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e5, vmax = 1.0e21),
            shading = "gouraud", # shading = "flat"
            rasterized = true
        )
        cbar = Plotter.colorbar(orientation = "vertical", extend = "both")

        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Electron Density (Equilibrium)", fontsize = 16)

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

        # hole density
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution_eq[iphip, :] .- solution_eq[ipsi, :]) .+ Ev) ./ (k_B * T))
        np_cm = np * 1.0e-6 # Convert to cm, so it is better comparable to Tesca

        # plot hole density
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
        Plotter.title("Hole Density (Equilibrium)", fontsize = 16)

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
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    biasValues_gate = range(0.0, stop = 5.0, length = 4) # 4 Steps in Tesca
    biasValues_drain = range(0.0, stop = 5.0, length = 43) # 12 steps in Tesca, bis length = 42 bricht es nach zwei Schritten ab, war 43

    IV_drain = zeros(0)

    inival = copy(solution_eq)
    solution = copy(solution_eq)

    # Bias Loop Gate
    for Δu_gate in biasValues_gate

        println("bias value at gate: Δu = ", Δu_gate, " V")

        set_contact!(ctsys, bregion_gate, Δu = Δu_gate)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

    end

    # Bias Loop Drain
    # Needed to calculate IV values, placed before loop to avoid allocations
    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
    tf = VoronoiFVM.testfunction(factory, [bregion_drain], [bregion_source, bregion_bulk])

    for Δu_drain in biasValues_drain

        println("bias value at drain: Δu = ", Δu_drain, " V")

        set_contact!(ctsys, bregion_drain, Δu = Δu_drain)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)

        current = IEdge[iphin] + IEdge[iphip]

        push!(IV_drain, abs.(zaus * current)) # zaus=wide of device, total current, zaus ist hier in m??

    end

    if test == false
        println("*** done\n")
    end

    if Plotter !== nothing
        ################################################################################
        # Solution with PythonPlot
        ################################################################################
        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        fig = Plotter.figure()
        ax = fig.add_subplot(111, projection = "3d")
        ax.plot_trisurf(XX[:], YY[:], solution[ipsi, :])

        Plotter.title("Electrostatic Potential", fontsize = 16)
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

        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)

        ################################################################################
        # Density plots with PythonPlot
        ################################################################################
        # electron density
        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution[iphin, :] .- solution[ipsi, :]) .+ Ec) ./ (k_B * T))
        nn_cm = nn * 1.0e-6 # Convert to cm, so it is better comparable to Tesca

        # electron density plot
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., vcat(nn_cm...),
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e4, vmax = 1.0e21),
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

        # hole density
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution[iphip, :] .- solution[ipsi, :]) .+ Ev) ./ (k_B * T))
        np_cm = np * 1.0e-6 # Convert to cm, so it is better comparable to Tesca
        # hole density plot
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., vcat(np_cm...),
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0, vmax = 1.0e16),
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

        ################################################################################
        # Current Density with PythonPlot
        ################################################################################
        nf = -VoronoiFVM.nodeflux(ctsys.fvmsys, solution)

        j = nf[:, iphin, :] ./ (1.0e4 * zaus) # damit A/ m^2 (zaus ist ja hier in m)
        jx = j[1, :]; jy = j[2, :]
        jabs = sqrt.(jx .^ 2 .+ jy .^ 2)

        # For the arrows
        coord = grid[Coordinates]
        x = coord[1, :]; y = coord[2, :]

        L = 0.03 * (maximum(y) - minimum(y))     # length arrows
        eps = 1.0e-30
        ux = jx ./ (jabs .+ eps);  uy = jy ./ (jabs .+ eps)
        idx = 1:10:length(x)

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
            x[idx], y[idx], (L * ux)[idx], (L * uy)[idx],
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

        ################################################################################
        # IV Curve with PythonPlot
        ################################################################################
        Plotter.figure()
        Plotter.plot(biasValues_drain, IV_drain, marker = "o", markersize = 3) # für Logarithmische Skala: semilogy
        Plotter.grid()

        Plotter.title("IV Curve", fontsize = 16)
        Plotter.xlabel("bias [V]", fontsize = 12)
        Plotter.ylabel("total current [A]", fontsize = 12)

        ax = Plotter.gca()
        ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(1.0))
        ax.tick_params(axis = "both", which = "major", labelsize = 10)

        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)

    end

    return
end # function main

end # module
