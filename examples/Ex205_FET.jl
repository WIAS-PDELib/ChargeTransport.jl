#=
# Field effect transistor.
([source code](@__SOURCE_URL__))

We consider an n-channel Metal-Oxide-Semiconductor (MOS) field effect transistor.
The material is silicon.
=#

module Ex205_FET

using ChargeTransport
using VoronoiFVM # for IV Curve
using LinearAlgebra # for current density
using GridVisualize # for current density

using ExtendableGrids

# This function is used for the plotting lateron
function tridata(grid::ExtendableGrid)
    coord = grid[Coordinates]
    cellnodes = Matrix(grid[CellNodes])
    return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
end

# thx https://discourse.julialang.org/t/meshgrid-function-in-julia/48679/4?u=j-fu
function meshgrid(rc)
    nx = length(rc[1])
    ny = length(rc[2])
    xout = zeros(ny, nx)
    yout = zeros(ny, nx)
    for jx in 1:nx
        for ix in 1:ny
            xout[ix, jx] = rc[1][jx]
            yout[ix, jx] = rc[2][ix]
        end
    end
    return xout, yout
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

    zaus = 15.0e-4 * cm                                 # depth of the device  [cm]
    thickness_ox = 0.044e-4 * cm                        # oxide thickness on gate [cm]

    ########## physical values of Si at room temperature ##########

    Ec = 0.562 * eV                                     # conduction band-edge energy, Tesca Manal S. 83
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
    Auger_n = 2.8e-31 #* cm^6 / s              # Tesca Manual S. 127: 2.8d-31
    Auger_p = 9.9e-32 #* cm^6 / s              # Tesca Manual S. 127: 9.9d-32
    SRH_TrapDensity_n = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10
    SRH_TrapDensity_p = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10
    SRH_LifeTime_n = 2.0e-4 * s               # Tesca Manual S. 128, tau_n0 = 2e-4 s
    SRH_LifeTime_p = 2.0e-6 * s               # Tesca Manual S. 128, trau_p0 = 2e-6 s

    # Doping (Vereinfacht, nur eine Dotierung pro Region)
    Na_gate = 1.0e16 / cm^3 # 1.0e16 = 1.0e22 * 1.0e-6
    Nd_drain = 1.0e20 / cm^3 # 1.0e20 = 1.0e26 *1.0e-6
    Nd_source = 1.0e20 / cm^3 # 1.0e20 = 1.0e26 * 1.0e-6
    Na_bulk = 1.0e15 / cm^3 # 1.0e15 = 1.0e21 * 1.0e-6

    # Voltage information
    U_add_gate = 0.55 * V                   # contact voltage on gate
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
    x3 = -0.8 * μm              # end n-doped zone drain
    x4 = 0.8 * μm               # beginning n-doped zone source
    x5 = 0.9 * μm               # end gate contact
    x6 = 1.0 * μm               # beginning right contact (source)
    x7 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y1 = -1.0 * μm              # middle of the device (bulk)
    ytest = -0.4 * μm # test um das Gitter gleichmäßiger zu machen
    y2 = -0.2 * μm              # beginning n-channel (gate region)
    y3 = 0.0 * μm               # top of the device

    # Refinement x-direction
    # X1 = geomspace(x0, x1, 2.0e-7, 3.0e-8)
    # X2 = collect(range(x1, x2, length = 4))
    # X3 = collect(range(x2, x3, length = 4))
    # X4 = collect(range(x3, x4, length = 20))
    # X5 = collect(range(x4, x5, length = 4))
    # X6 = collect(range(x5, x6, length = 4))
    # X7 = geomspace(x6, x7, 3.0e-8, 2.0e-7)
    # X_temp = glue(X1, X2)
    # X_temp = glue(X_temp, X3)
    # X_temp = glue(X_temp, X4)
    # X_temp = glue(X_temp, X5)
    # X_temp = glue(X_temp, X6)
    # X = glue(X_temp, X7)

    # oder vereinfacht, IV Curve bleibt gleich
    X = collect(range(x0, x7, length = 31))

    # Refinement y-direction
    Y1 = collect(range(y0, y1, length = 8))
    Y2 = collect(range(y1, y2, length = 11)) # bei length = 8 wäre ein Abstand: 8.57e-2, bei length = 7: 0.1, *10⁻6 ist dann Meter...
    Y3 = geomspace(y2, y3, 8.0e-8, 1.0e-10) # Achtung: Grenzen auch in μm
    #Y2 = collect(range(y1, ytest, length = 7)) # bei length = 8 wäre ein Abstand: 8.57e-2, bei length = 7: 0.1, *10⁻6 ist dann Meter...
    #Y3 = geomspace(ytest, y3, 1.0e-7, 1.0e-9) # Achtung: Grenzen auch in μm
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
        Plotter.pygui(true) # aktiviert/ deaktiviert Abbildungen als Pop-up
        ChargeTransport.gridplot(grid, Plotter = Plotter, title = "Grid", xlabel = "Width [m]", ylabel = "Height [m]")
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
    data.F .= Boltzmann

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
        # Surface plot equlibrium with PythonPlot
        ################################################################################

        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution_eq[ipsi, :])
        Plotter.title("Electrostatic Potential (Equilibrium)", fontsize = 12)
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.zlabel("potential [V]", fontsize = 12)
        ax = Plotter.gca()
        ax.tick_params(axis = "x", labelsize = 8)
        ax.tick_params(axis = "y", labelsize = 8)
        ax.tick_params(axis = "z", labelsize = 8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        ax.zaxis.get_offset_text().set_fontsize(8)
        current_figure = Plotter.gcf()
        display(current_figure)

        ################################################################################
        # Density plots equilibrium with PythonPlot
        ################################################################################

        # electron density
        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution_eq[iphin, :] .- solution_eq[ipsi, :]) .+ Ec) ./ (k_B * T))

        # plot electron density
        # Stetige FEM Approximation der Dichte
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., nn,
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e11, vmax = 1.0e27),
            shading = "gouraud", # shading = "flat"
            rasterized = true
        )
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Electron Density", fontsize = 12)
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        ax = Plotter.gca()
        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)

        # hole density
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution_eq[iphip, :] .- solution_eq[ipsi, :]) .+ Ev) ./ (k_B * T))

        # plot hole density
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., np,
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e9, vmax = 1.0e22),
            shading = "gouraud",
            rasterized = true
        )
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Hole Density", fontsize = 12)
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
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
    tf = VoronoiFVM.testfunction(factory, [2], [3]) # 2 drain, 3 source

    for Δu_drain in biasValues_drain

        println("bias value at drain: Δu = ", Δu_drain, " V")

        set_contact!(ctsys, bregion_drain, Δu = Δu_drain)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)

        current = sum(@view IEdge[1:(end - 1)])

        # current = 0.0
        # # no displacement as we have steady state, this way last one is taken out as it corresponds to electric potential
        # for ii in 1:(length(IEdge) - 1)
        #     current = current + IEdge[ii]
        # end

        push!(IV_drain, abs.(zaus * current)) # zaus=wide of device, total current

    end

    if test == false
        println("*** done\n")
    end

    if Plotter !== nothing

        ################################################################################
        # Solution with PythonPlot
        ################################################################################
        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution[ipsi, :])
        Plotter.title("Electrostatic Potential", fontsize = 12)
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.zlabel("potential [V]", fontsize = 12)
        ax = Plotter.gca()
        ax.tick_params(axis = "x", labelsize = 8)
        ax.tick_params(axis = "y", labelsize = 8)
        ax.tick_params(axis = "z", labelsize = 8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        ax.zaxis.get_offset_text().set_fontsize(8)
        current_figure = Plotter.gcf()
        display(current_figure)


        ################################################################################
        # Density plots with PythonPlot
        ################################################################################

        # electron density
        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution[iphin, :] .- solution[ipsi, :]) .+ Ec) ./ (k_B * T))

        # electron density plot
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., vcat(nn...),
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e11, vmax = 1.0e27),
            shading = "gouraud",
            rasterized = true
        )
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Electron Density", fontsize = 12)
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)


        # hole density
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution[iphip, :] .- solution[ipsi, :]) .+ Ev) ./ (k_B * T))

        # hole density plot
        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., vcat(np...),
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e5, vmax = 1.0e22),
            shading = "gouraud",
            rasterized = true
        )
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Hole Density", fontsize = 12)
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)

        ################################################################################
        # IV Curve with PythonPlot
        ################################################################################
        Plotter.figure()
        Plotter.plot(biasValues_drain, IV_drain, marker = "o", markersize = 3) # für Logarithmische Skala: semilogy
        Plotter.grid()
        Plotter.title("IV Curve", fontsize = 12)
        Plotter.xlabel("bias [V]", fontsize = 12)
        Plotter.ylabel("total current [A]", fontsize = 12)
        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)


        ################################################################################
        # Current Density with PythonPlot
        ################################################################################
        nodes = 1:num_nodes(grid) # weil ich habe ein SimplexGrid ohne Subgrids
        nft = VoronoiFVM.nodeflux(ctsys.fvmsys, solution)
        jPsi = nft[:, ipsi, nodes] ./ (ε_0 * εr)
        jPsiAbs = norm.(eachcol(jPsi))
        #jPsiAbs = vec(sqrt.(sum(abs2, jPsi; dims = 1))) # Alternative, noch nict getestet??

        # für die Pfeile?
        #rasterpoints = 5
        #rcE, rvE = vectorsample(grid, jPsi, rasterpoints = rasterpoints) # aus GridVisualize
        #scaleLines = 1.0e-10
        #X_mesh, Y_mesh = meshgrid(rcE)

        Plotter.figure()
        Plotter.tripcolor(
            tridata(grid)..., vcat(jPsiAbs...),
            norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e-2, vmax = 1.0e10),
            shading = "gouraud",
            rasterized = true
        )
        Plotter.xlabel("Width [m]", fontsize = 12)
        Plotter.ylabel("Height [m]", fontsize = 12)
        Plotter.title("Current Density", fontsize = 12)
        Plotter.colorbar(orientation = "vertical", label = "El. field strength [\$\\mathrm{V} \\mathrm{m}^{-1}\$]", extend = "both")
        #cbar = colorbar(orientation = "vertical", label = "El. field strength [\$\\mathrm{V} \\mathrm{m}^{-1}\$]", extend = "both")
        #Plotter.streamplot(
        #     X_mesh, Y_mesh,
        #     scaleLines .* rvE[1, :, :, 1]', scaleLines .* rvE[2, :, :, 1]',
        #     color = "black", density = (0.6, 0.8)
        # )
        Plotter.tight_layout()
        current_figure = Plotter.gcf()
        display(current_figure)

        # Hier ein Test mit GridVisualize
        # ispec = ipsi
        # vis = GridVisualizer(Plotter = Plotter)
        # scalarplot!(vis, grid, solution[ispec, :], clear = true, colormap = :summer)
        # vectorplot!(vis, grid, nf[:, ispec, :], clear = false)
        # reveal(vis)

    end

    return
end # function main

end # module
