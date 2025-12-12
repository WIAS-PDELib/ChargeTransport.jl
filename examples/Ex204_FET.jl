#=
# Field effect transistor.
([source code](@__SOURCE_URL__))

We consider an n-channel Metal-Oxide-Semiconductor (MOS) field effect transistor.
The material is silicon.
=#

module Ex204_FET

using ChargeTransport
using VoronoiFVM # for IV Curve

using ExtendableGrids
using PythonPlot

#import AMGCLWrap  # Algebraic multigrid solver
#using AlgebraicMultigrid  # Native Julia AMG

function main(; plotting = true, Plotter = PythonPlot, test = false)

    if plotting
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
    Nd_drain = 1.0e19 / cm^3 # eventuell 19
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
    x2 = -0.8 * μm              # beginning gate contact
    x3 = 0.8 * μm               # end gate contact
    x4 = 1.0 * μm               # beginning right contact (source)
    x5 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y_test = -1.2 * μm
    y1 = -0.2 * μm              # beginning n-channel (gate region)
    y2 = 0.0 * μm               # top of the device

    # Refinement x-direction
    X1 = geomspace(x0, x1, 2.0e-7, 3.0e-8)
    X2 = collect(range(x1, x2, length = 5))
    X3 = collect(range(x2, x3, length = 25))
    X4 = collect(range(x3, x4, length = 5))
    X5 = geomspace(x4, x5, 3.0e-8, 2.0e-7)
    X_temp = glue(X1, X2)
    X_temp = glue(X_temp, X3)
    X_temp = glue(X_temp, X4)
    X = glue(X_temp, X5)

    # Refinement y-direction
    Y1 = collect(y0:(0.2 * μm):y_test)
    Y2 = collect(y_test:(0.1 * μm):y1)
    Y3 = geomspace(y1, y2, 4.0e-8, 1.0e-8)
    Y12 = glue(Y1, Y2)
    Y = glue(Y12, Y3)

    grid = simplexgrid(X, Y)

    # cell regions
    cellmask!(grid, [x2, y1], [x3, y2], region_gate)
    cellmask!(grid, [x0, y1], [-0.8 * μm, y2], region_drain)
    cellmask!(grid, [0.8 * μm, y1], [x5, y2], region_source)
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
        gridplot(grid, Plotter = Plotter, title = "Grid")
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

        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg] = Auger
        params.recombinationAuger[iphip, ireg] = Auger
    end

    for ibreg in 1:grid[NumBFaceRegions] #boundary region data
        params.dielectricConstantOxide[ibreg] = εr_ox * ε_0
        params.thicknessOxide[ibreg] = thickness_ox
        params.surfacechargeDensityGate[ibreg] = qss
        params.additionalVoltageGate[ibreg] = U_add_gate
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
    data.boundaryType[bregion_gate] = GateContact # OhmicContact
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

    control = ChargeTransport.SolverControl()
    control.verbose = true
    control.maxiters = 70 # bei 69 Convergence Error im 5ten Schritt Equilibrium, ohne den damp_growth
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.damp_initial = 0.5
    control.max_round = 3
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

    if plotting
        pygui(true) # this opens plots in a different window (at least on my computer, windows ;))

        ################################################################################
        # Surface plot equlibrium with PythonPlot
        ################################################################################

        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution_eq[ipsi, :])
        Plotter.title("Electrostatic potential (Equilibrium)")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        display(gcf()) # I need this to open the plots (might also be an windows problem)

        ################################################################################
        # Density plots equilibrium with PythonPlot
        ################################################################################

        function tridata(grid::ExtendableGrid)
            coord = grid[Coordinates]
            cellnodes = Matrix(grid[CellNodes])
            return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
        end

        vmin = 1.0e15
        vmax = 1.0e25

        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution_eq[iphin, :] .- solution_eq[ipsi, :]) .+ Ec) ./ (k_B * T))
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution_eq[iphip, :] .- solution_eq[ipsi, :]) .+ Ev) ./ (k_B * T))

        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(nn...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        Plotter.xlabel(" \$x\$ [nm]", fontsize = 12)
        Plotter.ylabel(" \$y\$ [nm]", fontsize = 12)
        Plotter.title("electron density")
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        ax = Plotter.gca()
        ax.set_aspect("auto")
        Plotter.tight_layout()
        display(gcf())

        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(np...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        Plotter.xlabel(" \$x\$ [nm]", fontsize = 12)
        Plotter.ylabel(" \$y\$ [nm]", fontsize = 12)
        Plotter.title("hole density")
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        Plotter.tight_layout()
        display(gcf())
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    biasValues_gate = range(0.0, stop = 5.0, length = 4) # 4 Steps in Tesca
    biasValues_drain = range(0.0, stop = 5.0, length = 43) # 12 steps in Tesca, bis length = 42 bricht es nach zwei Schritten ab

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
    for Δu_drain in biasValues_drain

        println("bias value at drain: Δu = ", Δu_drain, " V")

        #set_contact!(ctsys, bregion_gate, Δu = 5.0)
        set_contact!(ctsys, bregion_drain, Δu = Δu_drain)

        solution = ChargeTransport.solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get I-V data
        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

        # 2 drain, 3 source
        tf = VoronoiFVM.testfunction(factory, [2], [3])

        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)

        current = 0.0
        # no displacement as we have steady state, this way last one is taken out as it corresponds to electric potential
        for ii in 1:(length(IEdge) - 1)
            current = current + IEdge[ii]
        end

        push!(IV_drain, abs.(zaus * current)) # zaus=wide of device, total current

    end

    if test == false
        println("*** done\n")
    end

    if plotting

        ################################################################################
        # Solution with PythonPlot
        ################################################################################
        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution[ipsi, :])
        Plotter.title("Electrostatic potential")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        display(gcf())

        ################################################################################
        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution[iphin, :])
        Plotter.title("Quasi Fermi potential for electrons")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        display(gcf())

        ################################################################################
        Plotter.figure()
        Plotter.surf(XX[:], YY[:], solution[iphip, :])
        Plotter.title("Quasi Fermi potential for holes")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        display(gcf())

        ################################################################################
        # Density plots with PythonPlot
        ################################################################################
        function tridata(grid::ExtendableGrid)
            coord = grid[Coordinates]
            cellnodes = Matrix(grid[CellNodes])
            return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
        end

        vmin = 1.0e9
        vmax = 1.0e28

        nn = Nc .* exp.(params.chargeNumbers[iphin] * (q * (solution[iphin, :] .- solution[ipsi, :]) .+ Ec) ./ (k_B * T))
        np = Nv .* exp.(params.chargeNumbers[iphip] * (q * (solution[iphip, :] .- solution[ipsi, :]) .+ Ev) ./ (k_B * T))

        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(nn...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        Plotter.xlabel(" \$x\$ [nm]", fontsize = 12)
        Plotter.ylabel(" \$y\$ [nm]", fontsize = 12)
        Plotter.title("electron density")
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        Plotter.tight_layout()
        display(gcf())

        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(np...), norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
        Plotter.xlabel(" \$x\$ [nm]", fontsize = 12)
        Plotter.ylabel(" \$y\$ [nm]", fontsize = 12)
        Plotter.title("hole density")
        Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
        Plotter.tight_layout()
        display(gcf())

        ################################################################################
        # IV Curve with PythonPlot
        ################################################################################
        Plotter.figure()
        Plotter.semilogy(biasValues_drain, IV_drain, marker = "o", markersize = 3)
        Plotter.grid()
        Plotter.title("IV Curve")
        Plotter.xlabel("bias [V]")
        Plotter.ylabel("total current [A]")
        Plotter.tight_layout()
        display(gcf())

    end

    return
end # function main

end # module
