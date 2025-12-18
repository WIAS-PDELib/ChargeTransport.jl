#=
# Simple laser structure with 5 layers.
([source code](@__SOURCE_URL__))

Simulating a simple laser structure with 5 layers.
The layers are defined by their material properties and thicknesses.
The simulation will solve the charge transport equations across the layers,
taking into account the stimulated recombination in the active region of the laser structure.
=#

module Ex204_Laser_simple

using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

###########################################################################
numberOfColoumns = Dict(
    "ref1" => [2, 4],
    "ref2" => [4, 8],
    "ref3" => [8, 16],
    "ref4" => [16, 32],
    "ref5" => [32, 64]
)

numberOfRows = Dict(
    "ref1" => [4, 8, 2, 8, 4],
    "ref2" => [8, 16, 4, 16, 8],
    "ref3" => [16, 32, 8, 32, 16],
    "ref4" => [32, 64, 16, 64, 32],
    "ref5" => [64, 128, 32, 128, 64]
)
###################################################################

""" Initializing X and Y coords for the tesca grid"""
function tesca_grid(; refinement = 1, showplot = false, airbox = false)

    @local_unitfactors μm

    ncol = numberOfColoumns["ref$(refinement)"]
    nrow = numberOfRows["ref$(refinement)"]

    widths_columns = [1.0  10.0] * μm
    heights_rows = [1.0  0.5  0.05  0.5  1.0] * μm

    coord_x1 = collect(range(0.0, widths_columns[1], length = ncol[1] + 1))
    coord_x2 = collect(range(widths_columns[1], sum(widths_columns[1:2]), length = ncol[2] + 1))
    X = glue(coord_x1, coord_x2)
    X = glue(reverse(-X), X)

    coord_y1 = collect(range(0.0, heights_rows[1], length = nrow[1] + 1))
    coord_y2 = collect(range(heights_rows[1], sum(heights_rows[1:2]), length = nrow[2] + 1))
    coord_y3 = collect(range(sum(heights_rows[1:2]), sum(heights_rows[1:3]), length = nrow[3] + 1))
    coord_y4 = collect(range(sum(heights_rows[1:3]), sum(heights_rows[1:4]), length = nrow[4] + 1))
    coord_y5 = collect(range(sum(heights_rows[1:4]), sum(heights_rows[1:5]), length = nrow[5] + 1))
    Y = glue(glue(glue(glue(coord_y1, coord_y2), coord_y3), coord_y4), coord_y5)

    grid = simplexgrid(X, Y)

    cellmask!(grid, [X[1], 0.0], [sum(widths_columns), heights_rows[1]], 1)
    cellmask!(grid, [X[1], heights_rows[1]], [sum(widths_columns), sum(heights_rows[1:2])], 2)
    cellmask!(grid, [X[1], sum(heights_rows[1:2])], [sum(widths_columns), sum(heights_rows[1:3])], 3)
    cellmask!(grid, [X[1], sum(heights_rows[1:3])], [sum(widths_columns), sum(heights_rows[1:4])], 4)
    cellmask!(grid, [X[1], sum(heights_rows[1:4])], [sum(widths_columns), sum(heights_rows[1:5])], 5)

    # AIR
    cellmask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [sum(widths_columns), sum(heights_rows[1:5])], 6)
    cellmask!(grid, [X[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows[1:5])], 6)

    bregionDonor1 = 1      # bottom boundary
    bregionAcceptor2 = 2   # top boundary
    bregionNoFlux = 3
    bregionAirBox = 4
    bfacemask!(grid, [-widths_columns[1], sum(heights_rows)], [widths_columns[1], sum(heights_rows)], bregionAcceptor2)
    bfacemask!(grid, [X[1], 0.0], [X[end], 0.0], bregionDonor1)

    bfacemask!(grid, [X[1], 0.0], [X[1], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [X[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [-widths_columns[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows)], bregionNoFlux)

    bfacemask!(grid, [X[end], 0.0], [X[end], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [X[end], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [widths_columns[1], sum(heights_rows)], bregionNoFlux)

    bfacemask!(grid, [X[1], sum(heights_rows[1:4])], [X[1], Y[end]], bregionAirBox)
    bfacemask!(grid, [X[1], Y[end]], [-widths_columns[1], Y[end]], bregionAirBox)
    bfacemask!(grid, [widths_columns[1], Y[end]], [X[end], Y[end]], bregionAirBox)
    bfacemask!(grid, [X[end], sum(heights_rows[1:4])], [X[end], Y[end]], bregionAirBox)

    if airbox == false
        grid = subgrid(grid, [1, 2, 3, 4, 5])
    end

    if showplot == true
        GridVisualize.gridplot(
            grid, Plotter = PyPlot, linewidth = 1, fontsize = 35, size = (1200, 900),
            legend = :best, show = true, aspect = 4, colorbar = false, title = "Device Geometry, values in [m]", xlabel = "x-coordinates", ylabel = "y-coordinates"
        )
    end

    return grid
end

function main(;
        refinement = 1, plotting = false, verbose = false, test = false, unknown_storage = :sparse,
        numberOfEigenvalues = 1,
        parameter_set = Params_Laser_simple # choose the parameter set
    )

    # parameter
    p = parameter_set()

    ################################################################################
    if test == false
        println("Set up grid.")
    end
    ################################################################################

    grid = tesca_grid(refinement = refinement, showplot = plotting, airbox = false)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in data
    data = Data(grid, p.numberOfCarriers, numberOfEigenvalues = numberOfEigenvalues)

    ## Possible choices: Stationary, Transient
    data.modelType = Stationary

    ## The default for electrons and holes is Boltzmann. Here, we set it to a more general statistics function
    data.F .= FermiDiracOneHalfTeSCA

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    data.boundaryType[p.bregionAcceptor2] = OhmicContact  # top boundary Dirichlet condition
    data.boundaryType[p.bregionDonor1] = OhmicContact     # bottom boundary Dirichlet condition
    #                                                     # rest is set to Neumann by default


    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    """ Data from Params_Laser_simple.jl: temperature T, band edge energies EC, EV, effective density of states NC, NV
        mobilities μn, μp, dielectricConstant εs, radiative recombination r0, life times τn, τp,
        Auger recombination coefficients Auger_Cn, Auger_Cp
        doping doping (or vcat(Nd,Na) = doping).
    """
    paramsoptical = ParamsOptical(grid, p.numberOfCarriers, numberOfEigenvalues)
    paramsoptical.laserWavelength = p.λ

    paramsoptical.absorption_0[:] = p.α0
    paramsoptical.gain_0[:] = p.gain0
    paramsoptical.refractiveIndex_0[:] = p.nTilde
    paramsoptical.refractiveIndex_d[:] = p.nTilde_d
    paramsoptical.refractiveIndex_γ[:] = p.γn
    paramsoptical.absorptionFreeCarriers[p.iphin, :] = p.fcnalf
    paramsoptical.absorptionFreeCarriers[p.iphip, :] = p.fcpalf

    paramsoptical.eigenvalues .= 1 + 1 * im   # dummy value for initializing

    data.params = Params(p)
    data.paramsoptical = paramsoptical

    ctsys = System(grid, data, unknown_storage = unknown_storage)

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
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.max_round = 3
    control.damp_initial = 0.8   # < 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    ## calculate equilibrium solution and set as initial guess
    psi0Vector = electroNeutralSolution(ctsys)

    inival = unknowns(ctsys)
    inival[1, :] = inival[2, :] .= 0.0
    inival[3, :] = psi0Vector

    solution = equilibrium_solve!(ctsys, inival = inival, control = control, nonlinear_steps = 20.0)
    inival = solution

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    maxBias = p.U[end]  # = 1.81 = topVoltageAcceptor2 # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 40)

    for Δu in biasValues

        if test == false
            println("bias value: Δu = ", Δu, " V")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, p.bregionAcceptor2, Δu = Δu)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

    end # bias loop

    if test == false
        println("*** done\n")
    end

    #########################################################

    ctsys.data.paramsoptical.oldSolution = solution
    currentSolution = solution
    inival = solution

    ############################
    ctsys.data.paramsoptical.eigenvalues = p.λ1
    ctsys.data.paramsoptical.eigenvectors = reshape(p.v1, length(p.v1), 1)  # reshaping because in system it must be a 2D array
    ctsys.data.paramsoptical.power = p.P[2]

    solution = solve(ctsys; inival = inival, control = control)
    ctsys.data.paramsoptical.oldSolution = solution
    currentSolution = solution
    inival = solution

    if plotting
        vis = GridVisualizer(; Plotter = PyPlot, fignumber = 2, resolution = (1200, 900))

        scalarplot!(
            vis, grid, solution[1, :], Plotter = PyPlot, legend = :best, clear = false, title = "Applied voltage Δu = $maxBias V",
            xlabel = "cross section space along \$x=0\$ [m]", ylabel = "potentials [V]", fontsize = 55, linewidth = 5,
            slice = :x => 0, label = "\$ \\varphi_n \$", color = "blue"
        )

        scalarplot!(
            vis, grid, solution[2, :], Plotter = PyPlot, linewidth = 5,
            slice = :x => 0, label = "\$ \\varphi_p \$", clear = false, color = "mediumvioletred"
        )

        scalarplot!(
            vis, grid, solution[3, :], Plotter = PyPlot, linewidth = 5,
            slice = :x => 0, label = "\$ \\psi \$", clear = false, color = "darkorange"
        )
    end


    testval = sum(solution) / length(solution)
    return testval

end # main

function test()
    testval = 0.5122451923673309
    return main(test = true) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
