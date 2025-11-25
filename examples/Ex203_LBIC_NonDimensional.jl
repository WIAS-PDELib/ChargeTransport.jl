#=
# 2D nondimensionalized LBIC measurement.
([source code](@__SOURCE_URL__))

Simulation of LBIC measurement technique with a focused laser in two dimensions, accompanying the manuscript

``Existence of solutions and uniform bounds for the stationary semiconductor equations with generation and ionic carriers``.

=#

module Ex203_LBIC_NonDimensional

using ChargeTransport
using VoronoiFVM
using ExtendableGrids  # grid initializer
using PyCall
# using PyPlot           # solution visualizer

# https://stackoverflow.com/questions/29443369/how-to-make-a-custom-colormap-using-pyplot-not-matplotlib-proper
matcolors = pyimport("matplotlib.colors")

# https://github.com/BIDS/colormap/blob/master/parula.py
cm_data = [
    [0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
    [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
    [0.1959047619, 0.2644571429, 0.7279], [
        0.1707285714, 0.2919380952,
        0.779247619,
    ], [0.1252714286, 0.3242428571, 0.8302714286],
    [0.0591333333, 0.3598333333, 0.8683333333], [
        0.0116952381, 0.3875095238,
        0.8819571429,
    ], [0.0059571429, 0.4086142857, 0.8828428571],
    [0.0165142857, 0.4266, 0.8786333333], [
        0.032852381, 0.4430428571,
        0.8719571429,
    ], [0.0498142857, 0.4585714286, 0.8640571429],
    [0.0629333333, 0.4736904762, 0.8554380952], [
        0.0722666667, 0.4886666667,
        0.8467,
    ], [0.0779428571, 0.5039857143, 0.8383714286],
    [0.079347619, 0.5200238095, 0.8311809524], [
        0.0749428571, 0.5375428571,
        0.8262714286,
    ], [0.0640571429, 0.5569857143, 0.8239571429],
    [0.0487714286, 0.5772238095, 0.8228285714], [
        0.0343428571, 0.5965809524,
        0.819852381,
    ], [0.0265, 0.6137, 0.8135], [
        0.0238904762, 0.6286619048,
        0.8037619048,
    ], [0.0230904762, 0.6417857143, 0.7912666667],
    [0.0227714286, 0.6534857143, 0.7767571429], [
        0.0266619048, 0.6641952381,
        0.7607190476,
    ], [0.0383714286, 0.6742714286, 0.743552381],
    [0.0589714286, 0.6837571429, 0.7253857143],
    [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
    [0.1452714286, 0.7097571429, 0.6646285714], [
        0.1801333333, 0.7176571429,
        0.6424333333,
    ], [0.2178285714, 0.7250428571, 0.6192619048],
    [0.2586428571, 0.7317142857, 0.5954285714], [
        0.3021714286, 0.7376047619,
        0.5711857143,
    ], [0.3481666667, 0.7424333333, 0.5472666667],
    [0.3952571429, 0.7459, 0.5244428571], [
        0.4420095238, 0.7480809524,
        0.5033142857,
    ], [0.4871238095, 0.7490619048, 0.4839761905],
    [0.5300285714, 0.7491142857, 0.4661142857], [
        0.5708571429, 0.7485190476,
        0.4493904762,
    ], [0.609852381, 0.7473142857, 0.4336857143],
    [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
    [0.7184095238, 0.7411333333, 0.3904761905],
    [0.7524857143, 0.7384, 0.3768142857], [
        0.7858428571, 0.7355666667,
        0.3632714286,
    ], [0.8185047619, 0.7327333333, 0.3497904762],
    [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
    [0.9139333333, 0.7257857143, 0.3062761905], [
        0.9449571429, 0.7261142857,
        0.2886428571,
    ], [0.9738952381, 0.7313952381, 0.266647619],
    [0.9937714286, 0.7454571429, 0.240347619], [
        0.9990428571, 0.7653142857,
        0.2164142857,
    ], [0.9955333333, 0.7860571429, 0.196652381],
    [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
    [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
    [0.9588714286, 0.8949, 0.1132428571], [
        0.9598238095, 0.9218333333,
        0.0948380952,
    ], [0.9661, 0.9514428571, 0.0755333333],
    [0.9763, 0.9831, 0.0538],
]

parula_map = matcolors.LinearSegmentedColormap.from_list("parula", cm_data)

# grid information
length_x = 8.0; length_y = 4.0
box_cx = length_x / 2;  box_cy = length_y / 2
box_x = 4.0;  box_y = 2.0

numberOfRegions = 2
bregion1 = 1; bregion2 = 2

## set indices of the quasi Fermi potentials
iphin = 1              # electron quasi Fermi potential
iphip = 2              # hole quasi Fermi potential
ipsi = 3               # electric potential
numberOfCarriers = 2

# We define the physical data.
zn = -1; zp = 1
En = 0.0; Ep = 0.0     # set the energies to 0
Nn = 1.0; Np = 1.0     # set the effective DOS to 1
μn = 1.0; μp = 1.0     # set the mobilities to 1
T = 1.0                # set temperature to 1

## recombination parameters
SRH_TrapDensity = 0.0
SRH_LifeTime = 1.0

# you can set verbose also to true to display some solver information
function main(;
        hmin = 0.1, h = 0.25, hmax = 0.25, # grid information
        Cn = 1.0e1, Cp = 1.0e1,            # doping
        λ = 1.0,                           # Debye length
        G0 = 1.0, x0 = 1.0, y0 = 2.0,      # laser information
        ####################################
        parameterStudy = false,
        ####################################
        Plotter = nothing, plotting = false,
        verbose = false, test = false
    )

    if plotting
        if isnothing(Plotter)
            @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
            plotting = false
        else
            if nameof(Plotter) != :PyPlot
                @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
                plotting = false
            end
        end
    end

    if plotting
        Plotter.rc("font", family = "sans-serif", size = 14)
        Plotter.rc("mathtext", fontset = "dejavusans")
        Plotter.close("all")
    end

    # constants
    # Here, we set the constants to unity_constants.
    # In particular, we set: q = k_B = ε_0 = 1
    # When defining `ChargeTransport.constants" the constants based on CODATA2022 are used.
    constants = ChargeTransport.unity_constants

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    X1 = geomspace(0.0, box_cx - box_x / 2, hmin, hmax)
    X2 = geomspace(box_cx - box_x / 2, box_cx, hmax, hmin)
    X3 = geomspace(box_cx, box_cx + box_x / 2, hmin, hmax)
    X4 = geomspace(box_cx + box_x / 2, length_x, hmax, hmin)
    X = glue(X1, X2); X = glue(X, X3); X = glue(X, X4)
    ##################
    Y1 = collect(0.0:h:(box_cy - box_y / 2))
    Y2 = geomspace(box_cy - box_y / 2, box_cy, hmax, hmin)
    Y3 = geomspace(box_cy, box_cy + box_y / 2, hmin, hmax)
    Y4 = collect((box_cy + box_y / 2):h:length_y)
    Y = glue(Y1, Y2); Y = glue(Y, Y3); Y = glue(Y, Y4)
    ##################
    grid = simplexgrid(X, Y)

    rect!(grid, [box_cx - box_x / 2, box_cy - box_y / 2], [box_cx + box_x / 2, box_cy + box_y / 2], region = 2, bregion = 0, tol = 1.0e-1)

    ## bfacemask! for setting different boundary regions.
    bfacemask!(grid, [0.0, 0.0], [0.0, length_y], bregion1)
    bfacemask!(grid, [length_x, 0.0], [length_x, length_y], bregion2)
    bfacemask!(grid, [0.0, 0.0], [length_x, 0.0], 0)
    bfacemask!(grid, [0.0, length_y], [length_x, length_y], 0)

    if plotting
        gridplot(grid, Plotter = Plotter)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    ## Optical generation
    sigma = 0.5
    G(x, y) = G0 * exp(- ((x - x0)^2 + (y - y0)^2) / (2 * sigma^2))

    genData = zeros(length(grid[Coordinates][1, :]))
    coord = grid[Coordinates]
    subg = subgrid(grid, [1, 2])
    iNode = subg[NodeParents] # to receive array with correct ordering of nodes

    for inode in iNode
        x = coord[1, inode]
        y = coord[2, inode]
        genData[inode] = G(x, y)
    end

    if plotting
        XX = grid[Coordinates][1, :]
        YY = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.tricontourf(XX, YY, genData, levels = 40)
        Plotter.colorbar(orientation = "vertical", label = "Generation \$ G \$ ")
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("\$ y \$", fontsize = 17)
        Plotter.axis("equal")
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # We initialize the Data instance and fill in predefined data.
    data = Data(grid, numberOfCarriers, generationData = genData, constants = constants)

    ## Following variable declares, if we want to solve stationary or transient problem
    data.modelType = Stationary

    ## Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk,
    ## FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F .= FermiDiracOneHalfTeSCA

    ## The desired recombination processes can be chosen here.
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true
    )

    # generation model
    data.generationModel = GenerationUserDefined

    ## flux discretization scheme
    data.fluxApproximation .= ExcessChemicalPotential

    ## Define the unity constants also in the discrete counterpart of the model
    data.constants = constants

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    # Define the Params struct
    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

    params.temperature = T
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp

    for ireg in 1:numberOfRegions # region data

        params.dielectricConstant[ireg] = λ^2 # Debye length, entering model as `DielectricConstant`, i.e., prefactor in the displacement flux

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn
        params.densityOfStates[iphip, ireg] = Np
        params.bandEdgeEnergy[iphin, ireg] = En
        params.bandEdgeEnergy[iphip, ireg] = Ep
        params.mobility[iphin, ireg] = μn
        params.mobility[iphip, ireg] = μp

        ## recombination parameters
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity

    end

    ## doping
    params.doping[iphin, 1] = Cn
    params.doping[iphip, 2] = Cp

    # Region dependent params is now a substruct of data which is again a substruct of the
    # system and will be parsed in next step.
    data.params = params

    # In the last step, we initialize our system with previous data which is likewise
    # dependent on the parameters.
    ctsys = System(grid, data, unknown_storage = :sparse)

    ## boundary model
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphin, bregion1, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphip, bregion1, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, ipsi, bregion1, 0.0)

    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphin, bregion2, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphip, bregion2, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, ipsi, bregion2, 0.0)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = VoronoiFVM.SolverControl()
    control.verbose = verbose
    control.maxiters = 50
    control.abstol = 1.0e-14
    control.reltol = 1.0e-14
    control.tol_round = 1.0e-8
    control.damp_initial = 0.1
    control.damp_growth = 1.21
    control.max_round = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Solving the nonlinear system of equations")
    end
    ################################################################################

    inival = ChargeTransport.unknowns(ctsys)
    inival .= 0.0

    sol = ChargeTransport.solve(ctsys, inival = inival, control = control)

    if parameterStudy
        return sol, ctsys
    end

    ###################################
    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
    tf = VoronoiFVM.testfunction(factory, [bregion1], [bregion2])

    I = VoronoiFVM.integrate(ctsys.fvmsys, tf, sol)

    current = I[1] + I[2]

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false && plotting
        println("Plotting")
    end
    ################################################################################

    if plotting
        # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
        function tridata(grid::ExtendableGrid)
            coord = grid[Coordinates]
            cellnodes = Matrix(grid[CellNodes])
            return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
        end

        XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

        nn = data.F[iphin].(zn * (sol[iphin, :] - sol[ipsi, :]))
        np = data.F[iphip].(zp * (sol[iphip, :] - sol[ipsi, :]))

        vmin = 0.005; vmax = 7.0

        Blues = Plotter.get_cmap(:Blues); Oranges = Plotter.get_cmap(:Oranges)
        Greens = Plotter.get_cmap(:Greens)

        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(nn...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", cmap = parula_map, rasterized = true)
        Plotter.title("\$ n_{\\mathrm{n}} \$")
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("\$ y \$", fontsize = 17)
        Plotter.colorbar(orientation = "vertical", label = "density")
        Plotter.axis("equal")
        Plotter.tight_layout()
        ########################################################
        Plotter.figure()
        Plotter.tripcolor(tridata(grid)..., vcat(np...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", cmap = parula_map, rasterized = true)
        Plotter.title("\$ n_{\\mathrm{p}} \$")
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("\$ y \$", fontsize = 17)
        Plotter.colorbar(orientation = "vertical", label = "density")
        Plotter.axis("equal")
        Plotter.tight_layout()
        ########################################################
        Plotter.figure()
        Plotter.surf(XX[:], YY[:], sol[ipsi, :], color = Blues(201))
        Plotter.title(" \$ \\psi \$")
        Plotter.xlabel("\$ x \$")
        Plotter.ylabel("\$ y \$")
        Plotter.zlabel("potential")
        # Set aspect ratio based on data ranges
        Plotter.gca().set_xlim3d(0.0, 8.0)
        Plotter.gca().set_ylim3d(0.0, 4.0)
        Plotter.gca().set_zlim3d(-3.0, 4.5)

        xlim = Plotter.gca().get_xlim3d()
        ylim = Plotter.gca().get_ylim3d()
        zlim = Plotter.gca().get_zlim3d()

        xr = xlim[2] - xlim[1]
        yr = ylim[2] - ylim[1]
        zr = zlim[2] - zlim[1]
        maxr = maximum([xr, yr, zr])
        Plotter.gca().set_box_aspect((xr / maxr, yr / maxr, zr / maxr))
        Plotter.tight_layout()
        ########################################################
        Plotter.figure()
        Plotter.surf(XX[:], YY[:], zn .* sol[iphin, :], color = Greens(201))
        Plotter.title("\$ v_{\\mathrm{n}} \$")
        Plotter.xlabel("\$ x \$")
        Plotter.ylabel("\$ y \$")
        Plotter.zlabel("potential")
        Plotter.gca().set_xlim3d(0.0, 8.0)
        Plotter.gca().set_ylim3d(0.0, 4.0)
        Plotter.gca().set_zlim3d(-1.0, 1.0)
        Plotter.gca().set_box_aspect((xr / maxr, yr / maxr, zr / maxr))
        Plotter.tight_layout()
        ########################################################
        Plotter.figure()
        Plotter.surf(XX[:], YY[:], zp .* sol[iphip, :], color = Oranges(201))
        Plotter.title("\$ v_{\\mathrm{p}} \$")
        Plotter.xlabel("\$ x \$")
        Plotter.ylabel("\$ y \$")
        Plotter.zlabel("potential")
        Plotter.gca().set_xlim3d(0.0, 8.0)
        Plotter.gca().set_ylim3d(0.0, 4.0)
        Plotter.gca().set_zlim3d(-1.0, 1.0)
        Plotter.gca().set_box_aspect((xr / maxr, yr / maxr, zr / maxr))
        Plotter.tight_layout()

    end

    if test == false && plotting
        println("*** done\n")
    end

    return current

end #  main


## This is one of the functions, with which the current data sets for the parameter study can be generated
function ParameterStudy1D(;
        hmin = 0.1, h = 0.25, hmax = 0.25, # grid information
        Cn = 1.0e1, Cp = 1.0e1,            # doping
        λ = 1.0,                           # Debye length
        G0 = 1.0,                          # laser information
        Plotter = nothing, plotting = false
    )

    if plotting
        if isnothing(Plotter)
            @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
            plotting = false
        else
            if nameof(Plotter) != :PyPlot
                @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
                plotting = false
            end
        end
    end

    if plotting
        Plotter.rc("font", family = "sans-serif", size = 14)
        Plotter.rc("mathtext", fontset = "dejavusans")
        Plotter.close("all")
    end

    X1 = geomspace(0.0, box_cx - box_x / 2, hmin, hmax)
    X2 = geomspace(box_cx - box_x / 2, box_cx, hmax, hmin)
    X3 = geomspace(box_cx, box_cx + box_x / 2, hmin, hmax)
    X4 = geomspace(box_cx + box_x / 2, length_x, hmax, hmin)
    X = glue(X1, X2); X = glue(X, X3); X = glue(X, X4)
    ##################
    Y1 = collect(0.0:h:(box_cy - box_y / 2))
    Y2 = geomspace(box_cy - box_y / 2, box_cy, hmax, hmin)
    Y3 = geomspace(box_cy, box_cy + box_y / 2, hmin, hmax)
    Y4 = collect((box_cy + box_y / 2):h:length_y)
    Y = glue(Y1, Y2); Y = glue(Y, Y3); Y = glue(Y, Y4)
    ##################
    grid = simplexgrid(X, Y)

    rect!(grid, [box_cx - box_x / 2, box_cy - box_y / 2], [box_cx + box_x / 2, box_cy + box_y / 2], region = 2, bregion = 0, tol = 1.0e-1)

    ## bfacemask! for setting different boundary regions.
    bfacemask!(grid, [0.0, 0.0], [0.0, length_y], bregion1)
    bfacemask!(grid, [length_x, 0.0], [length_x, length_y], bregion2)
    bfacemask!(grid, [0.0, 0.0], [length_x, 0.0], 0)
    bfacemask!(grid, [0.0, length_y], [length_x, length_y], 0)

    IVec = zeros(0)
    i = 1

    y0 = box_cy
    for xx in X

        println(i)
        i = i + 1
        x0 = xx

        sol, ctsys = Ex203_LBIC_NonDimensional.main(Cn = Cn, Cp = Cp, λ = λ, G0 = G0, x0 = x0, y0 = y0, test = true, plotting = false, parameterStudy = true)

        #########################################################
        #### calculate total current
        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
        tf = VoronoiFVM.testfunction(factory, [bregion1], [bregion2])

        I = VoronoiFVM.integrate(ctsys.fvmsys, tf, sol)
        current = I[1] + I[2]

        push!(IVec, current)
        #########################################################
    end

    if plotting
        Plotter.plot(X, IVec, linewidth = 4, marker = "o", markersize = 12)
        Plotter.axvspan(0.0, 2.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(2.0, 6.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.axvspan(6.0, 8.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.xlim(0.0, 8.0)
        Plotter.ylim(-0.8, 0.8)
        Plotter.xlabel("Laser position \$ x \$", fontsize = 17)
        Plotter.ylabel("Current", fontsize = 17)
        Plotter.tight_layout()
    end

    return nothing

end

## This is one of the functions, with which the data sets for the parameter study can be generated
function ParameterStudy2D(;
        hmin = 0.1, h = 0.25, hmax = 0.25, # grid information
        Cn = 1.0e1, Cp = 1.0e1,            # doping
        λ = 1.0,                           # Debye length
        G0 = 1.0,                          # laser information
        Plotter = nothing, plotting = false,
    )

    if plotting
        if isnothing(Plotter)
            @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
            plotting = false
        else
            if nameof(Plotter) != :PyPlot
                @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global enviroment via the package manager and choose `Plotter = PyPlot`."
                plotting = false
            end
        end
    end

    if plotting
        Plotter.rc("font", family = "sans-serif", size = 14)
        Plotter.rc("mathtext", fontset = "dejavusans")
        Plotter.close("all")
    end

    X1 = geomspace(0.0, box_cx - box_x / 2, hmin, hmax)
    X2 = geomspace(box_cx - box_x / 2, box_cx, hmax, hmin)
    X3 = geomspace(box_cx, box_cx + box_x / 2, hmin, hmax)
    X4 = geomspace(box_cx + box_x / 2, length_x, hmax, hmin)
    X = glue(X1, X2); X = glue(X, X3); X = glue(X, X4)
    ##################
    Y1 = collect(0.0:h:(box_cy - box_y / 2))
    Y2 = geomspace(box_cy - box_y / 2, box_cy, hmax, hmin)
    Y3 = geomspace(box_cy, box_cy + box_y / 2, hmin, hmax)
    Y4 = collect((box_cy + box_y / 2):h:length_y)
    Y = glue(Y1, Y2); Y = glue(Y, Y3); Y = glue(Y, Y4)
    ##################
    grid = simplexgrid(X, Y)

    rect!(grid, [box_cx - box_x / 2, box_cy - box_y / 2], [box_cx + box_x / 2, box_cy + box_y / 2], region = 2, bregion = 0, tol = 1.0e-1)

    ## bfacemask! for setting different boundary regions.
    bfacemask!(grid, [0.0, 0.0], [0.0, length_y], bregion1)
    bfacemask!(grid, [length_x, 0.0], [length_x, length_y], bregion2)
    bfacemask!(grid, [0.0, 0.0], [length_x, 0.0], 0)
    bfacemask!(grid, [0.0, length_y], [length_x, length_y], 0)

    IVec = zeros(0)
    i = 1

    for coord in eachcol(grid[Coordinates])

        println(i)
        i = i + 1
        x0 = coord[1]
        y0 = coord[2]

        sol, ctsys = Ex203_LBIC_NonDimensional.main(Cn = Cn, Cp = Cp, λ = λ, G0 = G0, x0 = x0, y0 = y0, test = true, plotting = false, parameterStudy = true)

        #########################################################
        #### calculate total current
        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
        tf = VoronoiFVM.testfunction(factory, [bregion1], [bregion2])

        I = VoronoiFVM.integrate(ctsys.fvmsys, tf, sol)
        current = I[1] + I[2]

        push!(IVec, current)

    end

    if plotting
        X = grid[Coordinates][1, :]; Y = grid[Coordinates][2, :]
        Plotter.surf(X[:], Y[:], IVec)
        Plotter.title(" \$ I \$")
        Plotter.xlabel("\$ x \$")
        Plotter.ylabel("\$ y \$")
        Plotter.zlabel("current ")
        Plotter.tight_layout()

        return nothing
    end

end

function test()
    testval = -0.052620329648232746
    return main(test = true) ≈ testval
end


end # module
