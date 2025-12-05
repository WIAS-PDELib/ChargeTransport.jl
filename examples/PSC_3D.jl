#=
# Illustrative example of a three dimensional simulation.
([source code](SOURCE_URL))

This code shows the capability of 3D simulations with ChargeTransport.jl.
For the sake of performance, we only do equilibrium calculations.

Here, a one-dimensional and a three-dimensional simulation of the same device are performed.

The parameters are based on the default parameter set of Ionmonger (with minor adjustments).
=#


module PSC_3D

using ChargeTransport
using ExtendableGrids
using GridVisualize

using PyPlot

# We strongly emphasize to use Plotter = GLMakie for the visualization here.
function main(;
        Plotter = PyPlot, plotting = false, test = false, verbose = false,
        parameter_set = Params_PSC_TiO2_MAPI_spiro, # choose the parameter set
    )

    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    @local_unitfactors μm cm s ns V K ps Hz W

    # parameter with variation
    p = parameter_set(numberOfBoundaryRegions = 5)
    bregionNoFlux = 5

    height = 2.0e-5 * cm
    width = 3.0e-5 * cm

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions for 1D and 3D")
    end
    ################################################################################

    ## 1D Grid
    n = 10
    coord_ndoping = collect(range(0.0, stop = p.h_ndoping, length = n))
    coord_intrinsic = collect(range(p.h_ndoping, stop = (p.h_ndoping + p.h_intrinsic), length = 2 * n))
    coord_pdoping = collect(range((p.h_ndoping + p.h_intrinsic), stop = (p.h_total), length = n))
    coord = glue(coord_ndoping, coord_intrinsic)
    coord = glue(coord, coord_pdoping)
    grid1D = simplexgrid(coord)

    cellmask!(grid1D, [0.0 * μm], [p.h_ndoping], p.regionDonor, tol = 1.0e-18)
    cellmask!(grid1D, [p.h_ndoping], [p.h_ndoping + p.h_intrinsic], p.regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid1D, [p.h_ndoping + p.h_intrinsic], [p.h_total], p.regionAcceptor, tol = 1.0e-18)

    bfacemask!(grid1D, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1) # first  inner interface
    bfacemask!(grid1D, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2) # second inner interface

    ## 3D Grid
    coord_height = collect(range(0.0, stop = height, length = n))
    coord_width = collect(range(0.0, stop = width, length = n))
    grid3D = simplexgrid(coord, coord_height, coord_width)

    cellmask!(grid3D, [0.0, 0.0, 0.0], [p.h_ndoping, height, width], p.regionDonor, tol = 1.0e-18)
    cellmask!(grid3D, [p.h_ndoping, 0.0, 0.0], [p.h_ndoping + p.h_intrinsic, height, width], p.regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid3D, [p.h_ndoping + p.h_intrinsic, 0.0, 0.0], [p.h_total, height, width], p.regionAcceptor, tol = 1.0e-18)

    ## metal interfaces [xmin, ymin, zmin], [xmax, ymax, zmax]
    bfacemask!(grid3D, [0.0, 0.0, 0.0], [0.0, height, width], p.bregionDonor) # BregionNumber = 1
    bfacemask!(grid3D, [p.h_total, 0.0, 0.0], [p.h_total, height, width], p.bregionAcceptor) # BregionNumber = 2

    ## interior interfaces
    bfacemask!(grid3D, [p.heightLayers[1], 0.0, 0.0], [p.heightLayers[1], height, width], p.bregionJ1) # first  inner interface
    bfacemask!(grid3D, [p.heightLayers[2], 0.0, 0.0], [p.heightLayers[2], height, width], p.bregionJ2) # second inner interface

    ## outer no flux interfaces
    bfacemask!(grid3D, [0.0, 0.0, 0.0], [p.h_total, 0.0, width], bregionNoFlux)
    bfacemask!(grid3D, [0.0, height, 0.0], [p.h_total, height, width], bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, 0.0], [p.h_total, height, 0.0], bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, width], [p.h_total, height, width], bregionNoFlux)

    vis = GridVisualizer(Plotter = Plotter, resolution = (1500, 1500), layout = (3, 2))

    if plotting == true && Plotter !== nothing # plotting is currently only tested with GLMakie and PyPlot
        gridplot!(vis[1, 1], grid1D)
        if nameof(Plotter) == :PyPlot
            gridplot!(vis[1, 2], grid3D, linewidth = 0.5, xplanes = [5.5e-7], zplanes = [1.5e-7])
        elseif nameof(Plotter) == :GLMakie
            gridplot!(vis[1, 2], grid3D, zplane = 1.0e-7, azim = 20, elev = 60, linewidth = 0.5, scene3d = :Axis3, legend = :none)
        end
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    ## Note that we define the data struct with respect to the three-dimensional grid, since we also defined there the outer no flux boundary conditions.
    data1D = Data(grid1D, p.numberOfCarriers)
    data1D.modelType = Transient
    data1.F[p.iphin] = FermiDiracOneHalfTeSCA
    data1.F[p.iphip] = FermiDiracOneHalfTeSCA
    data1D.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )
    data1D.boundaryType[p.bregionDonor] = OhmicContact
    data1D.boundaryType[p.bregionAcceptor] = OhmicContact

    enable_ionic_carrier!(data1D, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    ########################
    data3D = Data(grid3D, p.numberOfCarriers)
    data3D.modelType = Transient
    data3.F[p.iphin] = FermiDiracOneHalfTeSCA
    data3.F[p.iphip] = FermiDiracOneHalfTeSCA
    data3D.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )
    data3D.boundaryType[p.bregionDonor] = OhmicContact
    data3D.boundaryType[p.bregionAcceptor] = OhmicContact

    enable_ionic_carrier!(data3D, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    data1D.params = Params(p)
    data3D.params = Params(p)
    ctsys1D = System(grid1D, data1D, unknown_storage = :sparse)
    ctsys3D = System(grid3D, data3D, unknown_storage = :sparse)

    ipsi = data1D.index_psi

    if test == false
        show_params(ctsys1D)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 300
    control.abstol = 1.0e-10
    control.reltol = 1.0e-10
    control.tol_round = 1.0e-10
    control.max_round = 5
    control.damp_initial = 0.5
    control.damp_growth = 1.61 # >= 1

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    sol1D = equilibrium_solve!(ctsys1D, control = control, vacancyEnergyCalculation = false)

    # We already did the calculations to find the correct energy for the vacancies. Due to that, we just copied the value here.
    ctsys3D.data.params.bandEdgeEnergy[p.iphia, p.regionIntrinsic] = -4.461 * data3D.constants.q
    # In case you want that the solver also calculates this value, simply set vacancyEnergyCalculation = true.
    sol3D = equilibrium_solve!(ctsys3D, control = control, vacancyEnergyCalculation = false)

    if plotting == true
        #################################################
        scalarplot!(vis[2, 1], grid1D, sol1D[ipsi, :]; color = :blue, linewidth = 5, xlabel = "space [m]", ylabel = "potential [V]", title = "Electric potential (1D)")
        scalarplot!(vis[2, 2], grid3D, sol3D[ipsi, :]; scene3d = :Axis3, levels = 4, levelalpha = 0.9, outlinealpha = 0.0, xplanes = collect(range(0.0, stop = p.h_total, length = 100)), title = "Electric potential (3D)")

        grids1D = Array{typeof(grid1D), 1}(undef, p.numberOfRegions)
        densityn1D = Array{typeof(sol1D[p.iphin, :]), 1}(undef, p.numberOfRegions)

        grids3D = Array{typeof(grid3D), 1}(undef, p.numberOfRegions)
        densityn3D = Array{typeof(sol3D[p.iphin, :]), 1}(undef, p.numberOfRegions)
        logDens3D = Array{typeof(sol3D[p.iphin, :]), 1}(undef, p.numberOfRegions)

        for ireg in 1:p.numberOfRegions
            grids1D[ireg] = subgrid(grid1D, [ireg])
            densityn1D[ireg] = get_density(sol1D, ireg, ctsys1D, p.iphin)
            #############################################################
            grids3D[ireg] = subgrid(grid3D, [ireg])
            densityn3D[ireg] = get_density(sol3D, ireg, ctsys3D, p.iphin)
            logDens3D[ireg] = log.(densityn3D[ireg])
        end

        scalarplot!(vis[3, 1], grids1D, grid1D, densityn1D; color = :green, linewidth = 5, yscale = :log, xlabel = "space [m]", ylabel = "density [\$\\frac{1}{m^3}\$]", title = "Electron concentration (1D)")
        scalarplot!(vis[3, 2], grids3D, grid3D, densityn3D; scene3d = :Axis3, levels = 4, levelalpha = 0.9, outlinealpha = 0.0, xplanes = collect(range(0.0, stop = p.h_total, length = 100)), title = "Electron concentration (3D)")
    end

    if test == false
        println("*** done\n")
    end

    if test == false
        integral1D = integrated_density(ctsys1D, sol = sol1D, icc = p.iphia, ireg = p.regionIntrinsic)
        integral3D = integrated_density(ctsys3D, sol = sol3D, icc = p.iphia, ireg = p.regionIntrinsic)

        println("Calculated average vacancy density (1D) is: ", integral1D / ctsys1D.data.regionVolumes[p.regionIntrinsic])
        println("Calculated average vacancy density (3D) is: ", integral3D / ctsys3D.data.regionVolumes[p.regionIntrinsic])
        println(" ")
    end

    if test == false
        println("Value for vacancy energy (1D) is: ", ctsys1D.data.params.bandEdgeEnergy[p.iphia, p.regionIntrinsic] / data1D.constants.q, " eV. ")
        println("Value for vacancy energy (3D) is: ", ctsys3D.data.params.bandEdgeEnergy[p.iphia, p.regionIntrinsic] / data3D.constants.q, " eV. ")
        println("Save these values for later use. We recommend to calculate them on a fine grid.")
        println(" ")
    end

    testval = sum(filter(!isnan, sol1D)) / length(sol1D) + sum(filter(!isnan, sol3D)) / length(sol3D) # when using sparse storage, we get NaN values in solution

    return testval

end # main


function test()
    testval = -2.224244519821803
    return main(test = true) ≈ testval
end

end # module
