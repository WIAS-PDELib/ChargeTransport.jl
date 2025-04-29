# Illustrative example of a three dimensional simulation.
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examples/PSC_3D.jl))

This code shows the capability of 3D simulations with ChargeTransport.jl.
For the sake of performance, we only do equilibrium calculations.

Here, a one-dimensional and a three-dimensional simulation of the same device are performed.

The parameters are based on the default parameter set of Ionmonger (with minor adjustments).

````julia
module PSC_3D

using ChargeTransport
using ExtendableGrids
using GridVisualize

using GLMakie
using PyPlot
````

We strongly emphasize to use GLMakie for the visualization here.

````julia
function main(;Plotter = GLMakie, plotting = false, test = false, verbose = false,
              parameter_file = "../parameter_files/Params_PSC_TiO2_MAPI_spiro.jl", # choose the parameter file
             )

    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file) # include the parameter file we specified

    bregionNoFlux = 5
    height        = 2.00e-5 * cm
    width         = 3.00e-5 * cm

    if test == false
        println("*** done\n")
    end

    ################################################################################
     if test == false
        println("Set up grid and regions for 1D and 3D")
    end
    ################################################################################

    # 1D Grid
    n                = 10
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = n))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 2 * n))
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic), stop = (h_total), length = n))
    coord            = glue(coord_ndoping, coord_intrinsic)
    coord            = glue(coord, coord_pdoping)
    grid1D           = simplexgrid(coord)

    cellmask!(grid1D, [0.0 * μm],                 [h_ndoping],               regionDonor, tol = 1.0e-18)
    cellmask!(grid1D, [h_ndoping],                [h_ndoping + h_intrinsic], regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid1D, [h_ndoping + h_intrinsic],  [h_total],                 regionAcceptor, tol = 1.0e-18)

    bfacemask!(grid1D, [heightLayers[1]], [heightLayers[1]], bregionJ1) # first  inner interface
    bfacemask!(grid1D, [heightLayers[2]], [heightLayers[2]], bregionJ2) # second inner interface

    # 3D Grid
    coord_height     = collect(range(0.0, stop = height, length = n))
    coord_width      = collect(range(0.0, stop = width, length =  n))
    grid3D           = simplexgrid(coord, coord_height, coord_width)

    cellmask!(grid3D, [0.0, 0.0, 0.0],                   [h_ndoping, height, width],               regionDonor, tol = 1.0e-18)
    cellmask!(grid3D, [h_ndoping, 0.0, 0.0],             [h_ndoping + h_intrinsic, height, width], regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid3D, [h_ndoping+h_intrinsic, 0.0, 0.0], [h_total, height, width],                 regionAcceptor, tol = 1.0e-18)

    # metal interfaces [xmin, ymin, zmin], [xmax, ymax, zmax]
    bfacemask!(grid3D, [0.0, 0.0, 0.0],     [0.0, height, width],     bregionDonor) # BregionNumber = 1
    bfacemask!(grid3D, [h_total, 0.0, 0.0], [h_total, height, width], bregionAcceptor) # BregionNumber = 2

    # interior interfaces
    bfacemask!(grid3D, [heightLayers[1], 0.0, 0.0], [heightLayers[1], height, width], bregionJ1) # first  inner interface
    bfacemask!(grid3D, [heightLayers[2], 0.0, 0.0], [heightLayers[2], height, width], bregionJ2) # second inner interface

    # outer no flux interfaces
    bfacemask!(grid3D, [0.0, 0.0, 0.0],    [h_total, 0.0, width],    bregionNoFlux)
    bfacemask!(grid3D, [0.0, height, 0.0], [h_total, height, width], bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, 0.0],    [h_total, height, 0.0],   bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, width],  [h_total, height, width], bregionNoFlux)

    if plotting == true # plotting is currently only tested with GLMakie and PyPlot
        vis    = GridVisualizer(Plotter = Plotter, resolution=(1500,1500), layout=(3,2))
        gridplot!(vis[1,1], grid1D)
        if Plotter == PyPlot
            gridplot!(vis[1,2], grid3D, linewidth=0.5, xplanes=[5.5e-7], zplanes=[1.5e-7])
        elseif Plotter == GLMakie
            gridplot!(vis[1,2], grid3D, zplane=1.0e-7,azim=20,elev=60,linewidth=0.5, scene3d=:Axis3, legend=:none)
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

    # Initialize Data instance and fill in predefined data
    # Note that we define the data struct with respect to the three-dimensional grid, since we also defined there the outer no flux boundary conditions.
    data                               = Data(grid3D, numberOfCarriers)
    data.modelType                     = Transient
    data.F                             = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)
    data.boundaryType[bregionDonor]    = OhmicContact
    data.boundaryType[bregionAcceptor] = OhmicContact

    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    data.fluxApproximation            .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

      ################################################################################
      if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                              = Params(grid3D, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = zn
    params.chargeNumbers[iphip]                         = zp
    params.chargeNumbers[iphia]                         = za

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        # effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nn[ireg]
        params.densityOfStates[iphip, ireg]             = Np[ireg]
        params.densityOfStates[iphia, ireg]             = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = Ea[ireg]

        # recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, params, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, params, EI[ireg])
    end

    # interior doping
    params.doping[iphin, regionDonor]                   = Cn
    params.doping[iphia, regionIntrinsic]               = Ca
    params.doping[iphip, regionAcceptor]                = Cp

    data.params                                         = params
    ctsys1D                                             = System(grid1D, data, unknown_storage=:sparse)

    data.params                                         = params
    ctsys3D                                             = System(grid3D, data, unknown_storage=:sparse)

    ipsi = data.index_psi

    if test == false
        show_params(ctsys1D)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control              = SolverControl()
    control.verbose      = verbose
    control.maxiters     = 300
    control.abstol       = 1.0e-10
    control.reltol       = 1.0e-10
    control.tol_round    = 1.0e-10
    control.max_round    = 5
    control.damp_initial = 0.5
    control.damp_growth  = 1.61 # >= 1

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    sol1D = equilibrium_solve!(ctsys1D, control = control, nonlinear_steps = 20)
    sol3D = equilibrium_solve!(ctsys3D, control = control, nonlinear_steps = 20)

    if plotting == true
        #################################################
        scalarplot!(vis[2,1], grid1D, sol1D[ipsi, :]; color=:blue, linewidth = 5, xlabel = "space [m]", ylabel = "potential [V]", title = "Electric potential (1D)")
        scalarplot!(vis[2,2], grid3D, sol3D[ipsi, :]; scene3d=:Axis3, levels = 4, levelalpha = 0.9, outlinealpha = 0.00, xplanes = collect(range(0.0, stop = h_total, length = 100)), title = "Electric potential (3D)")

        grids1D    = Array{typeof(grid1D), 1}(undef, numberOfRegions)
        densityn1D = Array{typeof(sol1D[iphin, :]), 1}(undef, numberOfRegions)

        grids3D    = Array{typeof(grid3D), 1}(undef, numberOfRegions)
        densityn3D = Array{typeof(sol3D[iphin, :]), 1}(undef, numberOfRegions)
        logDens3D  = Array{typeof(sol3D[iphin, :]), 1}(undef, numberOfRegions)

        for ireg in 1:numberOfRegions
            grids1D[ireg]    = subgrid(grid1D, [ireg])
            densityn1D[ireg] = get_density(sol1D, ireg, ctsys1D, iphin)
            #############################################################
            grids3D[ireg]    = subgrid(grid3D, [ireg])
            densityn3D[ireg] = get_density(sol3D, ireg, ctsys3D, iphin)
            logDens3D[ireg]  = log.(densityn3D[ireg])
        end

        scalarplot!(vis[3,1], grids1D, grid1D, densityn1D; color=:green, linewidth = 5, yscale=:log, xlabel = "space [m]", ylabel = "density [\$\\frac{1}{m^3}\$]", title = "Electron concentration (1D)")
        scalarplot!(vis[3,2], grids3D, grid3D, densityn3D; scene3d=:Axis3, levels = 4, levelalpha = 0.9, outlinealpha = 0.00, xplanes = collect(range(0.0, stop = h_total, length = 100)), title = "Electron concentration (3D)")
    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, sol1D))/length(sol1D) + sum(filter(!isnan, sol3D))/length(sol3D) # when using sparse storage, we get NaN values in solution

    return testval

end # main


function test()
    testval = -2.2213072819274533
    main(test = true) ≈ testval
end

end # module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
