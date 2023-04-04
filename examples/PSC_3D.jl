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

using GLMakie
using PyPlot

# We strongly emphasize to use GLMakie for the visualization here.
function main(;Plotter = GLMakie, plotting = false, test = false, verbose = false)

     ################################################################################
     if test == false
        println("Set up grid and regions for 1D and 3D")
    end
    ################################################################################

    ## region numbers
    regionDonor      = 1                           # n doped region
    regionIntrinsic  = 2                           # intrinsic region
    regionAcceptor   = 3                           # p doped region
    regions          = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions  = length(regions)

    ## boundary region numbers
    bregionDonor     = 1
    bregionAcceptor  = 2
    bregionJunction1 = 3
    bregionJunction2 = 4
    bregionNoFlux    = 5

    ## grid
    h_ndoping        = 1.00e-5 * cm
    h_intrinsic      = 4.00e-5 * cm
    h_pdoping        = 2.00e-5 * cm
    h_total          = h_ndoping + h_intrinsic + h_pdoping
    heightLayers     = [h_ndoping, h_ndoping + h_intrinsic, h_ndoping + h_intrinsic + h_pdoping]
    height           = 2.00e-5 * cm
    width            = 3.00e-5 * cm

    ## 1D Grid
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

    bfacemask!(grid1D, [heightLayers[1]], [heightLayers[1]], bregionJunction1) # first  inner interface
    bfacemask!(grid1D, [heightLayers[2]], [heightLayers[2]], bregionJunction2) # second inner interface

    ## 3D Grid
    coord_height     = collect(range(0.0, stop = height, length = n))
    coord_width      = collect(range(0.0, stop = width, length =  n))
    grid3D           = simplexgrid(coord, coord_height, coord_width)

    cellmask!(grid3D, [0.0, 0.0, 0.0],                   [h_ndoping, height, width],               regionDonor, tol = 1.0e-18)
    cellmask!(grid3D, [h_ndoping, 0.0, 0.0],             [h_ndoping + h_intrinsic, height, width], regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid3D, [h_ndoping+h_intrinsic, 0.0, 0.0], [h_total, height, width],                 regionAcceptor, tol = 1.0e-18)

    ## metal interfaces [xmin, ymin, zmin], [xmax, ymax, zmax]
    bfacemask!(grid3D, [0.0, 0.0, 0.0],     [0.0, height, width],     bregionDonor) # BregionNumber = 1
    bfacemask!(grid3D, [h_total, 0.0, 0.0], [h_total, height, width], bregionAcceptor) # BregionNumber = 2

    ## interior interfaces
    bfacemask!(grid3D, [heightLayers[1], 0.0, 0.0], [heightLayers[1], height, width], bregionJunction1) # first  inner interface
    bfacemask!(grid3D, [heightLayers[2], 0.0, 0.0], [heightLayers[2], height, width], bregionJunction2) # second inner interface

    ## outer no flux interfaces
    bfacemask!(grid3D, [0.0, 0.0, 0.0],    [h_total, 0.0, width],    bregionNoFlux)
    bfacemask!(grid3D, [0.0, height, 0.0], [h_total, height, width], bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, 0.0],    [h_total, height, 0.0],   bregionNoFlux)
    bfacemask!(grid3D, [0.0, 0.0, width],  [h_total, height, width], bregionNoFlux)

    if plotting == true # plotting is currently only tested with GLMakie and PyPlot
        vis    = GridVisualizer(Plotter = Plotter, resolution=(1500,1500), layout=(2,2))
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
        println("Define physical parameters and model")
    end
    ################################################################################

    ## set indices of the quasi Fermi potentials
    iphin            = 1 # electron quasi Fermi potential
    iphip            = 2 # hole quasi Fermi potential
    iphia            = 3 # anion vacancy quasi Fermi potential
    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ipsi             = 4

    ## temperature
    T                = 300.0                 *  K

    ## band edge energies
    Ec_d             = -4.0                  *  eV
    Ev_d             = -5.8                  *  eV

    Ec_i             = -3.7                  *  eV
    Ev_i             = -5.4                  *  eV

    Ec_a             = -3.4                  *  eV
    Ev_a             = -5.1                  *  eV
    Ea_i             = -4.45                 *  eV

    EC               = [Ec_d, Ec_i, Ec_a]
    EV               = [Ev_d, Ev_i, Ev_a]
    EA               = [0.0,  Ea_i,  0.0]

    ## effective densities of state
    Nc_d             = 5.0e19                / (cm^3)
    Nv_d             = 5.0e19                / (cm^3)

    Nc_i             = 8.1e18                / (cm^3)
    Nv_i             = 5.8e18                / (cm^3)
    Nanion           = 1.0e21                / (cm^3)

    Nc_a             = 5.0e19                / (cm^3)
    Nv_a             = 5.0e19                / (cm^3)

    NC               = [Nc_d, Nc_i,  Nc_a]
    NV               = [Nv_d, Nv_i,  Nv_a]
    NAnion           = [0.0,  Nanion, 0.0]

    ## relative dielectric permittivity
    ε_d              = 10.0                  *  1.0
    ε_i              = 24.1                  *  1.0
    ε_a              = 3.0                   *  1.0

    ε                = [ε_d, ε_i, ε_a]

    ## radiative recombination
    r0_d             = 0.0e+0               * cm^3 / s
    r0_i             = 1.0e-12              * cm^3 / s
    r0_a             = 0.0e+0               * cm^3 / s

    r0               = [r0_d, r0_i, r0_a]

    ## life times and trap densities
    τn_d             = 1.0e100              * s
    τp_d             = 1.0e100              * s

    τn_i             = 3.0e-10              * s
    τp_i             = 3.0e-8               * s
    τn_a             = τn_d
    τp_a             = τp_d

    τn               = [τn_d, τn_i, τn_a]
    τp               = [τp_d, τp_i, τp_a]

    ## SRH trap energies
    Ei_d             = -5.0                 * eV
    Ei_i             = -4.55                * eV
    Ei_a             = -4.1                 * eV

    EI               = [Ei_d, Ei_i, Ei_a]

    ## doping
    Nd               = 1.03e18              / (cm^3)
    Na               = 1.03e18              / (cm^3)
    C0               = 1.6e19               / (cm^3)

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
    data                               = Data(grid3D, numberOfCarriers)
    data.modelType                     = Transient
    data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]
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
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    params.chargeNumbers[iphia]                         =  1

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]
        params.densityOfStates[iphia, ireg]             = NAnion[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = EA[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
    end

    ## boundary region data
    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    ## interior doping
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphia, regionIntrinsic]               = C0
    params.doping[iphip, regionAcceptor]                = Na

    ## boundary doping
    params.bDoping[iphip, bregionAcceptor]              = Na
    params.bDoping[iphin, bregionDonor]                 = Nd

    data.params                                         = params
    ctsys1D                                             = System(grid1D, data, unknown_storage=:sparse)

    data.params                                         = params
    ctsys3D                                             = System(grid3D, data, unknown_storage=:sparse)

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

        ## DA: We have a problem that logarithmic scaling of the colorbar is not working in 3D ....
        # grids1D    = Array{typeof(grid1D), 1}(undef, numberOfRegions)
        # densityn1D = Array{typeof(sol1D[iphin, :]), 1}(undef, numberOfRegions)

        # grids3D    = Array{typeof(grid3D), 1}(undef, numberOfRegions)
        # densityn3D = Array{typeof(sol3D[iphin, :]), 1}(undef, numberOfRegions)
        # logDens3D  = Array{typeof(sol3D[iphin, :]), 1}(undef, numberOfRegions)

        # for ireg in 1:numberOfRegions
        #     grids1D[ireg]    = subgrid(grid1D, [ireg])
        #     densityn1D[ireg] = get_density(sol1D, ireg, ctsys1D, iphin)
        #     #############################################################
        #     grids3D[ireg]    = subgrid(grid3D, [ireg])
        #     densityn3D[ireg] = get_density(sol3D, ireg, ctsys3D, iphin)
        #     logDens3D[ireg]  = log.(densityn3D[ireg])
        # end

        # scalarplot!(vis[2,1], grids1D, grid1D, densityn1D; color=:green, linewidth = 5, yscale=:log, xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", title = "Electron concentration (1D)")
        # scalarplot!(vis[2,2], grids3D, grid3D, densityn3D; scene3d=:Axis3, levels = 4, levelalpha = 0.9, outlinealpha = 0.00, xplanes = collect(range(0.0, stop = h_total, length = 100)), title = "Electron concentration (3D)")
    end

    return vis, grids3D, grid3D, densityn3D

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, sol1D))/length(sol1D) + sum(filter(!isnan, sol3D))/length(sol3D) # when using sparse storage, we get NaN values in solution

    return testval

end # main


function test()
    testval = -2.221100737029875
    main(test = true) ≈ testval
end

end # module