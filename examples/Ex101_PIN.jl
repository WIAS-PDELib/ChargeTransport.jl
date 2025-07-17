#=
# GaAs diode (1D).
([source code](SOURCE_URL))

We simulate charge transport in a GaAs pin diode, where we use the van Roosbroeck
system of equations as charge transport model. The unknowns are given by the quasi Fermi
potentials of electrons and holes $\varphi_n$, $\varphi_p$ and the electric potential $\psi$.
The simulations are performed out of equilibrium and for the stationary problem.
=#

module Ex101_PIN

using ChargeTransport  # drift-diffusion solver
using ExtendableGrids  # grid initializer
using PyPlot           # solution visualizer

## This function is used to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
    coord = glue(coord_ndoping, coord_intrinsic)
    coord = glue(coord, coord_pdoping)

    return coord
end

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(; n = 3, Plotter = PyPlot, plotting = false, verbose = "", test = false, unknown_storage = :sparse)

    if plotting
        Plotter.close("all")
    end
    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor = 1          # p doped region
    regionIntrinsic = 2          # intrinsic region
    regionDonor = 3          # n doped region
    regions = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions = length(regions)

    ## boundary region numbers
    ## Note that by convention we have 1 for the left boundary and 2 for the right boundary. If
    ## adding additional interior boundaries, continue with 3, 4, ...
    bregionAcceptor = 1
    bregionDonor = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    ## grid
    refinementfactor = 2^(n - 1)
    h_pdoping = 2.0 * μm
    h_intrinsic = 2.0 * μm
    h_ndoping = 2.0 * μm
    h_total = h_pdoping + h_intrinsic + h_ndoping
    w_device = 0.5 * μm  # width of device
    z_device = 1.0e-4 * cm  # depth of device
    coord = initialize_pin_grid(
        refinementfactor,
        h_pdoping,
        h_intrinsic,
        h_ndoping
    )

    grid = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)  # p-doped region = 1
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region = 3

    ## bfacemask! for setting different boundary regions. At exterior boundaries they are
    ## automatically set by ExtendableGridsjl. Thus, there the following two lines are actually
    ## unneccesarry, but are only written for completeness.
    bfacemask!(grid, [0.0], [0.0], bregionAcceptor)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionDonor)  # outer right boundary
    bfacemask!(grid, [h_pdoping], [h_pdoping], bregionJunction1) # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], bregionJunction2) # second inner interface

    if plotting
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
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
    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    numberOfCarriers = 2

    # We define the physical data.
    Ec = 1.424 * eV
    Ev = 0.0 * eV
    Nc = 4.35195989587969e17 / (cm^3)
    Nv = 9.139615903601645e18 / (cm^3)
    mun = 8500.0 * (cm^2) / (V * s)
    mup = 400.0 * (cm^2) / (V * s)
    εr = 12.9 * 1.0              # relative dielectric permittivity of GAs
    T = 300.0 * K

    ## recombination parameters
    Auger = 1.0e-29 * cm^6 / s
    SRH_TrapDensity = 1.0e10 / cm^3
    SRH_LifeTime = 1.0 * ns
    Radiative = 1.0e-10 * cm^3 / s

    ## doping
    dopingFactorNd = 1.0
    dopingFactorNa = 0.46
    Nd = dopingFactorNd * Nc
    Na = dopingFactorNa * Nv

    ## intrinsic concentration
    ni = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))

    ## contact voltage: we impose an applied voltage only on one boundary.
    ## At the other boundary the applied voltage is zero.
    voltageAcceptor = 1.5 * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # We initialize the Data instance and fill in predefined data.
    data = Data(grid, numberOfCarriers)

    ## Following variable declares, if we want to solve stationary or transient problem
    data.modelType = Stationary

    ## Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk,
    ## FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F .= Boltzmann

    ## Here, we need to specify which numbers are associated with electron and hole quasi
    ## Fermi potential. Further, the desired recombination processes can be chosen here.
    ## Note that, if you choose a SRH recombination you can further specify a transient SRH
    ## recombination by the method enable_trap_carrier! and adjusting the modelType. Otherwise, by
    ## default we use the stationary model for this type of recombination.
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Following choices are possible for boundary model: For contacts currently only
    ## OhmicContact and SchottkyContact are possible. For inner boundaries we have
    ## InterfaceNone, InterfaceRecombination.
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor] = OhmicContact

    ## Following choices are possible for the flux discretization scheme: ScharfetterGummel,
    ## ScharfetterGummelGraded, ExcessChemicalPotential, ExcessChemicalPotentialGraded,
    ## DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    # Define the Params struct. Params contains all necessary physical parameters. If one
    # wants to simulate space-dependent variables, one additionally needs to generate a
    # ParamsNodal struct, see Ex102.
    params = Params(grid, numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ireg in 1:numberOfRegions # region data

        params.dielectricConstant[ireg] = εr * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = mun
        params.mobility[iphip, ireg] = mup

        ## recombination parameters
        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg] = Auger
        params.recombinationAuger[iphip, ireg] = Auger

    end

    ## doping
    params.doping[iphin, regionDonor] = Nd     # data.doping   = [0.0  Na;
    params.doping[iphin, regionIntrinsic] = ni     #                  ni  0.0;
    params.doping[iphip, regionAcceptor] = Na     #                  Nd  0.0]

    # Region dependent params is now a substruct of data which is again a substruct of the
    # system and will be parsed in next step.
    data.params = params

    # In the last step, we initialize our system with previous data which is likewise
    # dependent on the parameters. It is important that this is in the end, otherwise our
    # VoronoiFVMSys is not dependent on the data we initialized but rather on default data.
    ctsys = System(grid, data, unknown_storage = unknown_storage)

    if test == false
        ## Here we can show region dependent physical parameters. show_params() only supports
        ## region dependent parameters, but, if one wishes to print nodal dependent parameters,
        ## currently this is possible with println(ctsys.data.paramsnodal). We neglected here,
        ## since in most applications where the numberOfNodes is >> 10 this would results in a
        ## large output in the terminal.
        show_params(ctsys)
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        ## set legend for plotting routines. Either you can use the predefined labels or write your own.
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        psi0 = electroNeutralSolution(ctsys)
        Plotter.figure()
        plot_energies(Plotter, ctsys, label_BEE)
        Plotter.figure()
        plot_doping(Plotter, ctsys, label_density)
        Plotter.figure()
        plot_electroNeutralSolutionBoltzmann(Plotter, grid, psi0, ; plotGridpoints = true)
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 50
    control.abstol = 1.0e-14
    control.reltol = 1.0e-14
    control.tol_round = 1.0e-8
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

    ## calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    maxBias = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 32)
    IV = zeros(0)

    for Δu in biasValues

        if test == false
            println("bias value: Δu = ", Δu, " V")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## Note that the old way of solving will be soon removed (see current API changes in VoronoiFVM)
        # solve!(solution, inival, ctsys, control = control, tstep = Inf)
        # inival .= solution

        ## get I-V data
        current = get_current_val(ctsys, solution)

        push!(IV, abs.(w_device * z_device * (current)))

    end # bias loop

    if test == false
        println("*** done\n")
    end

    ## plot solution and IV curve
    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_energy, plotGridpoints = false)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_solution, plotGridpoints = true)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_density, plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues, IV, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
    end

    testval = solution[15]
    return testval

    if test == false
        println("*** done\n")
    end

end #  main

function test()
    testval = 1.5068426833371802
    return main(test = true, unknown_storage = :dense) ≈ testval && main(test = true, unknown_storage = :sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
