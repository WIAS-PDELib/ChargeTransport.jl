# GaAs diode: transient with traps (1D).
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examples/Ex109_Traps.jl))

Simulating transient charge transport in a GaAs p-i-n diode with an electron trap.

````julia
module Ex109_Traps

using ChargeTransport
using ExtendableGrids
using PyPlot

# function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping   = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping   = collect(range((h_ndoping + h_intrinsic),
                                     stop = (h_ndoping + h_intrinsic + h_pdoping),
                                     length = 3 * refinementfactor))
    coord           = glue(coord_ndoping, coord_intrinsic)
    coord           = glue(coord, coord_pdoping)

    return coord
end

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionAcceptor   = 1                           # p doped region
    regionIntrinsic  = 2                           # intrinsic region
    regionDonor      = 3                           # n doped region
    regions          = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions  = length(regions)

    # boundary region numbers
    bregionAcceptor  = 1
    bregionDonor     = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    # grid
    refinementfactor = 2^(n-1)
    h_pdoping        = 2.0    * μm
    h_intrinsic      = 2.0    * μm
    h_ndoping        = 2.0    * μm
    h_total          = h_pdoping + h_intrinsic + h_ndoping
    w_device         = 0.5    * μm  # width of device
    z_device         = 1.0e-4 * cm  # depth of device

    coord            = initialize_pin_grid(refinementfactor,
                                           h_pdoping,
                                           h_intrinsic,
                                           h_ndoping)

    grid             = simplexgrid(coord)

    # set different regions in grid
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],               regionAcceptor)  # p-doped
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic], regionIntrinsic) # intrinsic
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_total],                 regionDonor)     # n-doped

    bfacemask!(grid, [h_pdoping],               [h_pdoping],               bregionJunction1, tol = 1.0e-18)
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], bregionJunction2, tol = 1.0e-18)

    if plotting
        gridplot(grid, Plotter = Plotter, legend=:lt)
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

    iphin            = 1 # index electron quasi Fermi potential
    iphip            = 2 # index hole quasi Fermi potential
    iphit            = 3 # index trap quasi Fermi potential
    numberOfCarriers = 3 # electrons, holes and traps

    # physical data
    Ec               = 1.424                             *  eV
    Ev               = 0.0                               *  eV
    Et               = 0.6                               *  eV
    Nc               = 4.351959895879690e17              / (cm^3)
    Nv               = 9.139615903601645e18              / (cm^3)
    Nt               = 1e16                              / (cm^3)
    mun              = 8500.0                            * (cm^2) / (V * s)
    mup              = 400.0                             * (cm^2) / (V * s)
    mut              = 0.0                               * (cm^2) / (V * s)  # such that there is no flux
    εr               = 12.9                              *  1.0              # relative dielectric permittivity of GAs
    T                = 300.0                             *  K

    # recombination parameters
    ni               = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))        # intrinsic concentration
    n0               = Nc * Boltzmann( (Et-Ec) / (kB*T) )                    # Boltzmann equilibrium concentration
    p0               = ni^2 / n0                                             # Boltzmann equilibrium concentration
    Auger            = 1.0e-29                           * cm^6 / s          # 1.0e-41
    SRH_LifeTime     = 1.0e-3                            * ns
    Radiative        = 1.0e-10                           * cm^3 / s          # 1.0e-16
    G                = 1.0e25                            / (cm^3 * s)

    # doping -- trap doping will not be set and thus automatically zero
    dopingFactorNd   = 1.0
    dopingFactorNa   = 0.46
    Nd               = dopingFactorNd * Nc
    Na               = dopingFactorNa * Nv

    # contact voltage
    voltageAcceptor  = 1.5                               * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # Initialize Data instance and fill in data
    data                               = Data(grid, numberOfCarriers)

    # Possible choices: Stationary, Transient
    data.modelType                     = Transient

    # Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    # FermiDiracMinusOne, Blakemore
    data.F                            .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    # Here, we enable the traps and parse the respective index and the regions where the trap is defined.
    enable_trap_carrier!(;data = data, trapCarrier = iphit, regions = regions)

    # Possible choices: GenerationNone, GenerationUniform
    data.generationModel               = GenerationUniform

    # Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    # InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

    # Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    # ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation            .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                              = Params(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    params.chargeNumbers[iphit]                         = -1 # trap charge number determines whether hole or electron trap is used

    for ireg in 1:numberOfRegions           # region data

        params.dielectricConstant[ireg]                 = εr * ε0

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.densityOfStates[iphit, ireg]             = Nt
        params.bandEdgeEnergy[iphin, ireg]              = Ec
        params.bandEdgeEnergy[iphip, ireg]              = Ev
        params.bandEdgeEnergy[iphit, ireg]              = Et
        params.mobility[iphin, ireg]                    = mun
        params.mobility[iphip, ireg]                    = mup
        params.mobility[iphit, ireg]                    = mut

        # recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = n0
        params.recombinationSRHTrapDensity[iphip, ireg] = p0
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger
        params.generationUniform[ireg]                  = G

    end

    # doping
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphin, regionIntrinsic]               = ni
    params.doping[iphip, regionIntrinsic]               = 0.0
    params.doping[iphip, regionAcceptor]                = Na

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control                = SolverControl()
    control.verbose        = verbose
    control.tol_round      = 1.0e-4
    control.damp_initial   = 0.5
    control.damp_growth    = 1.61
    control.max_round      = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # solve thermodynamic equilibrium and update initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival   = solution

    if test == false
        println("*** done\n")
    end

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        # add labels for traps
        label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
        label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"

        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
    end

    ################################################################################
    if test == false
        println("Loop for putting generation on")
    end
    ################################################################################

    # Scan rate and time steps
    scanrate      = 1.0 * V/s
    number_tsteps = 25
    endVoltage    = voltageAcceptor # bias goes until the given voltage at acceptor boundary

    # with fixed timestep sizes we can calculate the times
    # a priori
    tend          = endVoltage/scanrate
    tvalues       = range(0.0, stop = tend, length = number_tsteps)
    Δt            = tvalues[2] - tvalues[1]
````

these values are needed for putting the generation slightly on

````julia
    I             = collect(20:-1:0.0)
    LAMBDA        = 10 .^ (-I)

    IV            = zeros(0) # for IV values
    biasValues    = zeros(0) # for bias values

    for istep = 1:length(I)-1

        # turn slowly generation on
        data.λ2   = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival   = solution

    end # generation loop

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    for istep = 2:number_tsteps

        t  = tvalues[istep]          # Actual time
        Δu = t * scanrate            # Applied voltage
        Δt = t - tvalues[istep-1]    # Time step size

        # Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t) s")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)

        # get I-V data
        current  = get_current_val(ctsys, solution, inival, Δt)

        push!(IV, w_device * z_device * current)
        push!(biasValues, Δu)

        inival   = solution

    end # bias loop

    if test == false
        println("*** done\n")
    end

    # plot solution and IV curve
    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_solution)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, "bias \$\\Delta u\$ = $(biasValues[end])", plotGridpoints = true)
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

    if test == false
        println("*** done\n")
    end

end #  main

function test()
    testval = 0.9699245385329192
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
