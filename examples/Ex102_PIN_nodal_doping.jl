#=
# GaAs diode with spatially varying doping (1D).
([source code](@__SOURCE_URL__))

Simulating charge transport in a GaAs pin diode. This means the PDE problem corresponds to the
van Roosbroeck system of equations. The simulations are performed out of equilibrium and for
the stationary problem. A special feature here is that the doping is node-dependent.
=#

module Ex102_PIN_nodal_doping

using ChargeTransport
using ExtendableGrids
using GridVisualize
using LaTeXStrings

# supported Plotters are GLMakie and PythonPlot
# you can set verbose also to true to display some solver information
function main(; Plotter = nothing, verbose = false, test = false, unknown_storage = :sparse)

    # unit factors and constants
    @local_unitfactors μm cm s ns V K ps
    constants = ChargeTransport.constants

    eV = constants.q * V

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor = 1          # p doped region
    regionIntrinsic = 2         # intrinsic region
    regionDonor = 3             # n doped region
    regions = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions = length(regions)

    ## boundary region numbers
    ## Note that by convention we have 1 for the left boundary and 2 for the right boundary. If
    ## adding additional interior boundaries, continue with 3, 4, ...
    bregionAcceptor = 1
    bregionDonor = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    h_pdoping = 0.1 * μm
    h_intrinsic = 0.1 * μm
    h_ndoping = 0.1 * μm
    h_total = h_pdoping + h_intrinsic + h_ndoping
    w_device = 0.1 * μm  # width of device
    z_device = 1.0e-5 * cm  # depth of device

    coord = range(0.0, stop = h_ndoping + h_intrinsic + h_pdoping, length = 25)
    coord = collect(coord)
    grid = simplexgrid(coord)
    numberOfNodes = length(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor, tol = 1.0e-15)    # p-doped region = 1
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], regionIntrinsic, tol = 1.0e-15)    # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-15)    # n-doped region = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionAcceptor)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionDonor)  # outer right boundary
    bfacemask!(grid, [h_pdoping], [h_pdoping], bregionJunction1) # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], bregionJunction2) # second inner interface

    if Plotter !== nothing
        vis = GridVisualizer(; Plotter, layout = (3, 3), size = (1550, 800))
        gridplot!(vis[1, 1], grid; Plotter, legend = :lt, title = "Grid", xlabel = L"\text{space [m]}", show = true)
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

    ## Define the physical data.
    Ec = 1.424 * eV
    Ev = 0.0 * eV
    Nc = 4.35195989587969e17 / (cm^3)
    Nv = 9.139615903601645e18 / (cm^3)
    mun = 8500.0 * (cm^2) / (V * s)
    mup = 400.0 * (cm^2) / (V * s)
    εr = 12.9 * 1.0              # relative dielectric permittivity of GAs
    T = 300.0 * K

    ## recombination parameters
    SRH_TrapDensity_n = 4.760185435081902e5 / cm^3
    SRH_TrapDensity_p = 9.996936448738406e6 / cm^3
    SRH_LifeTime = 1.0 * ps

    ## contact voltage
    voltageAcceptor = 1.4 * V

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

    ## Possible choices: Stationary, Transient
    data.modelType = Stationary

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor] = OhmicContact

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    # Define the Params and ParamsNodal struct.
    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)
    paramsnodal = ParamsNodal(grid, numberOfCarriers)

    params.temperature = T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ireg in 1:numberOfRegions           # region data

        params.dielectricConstant[ireg] = εr * constants.ε_0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = mun
        params.mobility[iphip, ireg] = mup

        ## recombination parameters
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity_n
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity_p

    end

    ## initialize the space dependent doping (see FarrellPeschka2019, Computers & Mathematics with Applications, 2019).
    NDoping = 1.0e17 / cm^3
    κ = 500.0
    for icoord in 1:numberOfNodes
        paramsnodal.doping[icoord] = NDoping * 0.5 * (1.0 + tanh((0.1 - coord[icoord] / μm) * κ) - (1.0 + tanh((coord[icoord] / μm - 0.2) * κ)))
    end

    data.params = params
    data.paramsnodal = paramsnodal

    ctsys = System(grid, data, unknown_storage = unknown_storage)

    if test == false
        println("*** done\n")
    end

    if test == false
        show_params(ctsys)
    end

    if Plotter !== nothing
        ################################################################################
        println("Plot doping")
        ################################################################################
        plot_doping!(vis[1, 2], grid, paramsnodal)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.abstol = 1.0e-14
    control.reltol = 1.0e-14
    control.max_round = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    ## calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    if Plotter !== nothing
        ## set legend for plotting routines. Either you can use the predefined labels or write your own.
        label_solution, label_density, label_energy = set_plotting_labels(data)

        plot_energies!(vis[1, 3], ctsys, solution, "Equilibrium", label_energy)
        plot_densities!(vis[2, 1], ctsys, solution, "Equilibrium", label_density)
        plot_solution!(vis[2, 2], ctsys, solution, "Equilibrium", label_solution)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    maxBias = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 41)
    IV = zeros(0)

    for Δu in biasValues

        if test == false
            println("bias value: Δu  = ", Δu, " V")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

        ## get IV curve
        factory = TestFunctionFactory(ctsys)

        ## testfunction zero in bregionAcceptor and one in bregionDonor
        tf = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I = integrate(ctsys, tf, solution)

        push!(IV, abs.(w_device * z_device * (I[iphin] + I[iphip])))

    end # bias loop

    if test == false
        println("*** done\n")
    end

    if Plotter !== nothing # plot solution and IV curve

        plot_energies!(vis[2, 3], ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_energy)
        plot_densities!(vis[3, 2], ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_density, plotGridpoints = true)
        plot_solution!(vis[3, 1], ctsys, solution, "Applied voltage Δu = $(biasValues[end])", label_solution, plotGridpoints = true)
        plot_IV!(vis[3, 3], biasValues, IV, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)

        reveal(vis)
    end

    testval = solution[15]
    return testval

end #  main

function test()
    testval = 1.4676876548796856
    return main(test = true, unknown_storage = :dense) ≈ testval && main(test = true, unknown_storage = :sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module is successfully recompiled.")
end

end # module
