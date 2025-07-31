#=
# Nondimensionalized perovskite solar cell.
([source code](@__SOURCE_URL__))

Simulation of charge transport in a non-dimensionalized three-layer perovskite solar cell.
All physical parameters and constants are set to 1 besides the energies which are set to zero.
This example serves an academic purpose.
=#

module Ex109_PSC_NonDimensional

using ChargeTransport
using VoronoiFVM
using ExtendableGrids  # grid initializer
using PyPlot           # solution visualizer

# you can set verbose also to true to display some solver information
function main(;
        n = 3,               # for number of nodes in each layer
        Cn = 10, Cp = 10,    # doping
        λ = 1.0,             # Debye length
        DirichletVal = 1.0,  # Dirichlet value
        G0 = 1.0,            # photogeneration prefactor
        Plotter = PyPlot, plotting = false,
        verbose = false, test = false, unknown_storage = :sparse
    )

    if plotting
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # constants
    # Here, we set the constants to unity_constants.
    # In particular, we set: q = k_B = ε_0 = 1
    # When defining `ChargeTransport.constants" the constants based on CODATA2022 are used.
    constants = ChargeTransport.unity_constants

    # region numbers
    region1 = 1
    region2 = 2
    region3 = 3
    regions = [region1, region2, region3]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregion1 = 1
    bregion2 = 2

    ## grid
    h1 = 1.0; h2 = 4.0; h3 = 2.0
    h_total = h1 + h2 + h3

    coord1 = collect(range(0.0; stop = h1, length = n))
    coord2 = collect(range(h1; stop = h1 + h2, length = 4 * n))
    coord3 = collect(range(h1 + h2; stop = h_total, length = 2 * n))
    coord = glue(coord1, coord2)
    coord = glue(coord, coord3)

    grid = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0], [h1], region1)
    cellmask!(grid, [h1], [h1 + h2], region2)
    cellmask!(grid, [h1 + h2], [h_total], region3)

    ## bfacemask! for setting different boundary regions.
    bfacemask!(grid, [0.0], [0.0], bregion1)
    bfacemask!(grid, [h_total], [h_total], bregion2)

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
    λ = λ                  # Debye length, entering model as `DielectricConstant`, i.e., prefactor in the displacement flux.
    T = 1.0                # set temperature to 1

    ## recombination parameters
    SRH_TrapDensity = 0.0
    SRH_LifeTime = 1.0
    Radiative = 1.0

    ## doping
    Cn = Cn
    Cp = Cp
    Ca = 0.0

    # boundary value
    DirichletVal = DirichletVal

    # photogeneration
    G(x) = G0 .* exp.(- (x .- h1))
    genData = zeros(length(coord))
    genData[length(coord1):(length(coord1) + length(coord2) - 1)] = G.(coord2)

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
    data.F .= Boltzmann

    ## The desired recombination processes can be chosen here.
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
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

        params.dielectricConstant[ireg] = λ^2

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn
        params.densityOfStates[iphip, ireg] = Np
        params.bandEdgeEnergy[iphin, ireg] = En
        params.bandEdgeEnergy[iphip, ireg] = Ep
        params.mobility[iphin, ireg] = μn
        params.mobility[iphip, ireg] = μp

        ## recombination parameters
        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity

    end

    ## doping
    params.doping[iphin, region1] = Cn
    params.doping[iphin, region2] = Ca
    params.doping[iphip, region3] = Cp

    # Region dependent params is now a substruct of data which is again a substruct of the
    # system and will be parsed in next step.
    data.params = params

    # In the last step, we initialize our system with previous data which is likewise
    # dependent on the parameters.
    ctsys = System(grid, data, unknown_storage = unknown_storage)

    ## boundary model
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphin, bregion1, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphip, bregion1, 0.0)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, ipsi, bregion1, asinh(Cn / 2))

    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphin, bregion2, DirichletVal)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, iphip, bregion2, DirichletVal)
    VoronoiFVM.boundary_dirichlet!(ctsys.fvmsys, ipsi, bregion2, asinh(-Cp / 2) + DirichletVal)

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
    control.damp_initial = 0.5
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
    inival[iphin, :] = 0.0 .+ DirichletVal ./ h_total .* coord
    inival[iphip, :] = 0.0 .+ DirichletVal ./ h_total .* coord
    inival[ipsi, :] = asinh(Cn / 2) .+ ((asinh(-Cp / 2) + DirichletVal) - asinh(Cn / 2)) ./ h_total .* coord

    sol = ChargeTransport.solve(ctsys, inival = inival, control = control)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false && plotting
        println("Plotting")
    end
    ################################################################################

    nn = exp.(zn * (sol[iphin, :] - sol[ipsi, :]))
    np = exp.(zp * (sol[iphip, :] - sol[ipsi, :]))

    if plotting

        figure()
        PyPlot.plot(coord, sol[iphin, :], color = "green", linewidth = 5, label = "\$ \\varphi_{\\mathrm{n}}}\$")
        PyPlot.plot(coord, sol[iphip, :], color = "red", linewidth = 5, linestyle = "--", label = "\$ \\varphi_{\\mathrm{p}}}\$")
        PyPlot.plot(coord, sol[ipsi, :], color = "blue", linewidth = 5, label = "\$ \\psi\$")
        axvspan(0.0, 1.0, facecolor = [243 / 255 192 / 255 192 / 255])
        axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        axvspan(5.0, 7.0, facecolor = [211 / 255 232 / 255 208 / 255])
        xlabel("\$ x \$", fontsize = 17)
        ylabel("Potential", fontsize = 17)
        tight_layout()
        PyPlot.legend(loc = "center right", fontsize = 14)
        ########################################################

        figure()
        PyPlot.semilogy(coord, nn, color = "green", linewidth = 5, label = "\$ n_{\\mathrm{n}}}\$")
        PyPlot.semilogy(coord, np, color = "red", linewidth = 5, label = "\$ n_{\\mathrm{p}}}\$")
        PyPlot.legend(loc = "center right", fontsize = 14)
        axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        xlabel("\$ x \$", fontsize = 17)
        ylabel("Density", fontsize = 17)
        tight_layout()

    end

    if test == false && plotting
        println("*** done\n")
    end

    return sum(filter(!isnan, sol)) / length(sol)

end #  main

function test()
    testval = 0.4760637366166469
    return main(test = true, unknown_storage = :dense) ≈ testval && main(test = true, unknown_storage = :sparse) ≈ testval
end


end # module
