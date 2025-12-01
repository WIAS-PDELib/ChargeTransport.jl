#=
# Nondimensionalized perovskite solar cell.
([source code](@__SOURCE_URL__))

Simulation of charge transport in a non-dimensionalized three-layer perovskite solar cell.
All physical parameters and constants are set to 1 besides the energies which are set to zero. This example serves an academic purpose.
The main() method is doing the study for one specific prefactor G0 of the photogeneration, while
the GenerationStudy() method does the parameter study with respect to G0.
For the underlying manuscript

``Existence of solutions and uniform bounds for the stationary semiconductor equations with generation and ionic carriers``,
    by D. Abdel, A. Blaustein, M. Herda, C. Chainais-Hillairet, and J. Moatti.

the test cases are (for n = 80)
## A: enableIons = false, DirichletVal = 2.0
## B: enableIons = true,  DirichletVal = 1.0
=#

module Ex109_PSC_NonDimensional

using ChargeTransport
using VoronoiFVM
using ExtendableGrids  # grid initializer
# using PyPlot           # solution visualizer

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

iphin = 1              # electron quasi Fermi potential
iphip = 2              # hole quasi Fermi potential

# We define the physical data.
zn = -1; zp = 1; za = 1
En = 0.0; Ep = 0.0; Ea = 0.0  # set the energies to 0
Nn = 1.0; Np = 1.0; Na = 10.0 # set the effective DOS to 1
μn = 1.0; μp = 1.0; μa = 1.0  # set the mobilities to 1
T = 1.0                       # set temperature to 1

## recombination parameters
SRH_TrapDensity = 0.0
SRH_LifeTime = 1.0
Radiative = 1.0

## doping
Ca = 7.5

# you can set verbose also to true to display some solver information
function main(;
        n = 80,               # for number of nodes in each layer
        Cn = 10, Cp = 10,     # doping
        λ = 1.0,              # Debye length
        DirichletVal = 2.0,   # Dirichlet value
        G0 = 1.0,             # photogeneration prefactor
        enableIons = false,   # present vacancies or not
        #################################
        parameterStudy = false,
        parseInival = false, inival = Array{Float64, 2},
        #################################
        Plotter = nothing, plotting = false,
        verbose = false, test = false
    )

    if !isnothing(Plotter) && nameof(Plotter) != :PyPlot
        @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global environment via the package manager and choose `Plotter = PyPlot`."
        plotting = false
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

    if enableIons && DirichletVal != 1.0
        @warn "Caution, initial value for the ions is only correct for `DirichletVal = 1.0` as the average density need to be equal to Ca.
        We adjusted the Dirichlet value to 1.0"
        DirichletVal = 1.0
    end

    if enableIons
        iphia = 3          # vacancies
        ipsi = 4           # electric potential
        numberOfCarriers = 3
    else
        ipsi = 3           # electric potential
        numberOfCarriers = 2
    end

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
    # Define the unity constants also in the discrete counterpart of the model
    data = Data(grid, numberOfCarriers, generationData = genData, constants = constants)

    ## Following variable declares, if we want to solve stationary or transient problem
    data.modelType = Stationary

    ## Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk,
    ## FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F .= FermiDiracOneHalfTeSCA

    data.boundaryType[1] = OhmicContact
    data.boundaryType[2] = OhmicContact

    if enableIons
        data.F[iphia] = FermiDiracMinusOne
        enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [region2])
    end

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
    if enableIons
        params.chargeNumbers[iphia] = za
    end

    for ireg in 1:numberOfRegions # region data

        params.dielectricConstant[ireg] = λ^2

        # the effective density of states and mobilities are set by default to one and the band-edge by default to zero.
        # This is why they do not necessarily need to be parsed here.
        # ## effective DOS, band-edge energy and mobilities
        # params.densityOfStates[iphin, ireg] = Nn
        # params.densityOfStates[iphip, ireg] = Np
        # params.bandEdgeEnergy[iphin, ireg] = En
        # params.bandEdgeEnergy[iphip, ireg] = Ep
        # params.mobility[iphin, ireg] = μn
        # params.mobility[iphip, ireg] = μp

        ## recombination parameters
        params.recombinationRadiative[ireg] = Radiative
        params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity

    end

    ## doping
    params.doping[iphin, region1] = Cn
    params.doping[iphip, region3] = Cp

    ## vacancy parameters
    if enableIons
        params.densityOfStates[iphia, region2] = Na
        params.bandEdgeEnergy[iphia, region2] = Ea
        params.mobility[iphia, region2] = μa
        params.doping[iphia, region2] = Ca
    end

    # Region dependent params is now a substruct of data which is again a substruct of the
    # system and will be parsed in next step.
    data.params = params

    # In the last step, we initialize our system with previous data which is likewise
    # dependent on the parameters.
    ctsys = System(grid, data, unknown_storage = :sparse)

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
    control.abstol = 1.0e-6
    control.reltol = 1.0e-6
    control.tol_round = 1.0e-6
    control.damp_initial = 0.5
    control.damp_growth = 1.61
    control.max_round = 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Solving the nonlinear system of equations")
    end
    ################################################################################

    ## calculate equilibrium solution and as initial guess
    inival2 = equilibrium_solve!(ctsys, control = control)

    if parseInival
        inival2 = inival
    end

    ## we need to fix the average vacancy density
    if enableIons
        if G0 < 0.1
            inival2[iphia, :] .= 0.7
        elseif G0 == 0.1
            inival2[iphia, :] .= 0.705
        elseif G0 == 0.5
            inival2[iphia, :] .= 0.729
        elseif G0 == 1.0
            inival2[iphia, :] .= 0.757
        elseif G0 == 5.0e0
            inival2[iphia, :] .= 0.955
        elseif G0 == 1.0e1
            inival2[iphia, :] .= 1.143
        elseif G0 == 5.0e1
            inival2[iphia, :] .= 1.842
        elseif G0 == 1.0e2
            inival2[iphia, :] .= 2.209
        elseif G0 == 5.0e2
            inival2[iphia, :] .= 3.08
        elseif G0 == 1.0e3
            inival2[iphia, :] .= 3.435
        end

    end

    set_contact!(ctsys, 2, Δu = DirichletVal)
    sol = ChargeTransport.solve(ctsys, inival = inival2, control = control)

    if enableIons && test == false
        integral = integrated_density(ctsys, sol = sol, icc = iphia, ireg = region2)
        mOmega = data.regionVolumes[region2]

        println("Given vacancy density is: ", Ca)
        println("Calculated average vacancy density is: ", integral / mOmega)
    end

    if test == false
        println("*** done\n")
    end

    if parameterStudy
        return coord, sol
    end

    ################################################################################
    if test == false && plotting
        println("Plotting")
    end
    ################################################################################

    nn = Nn .* data.F[iphin].(zn * (sol[iphin, :] - sol[ipsi, :]))
    np = Np .* data.F[iphip].(zp * (sol[iphip, :] - sol[ipsi, :]))
    if enableIons
        na = Na .* data.F[iphia].(za * (sol[iphia, :] - sol[ipsi, :]))
    end

    if plotting

        Plotter.figure()
        Plotter.plot(coord, zn .* sol[iphin, :], color = "green", linewidth = 5, label = "\$ v_{\\mathrm{n}}}\$")
        Plotter.plot(coord, zp .* sol[iphip, :], color = "red", linewidth = 5, linestyle = "--", label = "\$ v_{\\mathrm{p}}}\$")
        Plotter.plot(coord, sol[ipsi, :], color = "blue", linestyle = ":", linewidth = 5, label = "\$ \\psi\$")
        if enableIons
            Plotter.plot(coord[n:5n], za .* sol[iphia, n:5n], color = "gold", linewidth = 5, linestyle = "--", label = "\$ v_{\\mathrm{a}}}\$")
        end
        Plotter.axvspan(0.0, 1.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.xlim(0.0, 7.0)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Potential", fontsize = 17)
        Plotter.tight_layout()
        Plotter.legend(loc = "center right", fontsize = 14)
        ########################################################

        Plotter.figure()
        Plotter.semilogy(coord, nn, color = "green", linewidth = 5, label = "\$ n_{\\mathrm{n}}}\$")
        Plotter.semilogy(coord, np, color = "red", linewidth = 5, label = "\$ n_{\\mathrm{p}}}\$")
        if enableIons
            Plotter.semilogy(coord[n:5n], na[n:5n], color = "gold", linewidth = 5, label = "\$ n_{\\mathrm{a}}}\$")
        end
        Plotter.legend(loc = "center right", fontsize = 14)
        Plotter.xlim(0.0, 7.0)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Density", fontsize = 17)
        Plotter.tight_layout()

    end

    return sum(filter(!isnan, sol)) / length(sol)

    if test == false && plotting
        println("*** done\n")
    end

end #  main


#######################################################################
#######################################################################

# For the underlying manuscript, the test cases are (for n = 80)
## A: enableIons = false, DirichletVal = 2.0
## B: enableIons = true,  DirichletVal = 1.0

function GenerationStudy(;
        n = 80,               # for number of nodes in each layer
        Cn = 10, Cp = 10,     # doping
        λ = 1.0,              # Debye length
        DirichletVal = 2.0,   # Dirichlet value
        enableIons = false,   # enabling ions
        Plotter = nothing, plotting = false,
        verbose = false
    )

    if plotting
        if isnothing(Plotter)
            @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global environment via the package manager and choose `Plotter = PyPlot`."
            plotting = false
        else
            if nameof(Plotter) != :PyPlot
                @warn "We need PyPlot as Plotter for this example. Please add PyPlot to your global environment via the package manager and choose `Plotter = PyPlot`."
                plotting = false
            end
        end
    end

    Plotter.rc("font", family = "sans-serif", size = 14)
    Plotter.rc("mathtext", fontset = "dejavusans")
    Plotter.close("all")

    if enableIons && DirichletVal != 1.0
        @warn "Caution, initial value for the ions is only correct for `DirichletVal = 1.0` as the average density need to be equal to Ca.
        We adjusted the Dirichlet value to 1.0"
        DirichletVal = 1.0
    end

    if enableIons
        iphia = 3          # vacancies
        ipsi = 4           # electric potential
        numberOfCarriers = 3
    else
        ipsi = 3           # electric potential
        numberOfCarriers = 2
    end

    G0Vec = [1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2, 1.0e-1, 5.0e-1, 1.0e0, 5.0e0, 1.0e1, 5.0e1, 1.0e2, 5.0e2, 1.0e3]

    soleM1 = Array{Float64, 2}; sole0 = Array{Float64, 2}
    sole1 = Array{Float64, 2}; sole2 = Array{Float64, 2}
    ####
    nnmaxVec = zeros(0); npmaxVec = zeros(0); namaxVec = zeros(0)
    phinmaxVec = zeros(0); phipmaxVec = zeros(0)
    phiamaxVec = zeros(0); psimaxVec = zeros(0)

    ##### first one
    G0 = G0Vec[1]
    println("G0 = ", G0)

    coord, sol = main(n = n, Cn = Cn, Cp = Cp, λ = λ, DirichletVal = DirichletVal, G0 = G0, enableIons = enableIons, parameterStudy = true, test = true, verbose = verbose)

    nn = FermiDiracOneHalfTeSCA.(zn * (sol[iphin, :] - sol[ipsi, :]))
    np = FermiDiracOneHalfTeSCA.(zp * (sol[iphip, :] - sol[ipsi, :]))

    ####
    nnmax = maximum(abs.(nn))
    npmax = maximum(abs.(np))
    phinmax = maximum(abs.(sol[iphin, :]))
    phipmax = maximum(abs.(sol[iphip, :]))
    psimax = maximum(abs.(sol[ipsi, :]))

    #####################################
    push!(nnmaxVec, nnmax)
    push!(npmaxVec, npmax)
    push!(phinmaxVec, phinmax)
    push!(phipmaxVec, phipmax)
    push!(psimaxVec, psimax)

    if enableIons
        na = filter(!isnan, Na .* FermiDiracMinusOne.(za * (sol[iphia, :] - sol[ipsi, :])))
        namax = maximum(abs.(na))
        phiamax = maximum(abs.(filter(!isnan, sol[iphia, :])))

        push!(namaxVec, namax)
        push!(phiamaxVec, phiamax)
    end

    inival = sol

    ######################################

    for G0 in G0Vec[2:end]

        println("G0 = ", G0)

        coord, sol = main(n = n, Cn = Cn, Cp = Cp, λ = λ, DirichletVal = DirichletVal, G0 = G0, enableIons = enableIons, parameterStudy = true, parseInival = true, inival = inival, test = true, verbose = verbose)

        nn = FermiDiracOneHalfTeSCA.(zn * (sol[iphin, :] - sol[ipsi, :]))
        np = FermiDiracOneHalfTeSCA.(zp * (sol[iphip, :] - sol[ipsi, :]))

        ####
        nnmax = maximum(abs.(nn))
        npmax = maximum(abs.(np))
        phinmax = maximum(abs.(sol[iphin, :]))
        phipmax = maximum(abs.(sol[iphip, :]))
        psimax = maximum(abs.(sol[ipsi, :]))

        #####################################
        push!(nnmaxVec, nnmax)
        push!(npmaxVec, npmax)
        push!(phinmaxVec, phinmax)
        push!(phipmaxVec, phipmax)
        push!(psimaxVec, psimax)

        if enableIons
            na = filter(!isnan, Na .* FermiDiracMinusOne.(za * (sol[iphia, :] - sol[ipsi, :])))
            namax = maximum(abs.(na))
            phiamax = maximum(abs.(filter(!isnan, sol[iphia, :])))

            push!(namaxVec, namax)
            push!(phiamaxVec, phiamax)
        end

        inival = sol

        if G0 == 1.0e-1
            soleM1 = copy(sol)
        elseif G0 == 1.0e0
            sole0 = copy(sol)
        elseif G0 == 1.0e1
            sole1 = copy(sol)
        elseif G0 == 1.0e2
            sole2 = copy(sol)
        end
    end

    ###################################
    if plotting
        size = 12
        if enableIons
            Plotter.loglog(G0Vec, Ca .* ones(length(G0Vec)), color = "gray", linestyle = ":", linewidth = 4, label = "\$  M_{\\mathrm{a}} \$ ")
        end
        Plotter.loglog(G0Vec, nnmaxVec, marker = "o", markersize = size, linewidth = 5, color = "darkgreen", label = "\$  || n_{\\mathrm{n}} ||_{\\infty} \$ ")
        #####
        Plotter.loglog(G0Vec, npmaxVec, marker = "o", linewidth = 5, linestyle = "--", markersize = size, color = "darkred", label = "\$  || n_{\\mathrm{p}} ||_{\\infty} \$ ")
        if enableIons
            Plotter.loglog(G0Vec, namaxVec, marker = "o", linewidth = 5, markersize = size, color = "gold", label = "\$  || n_{\\mathrm{a}} ||_{\\infty} \$ ")
        end
        Plotter.legend(fontsize = 15)
        Plotter.xlabel(" \$ G_0 \$ ", fontsize = 17)
        Plotter.ylabel("\$ L^{\\infty} \$ norm", fontsize = 17)
        Plotter.title("Cn = $Cn; Cp = $Cp; BC = $DirichletVal, \$ \\lambda = \$ $λ")
        Plotter.xlim(7.0e-4, 3.0e3)
        Plotter.ylim(6.0e-1, 30.0)
        Plotter.yticks([1.0e0, 1.0e1, 2.0e1])
        Plotter.tight_layout()

        ###################################
        Plotter.figure()
        Plotter.loglog(G0Vec, psimaxVec, marker = "o", markersize = size, color = "darkblue", linestyle = ":", linewidth = 5, label = "\$  ||\\psi||_{\\infty} \$ ")
        Plotter.loglog(G0Vec, phinmaxVec, marker = "o", markersize = size, color = "darkgreen", linewidth = 5, label = "\$  ||v_{\\mathrm{n}}||_{\\infty} \$ ")
        Plotter.loglog(G0Vec, phipmaxVec, marker = "o", markersize = size, linestyle = "--", color = "darkred", linewidth = 5, label = "\$  ||v_{\\mathrm{p}}||_{\\infty} \$ ")
        if enableIons
            Plotter.loglog(G0Vec, phiamaxVec, marker = "o", markersize = size, color = "gold", linewidth = 5, label = "\$  ||v_{\\mathrm{a}}||_{\\infty} \$ ")
        end

        Plotter.xlim(7.0e-4, 3.0e3)
        Plotter.ylim(6.0e-1, 30.0)
        Plotter.yticks([1.0e0, 1.0e1, 2.0e1])
        Plotter.legend(fontsize = 15)
        Plotter.xlabel(" \$ G_0 \$ ", fontsize = 17)
        Plotter.ylabel("\$ L^{\\infty} \$ norm", fontsize = 17)
        Plotter.title("Cn = $Cn; Cp = $Cp; BC = $DirichletVal, \$ \\lambda = \$ $λ")
        Plotter.tight_layout()

        ###################################
        Plotter.figure()
        Blues = Plotter.get_cmap(:Blues)
        Oranges = Plotter.get_cmap(:Oranges)
        Greens = Plotter.get_cmap(:Greens)
        Wistia = Plotter.get_cmap(:Wistia)

        if enableIons
            Plotter.plot(coord[n:5n]', za .* soleM1[iphia, n:5n], linewidth = 5, color = Wistia(201))
            Plotter.plot(coord[n:5n]', za .* sole0[iphia, n:5n], linewidth = 5, color = Wistia(171))
            Plotter.plot(coord[n:5n]', za .* sole1[iphia, n:5n], linewidth = 5, color = Wistia(131))
            Plotter.plot(coord[n:5n]', za .* sole2[iphia, n:5n], linewidth = 5, color = Wistia(101))
        end

        Plotter.plot(coord', zn .* soleM1[iphin, :], linewidth = 5, color = Greens(241))
        Plotter.plot(coord', zn .* sole0[iphin, :], linewidth = 5, color = Greens(201))
        Plotter.plot(coord', zn .* sole1[iphin, :], linewidth = 5, color = Greens(161))
        Plotter.plot(coord', zn .* sole2[iphin, :], linewidth = 5, color = Greens(121))


        Plotter.xlim(0.0, 7.0)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.ylim(-5.8, 7.5)
        Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Potential", fontsize = 17)
        Plotter.tight_layout()

        ###################################
        Plotter.figure()
        Plotter.plot(coord', zp .* soleM1[iphip, :], linewidth = 5, color = Oranges(241))
        Plotter.plot(coord', zp .* sole0[iphip, :], linewidth = 5, color = Oranges(201))
        Plotter.plot(coord', zp .* sole1[iphip, :], linewidth = 5, color = Oranges(161))
        Plotter.plot(coord', zp .* sole2[iphip, :], linewidth = 5, color = Oranges(121))
        Plotter.xlim(0.0, 7.0)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.ylim(-5.8, 7.5)
        Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Potential", fontsize = 17)
        Plotter.tight_layout()

        ###################################
        Plotter.figure()
        Plotter.plot(coord', soleM1[ipsi, :], linewidth = 5, color = Blues(241))
        Plotter.plot(coord', sole0[ipsi, :], linewidth = 5, color = Blues(201))
        Plotter.plot(coord', sole1[ipsi, :], linewidth = 5, color = Blues(161))
        Plotter.plot(coord', sole2[ipsi, :], linewidth = 5, color = Blues(121))
        Plotter.xlim(0.0, 7.0)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.ylim(-5.8, 7.5)
        Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Potential", fontsize = 17)
        Plotter.tight_layout()

        ###################################
        nneM1 = FermiDiracOneHalfTeSCA.(zn * (soleM1[iphin, :] - soleM1[ipsi, :]))
        npeM1 = FermiDiracOneHalfTeSCA.(zp * (soleM1[iphip, :] - soleM1[ipsi, :]))

        nne0 = FermiDiracOneHalfTeSCA.(zn * (sole0[iphin, :] - sole0[ipsi, :]))
        npe0 = FermiDiracOneHalfTeSCA.(zp * (sole0[iphip, :] - sole0[ipsi, :]))

        nne1 = FermiDiracOneHalfTeSCA.(zn * (sole1[iphin, :] - sole1[ipsi, :]))
        npe1 = FermiDiracOneHalfTeSCA.(zp * (sole1[iphip, :] - sole1[ipsi, :]))

        nne2 = FermiDiracOneHalfTeSCA.(zn * (sole2[iphin, :] - sole2[ipsi, :]))
        npe2 = FermiDiracOneHalfTeSCA.(zp * (sole2[iphip, :] - sole2[ipsi, :]))

        if enableIons
            naeM1 = Na .* FermiDiracMinusOne.(za * (soleM1[iphia, :] - soleM1[ipsi, :]))
            nae0 = Na .* FermiDiracMinusOne.(za * (sole0[iphia, :] - sole0[ipsi, :]))
            nae1 = Na .* FermiDiracMinusOne.(za * (sole1[iphia, :] - sole1[ipsi, :]))
            nae2 = Na .* FermiDiracMinusOne.(za * (sole2[iphia, :] - sole2[ipsi, :]))
        end

        Plotter.figure()
        Plotter.semilogy(coord', nneM1, linewidth = 5, color = Greens(241))
        Plotter.semilogy(coord', nne0, linewidth = 5, color = Greens(201))
        Plotter.semilogy(coord', nne1, linewidth = 5, color = Greens(161))
        Plotter.semilogy(coord', nne2, linewidth = 5, color = Greens(121))
        #################
        Plotter.semilogy(coord', npeM1, linewidth = 5, color = Oranges(241))
        Plotter.semilogy(coord', npe0, linewidth = 5, color = Oranges(201))
        Plotter.semilogy(coord', npe1, linewidth = 5, color = Oranges(161))
        Plotter.semilogy(coord', npe2, linewidth = 5, color = Oranges(121))
        ##########
        Plotter.xlim(0.0, 7.0)
        Plotter.ylim(3.0e-3, 1.5e1)
        Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
        Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
        Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
        Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
        Plotter.xlabel("\$ x \$", fontsize = 17)
        Plotter.ylabel("Density", fontsize = 17)
        Plotter.tight_layout()

        ###################################
        if enableIons
            Plotter.figure()
            Plotter.semilogy(coord', naeM1, linewidth = 5, color = Wistia(201))
            Plotter.semilogy(coord', nae0, linewidth = 5, color = Wistia(171))
            Plotter.semilogy(coord', nae1, linewidth = 5, color = Wistia(131))
            Plotter.semilogy(coord', nae2, linewidth = 5, color = Wistia(101))
            Plotter.xlim(0.0, 7.0)
            Plotter.yticks([1.0e0, 1.0e1])
            Plotter.xticks([0.0, 1.0, 3.0, 5.0, 7.0])
            Plotter.axvspan(0.0, 1.0, facecolor = [211 / 255 232 / 255 208 / 255])
            Plotter.axvspan(1.0, 5.0, facecolor = [233 / 255 226 / 255 215 / 255])
            Plotter.axvspan(5.0, 7.0, facecolor = [243 / 255 192 / 255 192 / 255])
            Plotter.xlabel("\$ x \$", fontsize = 17)
            Plotter.ylabel("Density", fontsize = 17)
            Plotter.tight_layout()
        end

    end


    return nothing

end

function test()
    testval = 0.9289261210695825
    return main(test = true) ≈ testval
end


end # module
