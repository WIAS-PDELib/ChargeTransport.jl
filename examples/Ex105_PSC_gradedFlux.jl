#=
# Three-layer PSC device with graded interfaces & Ohmic contacts (1D).
([source code](@__SOURCE_URL__))

Simulating a three layer PSC device Ti02| MAPI | spiro-OMeTAD without mobile ions. The simulations are
performed out of equilibrium, stationary and with two junctions between perovskite layer and
transport layers, to which we refer as graded interfaces. Hence, a graded flux discretization
with space dependent band-edge energies and density of states is tested here.
The difference here is that we adjusted the order of indexing the quasi Fermi potentials.

The parameters are based on the default parameter set of Ionmonger (with minor adjustments),
such that we can likewise compare with the software Driftfusion, see
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv
=#

module Ex105_PSC_gradedFlux

using ChargeTransport
using ExtendableGrids
using GridVisualize
using LaTeXStrings

## function for grading the physical parameters
function grading_parameter!(physicalParameter, coord, regionTransportLayers, regionJunctions, h, heightLayers, lengthLayers, values)
    for ireg in regionTransportLayers

        xcoord = lengthLayers[ireg]:lengthLayers[ireg + 1]
        physicalParameter[xcoord] .= values[ireg]

    end

    for ireg in regionJunctions

        xcoord = lengthLayers[ireg]:lengthLayers[ireg + 1]
        left = lengthLayers[ireg] - 3
        junction = h[ireg]
        right = lengthLayers[ireg + 2] - 3

        gradient = (physicalParameter[right] - physicalParameter[left]) / junction

        for index in xcoord
            physicalParameter[index] = physicalParameter[left] + (coord[index] - heightLayers[ireg - 1]) * gradient
        end

    end

    return physicalParameter
end

# supported Plotters are GLMakie and PythonPlot
# you can set verbose also to true to display some solver information
function main(; n = 2, Plotter = nothing, verbose = false, test = false, unknown_storage = :sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    @local_unitfactors μm cm s ns V K ps Hz W m

    constants = ChargeTransport.constants

    eV = constants.q * V

    ## region numbers
    regionDonor = 1          # n doped region
    regionJunction1 = 2
    regionIntrinsic = 3          # intrinsic region
    regionJunction2 = 4
    regionAcceptor = 5          # p doped region
    regions = [regionDonor, regionJunction1, regionIntrinsic, regionJunction2, regionAcceptor]
    regionTransportLayers = [regionDonor, regionIntrinsic, regionAcceptor]
    regionJunctions = [regionJunction1, regionJunction2]
    numberOfRegions = length(regions)

    ## boundary region numbers
    bregionDonor = 1
    bregionAcceptor = 2
    bregionDJ1 = 3
    bregionJ1I = 4
    bregionIJ2 = 5
    bregionJ2A = 6

    ## grid
    h_ndoping = 9.9e-6 * cm
    h_junction1 = 1.0e-7 * cm
    h_intrinsic = 4.0e-5 * cm
    h_junction2 = 1.0e-7 * cm
    h_pdoping = 1.99e-5 * cm
    h_total = h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping
    h = [h_ndoping, h_junction1, h_intrinsic, h_junction2, h_pdoping]
    heightLayers = [
        h_ndoping,
        h_ndoping + h_junction1,
        h_ndoping + h_junction1 + h_intrinsic,
        h_ndoping + h_junction1 + h_intrinsic + h_junction2,
        h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping,
    ]
    refinementfactor = 2^(n - 1)

    coord_ndoping = collect(range(0.0, stop = h_ndoping, length = 4 * refinementfactor))
    length_n = length(coord_ndoping)
    coord_junction1 = collect(
        range(
            h_ndoping,
            stop = h_ndoping + h_junction1,
            length = 3 * refinementfactor
        )
    )
    coord_intrinsic = collect(
        range(
            h_ndoping + h_junction1,
            stop = (h_ndoping + h_junction1 + h_intrinsic),
            length = 10 * refinementfactor
        )
    )
    coord_junction2 = collect(
        range(
            h_ndoping + h_junction1 + h_intrinsic,
            stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2),
            length = 3 * refinementfactor
        )
    )
    coord_pdoping = collect(
        range(
            (h_ndoping + h_junction1 + h_intrinsic + h_junction2),
            stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping),
            length = 4 * refinementfactor
        )
    )

    coord = glue(coord_ndoping, coord_junction1)
    length_j1 = length(coord)
    coord = glue(coord, coord_intrinsic)
    length_i = length(coord)
    coord = glue(coord, coord_junction2)
    length_j2 = length(coord)
    coord = glue(coord, coord_pdoping)

    grid = simplexgrid(coord)
    numberOfNodes = length(coord)
    lengthLayers = [1, length_n, length_j1, length_i, length_j2, numberOfNodes]

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionDonor)      # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionJunction1)  # first junction   = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionIntrinsic)  # intrinsic region = 3
    cellmask!(grid, [heightLayers[3]], [heightLayers[4]], regionJunction2)  # sec. junction    = 4
    cellmask!(grid, [heightLayers[4]], [heightLayers[5]], regionAcceptor)   # p-doped region   = 5

    # inner interfaces
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionDJ1)
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ1I)
    bfacemask!(grid, [heightLayers[3]], [heightLayers[3]], bregionIJ2)
    bfacemask!(grid, [heightLayers[4]], [heightLayers[4]], bregionJ2A)

    if Plotter !== nothing
        vis = GridVisualizer(; Plotter, layout = (4, 2), size = (1550, 800))
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
    iphin = 2 # electron quasi Fermi potential
    iphip = 1 # hole quasi Fermi potential
    numberOfCarriers = 2

    ## temperature
    T = 300.0 * K

    ## band edge energies
    Ec_d = -4.0 * eV
    Ev_d = -6.0 * eV

    Ec_i = -3.7 * eV
    Ev_i = -5.4 * eV

    Ec_a = -3.1 * eV
    Ev_a = -5.1 * eV

    ## these parameters at the junctions for E_\alpha and N_\alpha will be overwritten.
    Ec_j1 = Ec_d;     Ec_j2 = Ec_i
    Ev_j1 = Ev_d;     Ev_j2 = Ev_i

    EC = [Ec_d, Ec_j1, Ec_i, Ec_j2, Ec_a]
    EV = [Ev_d, Ev_j1, Ev_i, Ev_j2, Ev_a]

    ## effective densities of state
    Nc_d = 5.0e19 / (cm^3)
    Nv_d = 5.0e19 / (cm^3)

    Nc_i = 8.1e18 / (cm^3)
    Nv_i = 5.8e18 / (cm^3)

    Nc_a = 5.0e19 / (cm^3)
    Nv_a = 5.0e19 / (cm^3)

    Nc_j1 = Nc_d;     Nc_j2 = Nc_i
    Nv_j1 = Nv_d;     Nv_j2 = Nv_i

    NC = [Nc_d, Nc_j1, Nc_i, Nc_j2, Nc_a]
    NV = [Nv_d, Nv_j1, Nv_i, Nv_j2, Nv_a]

    ## mobilities
    μn_d = 3.89 * (cm^2) / (V * s)
    μp_d = 3.89 * (cm^2) / (V * s)

    μn_i = 6.62e1 * (cm^2) / (V * s)
    μp_i = 6.62e1 * (cm^2) / (V * s)

    μn_a = 3.89e-1 * (cm^2) / (V * s)
    μp_a = 3.89e-1 * (cm^2) / (V * s)

    μn_j1 = μn_d;     μn_j2 = μn_i
    μp_j1 = μp_d;     μp_j2 = μp_i

    μn = [μn_d, μn_j1, μn_i, μn_j2, μn_a]
    μp = [μp_d, μp_j1, μp_i, μp_j2, μp_a]

    ## relative dielectric permittivity
    ε_d = 10.0 * 1.0
    ε_i = 24.1 * 1.0
    ε_a = 3.0 * 1.0

    ε_j1 = ε_d;       ε_j2 = ε_a

    ε = [ε_d, ε_j1, ε_i, ε_j2, ε_a]

    ## radiative recombination
    r0_d = 0.0e+0 * cm^3 / s
    r0_i = 1.0e-12 * cm^3 / s
    r0_a = 0.0e+0 * cm^3 / s

    r0_j1 = r0_i;      r0_j2 = r0_i

    r0 = [r0_d, r0_j1, r0_i, r0_j2, r0_a]

    ## life times and trap densities
    τn_d = 1.0e100 * s
    τp_d = 1.0e100 * s

    τn_i = 3.0e-10 * s
    τp_i = 3.0e-8 * s
    τn_a = τn_d
    τp_a = τp_d

    τn_j1 = τn_i;     τn_j2 = τn_a
    τp_j1 = τp_i;     τp_j2 = τp_a

    τn = [τn_d, τn_j1, τn_i, τn_j2, τn_a]
    τp = [τp_d, τp_j1, τp_i, τp_j2, τp_a]

    ## SRH trap energies (needed for calculation of trap_density (SRH))
    # Ei_d             = -5.0                 * eV
    # Ei_i             = -4.55                * eV
    # Ei_a             = -4.1                 * eV

    # Ei_j1            = Ei_d;      Ei_j2     = Ei_i

    # EI               = [Ei_d, Ei_j1, Ei_i, Ei_j2, Ei_a]
    ## reference densities
    nτ_d = 7.94e8 / m^3
    pτ_d = 7.94e8 / m^3

    nτ_i = 4.26e10 / m^3
    pτ_i = 3.05e10 / m^3
    nτ_a = nτ_d
    pτ_a = pτ_d

    nτ_j1 = nτ_i;     nτ_j2 = nτ_a
    pτ_j1 = pτ_i;     pτ_j2 = pτ_a

    nτ = [nτ_d, nτ_j1, nτ_i, nτ_j2, nτ_a]
    pτ = [pτ_d, pτ_j1, pτ_i, pτ_j2, pτ_a]

    ## Auger recombination
    Auger = 0.0

    ## doping (doping values are from Driftfusion)
    Nd = 1.03e18 / (cm^3)
    Na = 1.03e18 / (cm^3)
    Ni_acceptor = 8.32e7 / (cm^3)

    ## contact voltage
    voltageAcceptor = 1.2 * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    data = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType = Stationary

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionDonor] = OhmicContact
    data.boundaryType[bregionAcceptor] = OhmicContact

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    ## Default is ExcessChemicalPotential. Here, we need a graded variant
    data.fluxApproximation .= ScharfetterGummelGraded

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    ## for region dependent parameters
    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)
    ## for space dependent parameters
    paramsnodal = ParamsNodal(grid, numberOfCarriers)

    params.temperature = T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    ## nodal band-edge energies
    paramsnodal.bandEdgeEnergy[iphin, :] = grading_parameter!(
        paramsnodal.bandEdgeEnergy[iphin, :],
        coord, regionTransportLayers, regionJunctions, h,
        heightLayers, lengthLayers, EC
    )
    paramsnodal.bandEdgeEnergy[iphip, :] = grading_parameter!(
        paramsnodal.bandEdgeEnergy[iphip, :],
        coord, regionTransportLayers, regionJunctions, h,
        heightLayers, lengthLayers, EV
    )
    ## nodal effective density of states
    paramsnodal.densityOfStates[iphin, :] = grading_parameter!(
        paramsnodal.densityOfStates[iphin, :],
        coord, regionTransportLayers, regionJunctions, h,
        heightLayers, lengthLayers, NC
    )
    paramsnodal.densityOfStates[iphip, :] = grading_parameter!(
        paramsnodal.densityOfStates[iphip, :],
        coord, regionTransportLayers, regionJunctions, h,
        heightLayers, lengthLayers, NV
    )

    for ireg in 1:numberOfRegions  ## region dependent data

        ## mobility
        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]

        params.dielectricConstant[ireg] = ε[ireg] * constants.ε_0
        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = nτ[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg] = pτ[ireg]
        params.recombinationAuger[iphin, ireg] = Auger
        params.recombinationAuger[iphip, ireg] = Auger

    end

    ## doping
    params.doping[iphin, regionDonor] = Nd
    params.doping[iphip, regionIntrinsic] = Ni_acceptor
    params.doping[iphip, regionAcceptor] = Na

    data.params = params
    data.paramsnodal = paramsnodal
    ctsys = System(grid, data, unknown_storage = unknown_storage)

    ## print data
    if test == false
        show_params(ctsys)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 200
    control.abstol = 1.0e-13
    control.reltol = 1.0e-13
    control.tol_round = 1.0e-13
    control.damp_initial = 0.5
    control.damp_growth = 1.61 # >= 1
    control.max_round = 5

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    if Plotter !== nothing
        label_solution, label_density, label_energy = set_plotting_labels(data)

        plot_energies!(vis[1, 2], ctsys, solution, "Equilibrium", label_energy)
        plot_densities!(vis[2, 1], ctsys, solution, "Equilibrium", label_density)
        # TODO MO: phip wird nicht angezeigt, auch mit dem ursprünglichen Plotting, stimmt das so?
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
    biasValues = range(0, stop = maxBias, length = 13)

    for Δu in biasValues
        if test == false
            println("Bias value: Δu = $(Δu) V")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # bias loop

    ## plotting
    if Plotter !== nothing
        plot_energies!(vis[3, 1], ctsys, solution, "Applied voltage Δu = $maxBias", label_energy)
        plot_densities!(vis[3, 2], ctsys, solution, "Applied voltage Δu = $maxBias", label_density)
        plot_solution!(vis[4, 1], ctsys, solution, "Applied voltage Δu = $maxBias", label_solution)

        reveal(vis)
    end

    if test == false
        println("*** done\n")
    end

    testval = solution[data.index_psi, 20]
    return testval

end #  main

function test()
    testval = -3.982748467515117
    return main(test = true, unknown_storage = :dense) ≈ testval && main(test = true, unknown_storage = :sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
