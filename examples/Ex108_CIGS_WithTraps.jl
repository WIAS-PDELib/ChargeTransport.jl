#=
# CIGS: stationary with traps and Schottky contacts.
([source code](SOURCE_URL))

Simulating stationary charge transport for CIGS with hole traps and mixed Schottky/Ohmic
contact conditions. Assume that SRH recombination only happens within a small regime.
=#

module Ex108_CIGS_WithTraps

using ChargeTransport
using ExtendableGrids
using PyPlot

## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_pdoping_left, h_pdoping_trap, h_pdoing_right)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 2 * refinementfactor))
    coord_pdoping_left  = collect(range(h_ndoping, stop = (h_ndoping + h_pdoping_left), length = 3 * refinementfactor))
    coord_pdoping_plus  = collect(range((h_ndoping + h_pdoping_left),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        length =  refinementfactor))
    coord_pdoping_right = collect(range((h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoing_right),
                                        length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_pdoping_left)
    coord            = glue(coord, coord_pdoping_plus)
    coord            = glue(coord, coord_pdoping_right)

    return coord
end

# you can also use other Plotters, if you add them to the example file
function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, AdditionalTrapSpecies = false)

    if plotting
        Plotter.close("all")
    end
    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionDonor          = 1                           # n doped region
    regionAcceptorLeft   = 2                           # p doped region
    regionAcceptorTrap   = 3                           # p doped region with trap
    regionAcceptorRight  = 4                           # p doped region
    regions              = [regionDonor, regionAcceptorLeft, regionAcceptorTrap, regionAcceptorRight]
    numberOfRegions      = length(regions)

    ## boundary region numbers
    bregionDonor         = 1
    bregionAcceptor      = 2
    bregionDALeft        = 3
    bregionALeftATrap    = 4
    bregionATrapARight   = 5

    ## grid
    refinementfactor     = 2^(n-1)
    h_ndoping            = 0.5    * μm
    h_pdoping_left       = 1.0    * μm
    h_pdoping_trap       = 0.1    * μm
    h_pdoing_right       = 1.0    * μm
    w_device             = 0.5    * μm  # width of device
    z_device             = 1.0e-4 * cm  # depth of device
    h_total              = h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoing_right
    coord                = initialize_pin_grid(refinementfactor,
                                               h_ndoping,
                                               h_pdoping_left,
                                               h_pdoping_trap,
                                               h_pdoing_right)

    grid                 = simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [h_ndoping], regionDonor) # n doped
    cellmask!(grid, [h_ndoping], [h_ndoping + h_pdoping_left], regionAcceptorLeft) # p doped
    cellmask!(grid, [h_ndoping + h_pdoping_left], [h_ndoping + h_pdoping_left + h_pdoping_trap], regionAcceptorTrap) # p doped with traps
    cellmask!(grid, [h_ndoping + h_pdoping_left + h_pdoping_trap], [h_total], regionAcceptorRight) # p doped

    bfacemask!(grid, [h_ndoping], [h_ndoping], bregionDALeft, tol = 1.0e-18)
    bfacemask!(grid, [h_ndoping + h_pdoping_left], [h_ndoping + h_pdoping_left], bregionALeftATrap, tol = 1.0e-18)
    bfacemask!(grid, [h_ndoping + h_pdoping_left + h_pdoping_trap], [h_ndoping + h_pdoping_left + h_pdoping_trap], bregionATrapARight, tol = 1.0e-18)

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

    iphin                = 1 # index electron quasi Fermi potential
    iphip                = 2 # index hole quasi Fermi potential
    if AdditionalTrapSpecies
        iphit            = 3 # index trap quasi Fermi potential
        numberOfCarriers = 3 # electrons, holes and traps
    else
        numberOfCarriers = 2 # electrons and holes
    end

    ## physical data
    T                 = 300.0                *  K

    ## band edge energies
    Ec_ZnO            = 3.4                  *  eV
    Ev_ZnO            = 0.0                  *  eV

    Ec_CIGS           = 3.4                  *  eV
    Ev_CIGS           = 2.3                  *  eV

    EC                = [Ec_ZnO, Ec_CIGS, Ec_CIGS, Ec_CIGS]
    EV                = [Ev_ZnO, Ev_CIGS, Ev_CIGS, Ev_CIGS]

    ## effective densities of states
    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)

    NC                = [Nc, Nc, Nc, Nc]
    NV                = [Nv, Nv, Nv, Nv]

    ## mobilities
    mun_ZnO           = 100                  * (cm^2) / (V * s)
    mup_ZnO           = 25                   * (cm^2) / (V * s)
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)

    μn                = [mun_ZnO, mun_CIGS, mun_CIGS, mun_CIGS]
    μp                = [mup_ZnO, mup_CIGS, mup_CIGS, mup_CIGS]

    ## relative dielectric permittivity
    εr_ZnO            = 9                    *  1.0
    εr_CIGS           = 13.6                 *  1.0

    ε                 = [εr_ZnO, εr_CIGS, εr_CIGS, εr_CIGS]

    ## trap information
    zt                = 1 # hole traps
    Et                = 2.8                  *  eV
    ET                = [0.0, 0.0, Et, 0.0]
    Nt                = 1.0e18               / (cm^3)
    NT                = [0, 0, Nt, 0]

    mu_t              = 0                    * (cm^2) / (V * s)
    μt                = [0.0, 0.0, mu_t, 0.0]

    ## recombination information parameters
    ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T))   # intrinsic concentration
    n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )                   # Boltzmann equilibrium concentration
    p0_ZnO            = ni_ZnO^2 / n0_ZnO                                        # Boltzmann equilibrium concentration
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )                  # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration

    p0                = [p0_ZnO, p0_CIGS, p0_CIGS, p0_CIGS]
    n0                = [n0_ZnO, n0_CIGS, n0_CIGS, n0_CIGS]

    # set the lifetime value high in all other regions, such that SRH recombination can be neglected there
    SRH_LifeTime      = [1.0e100, 1.0e100, 1.0e-3*ns, 1.0e100]

    Auger             = 1.0e-29  * cm^6 / s
    Radiative         = 1.0e-10  * cm^3 / s

    ## Schottky contact information
    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier           = 0.7 * eV

    ## doping information
    Nd                = 1.0e18 / (cm^3)
    Na                = 5.5e15 / (cm^3)

    ## we will impose this applied voltage on one boundary
    voltageAcceptor   = 1.0 * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)
    data.modelType                      = Stationary # R = Rn = Rp, since the model type is stationary
    if AdditionalTrapSpecies
        data.F                          = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
    else
        data.F                         .= FermiDiracOneHalfTeSCA
    end

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    data.boundaryType[bregionAcceptor]  = SchottkyContact
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ExcessChemicalPotential

    if AdditionalTrapSpecies
        ## Here, we enable the traps and parse the respective index and the regions where the trap is defined.
        enable_trap_carrier!(;data = data, trapCarrier = iphit, regions = [regionAcceptorTrap])
    else
        ## pass trap data in stationary setting since there is no separate trap species
        add_trap_density_Poisson!(data = data, zt = zt, Nt = NT)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    ## physical parameters
    params                                              = Params(grid, numberOfCarriers)
    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    if AdditionalTrapSpecies
        params.chargeNumbers[iphit]                     =  zt
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]
        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]
        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        if AdditionalTrapSpecies
            params.densityOfStates[iphit, ireg]         = NT[ireg]
            params.bandEdgeEnergy[iphit, ireg]          = ET[ireg]
            params.mobility[iphit, ireg]                = μt[ireg]
        end

        ## recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = n0[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg] = p0[ireg]
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphip, regionAcceptorLeft]            = Na
    params.doping[iphip, regionAcceptorTrap]            = Na
    params.doping[iphip, regionAcceptorRight]           = Na

    ## values for the schottky contacts
    params.SchottkyBarrier[bregionAcceptor]             = barrier
    params.bVelocity[iphin, bregionAcceptor]            = vn
    params.bVelocity[iphip, bregionAcceptor]            = vp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control              = SolverControl()
    control.verbose      = verbose
    control.tol_round    = 1.0e-7
    control.damp_initial = 0.5
    control.damp_growth  = 1.2
    control.maxiters     = 30
    control.max_round    = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## solve thermodynamic equilibrium and update initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival   = solution

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        if AdditionalTrapSpecies
            ## add labels for traps
            label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
            label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"
        end

        ## ##### set legend for plotting routines #####
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Stationary bias loop")
    end
    ################################################################################

    endVoltage      = voltageAcceptor       # final bias value
    biasValues      = collect(range(0, stop = endVoltage, length = 52))

    IV              = zeros(0)
    chargeDensities = zeros(0)

    for i in eachindex(biasValues)

        Δu = biasValues[i] # bias

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("bias: Δu = $(Δu) V")
        end

        ## solve time step problems with timestep Δt
        solution = solve(ctsys, inival = inival, control = control)
        inival   = solution

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store charge density in donor region (ZnO)
        push!(chargeDensities, charge_density(ctsys, solution)[regionDonor])


    end # bias loop

    ## compute static capacitance: check this is correctly computed
    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)

    ## plot solution and IV curve
    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage) V", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(endVoltage) V", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage) V", label_solution)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,chargeDensities, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.title("Charge density in donor region")
        Plotter.ylabel("Charge density [C]")
        Plotter.figure()
        plot_IV(Plotter, biasValues, staticCapacitance, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.title("Static capacitance in donor region")
        Plotter.ylabel("Static capacitance [C/V]")

    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval                  = 1.484831264268335
    testvalAdditionalSpecies = 1.1334257649339574

    main(test = true, AdditionalTrapSpecies = false) ≈ testval && main(test = true, AdditionalTrapSpecies = true) ≈ testvalAdditionalSpecies
end

if test == false
    println("This message should show when this module has successfully recompiled.")
end


end # module

