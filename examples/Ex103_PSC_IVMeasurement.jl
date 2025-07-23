#=
# PSC device with ions and different I-V scan protocols (1D).
([source code](@__SOURCE_URL__))

Simulating a three layer PSC device Ti02| MAPI | spiro-OMeTAD with mobile ions where
the ion vacancy accumulation is limited by the Fermi-Dirac integral of order -1.

The time-dependent simulations are performed with abrupt interfaces.
Two different I-V measurement protocols are included and the corresponding solution vectors
and the I-V curve after the scan can be depicted.


=#

module Ex103_PSC_IVMeasurement

using ChargeTransport
using ExtendableGrids
using PyPlot

# # you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(;
        n = 3, Plotter = PyPlot, plotting = false, verbose = "", test = false,
        parameter_set = Params_PSC_TiO2_MAPI_spiro, # choose the parameter set
        otherScanProtocol = false
    ) # you can choose between two scan protocols

    if plotting
        Plotter.close("all")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # parameters
    p = parameter_set()

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate = 1.0 * V / s
    number_tsteps = 31
    endVoltage = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = endVoltage / scanrate

    ## Define scan protocol function
    function linearScanProtocol(t)
        return if t == Inf
            0.0
        else
            scanrate * t
        end
    end

    ## Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, linearScanProtocol]

    ## Instead of a linear scan protocol, we can also apply other scan protocols which we
    ## define by our own and parse to the model generator via the struct Data
    if otherScanProtocol
        ## scan protocol parameter
        number_tsteps = 40
        frequency = 10.0 * Hz
        amplitude = 0.2 * V
        tend = 1 / frequency

        ## Define sinusoidal applied voltage
        function sinusoidalScanProtocol(t)
            return if t == Inf
                0.0
            else
                amplitude * sin(2.0 * pi * frequency * t)
            end
        end

        ## Apply zero voltage on left boundary and a linear scan protocol on right boundary
        contactVoltageFunction = [zeroVoltage, sinusoidalScanProtocol]
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u = collect(range(0.0, p.h_ndoping / 2, step = p.h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        p.h_ndoping / 2, p.h_ndoping,
        p.h_ndoping / (0.7 * δ), p.h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ndoping, p.h_ndoping + p.h_intrinsic / k,
        p.h_intrinsic / (2.8 * δ), p.h_intrinsic / (2.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ndoping + p.h_intrinsic / k, p.h_ndoping + p.h_intrinsic,
        p.h_intrinsic / (2.1 * δ), p.h_intrinsic / (2.8 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ndoping + p.h_intrinsic, p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2,
        p.h_pdoping / (1.6 * δ), p.h_pdoping / (1.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, p.h_ndoping + p.h_intrinsic + p.h_pdoping, step = p.h_pdoping / (1.3 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [p.heightLayers[1]], p.regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [p.heightLayers[1]], [p.heightLayers[2]], p.regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [p.heightLayers[2]], [p.heightLayers[3]], p.regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], p.bregionDonor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [p.h_total], [p.h_total], p.bregionAcceptor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2, tol = 1.0e-18) # second inner interface

    if plotting
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
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
    ## Currently, the way to go is to pass a contact voltage function exactly here.
    data = Data(grid, p.numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[p.bregionAcceptor] = OhmicContact
    data.boundaryType[p.bregionDonor] = OhmicContact

    ## With this method, the user enable the ionic carrier parsed to ionicCarrier and gives
    ## gives the information on which regions this ionic carrier is defined.
    ## In this application ion vacancies only live in active perovskite layer.
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params = Params(grid, p.numberOfCarriers)

    params.temperature = p.T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za

    for ireg in 1:p.numberOfRegions # region data

        params.dielectricConstant[ireg] = p.ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[p.iphin, ireg] = p.Nn[ireg]
        params.densityOfStates[p.iphip, ireg] = p.Np[ireg]
        params.densityOfStates[p.iphia, ireg] = p.Na[ireg]

        params.bandEdgeEnergy[p.iphin, ireg] = p.En[ireg]
        params.bandEdgeEnergy[p.iphip, ireg] = p.Ep[ireg]
        params.bandEdgeEnergy[p.iphia, ireg] = p.Ea[ireg]

        params.mobility[p.iphin, ireg] = p.μn[ireg]
        params.mobility[p.iphip, ireg] = p.μp[ireg]
        params.mobility[p.iphia, ireg] = p.μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = p.r0[ireg]
        params.recombinationSRHLifetime[p.iphin, ireg] = p.τn[ireg]
        params.recombinationSRHLifetime[p.iphip, ireg] = p.τp[ireg]
        params.recombinationSRHTrapDensity[p.iphin, ireg] = trap_density!(p.iphin, ireg, params, p.EI[ireg])
        params.recombinationSRHTrapDensity[p.iphip, ireg] = trap_density!(p.iphip, ireg, params, p.EI[ireg])
    end

    ## doping
    params.doping[p.iphin, p.regionDonor] = p.Cn
    params.doping[p.iphia, p.regionIntrinsic] = p.Ca
    params.doping[p.iphip, p.regionAcceptor] = p.Cp

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia] = "\$ n_a \$";      label_solution[iphia] = "\$ \\varphi_a\$"

    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.max_round = 5
    control.damp_initial = 0.1
    control.damp_growth = 1.21 # >= 1

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

    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    ## with fixed timestep sizes we can calculate the times a priori
    tvalues = range(0, stop = tend, length = number_tsteps)

    ## for saving I-V data
    IV = zeros(0)                   # for IV values
    ISRHn = zeros(0); ISRHp = zeros(0) # for SRH recombination current
    IRadn = zeros(0); IRadp = zeros(0) # for radiative recombination current

    for istep in 2:number_tsteps

        t = tvalues[istep]                                    # Actual time
        Δu = contactVoltageFunction[p.bregionAcceptor](t) # Applied voltage
        Δt = t - tvalues[istep - 1]                              # Time step size

        ## Apply new voltage by setting non equilibrium boundary conditions
        set_contact!(ctsys, p.bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t) s")
        end

        ## Solve time step problems with timestep Δt. inival plays the role of the solution
        ## from last timestep
        solution = solve(ctsys; inival = inival, control = control, tstep = Δt)
        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)
        IntSRH = integrate(ctsys, SRHRecombination!, solution)
        IntRad = integrate(ctsys, RadiativeRecombination!, solution)

        IntSRHnSum = 0.0; IntRadnSum = 0.0
        IntSRHpSum = 0.0; IntRadpSum = 0.0

        for ii in 1:p.numberOfRegions
            IntSRHnSum = IntSRHnSum - IntSRH[p.iphin, ii]
            IntRadnSum = IntRadnSum - IntRad[p.iphin, ii]

            IntSRHpSum = IntSRHpSum + IntSRH[p.iphip, ii]
            IntRadpSum = IntRadpSum + IntRad[p.iphip, ii]
        end

        push!(IV, current)
        push!(ISRHn, IntSRHnSum); push!(ISRHp, IntSRHpSum)
        push!(IRadn, IntRadnSum); push!(IRadp, IntRadpSum)

        inival = solution

    end # time loop

    if test == false
        println("*** done\n")
    end

    ## here in res the biasValues and the corresponding current are stored.
    ## res = [biasValues IV];

    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage)", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage)", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage)", label_solution)
    end

    biasValues = contactVoltageFunction[p.bregionAcceptor].(tvalues)

    if plotting
        Plotter.figure()
        Plotter.plot(tvalues, biasValues, marker = "o")
        Plotter.xlabel("time [s]")
        Plotter.ylabel("bias [V]")
        Plotter.figure()
        plot_IV(Plotter, biasValues[2:end], IV, "bias \$\\Delta u\$ = $(endVoltage)")
        ###############
        Plotter.figure()
        semilogy(biasValues[2:end], ISRHn .* (cm^2) .* 1.0e3, linewidth = 5, color = "darkblue", label = "SRH recombination")
        semilogy(biasValues[2:end], ISRHp .* (cm^2) .* 1.0e3, linewidth = 5, color = "lightblue", linestyle = ":")
        semilogy(biasValues[2:end], IRadn .* (cm^2) .* 1.0e3, linewidth = 5, color = "darkgreen", label = "Radiative recombination")
        semilogy(biasValues[2:end], IRadp .* (cm^2) .* 1.0e3, linewidth = 5, color = "lightgreen", linestyle = ":")

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]")
        PyPlot.ylabel("current density [mAcm\$^{-2} \$]")
    end

    testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
    return testval


end #  main

function test()
    testval = -0.6302819608784171; testvalOther = -1.123710261723505
    return main(test = true, otherScanProtocol = false) ≈ testval && main(test = true, otherScanProtocol = true) ≈ testvalOther
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
