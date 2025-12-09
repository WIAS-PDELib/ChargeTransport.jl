#=
# MoS2 with moving defects and immobile traps.
([source code](@__SOURCE_URL__))

Memristor simulation with additional moving positively charged defects and
immobile acceptor traps which are coupled to the conduction band.
=#

module Ex110_MoS2_withIons_withTraps

using ChargeTransport
using ExtendableGrids
using PyPlot

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(; Plotter = PyPlot, plotting = false, verbose = false, test = false, barrierLowering = true)

    

    ################################################################################
    if test == false
        println("Set up grid, regions and time mesh")
    end
    ################################################################################

    @local_unitfactors μm cm s ns V K ps Hz W m
    

    constants = ChargeTransport.constants
    (; q, k_B, ε_0, Planck_constant, m_e) = constants
    eV = q*V


    ## region numbers
    regionflake = 1

    ## boundary region numbers
    bregionLeft  = 1
    bregionRight = 2

    ## grid
    h_flake = 1.0 * μm # length of the conducting channel

    # non-uniform grid
    coord1 = geomspace(0.0, h_flake / 2, 5.0e-4 * h_flake, 2.0e-2 * h_flake)
    coord2 = geomspace(h_flake / 2, h_flake, 2.0e-2 * h_flake, 5.0e-4 * h_flake)
    coord  = glue(coord1, coord2)

    grid   = simplexgrid(coord)

    ## set region in grid
    cellmask!(grid, [0.0], [h_flake], regionflake, tol = 1.0e-18)

    if plotting
        # Plotter.close("all")
        gridplot(grid, Plotter = Plotter)
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

    ## set indices of unknowns
    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    iphit = 3 # Trap occupation level
    iphix = 4 # Vacancy quasi Fermi potential

    numberOfCarriers = 4 # electrons, holes, traps and ions

    # We define the physical data
    T = 300.0 * K
    εr =  9.0 * 1.0                   # relative dielectric permittivity
    εi =  1.0 * εr                    # image force dielectric permittivity

    Ec = - 4.0 * eV
    Ev = - 5.3 * eV
    Et = Ec - 0.5*eV                  # Trap level (0.5 eV below the conduction band)
    Ex = - 4.38 * eV

    Nc = 2 * (2 * pi * 0.55 * m_e * k_B * T / (Planck_constant^2))^(3 / 2) / m^3
    Nv = 2 * (2 * pi * 0.71 * m_e * k_B * T / (Planck_constant^2))^(3 / 2) / m^3
    
    Nt = 5.0e21 / (m^3)              # Trap density
    Nx = 1.0e28 / (m^3)              # Vacancy density

    μn = 1.00e-4 * (m^2) / (V * s)  # Electron mobility
    μp = 1.00e-4 * (m^2) / (V * s)  # Hole mobility
    μx = 0.8e-13 * (m^2) / (V * s)  # Vacancy mobility
                                    # Traps are frozen and have no flux expression
                                    # -> No mobility for trap species
    

    ## Schottky contact
    barrierLeft  = 0.225 * eV
    barrierRight = 0.215 * eV
    An = 4 * pi * q * 0.55 * m_e * k_B^2 / Planck_constant^3
    Ap = 4 * pi * q * 0.71 * m_e * k_B^2 / Planck_constant^3
    vn = An * T^2 / (q * Nc)
    vp = Ap * T^2 / (q * Nv)

    Nd = 1.0e10 / (cm^3) # doping

    Area = 2.1e-11 * m^2 # Area of electrode

    ## Trap capture and escape parameters
    vth = sqrt(3.0 * k_B * T / (0.55 * m_e) )       # Thermal velocity
    σₙ         = 1.0e-16 * (cm^2)                   # Scattering cross section
    sₙ         = σₙ * vth                           # Capture rate (escape calculated automatically from detailed balance)
    if test == false
        println("Capture rate: $(sₙ)") 
        # The capture rate assuming Boltzmann statistics can be written directly, 
        # however corrections due to non-Boltzmann statistics depends on the carrier
        # density in the band. 
        println("Approximate escape rate: $(sₙ*Nc*exp((Et-Ec)/(k_B*T)))") 
    end
    

    # Scan protocol information
    endTime   =  9.6 * s
    amplitude = 12.0 * V 
    scanrate  =  4.0 * amplitude / endTime

    ## Define scan protocol function
    function scanProtocol(t)

        if 0.0 <= t  && t <= endTime / 4
            biasVal = 0.0 + scanrate * t
        elseif t >= endTime / 4  && t <= 3 * endTime / 4
            biasVal = amplitude .- scanrate * (t - endTime / 4)
        elseif t >= 3 * endTime / 4 && t <= endTime
            biasVal = - amplitude .+ scanrate * (t - 3 * endTime / 4)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    # Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocol]

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    data = Data(grid, numberOfCarriers, contactVoltageFunction = contactVoltageFunction)
    data.modelType = Transient
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne, FermiDiracMinusOne]
    
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = false,
        bulk_recomb_trap=ChargeTransport.SingleStateTrap    # Set trap type to a singe-state trap
    )
    if barrierLowering
        data.boundaryType[bregionLeft] = SchottkyBarrierLowering
        data.boundaryType[bregionRight] = SchottkyBarrierLowering
    else
        data.boundaryType[bregionLeft] = SchottkyContact
        data.boundaryType[bregionRight] = SchottkyContact
    end

    data.fluxApproximation .= ExcessChemicalPotential

    # Populate trap and ionic carrier list from carrier list
    enable_trap_carrier!( data,  trapCarrier = iphit, regions = [regionflake])
    enable_ionic_carrier!(data, ionicCarrier = iphix, regions = [regionflake])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

    params.temperature = T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] =  1
    params.chargeNumbers[iphit] = -1    # Charge of trap in filled state: -1 -> Acceptor. +1 -> Donor
    params.chargeNumbers[iphix] =  2

    for ireg in 1:length([regionflake])           # region data

        params.dielectricConstant[ireg] = εr * ε_0
        params.dielectricConstantImageForce[ireg] = εi * ε_0

        ## effective DOS, band-edge energy and mobilities
        # Band carriers
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = μn
        params.mobility[iphip, ireg] = μp
        
        # Immobile traps density, depth and capture rate
        params.densityOfStates[iphit, ireg] = Nt
        params.bandEdgeEnergy[iphit, ireg] = Et
        params.recombinationTrapCaptureRates[iphit,iphin,ireg] = sₙ
        

        # Vacancies
        params.densityOfStates[iphix, ireg] = Nx
        params.bandEdgeEnergy[iphix, ireg] = Ex
        params.mobility[iphix, ireg] = μx

    end

    params.SchottkyBarrier[bregionLeft]  = barrierLeft
    params.SchottkyBarrier[bregionRight] = barrierRight
    params.bVelocity[iphin, bregionLeft] = vn
    params.bVelocity[iphin, bregionRight]= vn
    params.bVelocity[iphip, bregionLeft] = vp
    params.bVelocity[iphip, bregionRight]= vp

    ## interior doping
    params.doping[iphin, regionflake] = Nd

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

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
    control.damp_initial = 0.5
    control.damp_growth  = 1.61 # >= 1
    control.max_round    = 5

    control.abstol       = 1.0e-9
    control.reltol       = 1.0e-9
    control.tol_round    = 1.0e-9

    control.Δu_opt       = Inf
    control.Δt           = 1.0e-4
    control.Δt_min       = 1.0e-5
    control.Δt_max       = 5.0e-2
    control.Δt_grow      = 1.05

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

   
    ## initialize solution and starting vectors
    solEQ  = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 0)
    inival = copy(solEQ)
    

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)
        label_energy[1, iphit] = "\$E_t-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_t\$"
        label_density[iphit] = "\$ n_t\$";       label_solution[iphit] = "\$ \\varphi_t\$"
        label_energy[1, iphix] = "\$E_x-q\\psi\$"; label_energy[2, iphix] = "\$ - q \\varphi_x\$"
        label_density[iphix] = "\$ n_x\$";       label_solution[iphix] = "\$ \\varphi_x\$"

        Plotter.figure()
        plot_densities(Plotter, ctsys, solEQ, "Equilibrium", label_density)
        Plotter.legend()
        Plotter.figure()
        plot_energies(Plotter, ctsys, solEQ, "Equilibrium", label_energy)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    sol = solve(ctsys, inival = inival, times = (0.0,endTime), control = control)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    #########  IV curve calculation
    ################################################################################

    IV = zeros(0) # for saving I-V data

    tvalues = sol.t
    number_tsteps = length(tvalues)
    biasValues = scanProtocol.(tvalues)

    factory = TestFunctionFactory(ctsys)
    tf = testfunction(factory, [bregionLeft], [bregionRight])

    push!(IV, 0.0)
    for istep in 2:number_tsteps
        Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
        inival = sol.u[istep - 1]
        solution = sol.u[istep]

        I = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii in 1:(numberOfCarriers + 1)
            current = current + I[ii]
        end

        push!(IV, current)

    end

    if plotting
        Plotter.figure()
        Plotter.plot(tvalues, biasValues, marker = "x")
        Plotter.xlabel("time [s]")
        Plotter.ylabel("voltage [V]")
        Plotter.grid()
        tight_layout()

        Plotter.figure()
        Plotter.semilogy(biasValues, abs.(Area .* IV), linewidth = 5, color = "black")
        Plotter.grid()
        Plotter.xlabel("applied bias [V]")
        Plotter.ylabel("total current [A]")
        tight_layout()
    end
    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)
        label_energy[1, iphit] = "\$E_t-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_t\$"
        label_density[iphit] = "\$ n_t\$";       label_solution[iphit] = "\$ \\varphi_t\$"
        label_energy[1, iphix] = "\$E_x-q\\psi\$"; label_energy[2, iphix] = "\$ - q \\varphi_x\$"
        label_density[iphix] = "\$ n_x\$";       label_solution[iphix] = "\$ \\varphi_x\$"

        Plotter.figure()
        plot_densities(Plotter, ctsys, sol.u[end], "End of sweep", label_density)
        Plotter.legend()
        Plotter.figure()
        plot_energies(Plotter, ctsys, sol.u[end], "End of sweep", label_energy)
        Plotter.legend()

    end

    testval=sum(filter(!isnan, sol.u[end])) / length(sol.u[end])
    return testval

end #  main


function test()
    return main(;verbose="",test=true)≈-6894.342164365617
end

end # Module