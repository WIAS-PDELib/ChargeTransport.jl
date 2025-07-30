#=
# MoS2 with moving defects and Schottky Barrier Lowering.
([source code](@__SOURCE_URL__))

Memristor simulation with additional moving positively charged defects and
Schottky barrier lowering at the contacts.
=#


module Ex107_MoS2_withIons_BarrierLowering

using ChargeTransport
using ExtendableGrids
using PyPlot

# you can also use other Plotters, if you add them to the example file
# you can set verbose also to true to display some solver information
function main(; Plotter = PyPlot, plotting = false, verbose = "", test = false, barrierLowering = true)

    if plotting
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Set up grid, regions and time mesh")
    end
    ################################################################################

    @local_unitfactors μm cm eV s ns V K ps Hz W m

    constants = ChargeTransport.constants
    (; q, k_B, ε_0, Planck_constant, m_e) = constants


    ## region numbers
    regionflake = 1

    ## boundary region numbers
    bregionLeft = 1
    bregionRight = 2

    ## grid
    h_flake = 1.0 * μm # length of the conducting channel

    # non-uniform grid
    coord1 = geomspace(0.0, h_flake / 2, 5.0e-4 * h_flake, 2.0e-2 * h_flake)
    coord2 = geomspace(h_flake / 2, h_flake, 2.0e-2 * h_flake, 5.0e-4 * h_flake)
    coord = glue(coord1, coord2)

    grid = simplexgrid(coord)

    ## set region in grid
    cellmask!(grid, [0.0], [h_flake], regionflake, tol = 1.0e-18)

    if plotting
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
    iphix = 3

    numberOfCarriers = 3 # electrons, holes and ions

    # We define the physical data
    T = 300.0 * K
    εr = 9.0 * 1.0                   # relative dielectric permittivity
    εi = 1.0 * εr                                   # image force dielectric permittivity

    Ec = - 4.0 * eV
    Ev = - 5.3 * eV
    Ex = - 4.38 * eV

    Nc = 2 * (2 * pi * 0.55 * m_e * k_B * T / (Planck_constant^2))^(3 / 2) / m^3
    Nv = 2 * (2 * pi * 0.71 * m_e * k_B * T / (Planck_constant^2))^(3 / 2) / m^3
    Nx = 1.0e28 / (m^3)

    μn = 1.0e-4 * (m^2) / (V * s)
    μp = 1.0e-4 * (m^2) / (V * s)
    μx = 0.8e-13 * (m^2) / (V * s)

    ## Schottky contact
    barrierLeft = 0.225 * eV
    barrierRight = 0.215 * eV
    An = 4 * pi * q * 0.55 * m_e * k_B^2 / Planck_constant^3
    Ap = 4 * pi * q * 0.71 * m_e * k_B^2 / Planck_constant^3
    vn = An * T^2 / (q * Nc)
    vp = Ap * T^2 / (q * Nv)

    Nd = 1.0e17 / (m^3) # doping

    Area = 2.1e-11 * m^2                # Area of electrode

    # Scan protocol information
    endTime = 9.6 * s
    amplitude = 12.0 * V
    scanrate = 4 * amplitude / endTime

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
    data = Data(grid, numberOfCarriers, constants, contactVoltageFunction = contactVoltageFunction)
    data.modelType = Transient
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = false
    )
    if barrierLowering
        data.boundaryType[bregionLeft] = SchottkyBarrierLowering
        data.boundaryType[bregionRight] = SchottkyBarrierLowering
    else
        data.boundaryType[bregionLeft] = SchottkyContact
        data.boundaryType[bregionRight] = SchottkyContact
    end

    data.fluxApproximation .= ExcessChemicalPotential

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
    params.chargeNumbers[iphip] = 1
    params.chargeNumbers[iphix] = 2

    for ireg in 1:length([regionflake])           # region data

        params.dielectricConstant[ireg] = εr * ε_0
        params.dielectricConstantImageForce[ireg] = εi * ε_0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nc
        params.densityOfStates[iphip, ireg] = Nv
        params.bandEdgeEnergy[iphin, ireg] = Ec
        params.bandEdgeEnergy[iphip, ireg] = Ev
        params.mobility[iphin, ireg] = μn
        params.mobility[iphip, ireg] = μp
        params.densityOfStates[iphix, ireg] = Nx
        params.bandEdgeEnergy[iphix, ireg] = Ex
        params.mobility[iphix, ireg] = μx
    end

    params.SchottkyBarrier[bregionLeft] = barrierLeft
    params.SchottkyBarrier[bregionRight] = barrierRight
    params.bVelocity[iphin, bregionLeft] = vn
    params.bVelocity[iphin, bregionRight] = vn
    params.bVelocity[iphip, bregionLeft] = vp
    params.bVelocity[iphip, bregionRight] = vp

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
    if verbose == true
        control.verbose = verbose
    else
        control.verbose = "eda" # still print the time values
    end
    if test == true
        control.verbose = false # do not show time values in testing case
    end
    control.damp_initial = 0.9
    control.damp_growth = 1.61 # >= 1
    control.max_round = 5

    control.abstol = 1.0e-9
    control.reltol = 1.0e-9
    control.tol_round = 1.0e-9

    control.Δu_opt = Inf
    control.Δt = 1.0e-4
    control.Δt_min = 1.0e-5
    control.Δt_max = 5.0e-2
    control.Δt_grow = 1.05

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize solution and starting vectors
    solEQ = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 0)
    inival = solEQ

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)
        label_energy[1, iphix] = "\$E_x-q\\psi\$"; label_energy[2, iphix] = "\$ - q \\varphi_x\$"
        label_density[iphix] = "\$ n_x\$";       label_solution[iphix] = "\$ \\varphi_x\$"

        Plotter.figure()
        plot_densities(Plotter, ctsys, solEQ, "Equilibrium", label_density)
        Plotter.legend()
        Plotter.figure()
        plot_solution(Plotter, ctsys, solEQ, "Equilibrium", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    sol = solve(ctsys, inival = inival, times = (0.0, endTime), control = control)

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

        Plotter.figure()
        Plotter.semilogy(biasValues, abs.(Area .* IV), linewidth = 5, color = "black")
        Plotter.grid()
        Plotter.xlabel("applied bias [V]")
        Plotter.ylabel("total current [A]")
    end


    testval = sum(filter(!isnan, solEQ)) / length(solEQ)
    return testval

end #  main

function test()
    return main(test = true, barrierLowering = true) ≈ -1692.2303837883194
end


end # module
