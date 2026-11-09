#=
# PSC device with degradation/recovery reactions as an explicit customReaction! (1D).
([source code](@__SOURCE_URL__))

Three layer PSC device C60 | triple cation perovskite | PTAA with mobile anion vacancies
and a two-reaction degradation/recovery network in the perovskite layer:
holes are converted into anion vacancies (degradation, constant rate k₁) and back (recovery,
rate k₂ depending on the local vacancy and hole densities).

The reaction terms are written out as one function for this network inside `main` and
passed as `Data(...; customReaction! = ...)`.

Starting from the dark equilibrium state at zero bias, the light is switched on with a
smooth 2 ms ramp of the Beer-Lambert generation, written as the function `customGeneration`,
passed as `Data(...; customGeneration = ...)` and used with `GenerationCustom`, while the bias
is raised from 0 to 1.1 V with the same raised cosine; the device is then kept at 1.1 V
under light for 30 min.
Every 300 s a non-perturbing reverse J-V scan (1.2 V → 0 V at 0.11 V/s) is run from the
current state without changing the trajectory, and Jsc, Voc, FF and PCE are computed from it.
=#

module Ex111_PSC_Reaction_Custom

using ChargeTransport
using ExtendableGrids
using GridVisualize

# for convenience
parametersdir = ChargeTransport.parametersdir

# supported Plotters are GLMakie and PythonPlot
# you can set verbose also to true to display some solver information
function main(;
        n = 5,
        Plotter = nothing,
        verbose = false, test = false,
        ########################
        parameter_set = Params_PSC_C60_TripleCation_PTAA, # choose the parameter set
        ########################
        tend = 1800.0,                                   # duration of the light phase [s]
    )

    @local_unitfactors μm cm s ns V K ps Hz W m

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # parameters; vacancy density of states and temperature differ from the parameter file
    p = parameter_set(Na = [0.0, 1.0e26, 0.0] ./ m^3, T = 298.15 * K)

    ## applied bias V(t) on the right contact (contactVoltage stays 0). The type of
    ## contactVoltageFunction is fixed when Data is built, but the J-V scans need a different
    ## V(t) than the light phase, so V(t) is held in a container and swapped via bias[].
    bias = Ref{Function}(t -> 0.0)
    contactVoltageFunction = [zeroVoltage, t -> bias[](t)]

    ## light phase: generation and bias 0 → V_light with a 2 ms raised cosine, then held
    t_light_ramp = 2.0e-3 * s
    V_light = 1.1 * V
    light_ramp(t) = 0.5 * (1.0 - cos(pi * clamp(t, 0.0, t_light_ramp) / t_light_ramp))
    light_bias(t) = V_light * light_ramp(t)

    ## Beer-Lambert generation switched on with light_ramp; used with GenerationCustom
    function customGeneration(node, data)
        params = data.params
        ireg = node.region
        x = node.coord[node.index]
        G = params.generationIncidentPhotonFlux[ireg] * params.generationAbsorption[ireg] * exp(- params.invertedIllumination * params.generationAbsorption[ireg] * (x - params.generationPeak))

        return light_ramp(node.time) * G
    end

    ## non-perturbing reverse J-V scans
    t_pre = 0.2 * s            # raised-cosine transition V_light → V_start before the scan
    V_start = 1.2 * V
    V_end = 0.0 * V
    scan_rate = 0.11 * V / s
    P_in = 1000.0 * W / m^2    # incident light power, for the PCE
    ## output times of a scan: 20 points in the transition, then a 5 mV voltage grid
    n_pre = 20
    t_meas = abs(V_end - V_start) / scan_rate
    n_meas = round(Int, abs(V_end - V_start) / (5.0e-3 * V))
    direction = sign(V_end - V_start)
    scan_bias(t) = t <= t_pre ?
        V_light + (V_start - V_light) * 0.5 * (1.0 - cos(pi * t / t_pre)) :
        V_start + direction * scan_rate * min(t - t_pre, t_meas)
    tscan = [range(0.0, t_pre, length = n_pre + 1); (t_pre .+ range(0.0, t_meas, length = n_meas + 1))[2:end]]

    ## interface recombination velocities
    vETL = 2000.0 * cm / s
    vHTL = 500.0 * cm / s

    ############################################################################
    ## Reaction terms for holes ⇌ anion vacancies in the perovskite layer:
    ##   R₁ = k₁ n_p        degradation  p → a,  k₁ constant
    ##   R₂ = k₂ n_a        recovery     a → p,  k₂ = β_ion(n_a) · f_hole(n_p)
    ## the rate is added for the reactant and subtracted for the product, scaled by q z:
    ##   f[p] += q z_p (R₁ - R₂),   f[a] += q z_a (R₂ - R₁)
    function reaction_holes_vacancies!(f, u, node, data)
        ## same gate as the hard-coded reactions: off in the equilibrium solve and before
        ## generationComplete, perovskite only
        if !(data.calculationType == OutOfEquilibrium && data.generationComplete && node.region == 2)
            return nothing
        end

        ## recovery-rate constants
        β_min = 1.0e-5          # [1/s]  lower bound of β_ion
        β_max = 3.0e-5          # [1/s]  upper bound of β_ion
        log_n_ref = 16.5        # [log10(cm⁻³)]  midpoint of β_ion in n_a
        ion_steepness = 5.0     # [1]  slope of β_ion in log10(n_a)
        hole_min = 0.1          # [1]  lower bound of f_hole
        hole_max = 2.0          # [1]  upper bound of f_hole
        log_p_ref = 15.0        # [log10(cm⁻³)]  midpoint of f_hole in n_p
        hole_steepness = 5.0    # [1]  slope of f_hole in log10(n_p)

        n_p = ChargeTransport.get_density!(u, node, data, p.iphip) # [m⁻³]
        n_a = ChargeTransport.get_density!(u, node, data, p.iphia) # [m⁻³]

        k1 = 1.0e1                                  # [1/s]

        ## k₂ from the densities in cm⁻³:
        ##   β_ion  = β_min + (β_max - β_min) / (1 + exp(-s_a (log10 n_a - log_n_ref))), clamped
        ##   f_hole = h_min + (h_max - h_min) / (1 + exp(-s_p (log10 n_p - log_p_ref))), clamped
        n_a_cm = n_a * 1.0e-6                       # [cm⁻³]
        n_p_cm = n_p * 1.0e-6                       # [cm⁻³]
        if n_a_cm <= 0.0
            k2 = β_min * hole_min * one(n_a)       # same type as n_a (no allocation)
        else
            β_ion = β_min + (β_max - β_min) /
                (1.0 + exp(-ion_steepness * (log(n_a_cm) / log(10) - log_n_ref)))
            β_ion = clamp(β_ion, β_min, β_max)
            if n_p_cm <= 0.0
                k2 = β_ion * hole_min
            else
                f_hole = hole_min + (hole_max - hole_min) /
                    (1.0 + exp(-hole_steepness * (log(n_p_cm) / log(10) - log_p_ref)))
                f_hole = clamp(f_hole, hole_min, hole_max)
                k2 = β_ion * f_hole                 # [1/s]
            end
        end

        R1 = k1 * n_p                               # [m⁻³ s⁻¹]
        R2 = k2 * n_a                               # [m⁻³ s⁻¹]

        q = data.constants.q
        z_p = data.params.chargeNumbers[p.iphip]
        z_a = data.params.chargeNumbers[p.iphia]
        f[p.iphip] += q * z_p * (R1 - R2)
        f[p.iphia] += q * z_a * (R2 - R1)

        return nothing
    end
    ############################################################################

    ## scan at the end of the light ramp and every 300 s
    tscans = [t_light_ramp; 300.0:300.0:tend]
    tout = unique([0.0; tscans; tend])

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

    coord_n_u = collect(range(0.0, p.h_ETL / 2, step = p.h_ETL / (0.8 * δ)))
    coord_n_g = geomspace(
        p.h_ETL / 2,
        p.h_ETL,
        p.h_ETL / (0.7 * δ),
        p.h_ETL / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ETL,
        p.h_ETL + p.h_activePL / k,
        p.h_activePL / (2.8 * δ),
        p.h_activePL / (2.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ETL + p.h_activePL / k,
        p.h_ETL + p.h_activePL,
        p.h_activePL / (2.1 * δ),
        p.h_activePL / (2.8 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ETL + p.h_activePL,
        p.h_ETL + p.h_activePL + p.h_HTL / 2,
        p.h_HTL / (1.6 * δ),
        p.h_HTL / (1.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ETL + p.h_activePL + p.h_HTL / 2, p.h_ETL + p.h_activePL + p.h_HTL, step = p.h_HTL / (1.3 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [p.heightLayersPL[1]], p.regionETL, tol = 1.0e-18)                  # ETL         = 1
    cellmask!(grid, [p.heightLayersPL[1]], [p.heightLayersPL[2]], p.regionPero, tol = 1.0e-18)     # perovskite = 2
    cellmask!(grid, [p.heightLayersPL[2]], [p.heightLayersPL[3]], p.regionHTL, tol = 1.0e-18)      # HTL        = 3

    bfacemask!(grid, [p.heightLayersPL[1]], [p.heightLayersPL[1]], p.bregionJ1, tol = 1.0e-18)
    bfacemask!(grid, [p.heightLayersPL[2]], [p.heightLayersPL[2]], p.bregionJ2, tol = 1.0e-18)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    data = Data(
        grid, p.numberOfCarriers,
        contactVoltageFunction = contactVoltageFunction,
        customReaction! = reaction_holes_vacancies!,
        customGeneration = customGeneration
    )

    data.modelType = Transient

    ## The default for electrons and holes is Boltzmann. Here, we set it to a more general statistics function
    data.F[p.iphin] = FermiDiracOneHalfTeSCA
    data.F[p.iphip] = FermiDiracOneHalfTeSCA

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )
    data.boundaryType[p.bregionLeft] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionRight] = MixedOhmicSchottkyContact
    data.boundaryType[p.bregionJ1] = InterfaceRecombination
    data.boundaryType[p.bregionJ2] = InterfaceRecombination

    ## by default the statistics function is set to FermiDiracMinusOne to limit ion depletion
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionPero])

    ## dark during the equilibrium solve; the light is switched on afterwards
    data.generationModel = GenerationNone

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    data.params = Params(p)

    ## given vacancy energy level
    data.params.bandEdgeEnergy[p.iphia, p.regionPero] = p.Ea1D[p.regionPero]

    data.params.recombinationSRHvelocity[p.iphin, p.bregionJ1] = vETL
    data.params.recombinationSRHvelocity[p.iphip, p.bregionJ1] = vETL
    data.params.recombinationSRHvelocity[p.iphin, p.bregionJ2] = vHTL
    data.params.recombinationSRHvelocity[p.iphip, p.bregionJ2] = vHTL

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
    control.maxiters = 300
    control.max_round = 5
    control.damp_initial = 0.5
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
    solution = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = false)
    inival = solution
    solutionEQ = inival

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Light phase")
    end
    ################################################################################

    ctsys.data.generationModel = GenerationCustom
    bias[] = light_bias

    ## the ramp starts from the dark state, whose minority densities are so small that Newton
    ## needs tiny first steps; with Δu_opt = Inf the step grows by Δt_grow unless Newton fails
    control.Δt = 1.0e-10 * s
    control.Δt_min = 1.0e-12 * s
    control.Δt_max = 1.0 * s
    control.Δt_grow = 1.13
    control.Δu_opt = Inf

    sol = solve(ctsys, inival = inival, times = tout, control = control)
    ## VoronoiFVM only warns when a step fails; stop instead of continuing silently
    sol.t[end] ≈ tend || error("light phase stopped at t = $(sol.t[end]) s")

    ## average densities in the perovskite [m⁻³]
    average_density(icc) = [
        integrated_density(ctsys, sol = u, icc = icc, ireg = p.regionPero) / ctsys.data.regionVolumes[p.regionPero]
            for u in sol.u
    ]
    averages = (; time = sol.t, vacancies = average_density(p.iphia), holes = average_density(p.iphip))

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Non-perturbing J-V scans")
    end
    ################################################################################

    ## the scan's time restarts at 0, hence the time-independent full light (Beer-Lambert);
    ## ions are not frozen and the reactions stay on
    ctsys.data.generationModel = GenerationBeerLambert
    bias[] = scan_bias

    ## Δu_opt keeps its default here, unlike in the light phase
    controlScan = SolverControl()
    controlScan.verbose = verbose
    controlScan.maxiters = 300
    controlScan.max_round = 5
    controlScan.damp_initial = 0.5
    controlScan.damp_growth = 1.21 # >= 1
    controlScan.Δt = 1.0e-3 * t_pre
    controlScan.Δt_min = 1.0e-12 * V / scan_rate
    controlScan.Δt_max = 5.0e-3 * V / scan_rate
    controlScan.Δt_grow = 1.05

    Jsc = Float64[]; Voc = Float64[]; FF = Float64[]; PCE = Float64[]
    for tk in tscans
        inival = sol(tk)
        solscan = solve(ctsys; inival, times = tscan, control = controlScan)
        solscan.t[end] ≈ tscan[end] || error("scan at t = $tk s stopped at $(solscan.t[end]) s")

        ## current at each 5 mV point, from the solver step ending there; raw ChargeTransport sign
        voltage = scan_bias.(tscan[(n_pre + 2):end])
        current = map(tscan[(n_pre + 2):end]) do ti
            i = findfirst(==(ti), solscan.t)
            get_current_val(ctsys, solscan.u[i], solscan.u[i - 1], solscan.t[i] - solscan.t[i - 1])
        end

        v = reverse(voltage)
        j = -reverse(current)
        jsc = j[1]                         # the scan ends at V_end = 0
        i = findfirst(<(0), j)             # first point beyond Voc
        isnothing(i) && error("Voc above V_start in the scan at t = $tk s")
        voc = v[i - 1] - j[i - 1] * (v[i] - v[i - 1]) / (j[i] - j[i - 1])
        Pmax = maximum(v .* j)             # [W/m²]
        push!(Jsc, jsc)
        push!(Voc, voc)
        push!(FF, 100 * Pmax / (jsc * voc)) # [%]
        push!(PCE, 100 * Pmax / P_in)       # [%]
    end
    metrics = (; time = tscans, Jsc, Voc, FF, PCE)

    if Plotter !== nothing
        vis = GridVisualizer(; Plotter, layout = (2, 1), size = (700, 800))
        ## log time axis, so t = 0 is left out
        it = findall(>(0), averages.time)
        scalarplot!(
            vis[1, 1], averages.time[it], averages.vacancies[it];
            xscale = :log, yscale = :log, color = :red, label = "vacancies n_a", legend = :rc,
            xlabel = "time [s]", ylabel = "average density in the perovskite [m⁻³]"
        )
        scalarplot!(
            vis[1, 1], averages.time[it], averages.holes[it];
            xscale = :log, yscale = :log, color = :blue, label = "holes n_p", legend = :rc, clear = false
        )
        scalarplot!(
            vis[2, 1], metrics.time, metrics.PCE;
            markershape = :circle, xlabel = "time [s]", ylabel = "PCE [%]",
            title = "PCE of the non-perturbing J-V scans"
        )
        reveal(vis)
    end

    if test == false
        for (tk, pce) in zip(metrics.time, metrics.PCE)
            println("t = $tk s: PCE = $pce %")
        end
        println("*** done\n")
    end

    if test == false
        println("Average vacancy density at t = 0:    ", averages.vacancies[1])
        println("Average vacancy density at t = $(sol.t[end]) s: ", averages.vacancies[end])
        println(" ")
    end

    if test
        testval = sum(filter(!isnan, sol.u[end])) / length(sol.u[end]) # when using sparse storage, we get NaN values in solution
        return [testval, metrics.PCE[1]]
    end

    return (; sol, grid, ctsys, data = ctsys.data, averages, metrics)

end #  main

function test()
    ## [mean of the solution at t = 1 s, PCE [%] of the scan at the end of the light ramp]
    testval = [-0.6632438731028629, 15.211014883535208]
    return main(test = true, tend = 1.0) ≈ testval
end

end # module
