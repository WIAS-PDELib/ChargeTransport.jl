### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ b0138526-c79e-11ec-041a-156b0dfee367
# Load required packages.
begin
    import Pkg
    Pkg.activate(@__DIR__) # activate a temporary environment
    # use most current version of ChargeTransport.jl
    Pkg.add(path = "https://github.com/WIAS-PDELib/ChargeTransport.jl.git")
    using ChargeTransport
    using ExtendableGrids
    using LessUnitful: @unitfactors # this is needed for Pluto notebooks as @local_unitfactors do not work here.
    using PlutoVista
    using PlutoUI
    nothing
end

# ╔═╡ 39e1f60f-cd7a-49b5-b569-b3321f68c2ac
md"""
# Interactive 1D perovskite solar cell example.
"""

# ╔═╡ 8d25fad0-ff69-4bd5-ae01-d10b36fef92f
begin
    parameter_set = Params_PSC_TiO2_MAPI_spiro # choose the parameter set
    # parameters
    p = parameter_set()

    @unitfactors nm μm cm m s ns V K ps Hz W

    # constants
    constants = ChargeTransport.constants
    (; q, k_B, ε_0) = constants

    eV = q * V

    nothing
end

# ╔═╡ 6511e625-2af2-44c9-bc5c-d24e08109c3f
TableOfContents(; aside = false, depth = 5)

# ╔═╡ 997e8130-8e2d-45f6-a6b9-5ed78782d2b0
md"""
The purpose of this notebook is to provide a user-friendly and interactive introduction to investigating solar cell device behavior using Chargetransport.jl.
General information on the drift-diffusion charge transport model can be found in the [package documentation](https://wias-pdelib.github.io/ChargeTransport.jl/stable/) or in [this publication](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865?via%3Dihub).

It is possible to generalize to different
- dimensions
- material compositions and layer configurations
- contact and scan protocol types
- ...
"""

# ╔═╡ 1284fff2-af76-4d53-9444-233bde7cfaa9
md"""
# Device set-up
"""

# ╔═╡ c0f9f418-0c3b-428f-bd5b-e4d3d1ab9be2
md"""
For illustrative purposes we use a one-dimensional non-uniform grid. The used parameters are estimated, but represent the following set-up (n-i-p):

                            Ti02 | MAPI | spiro-OMeTAD.

We assume ohmic contacts on both sides of the boundary and the voltage is applied at the right acceptor boundary.

"""

# ╔═╡ 18823ee0-988d-459c-b340-7dfed6092b22
begin

    n = 8
    ################################################################################

    x0 = 0.0 * cm
    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u = collect(range(x0, p.h_ndoping / 2, step = p.h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        p.h_ndoping / 2,
        p.h_ndoping,
        p.h_ndoping / (0.7 * δ),
        p.h_ndoping / (1.0 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        p.h_ndoping,
        p.h_ndoping + p.h_intrinsic / k,
        p.h_intrinsic / (2.6 * δ),
        p.h_intrinsic / (1.8 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        p.h_ndoping + p.h_intrinsic / k,
        p.h_ndoping + p.h_intrinsic,
        p.h_intrinsic / (1.8 * δ),
        p.h_intrinsic / (2.6 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        p.h_ndoping + p.h_intrinsic,
        p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2,
        p.h_pdoping / (1.4 * δ),
        p.h_pdoping / (0.9 * δ),
        tol = t
    )
    coord_p_u = collect(range(p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, p.h_ndoping + p.h_intrinsic + p.h_pdoping, step = p.h_pdoping / (1.1 * δ)))

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

    bfacemask!(grid, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1, tol = 1.0e-18)
    bfacemask!(grid, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2, tol = 1.0e-18)
    nothing
end

# ╔═╡ bb701b88-6a76-4e77-ad81-a776ebfa5ec4
md"""
In total, we have $(length(coord)) number of nodes.

"""

# ╔═╡ fcaa2a1a-cc66-4c72-9dc8-6bf43b43830f
gridplot(grid, Plotter = PlutoVista, legend = :lt)

# ╔═╡ febcc48f-82ae-408a-b6a8-74cf3c76f61a
md"""
## Used parameter set
"""

# ╔═╡ 2d0bfda7-c23d-45b0-b3db-2000f1a291c3
md"""
We consider the device at a constant temperature T = 298 K with the following parameters for electrons, holes and anion vacancies. We assume that the device has an area of 0.1m x 0.1m.
"""

# ╔═╡ 3685ff11-b2f1-4fc1-bde6-723bf59ca57f
md"""
| physical quantity |   $\quad$symbol$\quad$   |    $\quad$ETL$\quad$    | $\quad$intrinsic$\quad$ |    $\quad$HTL$\quad$    |   $\quad$unit$\quad$   |
|-------------------|--------------------------|-------------------------|------------------------------------|-------------------------|------------------------|
|  layer thickness  |                 | $((p.h_ndoping)/nm) | $((p.h_intrinsic)/nm) |  $((p.h_pdoping)/nm) |  nm  |
|  conduction band-edge energy| $E_\text{c}$ | $(round(p.En[1]/q, sigdigits=3))     | $(round(p.En[2]/q, sigdigits=3))   | $(round(p.En[3]/q, sigdigits=3))     | eV  |
|  valence band-edge energy|    $E_\text{v}$ | $(round(p.Ep[1]/q, sigdigits=3))     | $(round(p.Ep[2]/q, sigdigits=3))   | $(round(p.Ep[3]/q, sigdigits=3))      | eV  |
|  conduction band-edge DOS|  $N_\text{c}$   | $(round(p.Nn[1]*cm^3, sigdigits=3))   | $(round(p.Nn[2]*cm^3, sigdigits=3))   | $(round(p.Nn[3]*cm^3, sigdigits=3))  | $1/\text{cm}^3$  |
|  valence band-edge DOS|    $N_\text{v}$    | $(round(p.Np[1]*cm^3, sigdigits=3))   | $(round(p.Np[2]*cm^3, sigdigits=3)) | $(round(p.Np[3]*cm^3, sigdigits=3))  | $1/\text{cm}^3$  |
|  max vacancy concentration|    $N_\text{a}$    | --  | $(round(p.Na_i*cm^3, sigdigits=3))  | --  | $1/\text{cm}^3$  |
|  electron mobility|  $\mu_\text{n}$   | $(round(p.μn[1]/cm^2, sigdigits = 3))   | $(p.μn[2]/cm^2)  | $(round(p.μn[3]/cm^2, sigdigits = 3))  | $\text{cm}^2/(Vs)$  |
|  hole mobility|    $\mu_\text{p}$    | $(p.μp[1]/cm^2)   | $(p.μp[2]/cm^2)   | $(round(p.μp[3]/cm^2, sigdigits = 3))   | $\text{cm}^2/(Vs)$  |
|  vacancy mobility|    $\mu_\text{a}$    | --  | $(p.μa[2]/cm^2)  | --  | $\text{cm}^2/(Vs)$  |
|  material permittivity|    $\varepsilon_r$ | $(p.ε[1])  | $(p.ε[2])  |$(p.ε[3])   |  |
| donor doping density|    $C_\text{n}$    | $(round(p.Cn*cm^3, sigdigits=3))  | --  | --  | $1/\text{cm}^3$ |
| acceptor doping density|    $C_\text{p}$    | --  | --  | $(round(p.Cp*cm^3, sigdigits=3))  | $1/\text{cm}^3$ |
| average vacancy density|    $C_\text{a}$    | --  | $(round(p.Ca*cm^3, sigdigits=3))  | --  | $1/\text{cm}^3$ |
"""


# ╔═╡ c72f9732-4fb3-485f-93e2-14999307d513
begin

    # solar cell area
    area = 0.1 * m * 0.1 * m
    nothing
end

# ╔═╡ c19b9329-1c9e-4ab3-8216-ff50dcb89e19
md"""
## Bulk recombination
"""

# ╔═╡ 2a43a7a8-b930-4fca-bf50-90220a3bb431
md"""
You have the option of turning the bulk recombination on or off.
"""

# ╔═╡ 5f0c398c-f389-4355-9f0e-6710c2b819a0
md"""
Turn bulk recombination $(@bind isBulkRecoOn Select(["true" => "on", "false" => "off"], default="true")).
"""

# ╔═╡ 004e4101-dbf3-4527-bfb8-47da01c70182
md"""
If turned on, we use the following predefined parameters.
"""

# ╔═╡ b0df6373-8d11-4581-9e42-62e3aee6c869
begin

    ## reference densities
    nτ = [7.93786e8, 6.15383e9, 6.01793e14] ./ (m^3)
    pτ = [1.81784e12  2.1087e11 1.14747e8] ./ (m^3)

    nothing

end

# ╔═╡ 66de390b-d0fd-48c4-93d4-50cb7bfc50cd
md"""
| physical quantity |   $\quad$symbol$\quad$   |    $\quad$ETL$\quad$    | $\quad$intrinsic$\quad$ |    $\quad$HTL$\quad$    |   $\quad$unit$\quad$   |
|-------------------|--------------------------|-------------------------|------------------------------------|-------------------------|------------------------|
| radiative recomb | $r_0$   | $(round(p.r0[1]/cm^3, sigdigits=3)) | $(round(p.r0[2]/cm^3, sigdigits=3)) |  $(round(p.r0[3]/cm^3, sigdigits=3)) |  $\text{cm}^3/s$  |
| lifetime, electron | $\tau_\text{n}$   | $(p.τn[1]) | $(p.τn[2])|  $(p.τn[3]) | s  |
| lifetime, hole | $\tau_\text{p}$   | $(p.τp[1]) | $(p.τp[2])|  $(p.τp[3]) | s  |
| reference dens, electron | $n_{\tau,\text{n}}$   | $(nτ[1]) | $(nτ[2])|  $(nτ[3]) | $1/\text{m}^3$   |
| reference dens, hole | $n_{\tau, \text{p}}$   | $(pτ[1]) | $(pτ[2])|  $(pτ[3]) | $1/\text{m}^3$   |
"""

# ╔═╡ 61155934-83a1-40ea-805e-2607fd8f9cd2
md"""
## Simulated scan protocol
"""

# ╔═╡ 380c5124-bd09-421f-9590-416400624374
md"""
For simplicity, we use a linear forward and backward scan protocol, where you can vary the scan rate. In general, you can use any function you want (preconditioning, sinusoidal signal, ...) as long as you do not have convergence issues of course.
"""

# ╔═╡ c0177ece-44e4-4431-abd6-4f28a8037d26
begin
    scanRateSlide = 0.2:0.1:0.8
    nothing
end

# ╔═╡ afcc0f3b-47aa-4418-aead-67202fc56119
md"""
Adjust the scanrate (in V/s) $(@bind scanrate  PlutoUI.Slider(scanRateSlide, default=0.4,show_value=true))

"""

# ╔═╡ 719df739-d611-4855-870c-64dc82444538
begin
    ## primary data for I-V scan protocol
    voltageAcceptor = 1.2 * V         # contact voltage
    endVoltage = voltageAcceptor # bias goes until the given voltage
    tend = endVoltage / scanrate

    ## Define scan protocol function
    function scanProtocol(t)

        if 0.0 <= t  && t <= tend
            biasVal = 0.0 + scanrate * t
        elseif t > tend  && t <= 2 * tend
            biasVal = scanrate * tend .+ scanrate * (tend - t)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    const contactVoltageFunction = [zeroVoltage, scanProtocol]

    nothing

end

# ╔═╡ 110f405e-7a3a-4385-82ef-9a1c3886bc98
begin
    tPlot = 0.0:0.01:(2 * tend)
    plot(tPlot, scanProtocol.(tPlot), xlabel = "time [s]", ylabel = "applied voltage [V]", xlimits = (0, 12))
end


# ╔═╡ c305fbd9-79ff-4ea3-8011-b07246f767c6
md"""
## Photogeneration based on Beer-Lambert
"""

# ╔═╡ bb4e474b-2028-4982-b202-3b22aa08c9d1
md"""
There are currently two methods implemented: a predefined user generation, a uniform photogeneration rate and a Beer-Lambert photogeneration rate. In this example we use the latter one.
"""

# ╔═╡ 29d920d6-ac56-4cb1-83a3-741ee6c876ae
begin
    function BeerLambert2(ctsys, ireg, node)
        params = ctsys.fvmsys.physics.data.params

        return params.generationIncidentPhotonFlux[ireg] .* params.generationAbsorption[ireg] .* exp.(- params.invertedIllumination .* params.generationAbsorption[ireg] .* (node' .- params.generationPeak))
    end
    nothing
end

# ╔═╡ 850bfa1a-4326-4aa3-97d2-c21d4ac2ba11
md"""
## Finite volume discretization
"""

# ╔═╡ 2b785f7b-4a96-4785-993b-ee7bbf6b0533
md"""
Using this predefined information, a finite volume discretization scheme is created for the model.

The solutions and the respective I-V characteristics are calculated with help of [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl).
"""

# ╔═╡ c00b33f1-8722-49e9-93b3-3703c5d0efb7
begin
    ## Initialize Data instance and fill in predefined data
    data = Data(grid, p.numberOfCarriers, contactVoltageFunction = contactVoltageFunction)
    data.modelType = Transient
    data.F = [Boltzmann, Boltzmann, FermiDiracMinusOne]
    if isBulkRecoOn == "true"
        data.bulkRecombination = set_bulk_recombination(;
            iphin = p.iphin, iphip = p.iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = true,
            bulk_recomb_SRH = true
        )
    elseif isBulkRecoOn == "false"
        data.bulkRecombination = set_bulk_recombination(;
            iphin = p.iphin, iphip = p.iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = false,
            bulk_recomb_SRH = false
        )
    end

    data.boundaryType[p.bregionAcceptor] = OhmicContact
    data.boundaryType[p.bregionDonor] = OhmicContact
    data.generationModel = GenerationBeerLambert
    data.fluxApproximation .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

    ######################################
    data.params = Params(p)

    data.params.generationIncidentPhotonFlux[2] = 1.5e19 / (m^2 * s)
    data.params.generationAbsorption[2] = 1.3e7 / m

    ctsys = System(grid, data, unknown_storage = :sparse)
    # otherwise there is a print in the terminal.
    nothing

end

# ╔═╡ 6f76b5cd-65a2-499c-a0be-0f3a0916e70c
begin
    vis = PlutoVistaPlot(resolution = (500, 300), titlefontsize = 20)
    #####
    ireg = 1
    subg = subgrid(grid, [ireg])
    plot!(vis, vec(subg[Coordinates]), vec(BeerLambert2(ctsys, ireg, subg[Coordinates])), label = "region $ireg")
    #####
    ireg = 2
    subg = subgrid(grid, [ireg])
    plot!(vis, vec(subg[Coordinates]), vec(BeerLambert2(ctsys, ireg, subg[Coordinates])), label = "region $ireg")
    #####
    ireg = 3
    subg = subgrid(grid, [ireg])
    plot!(vis, vec(subg[Coordinates]), vec(BeerLambert2(ctsys, ireg, subg[Coordinates])), label = "region $ireg", legend = :rt)

end

# ╔═╡ ab8d4426-9eda-4bab-a68c-1475042321db
begin
    control = SolverControl()
    control.maxiters = 300
    control.verbose = "" # due to allocations in julia 1.12
    control.max_round = 5
    control.damp_initial = 0.5
    control.damp_growth = 1.21 # >= 1
    control.Δt_max = 5.0e-2
    nothing
end

# ╔═╡ d6e4a543-d2e5-42a2-91af-1cf8b4d4632d
begin
    inival = unknowns(ctsys)
    ## calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution
    #############################################################
    # these values are needed for putting the generation slightly on
    I = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need to fix it in the bias loop,
    ## since we have no applied bias. Otherwise we get convergence errors
    ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop
    #############################################################
    ## put here back the homogeneous Neumann boundary conditions.
    ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 0.0
    ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

    sol = solve(ctsys, inival = inival, times = (0.0, tend), control = control)
    #############################################################
    inivalReverse = sol(tend)
    solReverse = solve(ctsys, inival = inivalReverse, times = (tend, 2 * tend), control = control)
    nothing

end

# ╔═╡ 9fa7dc02-4913-4b6e-a96d-0d73ccfee302
begin
    factory = TestFunctionFactory(ctsys)
    tf = testfunction(factory, [p.bregionDonor], [p.bregionAcceptor])

    tvalues = sol.t
    number_tsteps = length(tvalues)
    biasValues = scanProtocol.(tvalues)
    IV = zeros(0)

    for istep in 2:number_tsteps
        Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
        inival = sol[istep - 1]
        solution = sol[istep]

        I = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii in 1:(p.numberOfCarriers + 1)
            current = current + I[ii]
        end

        push!(IV, current)

    end

    tvaluesReverse = solReverse.t
    number_tstepsReverse = length(tvaluesReverse)
    biasValuesReverse = scanProtocol.(tvaluesReverse)
    IVReverse = zeros(0)

    for istep in 2:number_tstepsReverse
        Δt = tvaluesReverse[istep] - tvaluesReverse[istep - 1] # Time step size
        inival = solReverse[istep - 1]
        solution = solReverse[istep]

        I = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii in 1:(p.numberOfCarriers + 1)
            current = current + I[ii]
        end

        push!(IVReverse, current)

    end

end

# ╔═╡ 9b20814b-7f00-49da-ad83-63f1612d2f27
md"""
# Visualization of results
"""

# ╔═╡ 4a4ce5d2-c248-4d88-b4de-ba42f244e0e5
md"""
Once the finite volume system is solved, we can postprocess the solution. Here, the charge concentrations of electrons, holes and anion vacancies as well as the band-edges can be shown with respect to time. Further, forward and backward I-V characteristics are visualized.
"""

# ╔═╡ 09f362b4-a226-4b4e-8253-a6e87d979777
begin
    totalTime = glue(tvalues, tvaluesReverse)
    totalBias = scanProtocol.(totalTime)


    md"""
    Plot the concentrations and band-edge energies at a given time (in seconds)

    $(@bind printTime  Slider(totalTime, default=0.0,show_value=true))

    """
end

# ╔═╡ 3f40a78e-4b2e-4863-8404-75157d562459
begin
    indexPrintTime = findall(x -> x == printTime, totalTime)
    printBias = totalBias[indexPrintTime][1]
    nothing
end

# ╔═╡ 8dbd2c62-eae3-4857-8de7-8829554a847a
md"""
The corresponding voltage is Δu = $(printBias) V.
"""

# ╔═╡ 01fdd92e-4033-4701-a5d4-7012c7c6c063
begin
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

    ## add labels for anion vacancy
    label_energy[1, p.iphia] = "\$E_a-q\\psi\$"
    label_energy[2, p.iphia] = "\$ - q \\varphi_a\$"
    label_BEE[p.iphia] = "\$E_a\$"
    label_density[p.iphia] = "\$ n_a \$"
    label_solution[p.iphia] = "\$ \\varphi_a\$"
    nothing
end

# ╔═╡ 882f0ca5-e937-41c4-820b-23eae12e7c19
md"""
## IV curves
"""

# ╔═╡ 26255f99-1a36-43bf-b6f6-01cfa1e1c396
begin
    bias = biasValues[2:end]

    powerDensity = bias .* (-IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, -IV)

    IncidentLightPowerDens = 1000.0 * W / m^2 # for one sun

    efficiency = bias[indexPD] * -IV[indexPD] / (IncidentLightPowerDens * area) * 100
    fillfactor = (bias[indexPD] * -IV[indexPD]) / (-IV[1] * open_circuit) * 100

    visX = PlutoVistaPlot(resolution = (500, 300), titlefontsize = 20)

    plot!(visX, biasValues[2:end], -IV * (area * cm) * 1.0e3, linewidth = 5, label = "forward")
    plot!(visX, biasValuesReverse[2:end], -IVReverse * (area * cm) * 1.0e3, linestyle = ":", linewidth = 5, label = "reverse", legend = :rt, xlabel = "applied bias [V]", ylabel = "current density [Acm\$^{-2} \$]")

end

# ╔═╡ ce3423b6-8005-49da-b7c6-2e16e0675740
md"""
In case of 1 Sun and an solar cell area of 0.1m x 0.1m we receive the following values.

The open circuit voltage is $open_circuit V.

The fill factor is $fillfactor %.

The efficiency  is $efficiency %.

"""

# ╔═╡ 557b00c2-6a4b-4071-ba90-99b275dcadd8
begin

    function plot_densities2(Plotter, ctsys, solution, title, label_density, ; plotGridpoints = false)

        grid = ctsys.fvmsys.grid
        data = ctsys.fvmsys.physics.data
        numberOfRegions = grid[NumCellRegions]

        if plotGridpoints == true
            marker = "o"
        else
            marker = ""
        end

        params = data.params
        colors = ["green", "red", "gold", "purple", "orange"]

        vis2 = PlutoVistaPlot(resolution = (500, 300), titlefontsize = 20)

        ####### electrons
        icc = p.iphin
        label_icc = label_density[icc]
        ## ireg = 1
        ireg = 1
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 2
        ireg = 2
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 3
        ireg = 3
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")

        ####### holes
        icc = p.iphip
        label_icc = label_density[icc]
        ## ireg = 1
        ireg = 1
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 2
        ireg = 2
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 3
        ireg = 3
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]", title = title)

        #############################################################
        vis3 = PlutoVistaPlot(resolution = (500, 300), titlefontsize = 20)

        ####### vacancies
        icc = p.iphia
        label_icc = label_density[icc]
        ## ireg = 1
        ireg = 1
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 2
        ireg = 2
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]")
        ## ireg = 3
        ireg = 3
        subg = subgrid(grid, [ireg])
        ncc = get_density(solution, ireg, ctsys, icc)
        return plot!(vis2, vec(subg[Coordinates]), vec(1.0e-6 .* ncc), yscale = :log, marker = marker, label = label_icc, color = colors[icc], linewidth = 2, xlabel = "space [\$m\$]", title = title)

    end
    nothing
end

# ╔═╡ c13bb60b-b07c-4c56-b4f5-1fccbb194bb4
begin
    title = "Time value = $(printTime) s; bias value = $(printBias) V"
    if printTime <= tend
        printSol = sol(printTime)
    else
        printSol = solReverse(printTime)
    end
    plot_densities2(PlutoVista, ctsys, printSol, title, label_density, ; plotGridpoints = false)
end

# ╔═╡ 633ed076-9123-4989-b7e0-3ee078d1a7e0
md"""
# Parameter study example
"""

# ╔═╡ a6e14180-7eeb-48e7-afd2-8a147b32d870
md"""
Parameter studies can also be done with ChargeTransport.jl. For instance, we can find out how the band gap in the perovskite layer affects the fill factor, the open circuit voltage, and the power conversion efficiency.
"""

# ╔═╡ d03cb3b7-2d90-4ee3-842d-54d55e32db07
begin
    FillfactorVec = zeros(0)
    OpenCircuitVec = zeros(0)
    EfficiencyVec = zeros(0)
    EgTest = 1.3:0.05:2.1

    for Eg in EgTest

        Ev_iNew = p.En[p.regionIntrinsic] - Eg * eV
        ctsys.fvmsys.physics.data.params.bandEdgeEnergy[p.iphip, p.regionIntrinsic] = Ev_iNew

        ## calculate equilibrium solution and as initial guess
        local solution = equilibrium_solve!(ctsys, control = control)
        local inival = solution

        #######################################################################
        ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 1.0e30
        ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

        for istep in 1:(length(I) - 1)

            ## turn slowly generation on
            ctsys.data.λ2 = LAMBDA[istep + 1]
            solution = solve(ctsys, inival = inival, control = control)
            inival = solution

        end # generation loop

        ## put here back the homogeneous Neumann boundary conditions.
        ctsys.fvmsys.boundary_factors[p.iphia, p.bregionJ2] = 0.0
        ctsys.fvmsys.boundary_values[p.iphia, p.bregionJ2] = 0.0

        local sol = solve(ctsys, inival = inival, times = (0.0, tend), control = control)
        ###########################################################################

        local tvalues = sol.t
        local number_tsteps = length(tvalues)
        local biasValues = scanProtocol.(tvalues)
        local IV = zeros(0)

        for istep in 2:number_tsteps
            Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
            inival = sol[istep - 1]
            solution = sol[istep]

            local I = integrate(ctsys, tf, solution, inival, Δt)

            current = 0.0
            for ii in 1:(p.numberOfCarriers + 1)
                current = current + I[ii]
            end

            push!(IV, current)

        end
        ###########################################################################

        local bias = biasValues[2:end]
        local IV = -IV
        local powerDensity = bias .* (IV)           # power density function
        local MaxPD, indexPD = findmax(powerDensity)
        local open_circuit = compute_open_circuit_voltage(bias, IV)
        local IncLightPowerDens = 1000.0 * W / m^2 # for one sun

        local efficiency = bias[indexPD] * IV[indexPD] / (IncLightPowerDens * area) * 100
        local fillfactor = (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit) * 100

        push!(FillfactorVec, fillfactor)
        push!(OpenCircuitVec, open_circuit)
        push!(EfficiencyVec, efficiency)

    end # loop Eg
end

# ╔═╡ ccf30353-f057-47ad-8e69-e12ef4e01c00
begin
    plot(vec(EgTest), vec(EfficiencyVec), marker = "o", xlabel = "band gap (perovskite) [eV]", ylabel = "efficiency \$ \\eta \$ [%]")
end

# ╔═╡ 30ca5f4b-7303-4c43-bdfa-c5388603f091
begin
    plot(vec(EgTest), vec(OpenCircuitVec), marker = "o", xlabel = "band gap (perovskite) [eV]", ylabel = "open circuit voltage [V]")
end

# ╔═╡ 2c2d9521-5511-4766-9f27-cd524e7124f3
begin
    plot(vec(EgTest), vec(FillfactorVec), marker = "o", xlabel = "band gap (perovskite) [eV]", ylabel = "Fill factor [%]")
end

# ╔═╡ Cell order:
# ╟─39e1f60f-cd7a-49b5-b569-b3321f68c2ac
# ╠═b0138526-c79e-11ec-041a-156b0dfee367
# ╠═8d25fad0-ff69-4bd5-ae01-d10b36fef92f
# ╟─6511e625-2af2-44c9-bc5c-d24e08109c3f
# ╟─997e8130-8e2d-45f6-a6b9-5ed78782d2b0
# ╟─1284fff2-af76-4d53-9444-233bde7cfaa9
# ╟─c0f9f418-0c3b-428f-bd5b-e4d3d1ab9be2
# ╟─bb701b88-6a76-4e77-ad81-a776ebfa5ec4
# ╟─18823ee0-988d-459c-b340-7dfed6092b22
# ╟─fcaa2a1a-cc66-4c72-9dc8-6bf43b43830f
# ╟─febcc48f-82ae-408a-b6a8-74cf3c76f61a
# ╟─2d0bfda7-c23d-45b0-b3db-2000f1a291c3
# ╟─3685ff11-b2f1-4fc1-bde6-723bf59ca57f
# ╟─c72f9732-4fb3-485f-93e2-14999307d513
# ╟─c19b9329-1c9e-4ab3-8216-ff50dcb89e19
# ╟─2a43a7a8-b930-4fca-bf50-90220a3bb431
# ╟─5f0c398c-f389-4355-9f0e-6710c2b819a0
# ╟─004e4101-dbf3-4527-bfb8-47da01c70182
# ╟─b0df6373-8d11-4581-9e42-62e3aee6c869
# ╟─66de390b-d0fd-48c4-93d4-50cb7bfc50cd
# ╟─61155934-83a1-40ea-805e-2607fd8f9cd2
# ╟─380c5124-bd09-421f-9590-416400624374
# ╟─c0177ece-44e4-4431-abd6-4f28a8037d26
# ╟─afcc0f3b-47aa-4418-aead-67202fc56119
# ╟─110f405e-7a3a-4385-82ef-9a1c3886bc98
# ╟─719df739-d611-4855-870c-64dc82444538
# ╟─c305fbd9-79ff-4ea3-8011-b07246f767c6
# ╟─bb4e474b-2028-4982-b202-3b22aa08c9d1
# ╟─6f76b5cd-65a2-499c-a0be-0f3a0916e70c
# ╟─29d920d6-ac56-4cb1-83a3-741ee6c876ae
# ╟─850bfa1a-4326-4aa3-97d2-c21d4ac2ba11
# ╟─2b785f7b-4a96-4785-993b-ee7bbf6b0533
# ╟─c00b33f1-8722-49e9-93b3-3703c5d0efb7
# ╟─ab8d4426-9eda-4bab-a68c-1475042321db
# ╟─d6e4a543-d2e5-42a2-91af-1cf8b4d4632d
# ╟─9fa7dc02-4913-4b6e-a96d-0d73ccfee302
# ╟─9b20814b-7f00-49da-ad83-63f1612d2f27
# ╟─4a4ce5d2-c248-4d88-b4de-ba42f244e0e5
# ╠═3f40a78e-4b2e-4863-8404-75157d562459
# ╟─c13bb60b-b07c-4c56-b4f5-1fccbb194bb4
# ╟─09f362b4-a226-4b4e-8253-a6e87d979777
# ╟─8dbd2c62-eae3-4857-8de7-8829554a847a
# ╠═01fdd92e-4033-4701-a5d4-7012c7c6c063
# ╟─882f0ca5-e937-41c4-820b-23eae12e7c19
# ╟─26255f99-1a36-43bf-b6f6-01cfa1e1c396
# ╟─ce3423b6-8005-49da-b7c6-2e16e0675740
# ╟─557b00c2-6a4b-4071-ba90-99b275dcadd8
# ╟─633ed076-9123-4989-b7e0-3ee078d1a7e0
# ╟─a6e14180-7eeb-48e7-afd2-8a147b32d870
# ╟─ccf30353-f057-47ad-8e69-e12ef4e01c00
# ╟─30ca5f4b-7303-4c43-bdfa-c5388603f091
# ╟─2c2d9521-5511-4766-9f27-cd524e7124f3
# ╟─d03cb3b7-2d90-4ee3-842d-54d55e32db07
