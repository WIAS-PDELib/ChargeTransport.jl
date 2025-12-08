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
    # As default, we have for electrons and holes the Boltzmann statistics and for the ions the Fermi Dirac integral of order -1
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

The open circuit voltage is $(round(open_circuit, digits = 3)) V.

The fill factor is $(round(fillfactor, digits = 3)) %.

The efficiency  is $(round(efficiency, digits = 3)) %.

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ChargeTransport = "25c3eafe-d88c-11e9-3031-f396758f002a"
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
LessUnitful = "f29f6376-6e90-4d80-80c9-fb8ec61203d5"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"

[compat]
ChargeTransport = "~1.2.0"
ExtendableGrids = "~1.14.0"
LessUnitful = "~1.2.1"
PlutoUI = "~0.7.71"
PlutoVista = "~1.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "b7ed52be80051af150d0bcdc03a638596197b3f0"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "dbd8c3bbbdbb5c2778f85f4422c39960eac65a42"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.20.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "120e392af69350960b1d3b89d41dcc1d66543858"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "e35c672b239c5105f597963c33e740eeb46cf0ab"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.9.4"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChargeTransport]]
deps = ["Compat", "DocStringExtensions", "ExtendableGrids", "ForwardDiff", "GridVisualize", "Interpolations", "LessUnitful", "Printf", "Roots", "VoronoiFVM"]
git-tree-sha1 = "c17c01315ad87ca9dab60989874b31df98ab7102"
uuid = "25c3eafe-d88c-11e9-3031-f396758f002a"
version = "1.2.0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "cee1700673af54db57bd1c7fb834ad4ff31309a0"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.8"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "75e5697f521c9ab89816d3abeea806dfc5afb967"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.12"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Bijections", "Compat", "Dates", "DocStringExtensions", "ElasticArrays", "Graphs", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "StatsBase", "UUIDs", "WriteVTK"]
git-tree-sha1 = "815c12c33244415038b539a1ab1668a5f9033fd6"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "1.14.0"

    [deps.ExtendableGrids.extensions]
    ExtendableGridsGmshExt = "Gmsh"
    ExtendableGridsMetisExt = "Metis"
    ExtendableGridsTetGenExt = "TetGen"
    ExtendableGridsTriangulateExt = "Triangulate"

    [deps.ExtendableGrids.weakdeps]
    Gmsh = "705231aa-382f-11e9-3f0c-b7cb4346fdeb"
    Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b"
    TetGen = "c5d3f3f7-f850-59f6-8a2e-ffc6dc1317ea"
    Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "SparseArrays", "Sparspak", "StaticArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "c04d2c808f0b595dc3be5fa98d107213f81e77fb"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "1.7.1"

    [deps.ExtendableSparse.extensions]
    ExtendableSparseAMGCLWrapExt = "AMGCLWrap"
    ExtendableSparseAlgebraicMultigridExt = "AlgebraicMultigrid"
    ExtendableSparseIncompleteLUExt = "IncompleteLU"
    ExtendableSparseLinearSolveExt = "LinearSolve"
    ExtendableSparsePardisoExt = "Pardiso"

    [deps.ExtendableSparse.weakdeps]
    AMGCLWrap = "4f76b812-4ba5-496d-b042-d70715554288"
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
    IncompleteLU = "40713840-3770-5561-ab4c-a76e7d0d7895"
    LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "dc41303865a16274ecb8450c220021ce1e0cf05f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.2.1"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "Interpolations", "IntervalSets", "LinearAlgebra", "Observables", "OrderedCollections", "Printf", "StaticArrays"]
git-tree-sha1 = "1385a7536ecc73e955421e4f0a1e0d4a0a7ed8e0"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "1.15.4"

    [deps.GridVisualize.extensions]
    GridVisualizeMakieExt = "Makie"
    GridVisualizeMeshCatExt = "MeshCat"
    GridVisualizePlotsExt = "Plots"
    GridVisualizePlutoVistaExt = "PlutoVista"
    GridVisualizePyPlotExt = "PyPlot"
    GridVisualizeVTKViewExt = "VTKView"

    [deps.GridVisualize.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MeshCat = "283c5d60-a78f-5afe-a0af-af636b173e11"
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
    PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
    PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
    Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
    VTKView = "955f2c64-5fd0-11e9-0ad0-3332e913311a"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArrays", "StaticArraysCore"]
git-tree-sha1 = "7cfc079442c7bd2904bbfa32b76975054b06a639"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "3.0.2"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues", "TranscodingStreams"]
git-tree-sha1 = "d97791feefda45729613fafeccc4fbef3f539151"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.15"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "d1fc961038207e43982851e57ee257adc37be5e8"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.2"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "21057b6f4f5db1475e653735fda7d1de1c267b46"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.3"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LessUnitful]]
deps = ["PhysicalConstants", "Unitful"]
git-tree-sha1 = "6674eaca827a98225d140a6a6e5e2a0e5f46afec"
uuid = "f29f6376-6e90-4d80-80c9-fb8ec61203d5"
version = "1.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "d5d2e3abfb30ea9c2cff81d243e7235b51315ec2"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "6c22b14a5ea7fbcc140ea1f52f3cfe20d3da32e0"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.40.2"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "e5afce7eaf5b5ca0d444bcb4dc4fd78c54cbbac0"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.172"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "030f041d5502dbfa41f26f542aaac32bcbe89a64"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsMakieExt = "Makie"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "9be04ab1a9e54508348a2aacb5b5a4f04c85d397"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PlutoVista]]
deps = ["AbstractPlutoDingetjes", "ColorSchemes", "Colors", "DocStringExtensions", "GridVisualizeTools", "HypertextLiteral", "UUIDs"]
git-tree-sha1 = "2034ef3eb7e082c7eb46de790269d6b0e31f2718"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "1.2.1"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "96bef5b9ac123fff1b379acf0303cf914aaabdfd"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.37.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "SparseBandedMatrices", "StrideArraysCore", "TriangularSolve", "VectorizedRNG"]
git-tree-sha1 = "df8c2a4fa9bc79da18f7fbad5cc500127808bd9c"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.24"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "16fa030fb4bd4df373a677eca0460c3eee791ab2"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.120.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "024d829102878141aaee5cf8f8288bcabd2f57a0"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.7.2"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseBandedMatrices]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7bbeaac38a5b41674f49cb3161d92dbc88e66873"
uuid = "bd59d7e1-4699-4102-944e-d05209cb92aa"
version = "1.0.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "3c3a42a29f696f16273741ffe589b4003f539088"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.1.0"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerChainRulesCoreExt = "ChainRulesCore"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "9de43e0b9b976f1019bf7a879a686c4514520078"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.21"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "fcab7ea5354ffa3da57751c9a552fed0e3bcbda9"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.14"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "1e44e7b1dbb5249876d84c32466f8988a6b41bbb"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "617400a198bd433f921ca2a4e89999f835dd3fde"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.45"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TextWrap]]
git-tree-sha1 = "43044b737fa70bc12f6105061d3da38f881a3e3c"
uuid = "b718987f-49a8-5099-9789-dcd902bef87d"
version = "1.0.2"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "cec2df8cf14e0844a8c4d770d12347fda5931d72"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "d1d9a935a26c475ebffd54e9c7ad11627c43ea85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.72"

[[deps.VectorizedRNG]]
deps = ["Distributed", "Random", "SLEEFPirates", "UnPack", "VectorizationBase"]
git-tree-sha1 = "5ca83562ba95272d8709c6c91e31e23c3c4c9825"
uuid = "33b4df10-0173-11e9-2a0c-851a7edac40e"
version = "0.2.25"
weakdeps = ["Requires", "StaticArraysCore"]

    [deps.VectorizedRNG.extensions]
    VectorizedRNGStaticArraysExt = ["StaticArraysCore"]

[[deps.VoronoiFVM]]
deps = ["ADTypes", "BandedMatrices", "Colors", "CommonSolve", "DiffResults", "DifferentiationInterface", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "RecursiveFactorization", "SciMLBase", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "StaticArrays", "Statistics", "TextWrap"]
git-tree-sha1 = "f93d632518d6ecc4c463a9a1e976ddd7b9b34a80"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "2.12.2"

    [deps.VoronoiFVM.extensions]
    VoronoiFVMExtendableFEMBaseExt = "ExtendableFEMBase"

    [deps.VoronoiFVM.weakdeps]
    ExtendableFEMBase = "12fb9182-3d4c-4424-8fd1-727a0899810c"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "a329e0b6310244173690d6a4dfc6d1141f9b9370"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.21.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "c217bad8fccb3bbfef7d7902326eacfbd0d702ad"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.14.4+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
