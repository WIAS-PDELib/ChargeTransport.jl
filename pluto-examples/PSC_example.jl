### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    return quote
        local iv = try
            Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b0138526-c79e-11ec-041a-156b0dfee367
# Load required packages.
begin
    using Pkg
    #ENV["PYTHON"] = "" # if empty, create own Conda environment for Julia
    #using PyCall
    #using Conda
    #Pkg.build("PyCall");
    #Conda.add("matplotlib");
    #pyimport("matplotlib");
    using ChargeTransport
    using ExtendableGrids
    using PyPlot
    using PlutoUI
    #Pkg.resolve()
    # for convenience
    parametersdir = ChargeTransport.parametersdir
    include(parametersdir("Params_PSC_TiO2_MAPI_spiro.jl"))
    nothing
end

# ╔═╡ 39e1f60f-cd7a-49b5-b569-b3321f68c2ac
md"""
# Interactive 1D perovskite solar cell example.
"""

# ╔═╡ 6511e625-2af2-44c9-bc5c-d24e08109c3f
TableOfContents(; aside = false, depth = 5)

# ╔═╡ 997e8130-8e2d-45f6-a6b9-5ed78782d2b0
md"""
The purpose of this notebook is to provide a user-friendly and interactive introduction to investigating solar cell device behavior using Chargetransport.jl.
General information on the drift-diffusion charge transport model can be found in the [package documentation](https://patriciofarrell.github.io/ChargeTransport.jl/stable/PSC/) or in [this publication](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865?via%3Dihub).

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
    ## region numbers
    regionDonor = 1                           # n doped region
    regionIntrinsic = 2                           # intrinsic region
    regionAcceptor = 3                           # p doped region
    regions = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions = length(regions)

    ## boundary region numbers
    # By convention we have 1 for the left boundary and 2 for the right boundary. If
    # adding additional interior boundaries, continue with 3, 4, ...
    bregionDonor = 1
    bregionAcceptor = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    ## grid
    h_ndoping = 1.0e-5 * cm
    h_intrinsic = 4.0e-5 * cm
    h_pdoping = 2.0e-5 * cm
    heightLayers = [
        h_ndoping,
        h_ndoping + h_intrinsic,
        h_ndoping + h_intrinsic + h_pdoping,
    ]

    x0 = 0.0 * cm
    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u = collect(range(x0, h_ndoping / 2, step = h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        h_ndoping / 2,
        h_ndoping,
        h_ndoping / (0.7 * δ),
        h_ndoping / (1.0 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        h_ndoping,
        h_ndoping + h_intrinsic / k,
        h_intrinsic / (2.6 * δ),
        h_intrinsic / (1.8 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        h_ndoping + h_intrinsic / k,
        h_ndoping + h_intrinsic,
        h_intrinsic / (1.8 * δ),
        h_intrinsic / (2.6 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        h_ndoping + h_intrinsic,
        h_ndoping + h_intrinsic + h_pdoping / 2,
        h_pdoping / (1.4 * δ),
        h_pdoping / (0.9 * δ),
        tol = t
    )
    coord_p_u = collect(range(h_ndoping + h_intrinsic + h_pdoping / 2, h_ndoping + h_intrinsic + h_pdoping, step = h_pdoping / (1.1 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJunction1, tol = 1.0e-18)
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJunction2, tol = 1.0e-18)

    gridplot(grid, Plotter = PyPlot, legend = :lt)
end

# ╔═╡ bb701b88-6a76-4e77-ad81-a776ebfa5ec4
md"""
In total, we have $(length(coord)) number of nodes.

"""

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
|  layer thickness  |                 | $((h_ndoping)/nm) | $((h_intrinsic)/nm) |  $((h_pdoping)/nm) |  nm  |
|  conduction band-edge energy| $E_c$ | $(round(En[1]/q, sigdigits=3))     | $(round(En[2]/q, sigdigits=3))   | $(round(En[3]/q, sigdigits=3))     | eV  |
|  valence band-edge energy|    $E_v$ | $(round(Ep[1]/q, sigdigits=3))     | $(round(Ep[2]/q, sigdigits=3))   | $(round(Ep[3]/q, sigdigits=3))      | eV  |
|  conduction band-edge DOS|  $N_c$   | $(round(Nn[1]*cm^3, sigdigits=3))   | $(round(Nn[2]*cm^3, sigdigits=3))   | $(round(Nn[3]*cm^3, sigdigits=3))  | $1/\text{cm}^3$  |
|  valence band-edge DOS|    $N_v$    | $(round(Np[1]*cm^3, sigdigits=3))   | $(round(Np[2]*cm^3, sigdigits=3)) | $(round(Np[3]*cm^3, sigdigits=3))  | $1/\text{cm}^3$  |
|  max vacancy concentration|    $N_x$    | --  | $(round(Na_i*cm^3, sigdigits=3))  | --  | $1/\text{cm}^3$  |
|  electron mobility|  $\mu_n$   | $(μn[1]/cm^2)   | $(μn[2]/cm^2)  | $(μn[3]/cm^2)  | $\text{cm}^2/(Vs)$  |
|  hole mobility|    $\mu_p$    | $(μp[1]/cm^2)   | $(μp[2]/cm^2)   | $(μp[3]/cm^2)   | $\text{cm}^2/(Vs)$  |
|  vacancy mobility|    $\mu_a$    | --  | $(μa[2]/cm^2)  | --  | $\text{cm}^2/(Vs)$  |
|  electric permittivity|    $\varepsilon_r$ | $(ε[1])  | $(ε[2])  |$(ε[3])   |  |
| donor doping density|    $C_n$    | $(round(Cn*cm^3, sigdigits=3))  | --  | --  | $1/\text{cm}^3$ |
| acceptor doping density|    $C_p$    | --  | --  | $(round(Cp*cm^3, sigdigits=3))  | $1/\text{cm}^3$ |
| average vacancy density|    $C_a$    | --  | $(round(Ca*cm^3, sigdigits=3))  | --  | $1/\text{cm}^3$ |
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
| radiative recomb | $r_0$   | $(round(r0[1]/cm^3, sigdigits=3)) | $(round(r0[2]/cm^3, sigdigits=3)) |  $(round(r0[3]/cm^3, sigdigits=3)) |  $\text{cm}^3/s$  |
| lifetime, electron | $\tau_n$   | $(τn[1]) | $(τn[2])|  $(τn[3]) | s  |
| lifetime, hole | $\tau_p$   | $(τp[1]) | $(τp[2])|  $(τp[3]) | s  |
| reference dens, electron | $n_{\tau,n}$   | $(nτ[1]) | $(nτ[2])|  $(nτ[3]) | $1/\text{m}^3$   |
| reference dens, hole | $n_{\tau,p}$   | $(pτ[1]) | $(pτ[2])|  $(pτ[3]) | $1/\text{m}^3$   |
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
Adjust the scanrate (in V/s) $(@bind scanrate  Slider(scanRateSlide, default=0.4,show_value=true))

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
    PyPlot.clf()
    tPlot = 0.0:0.01:(2 * tend)
    PyPlot.plot(tPlot, scanProtocol.(tPlot))
    PyPlot.xlabel("time [s]")
    PyPlot.ylabel("applied voltage [V]")
    PyPlot.grid()
    PyPlot.gcf()
end


# ╔═╡ c305fbd9-79ff-4ea3-8011-b07246f767c6
md"""
## Photogeneration based on Beer-Lambert
"""

# ╔═╡ bb4e474b-2028-4982-b202-3b22aa08c9d1
md"""
There are currently two methods implemented: a uniform photogeneration rate and a Beer-Lambert photogeneration rate. In this example we use the latter one.

We plan to extend the software so that optical simulation rates can be read in. Get in touch with one of the authors if you're interested.
"""

# ╔═╡ 5573ea85-8bbb-411f-97b2-ef1e46981060
begin
    ## generation Beer Lambert
    photonflux_i = 1.5e19 / (m^2 * s)
    absorption_i = 1.3e7 / m
    incidentPhotonFlux = [0.0, photonflux_i, 0.0]
    absorption = [0.0, absorption_i, 0.0]
    generationPeak = h_ndoping

    ## generation uniform
    generation_a = 0.0
    generation_i = 8.2e20 / (cm^3 * s)
    generation_d = 0.0

    generation_uniform = [generation_a, generation_i, generation_d]

    nothing
end

# ╔═╡ 29d920d6-ac56-4cb1-83a3-741ee6c876ae
begin
    function BeerLambert(ctsys, ireg, node)
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
    data = Data(grid, numberOfCarriers, contactVoltageFunction = contactVoltageFunction)
    data.modelType = Transient
    data.F = [Boltzmann, Boltzmann, FermiDiracMinusOne]
    if isBulkRecoOn == "true"
        data.bulkRecombination = set_bulk_recombination(;
            iphin = iphin, iphip = iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = true,
            bulk_recomb_SRH = true
        )
    elseif isBulkRecoOn == "false"
        data.bulkRecombination = set_bulk_recombination(;
            iphin = iphin, iphip = iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = false,
            bulk_recomb_SRH = false
        )
    end

    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor] = OhmicContact
    data.generationModel = GenerationBeerLambert
    data.fluxApproximation .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])
    ###############################################################################
    ###############################################################################
    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp
    params.chargeNumbers[iphia] = za

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn[ireg]
        params.densityOfStates[iphip, ireg] = Np[ireg]
        params.densityOfStates[iphia, ireg] = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg] = En[ireg]
        params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg] = Ea[ireg]

        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]
        params.mobility[iphia, ireg] = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = nτ[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg] = pτ[ireg]

        ## generation parameters
        params.generationUniform[ireg] = generation_uniform[ireg]
        params.generationIncidentPhotonFlux[ireg] = incidentPhotonFlux[ireg]
        params.generationAbsorption[ireg] = absorption[ireg]

    end

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak = generationPeak

    ## interior doping
    params.doping[iphin, regionDonor] = Cn
    params.doping[iphia, regionIntrinsic] = Ca
    params.doping[iphip, regionAcceptor] = Cp

    #########################################################################

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)
    # otherwise there is a print in the terminal.
    nothing

end

# ╔═╡ 6f76b5cd-65a2-499c-a0be-0f3a0916e70c
begin
    PyPlot.clf()

    for ireg in 1:numberOfRegions
        subg = subgrid(grid, [ireg])
        PyPlot.plot(subg[Coordinates]', BeerLambert(ctsys, ireg, subg[Coordinates]), label = "region $ireg")
    end
    PyPlot.legend()
    PyPlot.grid()
    PyPlot.xlabel("space [\$m\$]")
    PyPlot.ylabel("photogeneration [\$\\frac{1}{cm^3s}\$]")
    PyPlot.gcf()
end

# ╔═╡ ab8d4426-9eda-4bab-a68c-1475042321db
begin
    control = SolverControl()
    control.maxiters = 300
    control.verbose = ""
    control.max_round = 5
    control.damp_initial = 0.5
    control.damp_growth = 1.21 # >= 1
    control.Δt_max = 5.0e-2
    nothing
end

# ╔═╡ d6e4a543-d2e5-42a2-91af-1cf8b4d4632d
begin
    ## calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution
    #############################################################
    # these values are needed for putting the generation slightly on
    I = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need to fix it in the bias loop,
    ## since we have no applied bias. Otherwise we get convergence errors
    ctsys.fvmsys.boundary_factors[iphia, bregionJunction2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJunction2] = 0.0

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop
    #############################################################
    ## put here back the homogeneous Neumann boundary conditions.
    ctsys.fvmsys.boundary_factors[iphia, bregionJunction2] = 0.0
    ctsys.fvmsys.boundary_values[iphia, bregionJunction2] = 0.0

    sol = solve(ctsys, inival = inival, times = (0.0, tend), control = control)
    #############################################################
    inivalReverse = sol(tend)
    solReverse = solve(ctsys, inival = inivalReverse, times = (tend, 2 * tend), control = control)
    nothing

end

# ╔═╡ 9fa7dc02-4913-4b6e-a96d-0d73ccfee302
begin
    factory = TestFunctionFactory(ctsys)
    tf = testfunction(factory, [bregionDonor], [bregionAcceptor])

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
        for ii in 1:(numberOfCarriers + 1)
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
        for ii in 1:(numberOfCarriers + 1)
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

    ## DA: Delete this once in main branch
    label_density[iphin] = "\$ n_n \$"
    label_density[iphip] = "\$ n_p \$"
    ## add labels for anion vacancy
    label_energy[1, iphia] = "\$E_a-q\\psi\$"
    label_energy[2, iphia] = "\$ - q \\varphi_a\$"
    label_BEE[iphia] = "\$E_a\$"
    label_density[iphia] = "\$ n_a \$"
    label_solution[iphia] = "\$ \\varphi_a\$"
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


    PyPlot.clf()
    PyPlot.plot(biasValues[2:end], -IV * (area * cm) * 1.0e3, linewidth = 5, label = "forward")
    PyPlot.plot(biasValuesReverse[2:end], -IVReverse * (area * cm) * 1.0e3, linestyle = ":", linewidth = 5, label = "reverse")

    #PyPlot.plot(biasValues[2:end], powerDensity*(area*cm)*1.0e3, linewidth = 3, label = "power density", linestyle = ":", color = "gray")
    PyPlot.grid()
    PyPlot.legend()
    PyPlot.xlabel("applied bias [V]")
    PyPlot.ylim(-0.2, 0.28)
    PyPlot.ylabel("current density [Acm\$^{-2} \$]")
    PyPlot.gcf()


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

        Plotter.clf()

        grid = ctsys.fvmsys.grid
        data = ctsys.fvmsys.physics.data
        numberOfRegions = grid[NumCellRegions]

        if dim_space(grid) > 1
            error("plot_densities is so far only tested in 1D")
        end

        if plotGridpoints == true
            marker = "o"
        else
            marker = ""
        end

        params = data.params
        colors = ["green", "red", "gold", "purple", "orange"]

        subplot(211)
        for icc in data.electricCarrierList

            # grids = Array{ExtendableGrid, 1}(undef, numberOfRegions)
            # nicc  = Array{Array{Float64, 1}, 1}(undef, numberOfRegions)

            ## first region for label
            label_icc = label_density[icc]
            subg = subgrid(grid, [1])
            ncc = get_density(solution, 1, ctsys, icc)

            Plotter.semilogy(subg[Coordinates]', 1.0e-6 .* ncc, marker = marker, label = label_icc, color = colors[icc], linewidth = 2)

            ## additional regions
            for ireg in 2:numberOfRegions
                subg = subgrid(grid, [ireg])
                ncc = get_density(solution, ireg, ctsys, icc)

                ## Note that this implies a 1D plot, for multidimensional plots, you may work with
                ## GridVisualize.jl or write your own code.
                Plotter.semilogy(subg[Coordinates]', 1.0e-6 .* ncc, marker = marker, color = colors[icc], linewidth = 2)
            end

        end

        Plotter.grid()
        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("density [\$\\frac{1}{cm^3}\$]")
        Plotter.legend(fancybox = true, loc = "best", fontsize = 11)
        PyPlot.ylim(2.0e-6, 1.0e19)
        PyPlot.tight_layout()
        Plotter.title(title)
        ####################################################################
        subplot(212)
        icc = iphia

        ## first region for label
        label_icc = label_density[icc]
        subg = subgrid(grid, [1])
        ncc = get_density(solution, 1, ctsys, icc)

        Plotter.semilogy(subg[Coordinates]', 1.0e-6 .* ncc, marker = marker, label = label_icc, color = colors[icc], linewidth = 2)

        ## additional regions
        for ireg in 2:numberOfRegions
            subg = subgrid(grid, [ireg])
            ncc = get_density(solution, ireg, ctsys, icc)

            ## Note that this implies a 1D plot, for multidimensional plots, you may work with
            ## GridVisualize.jl or write your own code.
            Plotter.semilogy(subg[Coordinates]', 1.0e-6 .* ncc, marker = marker, color = colors[icc], linewidth = 2)
        end
        Plotter.grid()
        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("density [\$\\frac{1}{cm^3}\$]")
        PyPlot.ylim(1.0e19, 2.5e19)
        PyPlot.tight_layout()
        return Plotter.legend(fancybox = true, loc = "best", fontsize = 11)

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
    PyPlot.clf()
    plot_densities2(PyPlot, ctsys, printSol, title, label_density, ; plotGridpoints = false)
    PyPlot.gcf()
end

# ╔═╡ f80b1946-73cb-448e-80a8-9e4770b50c79
begin
    function plot_energies2(Plotter, ctsys, solution, title, label_energy, ; plotGridpoints = false)

        Plotter.clf()

        grid = ctsys.fvmsys.grid
        data = ctsys.fvmsys.physics.data
        coord = grid[Coordinates]

        if length(coord[1]) != 1
            println("plot_energies is so far only implemented in 1D")
        end

        if plotGridpoints == true
            marker = "o"
        else
            marker = ""
        end

        colors = ["green", "red", "gold", "purple", "orange"]
        linestyles = ["-", ":", "--", "-.", "-"]

        for icc in data.electricCarrierList

            # grids = Array{ExtendableGrid, 1}(undef, numberOfRegions)
            # nicc  = Array{Array{Float64, 1}, 1}(undef, numberOfRegions)

            ## first region for label
            subg = subgrid(grid, [1])
            Ecc = get_BEE(icc, 1, ctsys)
            solpsi = view(solution[data.index_psi, :], subg)
            solcc = view(solution[icc, :], subg)

            Plotter.plot(subg[Coordinates]', Ecc ./ q .- solpsi, label = label_energy[1, icc], marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

            Plotter.plot(subg[Coordinates]', -solcc, label = label_energy[2, icc], marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[2])

            ## additional regions
            for ireg in 2:data.params.numberOfRegions
                subg = subgrid(grid, [ireg])
                Ecc = get_BEE(icc, ireg, ctsys)
                solpsi = view(solution[data.index_psi, :], subg)
                solcc = view(solution[icc, :], subg)

                ## Note that this implies a 1D plot, for multidimensional plots, you may work with
                ## GridVisualize.jl or write your own code.
                Plotter.plot(subg[Coordinates]', Ecc ./ q .- solpsi, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

                Plotter.plot(subg[Coordinates]', - solcc, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[2])
            end

        end

        for iicc in data.ionicCarrierList
            icc = iicc.ionicCarrier
            count = 0

            for ireg in 1:data.params.numberOfRegions
                if ireg ∈ iicc.regions
                    subg = subgrid(grid, [ireg])
                    Ecc = get_BEE(icc, ireg, ctsys)
                    solpsi = view(solution[data.index_psi, :], subg)
                    solcc = view(solution[icc, :], subg)

                    if count == 0
                        label1 = label_energy[1, icc]
                        label2 = label_energy[2, icc]
                    else
                        label1 = ""
                        label2 = ""
                    end
                    ## Note that this implies a 1D plot, for multidimensional plots, you may work with
                    ## GridVisualize.jl or write your own code.
                    Plotter.plot(subg[Coordinates]', Ecc ./ q .- solpsi, label = label1, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

                    Plotter.plot(subg[Coordinates]', - solcc, label = label2, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[2])

                    count = count + 1
                end
            end
        end

        Plotter.grid()
        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("energies [\$eV\$]")
        Plotter.legend(fancybox = true, loc = "best")
        return Plotter.title(title)
        #Plotter.ylim(-1.8, 1.7)

    end
    nothing

end

# ╔═╡ c5800a83-a57c-4eb0-9b6c-3c30266d0e3b
begin
    PyPlot.clf()
    plot_energies2(PyPlot, ctsys, printSol, "", label_energy)
    PyPlot.gcf()
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

        Ev_iNew = En[regionIntrinsic] - Eg * eV
        ctsys.fvmsys.physics.data.params.bandEdgeEnergy[iphip, regionIntrinsic] = Ev_iNew

        ## calculate equilibrium solution and as initial guess
        local solution = equilibrium_solve!(ctsys, control = control)
        local inival = solution

        #######################################################################
        ctsys.fvmsys.boundary_factors[iphia, bregionJunction2] = 1.0e30
        ctsys.fvmsys.boundary_values[iphia, bregionJunction2] = 0.0

        for istep in 1:(length(I) - 1)

            ## turn slowly generation on
            ctsys.data.λ2 = LAMBDA[istep + 1]
            solution = solve(ctsys, inival = inival, control = control)
            inival = solution

        end # generation loop

        ## put here back the homogeneous Neumann boundary conditions.
        ctsys.fvmsys.boundary_factors[iphia, bregionJunction2] = 0.0
        ctsys.fvmsys.boundary_values[iphia, bregionJunction2] = 0.0

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
            for ii in 1:(numberOfCarriers + 1)
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
    PyPlot.clf()

    PyPlot.subplot(311)
    PyPlot.plot(EgTest', EfficiencyVec, marker = "o")
    PyPlot.grid()
    PyPlot.xlabel("band gap (perovskite) [eV]")
    PyPlot.ylabel("efficiency \$ \\eta \$ [%]")
    PyPlot.tight_layout()

    PyPlot.subplot(312)
    PyPlot.plot(EgTest', OpenCircuitVec, marker = "o")
    PyPlot.grid()
    PyPlot.xlabel("band gap (perovskite) [eV]")
    PyPlot.ylabel("open circuit voltage [V]")

    PyPlot.subplot(313)
    PyPlot.plot(EgTest', FillfactorVec, marker = "o")
    PyPlot.grid()
    PyPlot.xlabel("band gap (perovskite) [eV]")
    PyPlot.ylabel("Fill factor [%]")
    PyPlot.gcf()


end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ChargeTransport = "25c3eafe-d88c-11e9-3031-f396758f002a"
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"

[compat]
ChargeTransport = "~0.2.4"
ExtendableGrids = "~0.9.17"
PlutoUI = "~0.7.50"
PyPlot = "~2.11.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "760fa544de395d94671bd1a702aeeb6c4ec561ef"

[[deps.ADTypes]]
git-tree-sha1 = "e6103228c92462a331003248fa31f00dcf41c577"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.1.1"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "a69dbe3b376ace7d9eebe2db43216e8b52ba6da9"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.29.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "6ef8fc1d77b60f41041d59ce61ef9eb41ed97a83"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.18"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChargeTransport]]
deps = ["DocStringExtensions", "ExtendableGrids", "ForwardDiff", "GridVisualize", "Interpolations", "Printf", "PyPlot", "Roots", "SparseArrays", "Test", "VoronoiFVM"]
git-tree-sha1 = "91da0a6f756a4e3f22ab27cb1726622fba4a6fe9"
uuid = "25c3eafe-d88c-11e9-3031-f396758f002a"
version = "0.2.4"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "e32a90da027ca45d84678b826fffd3110bb3fc90"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.8.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    IntervalSetsExt = "IntervalSets"
    StaticArraysExt = "StaticArrays"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "698124109da77b6914f64edd696be8dccf90229e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e1c40d78de68e9a2be565f0202693a158ec9ad85"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.11"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "2921bf0ffab4c8b7eda6a36c7b06a0dde6df0137"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.17"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "Sparspak", "SuiteSparse", "Test"]
git-tree-sha1 = "bbb16c582df45544612cd703fe1b8179339ccd2b"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "1.0.1"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "03fcb1c42ec905d15b305359603888ec3e65f886"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.19.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"
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
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "12165cfe9b04b67f0e349bf3bf370f08ea4cecbf"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "1.0.2"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "7c892c426f8d03a180366411566d0f6ac1790f6c"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "0.3.0"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "827f29c95676735719f8d6acbf0a3aaf73b3c9e5"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.3.2"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "d926e9c297ef4607866e8ef5df41cde1a642917f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.14"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "42c17b18ced77ff0be65957a591d34f4ed57c631"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.31"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "1d3e720d603557d697fedc036bd1af43fe7b3474"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.41.1"

    [deps.LinearSolve.extensions]
    LinearSolveHYPRE = "HYPRE"

    [deps.LinearSolve.weakdeps]
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

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
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "a282dbdbc2860134d6809acd951543ce359bcf15"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.155"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

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
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "0fe4e7c4d8ff4c70bfa507f0dd96fa161b115777"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "62f417f6ad727987c755549e9cd88c46578da562"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.95.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "92e7ca803b579b8b817f004e74b205a706d9a974"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "140cddd2c457e4ebb0cdc7c2fd14a7fbfbdf206e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.3"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "9088515ad915c99026beb5436d0a09cd8c18163e"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.18"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "b45deea4566988994ebb8fb80aa438a295995a6e"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.10"
weakdeps = ["ForwardDiff"]

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "f139e81a81e6c29c40f1971c9e5309b09c03f2c3"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SnoopPrecompile", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "392d3e28b05984496af37100ded94dc46fa6c8de"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.91.7"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "e61e48ef909375203092a6e83508c8416df55a83"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Reexport", "Requires", "SciMLOperators", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "VertexSafeGraphs"]
git-tree-sha1 = "aa5b879ce5fcd8adb0c069d93fa2567d9b68b448"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.0.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsZygote = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "08be5ee09a7632c32695d954a602df96a877bf0d"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.6"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "fd5f417fd7e103c121b0a0b4a6902f03991111f4"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.3.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f859ab67ca232b777a03a6cee588c1c15f7ec40a"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.9"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "bfbd444c209b41c7b2fef36b6e146a66da0be9f1"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.4"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "7ecd651e3829d2957478516e92f693f12d5b4781"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.2.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

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

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["BandedMatrices", "CommonSolve", "DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "33974fbb5978a37188e13f793a3e9740b9003f0d"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "1.2.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams"]
git-tree-sha1 = "49353f30da65f377cff0f934bb9f562a2c0441b9"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.17.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─39e1f60f-cd7a-49b5-b569-b3321f68c2ac
# ╠═b0138526-c79e-11ec-041a-156b0dfee367
# ╟─6511e625-2af2-44c9-bc5c-d24e08109c3f
# ╟─997e8130-8e2d-45f6-a6b9-5ed78782d2b0
# ╟─1284fff2-af76-4d53-9444-233bde7cfaa9
# ╟─c0f9f418-0c3b-428f-bd5b-e4d3d1ab9be2
# ╟─bb701b88-6a76-4e77-ad81-a776ebfa5ec4
# ╟─18823ee0-988d-459c-b340-7dfed6092b22
# ╟─febcc48f-82ae-408a-b6a8-74cf3c76f61a
# ╟─2d0bfda7-c23d-45b0-b3db-2000f1a291c3
# ╟─3685ff11-b2f1-4fc1-bde6-723bf59ca57f
# ╟─c72f9732-4fb3-485f-93e2-14999307d513
# ╟─c19b9329-1c9e-4ab3-8216-ff50dcb89e19
# ╟─2a43a7a8-b930-4fca-bf50-90220a3bb431
# ╟─5f0c398c-f389-4355-9f0e-6710c2b819a0
# ╟─004e4101-dbf3-4527-bfb8-47da01c70182
# ╟─66de390b-d0fd-48c4-93d4-50cb7bfc50cd
# ╟─b0df6373-8d11-4581-9e42-62e3aee6c869
# ╠═61155934-83a1-40ea-805e-2607fd8f9cd2
# ╟─380c5124-bd09-421f-9590-416400624374
# ╟─c0177ece-44e4-4431-abd6-4f28a8037d26
# ╟─afcc0f3b-47aa-4418-aead-67202fc56119
# ╟─110f405e-7a3a-4385-82ef-9a1c3886bc98
# ╟─719df739-d611-4855-870c-64dc82444538
# ╟─c305fbd9-79ff-4ea3-8011-b07246f767c6
# ╟─bb4e474b-2028-4982-b202-3b22aa08c9d1
# ╠═6f76b5cd-65a2-499c-a0be-0f3a0916e70c
# ╟─5573ea85-8bbb-411f-97b2-ef1e46981060
# ╠═29d920d6-ac56-4cb1-83a3-741ee6c876ae
# ╟─850bfa1a-4326-4aa3-97d2-c21d4ac2ba11
# ╟─2b785f7b-4a96-4785-993b-ee7bbf6b0533
# ╟─c00b33f1-8722-49e9-93b3-3703c5d0efb7
# ╠═ab8d4426-9eda-4bab-a68c-1475042321db
# ╠═d6e4a543-d2e5-42a2-91af-1cf8b4d4632d
# ╠═9fa7dc02-4913-4b6e-a96d-0d73ccfee302
# ╟─9b20814b-7f00-49da-ad83-63f1612d2f27
# ╟─4a4ce5d2-c248-4d88-b4de-ba42f244e0e5
# ╟─3f40a78e-4b2e-4863-8404-75157d562459
# ╟─c13bb60b-b07c-4c56-b4f5-1fccbb194bb4
# ╟─09f362b4-a226-4b4e-8253-a6e87d979777
# ╟─8dbd2c62-eae3-4857-8de7-8829554a847a
# ╟─01fdd92e-4033-4701-a5d4-7012c7c6c063
# ╟─c5800a83-a57c-4eb0-9b6c-3c30266d0e3b
# ╟─882f0ca5-e937-41c4-820b-23eae12e7c19
# ╟─26255f99-1a36-43bf-b6f6-01cfa1e1c396
# ╟─ce3423b6-8005-49da-b7c6-2e16e0675740
# ╟─557b00c2-6a4b-4071-ba90-99b275dcadd8
# ╟─f80b1946-73cb-448e-80a8-9e4770b50c79
# ╟─633ed076-9123-4989-b7e0-3ee078d1a7e0
# ╟─a6e14180-7eeb-48e7-afd2-8a147b32d870
# ╠═ccf30353-f057-47ad-8e69-e12ef4e01c00
# ╟─d03cb3b7-2d90-4ee3-842d-54d55e32db07
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
