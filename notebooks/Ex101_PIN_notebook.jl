### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ f34819f2-bb8b-4a1c-99ea-f0b34c4e4a5a
# Load Packages
begin
    using Pkg; Pkg.activate(".") # TODO remove this line before merging

    # for the numercis:
    using ChargeTransport                # Drift-diffusion solver
    using ExtendableGrids                # Grid manager
    using LessUnitful: @unitfactors      # Unit factors

    # for this notebook and plots:
    using PlutoUI                        # Fancy Notebook Buttons and Sliders
    using GridVisualize                  # Plotting Routines
    using CairoMakie                     # Plotting Backend
end

# ╔═╡ 79454036-46e1-47a9-952a-81773e17fb7a
TableOfContents()

# ╔═╡ ed10ce80-3786-49a9-a37b-987673d1b7c3
md"""
# GaAs PIN-Diode (1D Example)
We consider the example of a GaAs PIN-Diode with the following structure.
"""

# ╔═╡ c1597b55-f74b-46c9-9e49-f726ba7f6342
Resource(
    "https://upload.wikimedia.org/wikipedia/commons/a/ab/Pin-Diode.svg",
    :alt => "PIN diode",
    :width => 500
)

# ╔═╡ e36c2db4-911c-4190-8b13-b389ac506172
details(
    "Image source",
    md"""
    https://upload.wikimedia.org/wikipedia/commons/a/ab/Pin-Diode.svg
    """
)

# ╔═╡ 9658ee64-7467-4b3a-aa9a-deb1c0130921
md"""
# Introduction
"""

# ╔═╡ 1c22bea9-486d-4d63-a8bc-d40ce4a98c17
md"""
For an in-depth overview, see the handbook [Drift-Diffusion Models](https://www.taylorfrancis.com/chapters/edit/10.4324/9781315152318-25/drift-diffusion-models-patricio-farrell-nella-rotundo-duy-hai-doan-markus-kantner-j%C3%BCrgen-fuhrmann-thomas-koprucki) by Patricio Farrell, Nella Rotundo, Duy Hai Doan, Markus Kantner, Jürgen Fuhrmann, Thomas Koprucki

Preprint: [WIAS-Preprint](https://archive.wias-berlin.de/receive/wias_mods_00002331)
"""

# ╔═╡ e2cefaab-aa23-415b-b480-591c6d1715bc
md"""
## 1D Van Roosbroeck System
To introduce we shortly recab the stationary 1D van Roosbroeck system, a system to model the charge transport within 1D devices.

The stationary 1D van Roosbroeck system consists of three nonlinear ordinary differential equations for the unknown electrostatic potential $\psi(x)$, the quasi-Fermi potential for electrons $\varphi_n(x)$ and the quasi-Fermi potential for holes $\varphi_p(x)$.
It links Poisson’s equation for the electric field to the continuity equations for the carrier densities as follows:
```math
\begin{align}
	 -\frac{d}{dx}\left( \varepsilon_s \frac{d}{dx} \psi \right) &= q \Big( (p(\psi, \varphi_p) - C_p) - (n(\psi, \varphi_n) - C_n) \Big), \quad & (1)\\
	 \frac{d}{dx} j_n &= q R(\psi, \varphi_n, \varphi_p), \quad & (2)\\
	 \frac{d}{dx} j_p &= - q R(\psi, \varphi_n, \varphi_p). \quad & (3)
\end{align}
```

The fluxe are given by
```math
\begin{align}
	 j_n = - q \mu_n n(\psi, \varphi_n) \frac{d}{dx} \varphi_n 
	 \quad \text{and} \quad
	 j_p = - q \mu_p p(\psi, \varphi_p) \frac{d}{dx} \varphi_p.
\end{align}
```

The electron density $n$ and the hole density $p$ are related to the electric potential $\psi$ as well as the quasi-Fermi potentials of electrons and holes via
```math
\begin{align}
	 n(\psi, \varphi_n) = N_c \mathcal F \left( \frac{q(\psi - \varphi_n) - E_c}{k_B T} \right)
	 \quad \text{and} \quad
	 p(\psi, \varphi_p) = N_v \mathcal F \left( \frac{q(\varphi_p - \psi) + E_v}{k_B T} \right).
\end{align}
```
"""

# ╔═╡ 42878cd0-dd38-42bc-97c6-f356707eb2f9
md"""
## Model Choices for PIN Diode

**Statistics**\
In this example we choose the Bolzmann Approximation as statistical distribution:
```math
\mathcal{F} (\eta) \approx \exp(\eta), \quad \eta \leq -2.
```
In general more choices are possible:

- Fermi–Dirac integral of order $\frac{1}{2}$ with corresponding Blakemore or Bolzmann approximation
- Gauss–Fermi integral with corresponding Blakemore or Bolzmann approximation.
"""

# ╔═╡ 88f2372b-8000-4191-9612-00b57ddd59c4
md"""
View statistic: 
$(@bind statistic Select(["Fermi-Dirac integral of order 1/2", "Gauss-Fermi integral"]))
"""

# ╔═╡ 3cbc7c62-eb28-4b14-b265-aada673b2114
if statistic == "Fermi-Dirac integral of order 1/2"
    md"""
    ```math
    \mathcal{F}(\eta) := \mathcal{F}_{1/2}(\eta) = \frac{2}{\sqrt{\pi}} \int_{0}^{\infty} 
    \frac{\xi^{1/2}}{\exp(\xi - \eta) + 1} \, d\xi
    ```

    This statistic is suitable for three-dimensional bulk semiconductors.
    """
else
    md"""
    ```math
    \mathcal{F}(\eta) := G(\eta; \sigma)
    		= \frac{1}{\sigma} \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \exp\!\left( -\frac{\xi^{2}}{2\sigma^{2}} \right) \frac{1}{\exp(\xi - \eta) + 1} \, d\xi
    ```

    This statistic is suitable for organic semiconductors.
    """
end

# ╔═╡ e2b97c6a-827c-40d0-a15e-b75bc53c054a
md"""
View corresponding approximation: 
$(@bind approx Select(["Bolzmann", "Blakemore"]))
"""

# ╔═╡ 1e32f86e-e9f8-474d-b126-e27446d11a0f
if statistic == "Fermi-Dirac integral of order 1/2"
    if approx == "Blakemore"
        md"""
        ```math
        F(\eta) \approx \frac{1}{\exp(-\eta) + 0.27}, \qquad \eta \le 1.3
        ```
        """
    else
        md"""
        ```math
        \mathcal{F}(\eta) \approx \exp(\eta), \qquad \eta \le -2
        ```
        """
    end
else
    if approx == "Blakemore"
        md"""
        ```math
        F(\eta) \approx \frac{1}{\exp(-\eta) + 1}, \qquad \sigma \to 0
        ```
        """
    else
        md"""
        ```math
        F(\eta) \approx \exp\!\left(\frac{\sigma^{2}}{2}\right)\exp(\eta),
        \qquad \eta \le -\sigma^{2} - 4
        ```
        """
    end
end

# ╔═╡ 35c77f4b-7734-4d36-bc03-8ae3fddc9d36
md"""
**Recombination**\
The function $R(\psi, \varphi_n, \varphi_p)$ models the recombination and generation of electron-hole pairs. 
The following recombination processes are taken into account:
"""

# ╔═╡ dcad13fc-485c-4079-b66e-b421e1d543aa
Resource(
    "https://lab.wias-berlin.de/osewold/arise-chargetransport-jl-tutorial/-/raw/main/Recombination-processes.svg",
)

# ╔═╡ 4192ddae-7a73-4eb3-b9e1-3ccdb01552f7
md"""
**Doping**\
The doping profiles ``C_p = C_p(x)`` and ``C_n = C_n(x)`` describe the spatial distribution of ionized dopants. It is defined by the local concentrations of donors $N_D$ and acceptors $N_A$:
```math
\begin{align}
	\text{p-doping:}\quad C_p(x) \coloneqq N_A (x)\\
	\text{n-doping:}\quad C_n(x) \coloneqq N_D (x).
\end{align}
```
Signs are handled internally by ChargeTransport.jl and do not need to be taken into account when entering data.

**Boundary model**\
We choose OhmicContacts for both boundaries. For contacts currently the following are possible:

- OhmicContact
- SchottkyContact. 

For inner boundaries we have InterfaceNone, InterfaceRecombination.
"""

# ╔═╡ ca761626-09f1-4168-9a05-6a3890e3641e
md"""
## Discretization of the System

**Flux Discretization**\
We choose the Sedan Scheme. In general there are more choices possible:

- ScharfetterGummel, ScharfetterGummelGraded
- ExcessChemicalPotential, ExcessChemicalPotentialGraded (Sedan Scheme)
- DiffusionEnhanced
- GeneralizedSG.

**Finite Volume Method**\
We do a finite volume discretization of the Poisson equation (1) and the two continuity equations (2) and (3).
"""

# ╔═╡ 17087764-3a10-47c3-a5a0-40c91a9c281c
md"""
# Numerical Simulation
"""

# ╔═╡ 5fab15af-ea68-4456-8c8e-1db64e250176
md"""
👉 `CairoMakie` is a very stable backend but it takes __very__ long to compile on the first run.
You can also try `use PythonPlot` and change `Plotter = PythonPlot` in the next cell.
"""

# ╔═╡ e90bad28-9fe6-407f-b4f4-f248ea4609b7
Plotter = CairoMakie;

# ╔═╡ e2220ac5-ad1b-4016-96da-428d1b3542dc
md"""
### Factors and Constants
We globally define our factors and constants.
"""

# ╔═╡ d63cb9f8-e208-40fd-a723-1b4c83881477
begin
    @unitfactors μm cm s ns V K

    (; q, k_B, ε_0) = ChargeTransport.constants

    eV = q * V
end;

# ╔═╡ 66d27904-3bde-49df-9ab4-b871c3d9e191
md"""
### PIN-Parameters
The PIN-Paremeters are characterizing the PIN-Diode.
"""

# ╔═╡ 0bce7b85-f45a-4dd0-a9aa-2ec7ad13979f
@kwdef struct PIN_Parameter

    ########## charge carriers ##########
    iphin = 1                            # quasi Fermi potential for electrons
    iphip = 2                            # quasi Fermi potential for holes
    numberOfCarriers = 2

    ########## device geometry ##########
    # region numbers
    regionAcceptor = 1           # p doped region
    regionIntrinsic = 2          # intrinsic region
    regionDonor = 3              # n doped region
    regions = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionAcceptor = 1         # outer left boundary
    bregionDonor = 2            # outer right boundary
    bregionJunction1 = 3        # first inner interface
    bregionJunction2 = 4        # second inner interface

    h_pdoping = 2.0 * μm        # length p-doped layer
    h_intrinsic = 2.0 * μm      # length intrinsic layer (adjustable for this 																notebook, default = 2.0)
    h_ndoping = 2.0 * μm        # length n-doped layer
    w_device = 0.5 * μm         # width of device
    z_device = 1.0e-4 * cm      # depth of device

    ########## physical values ##########
    # band edge energies
    Ec = 1.424 * eV                      # conduction band-edge energy
    Ev = 0.0 * eV                        # valence band-edge energy
    Nc = 4.35195989587969e17 / (cm^3)    # conduction band density of states
    Nv = 9.139615903601645e18 / (cm^3)   # valence band densitiy of states

    # mobilities
    mun = 8500.0 * (cm^2) / (V * s)      # electron mobility
    mup = 400.0 * (cm^2) / (V * s)       # hole mobility
    εr = 12.9 * 1.0                      # relative dielectric permittivity of GaAs
    T = 300.0 * K                        # room temperature

    Ni = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * k_B * T)) # intrinsic concentration

    # Recombination parameters
    Auger = 1.0e-29 * cm^6 / s
    SRH_TrapDensity = 1.0e10 / cm^3
    SRH_LifeTime = 1.0 * ns
    Radiative = 1.0e-10 * cm^3 / s

    # doping
    dopingFactorNd = 1.0       # Donor-factor
    dopingFactorNa = 0.46      # Acceptor-factor
    Nd = dopingFactorNd * Nc   # Donor-concentration
    Na = dopingFactorNa * Nv   # Acceptor-concentration

    # contact voltage: We impose an applied voltage only on one boundary. At the 		other boundary the applied voltage is zero.
    voltageAcceptor = 1.5 * V
end;

# ╔═╡ 0373265a-902b-4d6f-a007-10ce8187fded
md"""
### Grid
We use a one-dimensional uniform grid.
"""

# ╔═╡ 63e484be-5bec-48e4-be8a-60344b3f9592
md"""
### Store PIN-Parameters in Params

`Pin_Parameter` was a user-defined struct. Now we need a strict format for `ChargeTransport`. 
This `Params` object contains all necessary physical parameters.
"""

# ╔═╡ 6594723e-a59e-45ab-9a4a-7773eabb8517
function make_params(p::PIN_Parameter, grid)

    (; iphin, iphip) = p

    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], p.numberOfCarriers)

    params.temperature = p.T
    params.chargeNumbers[iphin] = -1
    params.chargeNumbers[iphip] = 1

    for ireg in 1:grid[NumCellRegions] # region data

        params.dielectricConstant[ireg] = p.εr * ε_0

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg] = p.Nc
        params.densityOfStates[iphip, ireg] = p.Nv
        params.bandEdgeEnergy[iphin, ireg] = p.Ec
        params.bandEdgeEnergy[iphip, ireg] = p.Ev
        params.mobility[iphin, ireg] = p.mun
        params.mobility[iphip, ireg] = p.mup

        # recombination parameters
        params.recombinationRadiative[ireg] = p.Radiative
        params.recombinationSRHLifetime[iphin, ireg] = p.SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg] = p.SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = p.SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = p.SRH_TrapDensity
        params.recombinationAuger[iphin, ireg] = p.Auger
        params.recombinationAuger[iphip, ireg] = p.Auger

    end

    # doping                                                        #   n    p
    params.doping[iphin, p.regionDonor] = p.Nd                      # [Nd   0.0;
    params.doping[[iphin, iphip], p.regionIntrinsic] .= p.Ni        #  ni   ni;
    params.doping[iphip, p.regionAcceptor] = p.Na                   #  0.0  Na]

    return params
end;

# ╔═╡ 6a8b7115-6e2b-4869-8e22-38d2f4af556b
md"""
### Store Model Information in Data

Data contains all model choices, including the physical parameters defined in `Params`, and the `grid`.
"""

# ╔═╡ e9996891-cb52-493e-be0a-5e3e06c402fa
function make_data(p::PIN_Parameter, grid)

    data = Data(grid, p.numberOfCarriers)

    # stationary or transient problem
    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    # boundary model
    data.boundaryType[p.bregionAcceptor] = OhmicContact
    data.boundaryType[p.bregionDonor] = OhmicContact

    # flux discretization
    data.fluxApproximation .= ExcessChemicalPotential           # Sedan scheme

    data.params = make_params(p, grid)

    return data
end;

# ╔═╡ 3c7bfd68-29b4-4ecf-bd7d-f6ae0e069abd
md"""
### ChargeTransport System

The `System` implements a `VoronoiFVM` system. This is passed to the underlying solver.

We initialize our system with previously defined `data` which is likewise dependent on all parameters.
"""

# ╔═╡ 7e416f00-17ce-4c60-b11b-99e847284cb7
md"""
### Model Plots
"""

# ╔═╡ 822a4096-c5eb-4a76-962c-a28717aae9a6
md"""
**Band-edge energies**\
Here the conduction and valence band-edge energies $E_{\rm C}$ and $E_{\rm V}$ are shown. It is assumed that the band-edge energies are constant along the diode material.
"""

# ╔═╡ c5af3d6b-5840-4445-9851-cf89f9a8db9b
md"""
**Doping**\
Looking from the left boundary, ``2`` µm are p-doped with 
```math
	N_a = 4.204223315656756 \ 10^{18} \text{cm}^{-3}.
```

Looking from the right boundary ``2`` μm are n-doped with
```math
	N_d = 4.35195989587969 \ 10^{17} \text{cm}^{-3}.
```

The intrinsic layer is slightly n- and p-doped with 
```math
	N_i = 2.1813876738858657 \ 10^6 \text{cm}^{-3}.
```
"""

# ╔═╡ 03bd6c73-786b-4dc8-965d-ea3d763dba68
md"""
**Electroneutral potential**\
Local charge neutrality is characterized by a vanishing left-hand side in the Poisson equation (1) and thus ``\psi_0`` is the solution of

```math
0 = N_v \exp\left(\frac{E_v - q \psi_0}{k_B T} \right) - N_c \exp\left(\frac{q \psi_0 - E_c}{k_B T} \right) + C(x).
```
"""

# ╔═╡ b04e92a2-50dc-4f02-9fbf-23a7125338b2
md"""
### Control Parameters for Solver
These parameters are used during the Newton Step within the VoronoiFVM.jl package.

Overview options: [https://wias-pdelib.github.io/VoronoiFVM.jl/stable/solver/#Solver-control](https://wias-pdelib.github.io/VoronoiFVM.jl/stable/solver/#Solver-control)
"""

# ╔═╡ f25d998d-007e-4631-8734-27845c16d3d4
begin
    function make_solver_control()
        control = SolverControl()
        control.verbose = true # or false
        control.maxiters = 50
        control.abstol = 1.0e-14
        control.reltol = 1.0e-14
        control.tol_round = 1.0e-8
        control.damp_initial = 0.5
        control.max_round = 3
        return control
    end

    solver_control = make_solver_control()
end

# ╔═╡ 0d0a59e7-6f51-4889-9659-c76ccd1328b1
md"""
### Solution in Thermodynamic Equilibrium
The thermodynamic equilibrium is defined by vanishing currents $j_n$ and $j_p$. Thus the three differential equations (1), (2) and (3) reduce to a 1D nonlinear Possion equation in ``\Omega = [0,L]``:\
\

```math
- \frac{d}{dx} \left( \varepsilon_s \frac{d}{dx} \psi \right) = q \left(C + N_v \exp\left(\frac{E_v - q\psi}{k_B T}\right) - N_c \exp\left(\frac{q\psi - E_c}{k_B T}\right) \right)
```
\
with boundary conditions
```math
\begin{align}
	\psi(0) = \psi_0(0) \quad \text{and} \quad \psi(L) = \psi_0(L).
\end{align}
```

"""

# ╔═╡ f035cd86-16c4-4757-bfe1-f356923cf6e1
md"""
### Solution Plots
""" |> WideCell

# ╔═╡ 31b6a3eb-fa52-4d34-b442-267e246bbd39
WideCell(
    md"""
    $(@bind reset_to_defaults CounterButton("Reset to defaults"))
    """
)

# ╔═╡ 2d918065-bd3a-4177-8f43-8aa41f3b0d62
WideCell(
    begin
        reset_to_defaults
        @bind values PlutoUI.combine() do Child
            md"""
            Here you have the option to adjust parameters.

            *Keep in mind that this is a real numerical simulation, so some parameter combinations may lead to a crash. As with every other project, there may be hidden instabilities in the solver.* 😉 
            			
            | | | | |
            |:---|:---|:---|:---|
            | **Geometry parameter:** | | | | 
            |intrinsic layer [μm]: | $(Child("h_intrinsic", PlutoUI.Slider(0.1:0.1:5.0, default = 2, show_value = true))) | | |
            | | | | |
            | | | | |
            | | | | |
            | **Physical parameters:** | | | |
            |voltage acceptor [V]:| $(Child( "voltageAcceptor", PlutoUI.Slider(-3.5:0.1:5.0, default = 1.5, show_value = true))) | n-doping factor: | $(Child( "dopingFactorNd", PlutoUI.Slider(0.4:0.1:1.2, default = 1.0, show_value = true))) |
            |conduction band-edge energy [eV]:| $(Child( "Ec", PlutoUI.Slider(1.0:0.01:2.0, default = 1.424, show_value = true))) | electron mobility [cm^2/(V*s)]: | $(Child( "mun", PlutoUI.Slider(5000:100:12000, default = 8500, show_value = true))) |
            | | | | |
            | | | | |
            | | | | |
            | **Simulation parameters:** | | | |
            |grid refinement: | $(Child( "number_of_refinements", PlutoUI.Slider(1:5, default = 2, show_value = true))) | bias steps: | $(Child( "bias_steps", PlutoUI.Slider(11:100, default = 32, show_value = true))) |

            $(@bind go CounterButton("Start Simulation 😃"))
            """
        end
    end
)

# ╔═╡ df962a5e-f5ed-4c9a-b2c5-63f580545243
WideCell(
    md"""
    **Energies and carrier densities**
    """
)

# ╔═╡ a05e656f-2e2f-40f1-993d-f6515e95cffc
md"""
**Solution**

Recall:
- ``\psi(x)`` electrostatic potential
- ``\varphi_n(x)`` quasi-Fermi potential for electrons
- ``\varphi_p(x)`` quasi-Fermi potential for holes
"""

# ╔═╡ e03b802a-00de-4e37-b170-4d1c8ba420a0
md"""
**Current-Voltage curve**
"""

# ╔═╡ a9cc27c3-6f41-4dd4-974b-9bdf6e2346b7
@kwdef mutable struct SliderBoxData
    h_intrinsic = 2.0
    voltageAcceptor = 1.5
    dopingFactorNd = 1.0
    Ec = 1.424
    mun = 8500.0
    number_of_refinements = 3
    bias_steps = 32
end;

# ╔═╡ 5c4fcc36-a3e1-40a8-aa0b-2f14eeb084fe
slider_box_data_dummy = SliderBoxData();

# ╔═╡ c5a26175-b0a2-498d-bd4a-9e71a3e6ec29
for prop in propertynames(values)
    setproperty!(slider_box_data_dummy, prop, getproperty(values, prop))
end;

# ╔═╡ 1602d7ce-86be-40dd-8128-a014d876ecb0
begin
    go
    slider_box_data = slider_box_data_dummy
end;

# ╔═╡ 89dc012a-206b-4e5f-bef1-ab476168f028
function make_grid(p::PIN_Parameter)
    # This function is used to initialize the grid for a possible extension to other p-i-n devices.
    function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
        coord_ndoping = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
        coord_intrinsic = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
        coord_pdoping = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
        coord = glue(coord_ndoping, coord_intrinsic)
        coord = glue(coord, coord_pdoping)
        return coord
    end

    # grid (adjustable for this notebook, default = 3)
    refinementfactor = 2^(slider_box_data.number_of_refinements - 1)

    # extract some parameters from PIN_Parameters
    (; h_pdoping, h_intrinsic, h_ndoping) = p

    h_total = h_pdoping + h_intrinsic + h_ndoping   # length of device

    coord = initialize_pin_grid(
        refinementfactor,
        h_pdoping,
        h_intrinsic,
        h_ndoping
    )

    grid = simplexgrid(coord)

    # cellmask! to define subregions and assigning region number
    cellmask!(grid, [0.0 * μm], [h_pdoping], p.regionAcceptor)
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], p.regionIntrinsic)
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], p.regionDonor)

    # bfacemask! to set different boundary regions
    bfacemask!(grid, [0.0], [0.0], p.bregionAcceptor)
    bfacemask!(grid, [h_total], [h_total], p.bregionDonor)
    bfacemask!(grid, [h_pdoping], [h_pdoping], p.bregionJunction1)
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], p.bregionJunction2)

    return grid
end;

# ╔═╡ 9c0b686b-9973-4599-bead-6b9310bd56db
# tmp_string for double interpolation (LaTeX and Julia)
let
    tmp_string = """
    ### Bias Loop
    We increase the bias step by step until we reach the maximal voltage. The voltage is applied to the acceptor boundary, at the p-doped region with a maximal voltage of

    ```math
    	U_{\\rm end} = $(slider_box_data.voltageAcceptor) \\text{ V.}
    ```
    """
    Markdown.parse(tmp_string)
end

# ╔═╡ 9dc4cae5-b0bd-4baa-993f-1aa18ee42148
md"""
# Appendix to modify PIN parameters
"""

# ╔═╡ c9492318-eecd-4941-a25e-93bb43cf79b7
pin_parameter = PIN_Parameter(
    h_intrinsic = slider_box_data.h_intrinsic * μm,
    voltageAcceptor = slider_box_data.voltageAcceptor,
    dopingFactorNd = slider_box_data.dopingFactorNd,
    Ec = slider_box_data.Ec * eV,
    mun = slider_box_data.mun * (cm^2) / (V * s)
);

# ╔═╡ 0ae253cb-2e0d-4113-b835-ec0dc6c1fc91
grid = make_grid(pin_parameter)

# ╔═╡ e091f990-e99c-4ec3-bed9-44649eda06b9
let
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    GridVisualize.gridplot!(vis[1, 1], grid, Plotter = Plotter, legend = :none, title = "Grid", xlabel = "space [m]", size = (1300, 500))

    reveal(vis)
end

# ╔═╡ 13d8795c-b937-4c62-9064-cf6de3b9ba21
begin
    # Charge Transport System
    data = make_data(pin_parameter, grid)
    chargetransport_system = System(grid, data, unknown_storage = :sparse) # store sparse data types whenever possible

    # Printout Parameters for double-checking 😄
    show_params(chargetransport_system)
end

# ╔═╡ 1a2f9cca-3faf-4397-b2ff-e67c389d6b36
let
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    plot_energies!(vis[1, 1], chargetransport_system, label_BEE)
    reveal(vis)
end

# ╔═╡ 272dd019-6192-4935-8c79-7cf938b2babc
let
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    plot_doping!(vis[1, 1], chargetransport_system, label_density)
    reveal(vis)
end

# ╔═╡ 25e956c5-be10-4298-8ba0-61588cf73c13
begin
    psi0 = electroNeutralSolution(chargetransport_system)

    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    plot_electroNeutralSolutionBoltzmann!(vis[1, 1], grid, psi0; plotGridpoints = true)
    reveal(vis)
end

# ╔═╡ a5f82ebf-9b5d-482e-8c14-7d0da3f86f2d
solution_eq = equilibrium_solve!(chargetransport_system; control = solver_control)

# ╔═╡ 499b123d-dd70-4014-9210-a31375cc4b1c
begin
    maxBias = pin_parameter.voltageAcceptor
    biasValues = range(0, stop = maxBias, length = slider_box_data.bias_steps)
    IV = zeros(0)

    inival = copy(solution_eq)
    solution = copy(solution_eq)

    for Δu in biasValues

        println("bias value: Δu = ", Δu, " V")

        # set non equilibrium boundary conditions
        set_contact!(chargetransport_system, pin_parameter.bregionAcceptor, Δu = Δu)

        solution .= solve(chargetransport_system; inival, control = solver_control)
        inival .= solution

        # get I-V data
        current = get_current_val(chargetransport_system, solution)

        push!(IV, abs.(pin_parameter.w_device * pin_parameter.z_device * (current)))

    end
end

# ╔═╡ 46b13696-58e5-4f1c-b85f-bef15c499127
WideCell(
    let
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
        vis = GridVisualizer(; Plotter = Plotter, layout = (1, 2), size = (1500, 500))

        plot_energies!(vis[1, 1], chargetransport_system, solution, "Applied voltage Δu = $(biasValues[end])", label_energy, plotGridpoints = true)

        plot_densities!(vis[1, 2], chargetransport_system, solution, "Applied voltage Δu = $(biasValues[end])", label_density, plotGridpoints = true)

        reveal(vis)
    end
)

# ╔═╡ 645762b3-71ee-45ad-bc90-5423423a94b6
let
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    plot_solution!(vis[1, 1], chargetransport_system, solution, "Applied voltage Δu = $(biasValues[end])", label_solution, plotGridpoints = true)
    reveal(vis)
end

# ╔═╡ 3e981a6f-160a-471c-9bdf-2abd856094ec
begin
    # Print Solution
    φn = solution[1, :]
    φp = solution[2, :]
    ψ = solution[3, :]

    println(
        "\nQuasi Fermi potential electrons: \nφn = ", φn,

        "\n\nQuasi Fermi potential holes: \nφp = ", φp,

        "\n\nElectric potential: \nψ = ", ψ
    )
end

# ╔═╡ 84868da5-bf7b-4273-82f3-882a1da40152
let
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (900, 500))

    plot_IV!(vis[1, 1], biasValues, IV, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
    reveal(vis)
end

# ╔═╡ Cell order:
# ╟─79454036-46e1-47a9-952a-81773e17fb7a
# ╟─ed10ce80-3786-49a9-a37b-987673d1b7c3
# ╟─c1597b55-f74b-46c9-9e49-f726ba7f6342
# ╟─e36c2db4-911c-4190-8b13-b389ac506172
# ╟─9658ee64-7467-4b3a-aa9a-deb1c0130921
# ╟─1c22bea9-486d-4d63-a8bc-d40ce4a98c17
# ╟─e2cefaab-aa23-415b-b480-591c6d1715bc
# ╟─42878cd0-dd38-42bc-97c6-f356707eb2f9
# ╟─88f2372b-8000-4191-9612-00b57ddd59c4
# ╟─3cbc7c62-eb28-4b14-b265-aada673b2114
# ╟─e2b97c6a-827c-40d0-a15e-b75bc53c054a
# ╟─1e32f86e-e9f8-474d-b126-e27446d11a0f
# ╟─35c77f4b-7734-4d36-bc03-8ae3fddc9d36
# ╟─dcad13fc-485c-4079-b66e-b421e1d543aa
# ╟─4192ddae-7a73-4eb3-b9e1-3ccdb01552f7
# ╟─ca761626-09f1-4168-9a05-6a3890e3641e
# ╟─17087764-3a10-47c3-a5a0-40c91a9c281c
# ╠═f34819f2-bb8b-4a1c-99ea-f0b34c4e4a5a
# ╟─5fab15af-ea68-4456-8c8e-1db64e250176
# ╠═e90bad28-9fe6-407f-b4f4-f248ea4609b7
# ╟─e2220ac5-ad1b-4016-96da-428d1b3542dc
# ╠═d63cb9f8-e208-40fd-a723-1b4c83881477
# ╟─66d27904-3bde-49df-9ab4-b871c3d9e191
# ╠═0bce7b85-f45a-4dd0-a9aa-2ec7ad13979f
# ╟─0373265a-902b-4d6f-a007-10ce8187fded
# ╠═89dc012a-206b-4e5f-bef1-ab476168f028
# ╠═0ae253cb-2e0d-4113-b835-ec0dc6c1fc91
# ╟─e091f990-e99c-4ec3-bed9-44649eda06b9
# ╟─63e484be-5bec-48e4-be8a-60344b3f9592
# ╠═6594723e-a59e-45ab-9a4a-7773eabb8517
# ╟─6a8b7115-6e2b-4869-8e22-38d2f4af556b
# ╠═e9996891-cb52-493e-be0a-5e3e06c402fa
# ╟─3c7bfd68-29b4-4ecf-bd7d-f6ae0e069abd
# ╠═13d8795c-b937-4c62-9064-cf6de3b9ba21
# ╟─7e416f00-17ce-4c60-b11b-99e847284cb7
# ╟─822a4096-c5eb-4a76-962c-a28717aae9a6
# ╟─1a2f9cca-3faf-4397-b2ff-e67c389d6b36
# ╟─c5af3d6b-5840-4445-9851-cf89f9a8db9b
# ╟─272dd019-6192-4935-8c79-7cf938b2babc
# ╟─03bd6c73-786b-4dc8-965d-ea3d763dba68
# ╟─25e956c5-be10-4298-8ba0-61588cf73c13
# ╟─b04e92a2-50dc-4f02-9fbf-23a7125338b2
# ╠═f25d998d-007e-4631-8734-27845c16d3d4
# ╟─0d0a59e7-6f51-4889-9659-c76ccd1328b1
# ╠═a5f82ebf-9b5d-482e-8c14-7d0da3f86f2d
# ╟─9c0b686b-9973-4599-bead-6b9310bd56db
# ╠═499b123d-dd70-4014-9210-a31375cc4b1c
# ╟─f035cd86-16c4-4757-bfe1-f356923cf6e1
# ╟─2d918065-bd3a-4177-8f43-8aa41f3b0d62
# ╟─31b6a3eb-fa52-4d34-b442-267e246bbd39
# ╟─df962a5e-f5ed-4c9a-b2c5-63f580545243
# ╟─46b13696-58e5-4f1c-b85f-bef15c499127
# ╟─a05e656f-2e2f-40f1-993d-f6515e95cffc
# ╟─645762b3-71ee-45ad-bc90-5423423a94b6
# ╟─3e981a6f-160a-471c-9bdf-2abd856094ec
# ╟─e03b802a-00de-4e37-b170-4d1c8ba420a0
# ╟─84868da5-bf7b-4273-82f3-882a1da40152
# ╟─a9cc27c3-6f41-4dd4-974b-9bdf6e2346b7
# ╟─5c4fcc36-a3e1-40a8-aa0b-2f14eeb084fe
# ╟─c5a26175-b0a2-498d-bd4a-9e71a3e6ec29
# ╟─1602d7ce-86be-40dd-8128-a014d876ecb0
# ╟─9dc4cae5-b0bd-4baa-993f-1aa18ee42148
# ╠═c9492318-eecd-4941-a25e-93bb43cf79b7
