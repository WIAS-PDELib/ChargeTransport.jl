Perovskite solar cell
================================
We simulate charge transport in perovskite solar cells (PSCs), where we have apart from holes and electrons also ionic charge carriers. Here, we assume to have three domains, denoted by
$\mathbf{\Omega} = \mathbf{\Omega}_{\text{HTL}} \cup \mathbf{\Omega}_{\text{PVK}} \cup \mathbf{\Omega}_{\text{ETL}}  $.
The unknowns are the quasi Fermi potentials of electrons, holes and anion vacancies
$\varphi_n, \varphi_p, \varphi_a$
as well as the electric potential
$\psi$.
The underlying PDEs are given by
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (n_\text{p}(\psi, \varphi_\text{p}) - C_\text{p} ) - (n_\text{n}(\psi, \varphi_\text{n}) - C_\text{n}) \Big),\\
	q \partial_t n_\text{n}(\psi, \varphi_\text{n}) - \nabla \cdot \mathbf{j}_\text{n} &= q\Bigl(G(\mathbf{x}) - R(n_\text{n},n_\text{p}) \Bigr), \\
	q \partial_t n_\text{p}(\psi, \varphi_\text{p}) + \nabla \cdot \mathbf{j}_\text{p} &= \Bigl(G(\mathbf{x}) - R(n_\text{n},n_\text{p}) \Bigr),
\end{aligned}
```
for
$\mathbf{x} \in \mathbf{\Omega}_{\text{HTL}} \cup  \mathbf{\Omega}_{\text{ETL}} $, $t \in [0, t_F]$. In the middle, intrinsic region ($ \mathbf{x} \in \mathbf{\Omega}_{\text{PVK}} $), we have
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( n_\text{p}(\psi, \varphi_\text{p})  - n_\text{n}(\psi, \varphi_\text{n}) + n_\text{a}(\psi, \varphi_\text{a}) - C_\text{a} \Big),\\
q \partial_t n_\text{n}(\psi, \varphi_\text{n})	- \nabla \cdot \mathbf{j}_\text{n} &= \Bigl(G(\mathbf{x}) - R(n_\text{n},n_\text{p}) \Bigr), \\
	q \partial_t n_\text{p}(\psi, \varphi_\text{p}) + \nabla \cdot \mathbf{j}_\text{p} &= \Bigl(G(\mathbf{x}) - R(n_\text{n},n_\text{p}) \Bigr),\\
	q \partial_t n_\text{a}(\psi, \varphi_\text{a}) + \nabla \cdot \mathbf{j}_\text{a} &= 0,
\end{aligned}
```
see [Abdel2021](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865).

Differences to the previous example include
- an additional charge carrier (the anion vacancy)
- parameter jumps across heterojunctions
- the transient case
- a generation rate $G$
- higher dimensional problem.

# Simulating electronic-ionic charge transport

For example, for perovskite solar cells an average vacancy density is given in literature.
If the charge carrier densities are used as the set of unknowns, the initial condition for the anion vacancy density can be given by
```math
\begin{aligned}
	n_{\text{a}}^0 = C_{\text{a}},
\end{aligned}
```
where $C_{\text{a}}$ corresponds to the uniform density of cation vacancies and is set equal to the average anion vacancy density to ensure global charge neutrality.
With homogeneous no-flow Neumann boundary conditions around the perovskite layer for the anion vacancies, the total mass of anions is conserved at all times, i.e.
```math
\begin{aligned}
	\frac{1}{|\Omega_{\text{PVK}} |} \int_{\Omega_{\text{PVK}}} n_{\text{a}} (\mathbf{x}, t)\, d\mathbf{x}
    = \frac{1}{|\Omega_{\text{PVK}} |}  \int_{\Omega_{\text{PVK}}} n_{\text{a}}^0 (\mathbf{x})\, d\mathbf{x}
    = C_{\text{a}}, \quad \text{for all} \quad t \geq 0.
\end{aligned}
```
Since we do not use charge carrier densities directly but instead employ quasi Fermi potentials as the unknowns, due to mathematical, physical, and numerical advantages, we must find a workaround to properly fix the initial condition for the vacancy density.

## Quasi Fermi potential notation

The statistical relation between the vacancy density and the potentials (our chosen unknowns) reads
```math
\begin{aligned}
    n_{\text{a}}  = N_{\text{a}}  F_{-1} \Bigl(\eta_{\text{a}} (\psi, \varphi_{\text{a}} ) \Bigr), \quad \eta_{\text{a}}  = z_{\text{a}}  \frac{q (\varphi_{\text{a}}  - \psi) + E_{\text{a}} }{k_B T},
\end{aligned}
```
where $N_{\text{a}}$ denotes the maximum vacancy density allowed, $F_{-1} = (\exp(-x) +1)^{-1}$ is the Fermi-Dirac integral of order $-1$, and we refer to $E_{\text{a}}$ is the intrinsic vacancy energy level (somehow a model parameter).
In equilibrium, we set the vacancy quasi-Fermi potential to zero, i.e., $\varphi_{\text{a}} = 0$, when the applied voltage is $V = 0$.
Because the initial condition for the vacancy density is prescribed as $C_\text{a}$ and the electric potential $\psi$ is an unknown, the only remaining free parameter is the vacancy energy $E_{\text{a}}$.
> The value of $E_{\text{a}}$ must be chosen such that the average vacancy density remains conserved for all time steps.

## Current way of handling $E_\text{a}$

We implement a method that calculates the appropriate `Ea` values internally via the secant method.
To make use of this feature, you can add in the equilibrium solving the flag `vacancyEnergyCalculation=true`.

```julia
solution = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = true)
```
To check, if, indeed, the average vacancy density is maintained, you can calculate that value and print the chosen vacancy energy level.
```julia
integral = integrated_density(ctsys, sol = solution, icc = iphia, ireg = regionIntrinsic)

println("Calculated average vacancy density is: ", integral / data.regionVolumes[regionIntrinsic])

vacancyEnergy = data.params.bandEdgeEnergy[iphia, regionIntrinsic] / data.constants.q
println("Value for vacancy energy is: ", vacancyEnergy, " eV")
```
## Remarks

- For **1D simulations**, this approach is sufficient.
- For **multi-dimensional simulations**, however, we recommend precomputing the `Ea` values and storing them in case of multiple computations with the same parameter set.


Next, we give a quick survey on how to use `ChargeTransport.jl` to adjust the input parameters such that all mentioned features can be simulated will be given in the following.

# Example 1: Graded interfaces
By default, we assume abrupt inner interfaces. If one wishes to simulate graded interfaces, where for example the effective density of states and the band-edge energy may vary, we refer to [this](https://github.com/WIAS-PDELib/ChargeTransport.jl/blob/master/examples/Ex105_PSC_gradedFlux.jl) example.

We sketch the relevant parts here.
First, we need to import the constants and units.
```julia
# unit factors
@local_unitfactors μm cm s ns V K

# constants
constants = ChargeTransport.constants
(; q, k_B, ε_0) = constants

eV = q * V
```

Then, we need to define two additional thin interface layers

```julia

# region numbers
regionDonor = 1      # n doped region
regionJunction1 = 2
regionIntrinsic = 3  # intrinsic region
regionJunction2 = 4
regionAcceptor = 5   # p doped region
```
which need to be taken into account by the initialization of the grid.

Second, since we allow varying parameters within the thin interface layers, the flux discretization scheme needs to be chosen accordingly and we need to construct a nodally dependent parameter struct

```julia
data.fluxApproximation = ScharfetterGummelGraded

paramsnodal = ParamsNodal(grid, numberOfCarriers)
```

Finally, we introduce graded parameters. Currently, only a linear grading is implemented.

```julia
paramsnodal.bandEdgeEnergy[iphin, :] = grading_parameter!(
    paramsnodal.bandEdgeEnergy[iphin, :],
    coord, regionTransportLayers, regionJunctions, h,
    heightLayers, lengthLayers, EC
)
```

## Example 2: Linear IV scan protocol
Here, we summarize the main parts of [this](https://github.com/WIAS-PDELib/ChargeTransport.jl/blob/master/examples/Ex103_PSC_IVMeasurement.jl) example.
Define three charge carriers.
```julia
iphin = 1 # electrons
iphip = 2 # holes
iphia = 3 # anion vacancies
numberOfCarriers = 3
```
Consider the transient problem and enable the ionic charge carriers only in the active layer:
```julia
data.modelType = Transient
enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])
```

Following specification is needed for a linear I-V scan protocol.

```julia
scanrate = 1.0 * V / s
endVoltage = voltageAcceptor # bias goes until the given voltage at acceptor boundary
tend = endVoltage / scanrate
```

### Use internal time stepping
To make use of internal time stepping, the scan protocol need to be previously defined, e.g.

```julia
function linearScanProtocol(t)
    if t == Inf
        0.0
    else
        scanrate * t
    end
end

## Apply zero voltage on left boundary and a linear scan protocol on right boundary
contactVoltageFunction = [zeroVoltage, linearScanProtocol]
```
And then, need to be parsed into the data construction method
```julia
data = Data(grid, numberOfCarriers, contactVoltageFunction = contactVoltageFunction)
```
This makes it possible to use the internal time solving method

```julia
sol = solve(ctsys, inival = inival, times=(0.0, tend), control = control)
```

## Example 3: Illumination
Add uniform photogeneration to the previous code by setting

```julia
data.generationModel = GenerationUniform
```
and specify the uniform generation rate in each region, i.e.

```julia
for ireg in 1:numberOfRegions
    params.generationUniform[ireg] = generationUniform[ireg]
end
```
for given data stored in `generationUniform`.
If one wishes to use the Beer-Lambert generation, then the corresponding code would be
```julia
data.generationModel = GenerationBeerLambert

for ireg in 1:numberOfRegions
    params.generationIncidentPhotonFlux[ireg] = incidentPhotonFlux[ireg]
    params.generationAbsorption[ireg] = absorption[ireg]
end

params.generationPeak = generationPeak
```
If one wishes to invert the illumination, one needs to define
```julia
params.invertedIllumination = -1
```
where this value is by default set to one (for light entering from the left).
Furthermore, we recommend performing a time loop while increasing the generation rate and afterwards applying the scan protocol with a full generation due to numerical stability, see this [example](https://github.com/WIAS-PDELib/ChargeTransport.jl/blob/master/examples/Ex104_PSC_Photogeneration.jl).

## Example 4: Multi-dimensional problems
It is also possible to perform multi-dimensional simulations.

For a 2D mesh you may use a structured grid via [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl) or an unstructured mesh via the Julia wrapper [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) for Jonathan Richard Shewchuk's Triangle mesh generator.
Respective examples can be likewise found within this package.

Lastly, with help of the [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl) wrapper, three-dimensional tetrahedral meshes can be generated, see [this](https://github.com/WIAS-PDELib/ChargeTransport.jl/blob/master/examples/Grid_3D.jl) example.
