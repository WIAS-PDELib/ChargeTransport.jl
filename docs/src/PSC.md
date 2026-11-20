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
Further ionic carriers, e.g. cation vacancies, can be added in the same way.

Differences to the previous example include
- an additional charge carrier (the anion vacancy)
- parameter jumps across heterojunctions
- the transient case
- a generation rate $G$
- higher dimensional problem.

# Simulating electronic-ionic charge transport

For perovskite solar cells, an average vacancy density is typically given in literature.
We consider ionic carriers $\alpha$ (e.g. anion vacancies $\alpha = \text{a}$ and possibly cation vacancies $\alpha = \text{c}$) with charge numbers $z_\alpha$, each present in a region $\Omega_\alpha$ (here $\Omega_{\text{PVK}}$).
If the charge carrier densities are used as the set of unknowns, the initial condition can be given by
```math
\begin{aligned}
	n_\alpha^0 = \bar n_\alpha^{\,\text{target}},
\end{aligned}
```
where $\bar n_\alpha^{\,\text{target}}$ is either the background density $C_\alpha$ (e.g. anion vacancies compensated by immobile cation vacancies, to ensure global charge neutrality) or the value given via `averageDensity` in `enable_ionic_carrier!` (e.g. two mobile, compensating vacancy species).
With homogeneous no-flow Neumann boundary conditions around $\Omega_\alpha$, the total mass of each ionic carrier is conserved at all times, i.e.
```math
\begin{aligned}
	\frac{1}{|\Omega_\alpha|} \int_{\Omega_\alpha} n_\alpha (\mathbf{x}, t)\, d\mathbf{x}
    = \frac{1}{|\Omega_\alpha|}  \int_{\Omega_\alpha} n_\alpha^0 (\mathbf{x})\, d\mathbf{x}
    = \bar n_\alpha^{\,\text{target}}, \quad \text{for all} \quad t \geq 0.
\end{aligned}
```
Since we do not use charge carrier densities directly but instead employ quasi Fermi potentials as the unknowns, due to mathematical, physical, and numerical advantages, we must find a workaround to properly fix the initial condition for the vacancy densities.

## Quasi Fermi potential notation

The statistical relation between the vacancy densities and the potentials (our chosen unknowns) reads
```math
\begin{aligned}
    n_\alpha  = N_\alpha  F_{-1} \Bigl(\eta_\alpha (\psi, \varphi_\alpha ) \Bigr), \quad \eta_\alpha  = z_\alpha  \frac{q (\varphi_\alpha  - \psi) + E_\alpha }{k_B T},
\end{aligned}
```
where $N_\alpha$ denotes the maximum vacancy density allowed, $F_{-1} = (\exp(-x) +1)^{-1}$ is the Fermi-Dirac integral of order $-1$, and $E_\alpha$ is the intrinsic vacancy energy level (somehow a model parameter).
In equilibrium, we set the vacancy quasi Fermi potentials to zero, i.e., $\varphi_\alpha = 0$, when the applied voltage is $V = 0$.
Because the average vacancy densities are prescribed and the electric potential $\psi$ is an unknown, the only remaining free parameters are the vacancy energies $E_\alpha$.
> The values of $E_\alpha$ must be chosen such that the average vacancy densities remain conserved for all time steps.

## Equilibrium problem

In equilibrium, all quasi Fermi potentials vanish and the electric potential $\psi$ solves (here, we omit fixed background charges for the sake of simplicity)
```math
-\nabla\cdot\big(\varepsilon_s\,\nabla\psi\big)
= q\Big( n_\text{p}(\psi) - n_\text{n}(\psi) + \sum_{\alpha} z_\alpha\, n_\alpha(\psi, E_\alpha) \Big),
\quad
n_\alpha = N_\alpha\, F_{-1}\!\left( \frac{z_\alpha\,(E_\alpha - q\psi)}{k_BT} \right),
```
and the energies $E_\alpha$ have to satisfy
```math
\frac{1}{|\Omega_\alpha|}\int_{\Omega_\alpha} n_\alpha\, d\mathbf{x} = \bar n_\alpha^{\,\text{target}} .
```

### Are the energies unique?

Yes, as long as $\psi$ is fixed at the contacts (Dirichlet), and interestingly, we do not need charge neutrality for this.

We sketch the proof for anion and cation vacancies, $\alpha \in \{a, c\}$. For this, we take two solutions and set $w = \psi^* - \psi$, $\Delta E_\alpha = E_\alpha^* - E_\alpha$ and $\Delta n_\alpha = n_\alpha^* - n_\alpha$.

*Step 1: multiply by $w$ and integrate.* We subtract the two Poisson equations, multiply by $w$ and integrate. Since $w = 0$ at the contacts, we get
```math
\int_\Omega \varepsilon_s\,|\nabla w|^2 \, d\mathbf{x} = \underbrace{q\int_\Omega (\Delta n_\text{p} - \Delta n_\text{n})\,w \, d\mathbf{x}}_{\le 0} \;+\; \sum_{\alpha} q z_\alpha \int_{\Omega_\alpha} \Delta n_\alpha\, w \, d\mathbf{x}.
```

*Step 2: add a constant.* Since the average densities are fixed, we have $\int_{\Omega_\alpha} \Delta n_\alpha \, d\mathbf{x} = 0$ and we can add a constant to $w$ without changing the ion terms, i.e.
```math
q z_\alpha \int_{\Omega_\alpha} \Delta n_\alpha\, w \, d\mathbf{x} = -\int_{\Omega_\alpha} \Delta n_\alpha\, z_\alpha\big(\Delta E_\alpha - q w\big) \, d\mathbf{x}.
```

*Step 3: sign rule.* Since $F_{-1}$ is increasing, $\Delta n_\alpha$ and $z_\alpha(\Delta E_\alpha - q w)$ have the same sign, and hence every ion term is non-positive. The same holds for electrons and holes, whose energies are fixed anyway.

Thus, we have $\int_\Omega \varepsilon_s |\nabla w|^2 \, d\mathbf{x} \le 0$, which gives us $w = 0$ and then also $\Delta E_a = \Delta E_c = 0$.
The same argument works for an arbitrary number of ionic carriers. Note that with only Neumann boundary conditions, we could shift all energies together with $\psi$ and would get another solution, i.e. we would lose uniqueness.

## Calculation of the energies

**Root-finding problem.** For $m$ ionic carriers, find $\mathbf{E} \in \mathbb{R}^m$ such that $\mathbf{G}(\mathbf{E}) = \mathbf{0}$, where
```math
G_\alpha(\mathbf{E}) = F_{-1}^{-1}\!\left(\frac{\bar n_\alpha(\mathbf{E})}{N_\alpha}\right) - F_{-1}^{-1}\!\left(\frac{\bar n_\alpha^{\,\text{target}}}{N_\alpha}\right),
\qquad
\bar n_\alpha(\mathbf{E}) = \frac{1}{|\Omega_\alpha|}\int_{\Omega_\alpha} n_\alpha\big(\psi[\mathbf{E}], E_\alpha\big)\, d\mathbf{x} .
```
We apply $F_{-1}^{-1}$ here, since this makes $\mathbf{G}$ nearly linear in $\mathbf{E}$.

**Newton's method.** The classical Newton method approximates the $(k+1)$-th iterate via
```math
\mathbf{E}^{k+1} = \mathbf{E}^k - J(\mathbf{E}^k)^{-1}\, \mathbf{G}(\mathbf{E}^k),
\qquad
J_{\alpha\beta} = \frac{\partial G_\alpha}{\partial E_\beta}.
```
However, since $\mathbf{G}$ depends on $\mathbf{E}$ through $\psi[\mathbf{E}]$, we do not know the Jacobian explicitly.

**Broyden's method.** Instead, we approximate the Jacobian by finite differences only once at the start and afterwards update it from the last two iterates. With $\mathbf{s}^k = \mathbf{E}^{k+1} - \mathbf{E}^k$ and $\Delta\mathbf{G}^k = \mathbf{G}(\mathbf{E}^{k+1}) - \mathbf{G}(\mathbf{E}^k)$,
```math
J_{k+1} = J_k + \frac{\big(\Delta\mathbf{G}^k - J_k\, \mathbf{s}^k\big)\, (\mathbf{s}^k)^T}{\|\mathbf{s}^k\|^2},
```
we get $J_{k+1}\,\mathbf{s}^k = \Delta\mathbf{G}^k$, and thus each iteration costs only one PDE solve.
For one ionic carrier ($m = 1$), this reduces to the secant method:
```math
J_{k+1} = \frac{G(E^{k+1}) - G(E^k)}{E^{k+1} - E^k}.
```

**Stopping.** We round all energies to `energyDigits` decimal places (in eV), e.g. to 1 meV for `energyDigits = 3`. The iteration stops once the Newton correction is below this resolution, i.e. when rounding gives back the same energies, $\mathbf{E}^{k+1} = \mathbf{E}^k$.
A stopping criterion on the density would not be practical, since the reachable accuracy of the density is itself limited by the energy resolution and may be not reached.

## Usage

To make use of this feature, add the flag `vacancyEnergyCalculation = true` in the equilibrium solve.
`energyDigits` sets the resolution of the energy grid (default 4, i.e. 0.1 meV) and `verboseEnergy = true` prints the iterations.

```julia
solution = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = true)
```
For ionic carriers without background charge (e.g. two compensating vacancy species), prescribe the mean density directly:
```julia
enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic], averageDensity = [Ca])
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
- For **multi-dimensional simulations**, however, we recommend precomputing the energies $E_\alpha$ and storing them in case of multiple computations with the same parameter set.


Next, we give a quick survey on how to adjust the input parameters in `ChargeTransport.jl` such that all mentioned features can be simulated.

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
# The flux approximation is by default the excess chemical potential scheme.
# In this example, we need a graded scheme.
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
