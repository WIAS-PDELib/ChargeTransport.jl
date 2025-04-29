Perovskite solar cell
================================
We simulate charge transport in perovskite solar cells (PSCs), where we have apart from holes and electrons also ionic charge carriers. Here, we assume to have three domains, denoted by
$\mathbf{\Omega} = \mathbf{\Omega}_{\text{HTL}} \cup \mathbf{\Omega}_{\text{intr}} \cup \mathbf{\Omega}_{\text{ETL}}  $.
The unknowns are the quasi Fermi potentials of electrons, holes and anion vacancies
$\varphi_n, \varphi_p, \varphi_a$
as well as the electric potential
$\psi$.
The underlying PDEs are given by
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (p(\psi, \varphi_p) - C_p ) - (n(\psi, \varphi_n) - C_n) \Big),\\
	q \partial_t n(\psi, \varphi_n) - \nabla \cdot \mathbf{j}_n &= q\Bigl(G(\mathbf{x}) - R(n,p) \Bigr), \\
	q \partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr),
\end{aligned}
```
for
$\mathbf{x} \in \mathbf{\Omega}_{\text{HTL}} \cup  \mathbf{\Omega}_{\text{ETL}} $, $t \in [0, t_F]$. In the middle, intrinsic region ($ \mathbf{x} \in \mathbf{\Omega}_{\text{intr}} $), we have
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( p(\psi, \varphi_p)  - n(\psi, \varphi_n) + a(\psi, \varphi_a) - C_a \Big),\\
q \partial_t n(\psi, \varphi_n)	- \nabla \cdot \mathbf{j}_n &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr), \\
	q \partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr),\\
	q \partial_t a(\psi, \varphi_a) + \nabla \cdot \mathbf{j}_a &= 0,
\end{aligned}
```
see [Abdel2021](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865).

Differences to the previous example include
- an additional charge carrier (the anion vacancy)
- parameter jumps across heterojunctions
- the transient case
- a generation rate $G$
- higher dimensional problem.

A quick survey on how to use `ChargeTransport.jl` to adjust the input parameters such that these features can be simulated will be given in the following.

## Example 1: Graded interfaces
By default, we assume abrupt inner interfaces. If one wishes to simulate graded interfaces, where for example the effective density of states and the band-edge energy may vary, we refer to [this](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex104_PSC_gradedFlux_Schottky_contacts.jl) or [this](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex105_PSC_gradedFlux.jl) example.

We sketch the relevant parts here. First, we need to define two additional thin interface layers

```julia
# region numbers
regionDonor     = 1       # n doped region
regionJunction1 = 2
regionIntrinsic = 3       # intrinsic region
regionJunction2 = 4
regionAcceptor  = 5       # p doped region
```
which need to be taken into account by the initialization of the grid.

Second, since we allow varying parameters within the thin interface layers, the flux discretization scheme needs to be chosen accordingly and we need to construct a nodally dependent parameter struct

```julia
data.fluxApproximation = ScharfetterGummelGraded

paramsnodal            = ParamsNodal(grid, numberOfCarriers)
```

Finally, we introduce graded parameters. Currently, only a linear grading is implemented.

```julia
paramsnodal.bandEdgeEnergy[iphin, :] = grading_parameter!(paramsnodal.bandEdgeEnergy[iphin, :],
                                                         coord, regionTransportLayers, regionJunctions,
                                                         h, heightLayers, lengthLayers, EC)
```

## Example 2: Linear IV scan protocol
Here, we summarize the main parts of [this](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex106_PSC_withIons_IVMeasurement.jl) example.
Define three charge carriers.
```julia
iphin                    = 2 # electrons
iphip                    = 1 # holes
iphia                    = 3 # anion vacancies
numberOfCarriers         = 3
```
Consider the transient problem and enable the ionic charge carriers only in the active layer:
```julia
data.modelType           = Transient
enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])
```

Following specification is needed for a linear I-V scan protocol.

```julia
scanrate                 = 1.0 * V/s
number_tsteps            = 31
endVoltage               = voltageAcceptor # bias goes until the given voltage at acceptor boundary
tend                     = endVoltage/scanrate
```

### Variant A: Solve the transient problem manually
```julia
## with fixed timestep sizes we can calculate the times a priori
tvalues                    = range(0, stop = tend, length = number_tsteps)

for istep = 2:number_tsteps

    t  = tvalues[istep]                  # current time
    Δu = t * scanrate                    # applied voltage
    Δt = t - tvalues[istep-1]            # time step
    set_contact!(ctsys, bregionAcceptor, Δu = Δu)
    solution = solve(ctsys, inival = inival, control = control, tstep = Δt) # provide time step
    inival   = solution

end
```

### Variant B: Use internal time stepping
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
Add uniform illumination to the previous code by setting

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
data.generationModel                          = GenerationBeerLambert

for ireg in 1:numberOfRegions
    params.generationIncidentPhotonFlux[ireg] = incidentPhotonFlux[ireg]
    params.generationAbsorption[ireg]         = absorption[ireg]
end

params.generationPeak                         = generationPeak
```
If one wishes to invert the illumination, one needs to define
```julia
params.invertedIllumination                   = -1
```
where this value is by default set to one (for light entering from the left).
Furthermore, we recommend performing a time loop while increasing the generation rate and afterwards applying the scan protocol with a full generation due to numerical stability, see this [example](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex108_PSC_uniform_Generation.jl) or [this one](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex108_PSC_BeerLambert_Generation.jl).

## Example 4: Multi-dimensional problems
It is also possible to perform multi-dimensional simulations.

For a 2D mesh you may use a structured grid via [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl) or an unstructured mesh via the Julia wrapper [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) for Jonathan Richard Shewchuk's Triangle mesh generator.
Respective examples can be likewise found within this package.

Lastly, with help of the [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl) wrapper, three dimensional tetrahedral meshes can be generated, see [this](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Grid_3D.jl) example.
