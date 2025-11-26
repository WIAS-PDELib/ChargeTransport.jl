
Mathematical drift-diffusion models
================================
`ChargeTransport.jl` aims to discretize charge transport models based on drift-diffusion equations. The bipolar case is sometimes referred to as van Roosbroeck system. This nonlinear system of partial differential equations
couples Poisson's equation to several continuity equations. The precise type and amount will vary with the specific application.

In this section, we would like to
describe the mathematical theory a bit more in detail. We denote with $\alpha$ the charge
carrier, with $n_\alpha$ its corresponding density in a device region $\mathbf{\Omega}$ during a
finite time interval $[0, t_F]$.

## Poisson's equation

Poisson's equation for the electric potential $\psi$ is given by
```math
\begin{aligned}
- \nabla \cdot \Bigl(\varepsilon_s \nabla \psi(\mathbf{x}, t) \Bigr) &= q \sum_{\alpha} z_\alpha \Bigl( n_\alpha(\mathbf{x}, t) - C_\alpha(\mathbf{x}) \Bigr).
\end{aligned}
```
Here,
$\varepsilon_s$
denotes the dielectric permittivity and $ q $ the elementary charge. The right-hand side of Poisson's equation, the space charge density, is the sum of charge carrier densities
$n_\alpha$
multiplied by their respective charge numbers
$z_\alpha$
and some corresponding fixed charges, the doping $ C_\alpha $.

## Continuity equations

Poisson's equation is coupled to additional continuity equations for each charge carrier $\alpha$, which describe the motion of free charge carriers in an electric field
```math
\begin{aligned}
z_\alpha q \partial_t n_\alpha +  \nabla\cdot \mathbf{j}_\alpha
	&=
	z_\alpha q	r_\alpha.
\end{aligned}
```
Here, the flux
$\mathbf{j}_\alpha$
refers to the the carrier's current density and $r_\alpha$ to some production/reduction rates.
These rates may be chosen to represent different recombination or generation models such as Shockley-Read-Hall, Auger or direct recombination.

The amount and type of charge carriers will be dependent on the specific application. The standard semiconductor equations use electrons $\alpha=n$ and holes $\alpha=p$.

## Drift-diffusion fluxes
Our code uses as independent variables the electrostatic potential $\psi$ as well as the quasi Fermi
potentials $\varphi_\alpha$. The charge carrier densities $n_\alpha$ are linked to the corresponding quasi Fermi potentials via the state equations
```math
\begin{aligned}
n_\alpha = N_\alpha \mathcal{F}_\alpha \Bigl(\eta_\alpha(\psi, \varphi_\alpha) \Bigr), \quad \eta_\alpha = z_\alpha \frac{q (\varphi_\alpha - \psi) + E_\alpha}{k_B T},
\end{aligned}
```
where the physical parameters are defined [in the list of notations](@ref notation). With this definition we can formulate the carrier current given by
```math
\begin{aligned}
    \mathbf{j}_\alpha
	=
    - (z_\alpha)^2 q \mu_\alpha
    n_\alpha
    \nabla\varphi_\alpha
    ~
\end{aligned}
```
with the negative gradients of the quasi Fermi potentials as driving forces. Using the state equations one may rewrite these fluxes in a drift-diffusion form.

!!! note

    The unknowns in `ChargeTransport.jl` are always defined as the quasi Fermi potentials $ \varphi_\alpha$ and the electric potential $\psi$.

## Boundary conditions
Currently, ohmic contacts, Schottky contacts and Schottky barrier lowering boundary conditions are implemented. For further model information, please look closer to the types, constructors and methods section.

## Background literature

For a comprehensive overview of drift-diffusion models, semiconductor applications as well as the underlying numerical methods, we recommend the following sources:

1. P. Farrell, D. H. Doan, M. Kantner, J. Fuhrmann, T. Koprucki, and N. Rotundo. [“Drift-Diffusion Models”](https://www.taylorfrancis.com/chapters/edit/10.4324/9781315152318-25/drift-diffusion-models-patricio-farrell-nella-rotundo-duy-hai-doan-markus-kantner-j%C3%BCrgen-fuhrmann-thomas-koprucki). In: Optoelectronic Device Modeling and Simulation: Fundamentals, Materials, Nanostructures, LEDs, and Amplifiers. CRC Press Taylor & Francis Group, 2017, pp. 733–771.
2. S. Selberherr. [Analysis and Simulation of Semiconductor Devices](https://link.springer.com/book/10.1007/978-3-7091-8752-4). Springer-Verlag, 1984.
3. S. M. Sze and K. K. Ng. [Physics of Semiconductor Devices](https://onlinelibrary.wiley.com/doi/book/10.1002/0470068329). Wiley, 2006.


# [Notation](@id notation)

| **symbol** | **physical quantity** |   |   |   |   | **symbol** | **physical quantity** |
| :---:         |     :---:      |          :---: |          :---: |          :---: |          :---: |          :---: |          :---: |
| $ \alpha $   | mobile charge carrier     |      |      |      |      | $ n_\alpha $    | charge carrier density of $ \alpha $    |
| $\varepsilon_s$     | dielectric permittivity       |      |      |      |      | $ \psi $      | electrostatic potential      |
| $ q $     | elementary charge       |      |      |      |      | $ C_\alpha $      | doping/background charge      |
| $ z_\alpha $     | charge number for $ \alpha $       |      |      |      |      | $ r_\alpha $     | production/reaction rate for $ \alpha $       |      |      |      |      | $ \mathbf{j}_\alpha $      | current density for $ \alpha $      |
| $ N_\alpha $     | effective density of states for $ \alpha $       |      |      |      |      | $ \mathcal{F}_\alpha $      | statistics function      |
| $ \varphi_\alpha $     | quasi Fermi potential for $ \alpha $       |      |      |      |      | $ E_\alpha $      | band-edge energy for $ \alpha $      |
| $ k_B $     | Boltzmann constant       |      |      |      |      | $ T $      | temperature      |
| $ \mu_\alpha $     | mobility of carrier $ \alpha $      |      |      |      |      |        |        |
