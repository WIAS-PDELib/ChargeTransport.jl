"""
    Famous constants in SI units
"""

### Below is a first try with Unitful/PhysicalConstants. Problem: Units need to be erased.

# using Unitful
# using PhysicalConstants.CODATA2018

# kB               = BoltzmannConstant           # JK^{-1}  --- Boltzmann constant
# Planck_constant  = PlanckConstant              # Js       --- Planck constant
# mₑ               = ElectronMass                # kg       --- electron rest mass
# q                = ElementaryCharge            # C        --- elementary charge
# ε0               = VacuumElectricPermittivity  # C/(V*m)  --- absolute dielectric permittivity of classical vacuum

# Famous constants
global k_B = 1.380649e-23 # 1.380649e-23    # JK^{-1}     --- Boltzmann constant
global Planck_constant = 6.62607015e-34    # Js          --- Planck constant
global mₑ = 9.1093837015e-31  # kg          --- electron rest mass
global q = 1.602176634e-19   # C           --- elementary charge
global ε_0 = 8.8541878128e-12  # C/(V*m)     --- absolute dielectric permittivity of classical vacuum

# Numerical parameters
const tiny_penalty_value = 1.0e-10        # tiny penalty value

# Set this as temporary values from pdelib
function set_pdelib_constants()
    global k_B = 1.3806503e-23
    global q = 1.602176462e-19
    global ε_0 = 8.85418781762039e-12
    global eV = q * V
    return nothing
end

function set_TeSCA_constants()
    global k_B = 1.380662e-23
    global q = 1.6021e-19
    global ε_0 = 8.85419e-12
    global eV = q * V
    return nothing
end

# set unity for constants
function set_unity_constants()
    global k_B = 1.0
    global q = 1.0
    global ε_0 = 1.0
    global eV = 1.0
    return nothing
end
